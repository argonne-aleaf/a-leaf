# Copyright © 2026, UChicago Argonne, LLC. All Rights Reserved.
# SPDX-License-Identifier: BSD-3-Clause  (see LICENSE)

# Solver wrappers: solve JuMP models and collect primal/dual outputs.

const _ED_FAILURE_LP_DUMP_COUNTER = Threads.Atomic{Int}(0)
const _ED_FAILURE_LP_DUMP_LIMIT = 3

# CPLEX parameters that only apply to mixed-integer problems. Setting any of them on a
# pure LP makes CPLEX print "CPLEX Error 3003: Not a mixed-integer problem" to the console.
# Both the CPX_PARAM_* and modern CPXPARAM_* spellings are listed. Note CPX_PARAM_EPGAP
# (the relative MIP-gap tolerance) has no "MIP" in its name, so an explicit list is needed.
const _CPLEX_MIP_ONLY_PARAMS = Set([
    "CPX_PARAM_EPGAP", "CPXPARAM_MIP_Tolerances_MIPGap",
    "CPX_PARAM_EPAGAP", "CPXPARAM_MIP_Tolerances_AbsMIPGap",
    "CPX_PARAM_EPINT", "CPXPARAM_MIP_Tolerances_Integrality",
    "CPX_PARAM_MIPDISPLAY", "CPXPARAM_MIP_Display",
    "CPX_PARAM_MIPEMPHASIS", "CPXPARAM_Emphasis_MIP",
    "CPX_PARAM_VARSEL", "CPXPARAM_MIP_Strategy_VariableSelect",
    "CPX_PARAM_STARTALG", "CPXPARAM_MIP_Strategy_StartAlgorithm",
    "CPX_PARAM_NODESEL", "CPXPARAM_MIP_Strategy_NodeSelect",
    "CPX_PARAM_NODEFILEIND", "CPXPARAM_MIP_Strategy_File",
    "CPX_PARAM_FLOWCOVERS", "CPXPARAM_MIP_Cuts_FlowCovers",
    "CPX_PARAM_CUTSFACTOR", "CPXPARAM_MIP_Limits_CutsFactor",
])

"true if the given solver-setting dict directs CPLEX to solve as an LP relaxation,
i.e. an active CPXPARAM_SolutionType (a.k.a. CPX_PARAM_SOLUTIONTYPE) row with value 2
(no-crossover LP). When this is set, CPLEX solves even an integer-formulated model as a
continuous LP and then rejects MIP-only parameters like CPX_PARAM_EPGAP with
'CPLEX Error 3003: Not a mixed-integer problem'. The setting dict is the nested
row-keyed form ({\"Parameter\", \"Flag\", \"Value\"} per row) used by solve_model_GTEP!.
Returns false for any other shape or when the directive is absent/disabled."
function solver_setting_forces_lp_relaxation(solver_setting)
    solver_setting isa AbstractDict || return false
    for (_, row) in solver_setting
        row isa AbstractDict || continue
        param = get(row, "Parameter", nothing)
        (param isa AbstractString) || continue
        pu = uppercase(String(param))
        if pu == "CPXPARAM_SOLUTIONTYPE" || pu == "CPX_PARAM_SOLUTIONTYPE"
            get(row, "Flag", false) == true || continue
            val = get(row, "Value", nothing)
            if val isa Real && Int(round(val)) == 2
                return true
            end
        end
    end
    return false
end

"Set a solver attribute, but skip CPLEX MIP-only parameters when CPLEX will solve the
model as an LP, to avoid the harmless-but-noisy 'CPLEX Error 3003: Not a mixed-integer
problem' console message. A model is treated as LP-solved when it has no integer/binary
variables OR when `lp_relaxation` is true (the solver settings force a no-crossover LP
relaxation via CPXPARAM_SolutionType=2, which strips integrality at solve time even for a
MIP-formulated model). Non-CPLEX solvers and genuine MIP branch-and-bound solves are
unaffected: MIP-only params are still set for them."
function set_optimizer_attribute_guarded(model::JuMP.AbstractModel, param, value; lp_relaxation::Bool=false)
    if (param isa AbstractString) && (String(param) in _CPLEX_MIP_ONLY_PARAMS) &&
       (lp_relaxation || !model_has_integrality(model))
        return
    end
    JuMP.set_optimizer_attribute(model, param, value)
end

# CPLEX's error-channel hygiene lives in ext/ALEAFCplexExt.jl (no-op unless CPLEX is loaded).
silence_cplex_error_channel!(model::JuMP.AbstractModel) = cplex_silence_error_channel!(Val(:CPLEX), model)

function solve_model!(am::Abstract_ALEAF_Model)
    start_time = time()

    "set solver"
    solver_name = am.setting["Solver Setting"]["solver_name"]
    if solver_name  == "CPLEX"
        optimizer = cplex_optimizer_type(Val(:CPLEX))
        if am.setting["Solver Setting"]["1"]["Value"] == false  # CPLEX needs direct mode
            JuMP.set_optimizer(am.model, optimizer)
        end
    elseif solver_name == "GLPK"
        optimizer = GLPK.Optimizer
        if am.setting["Solver Setting"]["solver_direct_mode_flag"] == false
            JuMP.set_optimizer(am.model, optimizer)
        end
    elseif solver_name == "CBC"
        optimizer = Cbc.Optimizer
        if am.setting["Solver Setting"]["solver_direct_mode_flag"] == false
            JuMP.set_optimizer(am.model, optimizer)
        end
    elseif solver_name == "HiGHS"
        optimizer = HiGHS.Optimizer
        if am.setting["Solver Setting"]["solver_direct_mode_flag"] == false
            JuMP.set_optimizer(am.model, optimizer)
        end
    elseif solver_name == "Gurobi"
        optimizer = Gurobi.Optimizer
        if am.setting["Solver Setting"]["solver_direct_mode_flag"] == false
            JuMP.set_optimizer(am.model, optimizer)
        end
    elseif solver_name == "MadNLP"
        optimizer = gpu_optimizer(Val(:MadNLP), am.setting["Solver Setting"])
        if am.setting["Solver Setting"]["solver_direct_mode_flag"] == false
            JuMP.set_optimizer(am.model, optimizer)
        end
    elseif solver_name == "cuOpt"
        optimizer = gpu_optimizer(Val(:cuOpt))
        if am.setting["Solver Setting"]["solver_direct_mode_flag"] == false
            JuMP.set_optimizer(am.model, optimizer)
            for (k, v) in am.setting["Solver Setting"]
                if v isa AbstractDict && get(v, "Flag", false) == true
                    JuMP.set_optimizer_attribute(am.model, v["Parameter"], v["Value"])
                end
            end
        end
    end

    "set solver setting"
    if am.setting["Solver Setting"]["num_solver_setting"] != 0
        solver_setting_list_data = am.setting["Solver Setting"]["solver_setting_list"]
        split_key = split(solver_setting_list_data, ", ")
        solver_setting_list = []
        for idx in 1:am.setting["Solver Setting"]["num_solver_setting"]
            push!(solver_setting_list, split_key[idx])
        end
        lp_relaxation = solver_setting_forces_lp_relaxation(am.setting["Solver Setting"])
        for setting_keys in keys(am.setting["Solver Setting"])
            if setting_keys in solver_setting_list
                set_optimizer_attribute_guarded(am.model, setting_keys, am.setting["Solver Setting"][setting_keys]; lp_relaxation)
            end
        end
    end

    "solve"
    silence_cplex_error_channel!(am.model)
    _, solve_time, solve_bytes_alloc, sec_in_gc = @timed JuMP.optimize!(am.model, )

    "check solution"
    start_time = time()
    if JuMP.termination_status(am.model) == _MOI.OPTIMAL
        result = collect_result_multi_thread(am, solve_time)
        am.solution = result["solution"]
    elseif JuMP.termination_status(am.model) == _MOI.TIME_LIMIT && JuMP.has_values(am.model)
        @aleaf_info "Time limit reached. Suboptimal solution exists."        
        result = collect_result(am, solve_time)        
        am.solution = result["solution"]
    else
        @aleaf_info "**** No solution found! Termination Status:"

        JuMP.set_optimizer_attribute(am.model, "CPX_PARAM_SCRIND", 1)
        JuMP.set_optimizer_attribute(am.model, "CPX_PARAM_NUMERICALEMPHASIS", 1)
        JuMP.set_optimizer_attribute(am.model, "CPX_PARAM_EPMRK", 0.9)
        
        
        _, solve_time, solve_bytes_alloc, sec_in_gc = @timed JuMP.optimize!(am.model, )

        @aleaf_info "**** Resolved Termination Status:"
        result = collect_result(am, solve_time)
    end

    return result

end


function solve_model_GTEP!(JuMP_model::JuMP.AbstractModel, decomp_group::Int, solution_list::Dict{Symbol,<:Any}, solver_setting; iteration::Int=1, PH_flag::Bool=false, terminate_if_error::Bool=true)
    
    solver_name = solver_setting["solver_name"]
    sol_iteration = 1
    if solver_name  == "CPLEX"
        optimizer = cplex_optimizer_type(Val(:CPLEX))
        if solver_setting["1"]["Value"] == false  # CPLEX needs direct mode
            JuMP.set_optimizer(JuMP_model, optimizer)
        end
    elseif solver_name == "GLPK"
        optimizer = GLPK.Optimizer
        if solver_setting["1"]["Value"] == false
            JuMP.set_optimizer(JuMP_model, optimizer)
        end
    elseif solver_name == "CBC"
        optimizer = Cbc.Optimizer
        if solver_setting["1"]["Value"] == false
            JuMP.set_optimizer(JuMP_model, optimizer)
        end
    elseif solver_name == "HiGHS"
        optimizer = HiGHS.Optimizer
        if solver_setting["1"]["Value"] == false
            JuMP.set_optimizer(JuMP_model, optimizer)
        end        
    elseif solver_name == "Gurobi"
        optimizer = Gurobi.Optimizer
        if solver_setting["1"]["Value"] == false
            JuMP.set_optimizer(JuMP_model, optimizer)
        end
    elseif solver_name == "MadNLP"
        optimizer = gpu_optimizer(Val(:MadNLP), solver_setting)
        if solver_setting["1"]["Value"] == false
            JuMP.set_optimizer(JuMP_model, optimizer)
        end
    elseif solver_name == "cuOpt"
        optimizer = gpu_optimizer(Val(:cuOpt))
        if solver_setting["1"]["Value"] == false
            JuMP.set_optimizer(JuMP_model, optimizer)
        end
    end

    # Set solver setting. MadNLP options are ExaModels.Optimizer constructor arguments
    # (applied in gpu_optimizer), not MOI attributes, so skip this attribute loop for it.
    if solver_name != "MadNLP"
        lp_relaxation = solver_setting_forces_lp_relaxation(solver_setting)
        for setting_id in keys(solver_setting)
            if !(setting_id in ["optimizer", "solver_name", "1", "solution_type"])
                if solver_setting[setting_id]["Flag"] == true
                    set_optimizer_attribute_guarded(JuMP_model, solver_setting[setting_id]["Parameter"], solver_setting[setting_id]["Value"]; lp_relaxation)
                end
            end
        end
    end

    silence_cplex_error_channel!(JuMP_model)

    start_time = time()
    solver_exception = nothing
    try
        _, solve_time, solve_bytes_alloc, sec_in_gc = @timed JuMP.optimize!(JuMP_model)
    catch e
        solver_exception = e
        err_str = sprint(showerror, e)
        bt_str = sprint(Base.show_backtrace, catch_backtrace())
        if PH_flag
            println("--- [GTEP PH of decomposition group $decomp_group]\t ERROR; PH iteration: $iteration\n$err_str$bt_str")
        else
            println("[ALEAF LC_GTEP Expansion Model]:\tERROR\n$err_str$bt_str")
        end
        flush(stdout); flush(stderr)
    end
    solve_time = time() - start_time
    status = JuMP.termination_status(JuMP_model)

    if status in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED] || JuMP.has_values(JuMP_model)
        result = collect_result_distributed(JuMP_model, solution_list, solve_time)
        return result 
        
    elseif status == _MOI.TIME_LIMIT
        if JuMP.has_values(JuMP_model)
           
            result = collect_result_distributed(JuMP_model, solution_list, solve_time)
            try 
                gap = round(result["relative_gap"], digits=5)
            catch
                gap = 100.0
            end

            @aleaf_error("[ALEAF LC_GTEP Model]:\tTime Limit Reached. Sub-optimal solution exists; solution_time: $solve_time, gap: $gap")
            return result
        else
            @aleaf_error("[ALEAF LC_GTEP Model]:\tTime Limit Reached without solution; solution_time: $solve_time")
            
            if terminate_if_error == true
                terminate_with_error(; msg="Time Limit Reached without solution. Terminate the program")
            else
                return "No Solution"
            end
        end
    else
        if JuMP.has_values(JuMP_model)
            result = collect_result_distributed(JuMP_model, solution_list, solve_time)
            gap = 0.0
            try 
                gap = round(result["relative_gap"], digits=5)
            catch
            end
            @aleaf_error("[ALEAF LC_GTEP Model]:\tSub-optimal solution exists; solution_time: $solve_time, gap: $gap")

            JuMP.write_to_file(JuMP_model, string("subopt_LC_GTEP_model_",decomp_group,".lp"))
            @aleaf_error("[ALEAF LC_GTEP Model]:\tLP file has generated: $solve_time")

            return result
        else
            worker_id_local = try myid() catch; 0 end
            n_dumped_prev = Threads.atomic_add!(_ED_FAILURE_LP_DUMP_COUNTER, 1)
            lp_path = ""
            if n_dumped_prev < _ED_FAILURE_LP_DUMP_LIMIT
                try
                    lp_path = string("infeas_LC_GTEP_w", worker_id_local, "_dg", decomp_group, "_n", n_dumped_prev + 1, ".lp")
                    JuMP.write_to_file(JuMP_model, lp_path)
                    @info "[ALEAF LC_GTEP Model]:\tFailing LP written to $(lp_path)"
                    flush(stdout); flush(stderr)
                catch lp_err
                    @warn "[ALEAF LC_GTEP Model]:\tFailed to write LP file $(lp_path): $(sprint(showerror, lp_err))"
                    flush(stdout); flush(stderr)
                end
            end

            diag = collect_solver_diagnostics(JuMP_model)
            prim_st = safe_jump_get(() -> JuMP.primal_status(JuMP_model), nothing)
            ex_summary = solver_exception === nothing ? "none" : sprint(showerror, solver_exception)

            @aleaf_error("[ALEAF LC_GTEP Model]:\tFailed to find a solution; worker=$worker_id_local, decomp_group=$decomp_group, solution_time=$solve_time, termination_status=$status, primal_status=$prim_st, raw_status=\"$(diag["raw_status"])\", result_count=$(diag["result_count"]), lp_file=\"$lp_path\", solver_exception=$ex_summary")

            if terminate_if_error == true
                terminate_with_error(; msg="Failed to find a solution. Terminate the program")
            else
                return "No Solution"
            end
        end
    end
end


"use seperate function for SCIM for now. will merge later after testing"
function safe_moi_get(JuMP_model::JuMP.AbstractModel, attr, default)
    try
        return _MOI.get(JuMP_model, attr)
    catch
        return default
    end
end

function safe_jump_get(fn::Function, default)
    try
        return fn()
    catch
        return default
    end
end

"true if the model contains integer/binary variables (a MIP). MIP-only solver queries
(NodeCount, relative gap, objective bound) make CPLEX print 'CPLEX Error 3003: Not a
mixed-integer problem' when called on a pure LP, so callers guard on this."
function model_has_integrality(JuMP_model::JuMP.AbstractModel)
    for (_, S) in JuMP.list_of_constraint_types(JuMP_model)
        if S === _MOI.Integer || S === _MOI.ZeroOne
            return true
        end
    end
    return false
end

function collect_solver_diagnostics(JuMP_model::JuMP.AbstractModel)
    is_mip = model_has_integrality(JuMP_model)
    return Dict{String,Any}(
        "raw_status" => safe_jump_get(() -> JuMP.raw_status(JuMP_model), ""),
        "result_count" => safe_moi_get(JuMP_model, _MOI.ResultCount(), 0),
        "node_count" => is_mip ? safe_moi_get(JuMP_model, _MOI.NodeCount(), missing) : missing,
        "solve_time_sec" => safe_moi_get(JuMP_model, _MOI.SolveTimeSec(), missing),
        "simplex_iterations" => safe_moi_get(JuMP_model, _MOI.SimplexIterations(), missing),
        "barrier_iterations" => safe_moi_get(JuMP_model, _MOI.BarrierIterations(), missing),
        # MIP-only queries; skip on LP to avoid CPLEX Error 3003 console spam.
        "relative_gap" => is_mip ? safe_jump_get(() -> JuMP.relative_gap(JuMP_model), 0.0) : 0.0,
        "objective_bound" => is_mip ? safe_jump_get(() -> JuMP.objective_bound(JuMP_model), missing) : missing,
    )
end

function collect_result_distributed(JuMP_model::JuMP.AbstractModel, solution_list::Dict{Symbol,<:Any}, solve_time)

    result_count = _MOI.get(JuMP_model, _MOI.ResultCount())
    solution = Dict{String,Any}()
    solver_diagnostics = collect_solver_diagnostics(JuMP_model)
    if result_count > 0
        solution = collect_solution_distributed(solution_list)
        gap = solver_diagnostics["relative_gap"]

        result = Dict{String,Any}(
            "result_count" => result_count,
            "optimizer" => JuMP.solver_name(JuMP_model),
            "termination_status" => JuMP.termination_status(JuMP_model),
            "primal_status" => JuMP.primal_status(JuMP_model),
            "dual_status" => JuMP.dual_status(JuMP_model),
            "objective" => JuMP.objective_value(JuMP_model),
            "objective_lb" => (solver_diagnostics["objective_bound"] === missing ? JuMP.objective_value(JuMP_model) : solver_diagnostics["objective_bound"]),
            "solve_time" => solve_time,
            "solution" => solution,
            "relative_gap" => gap,
            "solver_diagnostics" => solver_diagnostics
            )        
    else
        @aleaf_warn "Model has no results, solution dict will be empty"

        result = Dict{String,Any}(
            "result_count" => result_count,
            "optimizer" => JuMP.solver_name(JuMP_model),
            "termination_status" => JuMP.termination_status(JuMP_model),
            "primal_status" => JuMP.primal_status(JuMP_model),
            "dual_status" => JuMP.dual_status(JuMP_model),
            "objective" => 0,
            "objective_lb" => 0,
            "solve_time" => solve_time,
            "solution" => 0,
            "relative_gap" => 0.0,
            "solver_diagnostics" => solver_diagnostics
            )     
    end

    return result
end


"this will collect solution for multiple networks"
function collect_solution_distributed(solution_list::Dict{Symbol,<:Any})

    sol = Dict{String, Any}()
    for solution_category in keys(solution_list)
        sol[string(solution_category)] = Dict{String, Any}()
        for var_idx in keys(solution_list[solution_category])
            sol[string(solution_category)][string(var_idx)] = Dict{String, Any}()
        end
    end
    
    for solution_category in keys(solution_list)
        Threads.@threads for var_idx in collect(keys(solution_list[solution_category]))
            for var in keys(solution_list[solution_category][var_idx])
                sol[string(solution_category)][string(var_idx)][string(var)] = ALEAF.JuMP.value(solution_list[solution_category][var_idx][var])
            end
        end
    end

    return sol
end


function collect_result(am::Abstract_ALEAF_Model, solve_time)

    result_count = _MOI.get(am.model, _MOI.ResultCount())
    solution = Dict{String,Any}()
    solver_diagnostics = collect_solver_diagnostics(am.model)
    if result_count > 0

        solution = collect_solution_multi_thread(am)

        result = Dict{String,Any}(
            "result_count" => result_count,
            "optimizer" => JuMP.solver_name(am.model),
            "termination_status" => JuMP.termination_status(am.model),
            "primal_status" => JuMP.primal_status(am.model),
            "dual_status" => JuMP.dual_status(am.model),
            "objective" => JuMP.objective_value(am.model),
            "objective_lb" => (solver_diagnostics["objective_bound"] === missing ? JuMP.objective_value(am.model) : solver_diagnostics["objective_bound"]),
            "solve_time" => solve_time,
            "solution" => solution,
            "relative_gap" => solver_diagnostics["relative_gap"],
            "solver_diagnostics" => solver_diagnostics
            )        
    else
        @aleaf_warn "Model has no results, solution dict will be empty"

        result = Dict{String,Any}(
            "result_count" => result_count,
            "optimizer" => JuMP.solver_name(am.model),
            "termination_status" => JuMP.termination_status(am.model),
            "primal_status" => JuMP.primal_status(am.model),
            "dual_status" => JuMP.dual_status(am.model),
            "objective" => 0,
            "objective_lb" => 0,
            "solve_time" => solve_time,
            "solution" => 0,
            "relative_gap" => 0.0,
            "solver_diagnostics" => solver_diagnostics
            )     
    end

    return result
end


function collect_result_multi_thread(am::Abstract_ALEAF_Model, solve_time)

    result_count = _MOI.get(am.model, _MOI.ResultCount())
    solution = Dict{String,Any}()
    solver_diagnostics = collect_solver_diagnostics(am.model)
    if result_count > 0

        solution = collect_solution_multi_thread(am)

        # objective_bound is `missing` on an LP solve (collect_solver_diagnostics skips the
        # MIP-only getter); fall back to the optimal objective, which is its own bound. This
        # matches collect_result / collect_result_distributed and reuses the already-collected
        # bound instead of re-querying (the previous `am.data["model_type"] == "LP"` check
        # compared a Dict to a String and was always false).
        objective_lb = solver_diagnostics["objective_bound"] === missing ? JuMP.objective_value(am.model) : solver_diagnostics["objective_bound"]

        result = Dict{String,Any}(
            "result_count" => result_count,
            "optimizer" => JuMP.solver_name(am.model),
            "termination_status" => JuMP.termination_status(am.model),
            "primal_status" => JuMP.primal_status(am.model),
            "dual_status" => JuMP.dual_status(am.model),
            "objective" => JuMP.objective_value(am.model),
            "objective_lb" => objective_lb,
            "solve_time" => solve_time,
            "solution" => solution,
            "relative_gap" => solver_diagnostics["relative_gap"],
            "solver_diagnostics" => solver_diagnostics
            )        
    else
        @aleaf_warn "Model has no results, solution dict will be empty"

        result = Dict{String,Any}(
            "result_count" => result_count,
            "optimizer" => JuMP.solver_name(am.model),
            "termination_status" => JuMP.termination_status(am.model),
            "primal_status" => JuMP.primal_status(am.model),
            "dual_status" => JuMP.dual_status(am.model),
            "objective" => 0,
            "objective_lb" => 0,
            "solve_time" => solve_time,
            "solution" => 0,
            "relative_gap" => 0.0,
            "solver_diagnostics" => solver_diagnostics
            )     
    end

    return result
end


"this will collect solution for multiple networks"
function collect_solution_multi(am::Abstract_ALEAF_Model, genco_id::Int, day_id::Int)
    sol = collect_solution_values(am.sol[:nw][genco_id][day_id])

    for (k,v) in sol
        sol[k] = v
    end

    return sol
end


function collect_solution(am::Abstract_ALEAF_Model)
    sol = collect_solution_values(am.sol)

    for (k,v) in sol["nw"]["$(am.cnw)"]
        sol[k] = v
    end
    delete!(sol, "nw")

    return sol
end


function collect_solution_multi_thread(am::Abstract_ALEAF_Model)
    
    sol = Dict{String, Any}()
    nw_idx = am.cnw
    for solution_category in keys(am.sol[:nw][nw_idx])
        sol[string(solution_category)] = Dict{String, Any}()
        for var_idx in keys(am.sol[:nw][nw_idx][solution_category])
            sol[string(solution_category)][string(var_idx)] = Dict{String, Any}()
        end
    end
    
    for solution_category in keys(am.sol[:nw][nw_idx])
        Threads.@threads for var_idx in collect(keys(am.sol[:nw][nw_idx][solution_category]))
            for var in keys(am.sol[:nw][nw_idx][solution_category][var_idx])
                sol[string(solution_category)][string(var_idx)][string(var)] = ALEAF.JuMP.value(am.sol[:nw][nw_idx][solution_category][var_idx][var])
            end
        end
    end

    return sol

end


