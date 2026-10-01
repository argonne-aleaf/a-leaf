# Copyright © 2026, UChicago Argonne, LLC. All Rights Reserved.
# SPDX-License-Identifier: BSD-3-Clause  (see LICENSE)

# LC_GTEP orchestration: preprocess, build, solve, and report expansion + operation runs.
# Uses shared variables/constraints from src/component.


using XLSX
using DataFrames
using Dates
using JSON
using DelimitedFiles
using Base.Threads: @threads, @spawn, @sync
using Distributed


# Dense-array backing for VRE zdt shapes: dim3 = tech. Replaces per-(y,b,d,h) untyped nested Dicts
# that ballooned memory at nodal scale (~500x smaller as one Array{Float64,3} per [year][bus]).
const VRE_ZDT_TECH_IDX = Dict("csp_shape"=>1,"pv_shape"=>2,"rtpv_shape"=>3,"wind_ons_shape"=>4,"wind_ofs_shape"=>5,"hydro_shape"=>6)

@inline function get_vre_zdt_shape(nw0, y, b, d, h, type_key)
    idx = nw0[:vre_zdt_index]
    return nw0[:vre_aggregated_data_zdt][y][b][idx.d[d], idx.h[h], VRE_ZDT_TECH_IDX[type_key]]
end
@inline function set_vre_zdt_shape!(nw0, y, b, d, h, type_key, val)
    idx = nw0[:vre_zdt_index]
    nw0[:vre_aggregated_data_zdt][y][b][idx.d[d], idx.h[h], VRE_ZDT_TECH_IDX[type_key]] = val
    return val
end


function build_run_LC_GTEP_model(ALEAF_setting::Dict{String,<:Any})
        
    if ALEAF_setting["Simulation Setting"]["run_GTEP_in_parallel_flag"] == true
        execution_status = execute_LC_GTEP_runs_parallel(ALEAF_setting)
    else
        execution_status = execute_LC_GTEP_runs(ALEAF_setting)
    end

    return execution_status

end

function execute_LC_GTEP_runs(ALEAF_setting::Dict{String,<:Any})

    # Pre-processing-------------------------------------------------------------------
    case_set = sort(collect(keys(filter(x->(x.second["Run_Flag"] == true), ALEAF_setting["Simulation Configuration"]))))

    # ALEAF_CASE_ID restricts the run to one Run_Flag case (one case per Julia process). Useful for GPU solvers:
    # only process exit reliably returns cuOpt's GPU/native memory pool to the OS.
    if haskey(ENV, "ALEAF_CASE_ID") && !isempty(strip(ENV["ALEAF_CASE_ID"]))
        target = strip(ENV["ALEAF_CASE_ID"])
        case_set = filter(c -> c == target, case_set)
        isempty(case_set) && @aleaf_warn "ALEAF_CASE_ID=$target is not a Run_Flag case; nothing to run."
    end

    @aleaf_info "Cases to run = $case_set"

    execution_status = "success"

    for case_id in case_set
        sim_status = execute_LC_GTEP_case(ALEAF_setting, case_id)
        if sim_status == "completed"
            @aleaf_info "Finished executing case ID: $case_id"
        else
            @aleaf_info "Failed to execute case ID: $case_id"
            execution_status = "failed"
        end
        # Reclaim the prior case's heap and run finalizers (frees the solver's native/GPU
        # memory) before the next case allocates, so memory doesn't accumulate across cases.
        GC.gc()
    end

    cleanup_multi_round_info_files(ALEAF_setting, case_set)

    return execution_status
end


# The checkpoint JSON must exist during the run (RA handoff/resume/reuse), so it is always written;
# remove it after the run only when its report flag is explicitly false (default-on: kept).
function cleanup_multi_round_info_files(ALEAF_setting::Dict{String,<:Any}, case_set)
    report_enabled(ALEAF_setting, "report_multi_round_summary_json_EXP_flag") && return
    for case_id in case_set
        file_path = joinpath(define_output_path(ALEAF_setting, case_id), "GTEP_multi_round_info.json")
        if isfile(file_path)
            rm(file_path; force=true)
            @aleaf_info "[ALEAF LC_GTEP]: Removed GTEP_multi_round_info.json for case $case_id (report_multi_round_summary_json_EXP_flag = false)"
        end
    end
end


function execute_LC_GTEP_runs_parallel(ALEAF_setting::Dict{String,<:Any})

    # Pre-processing-------------------------------------------------------------------
    case_set = sort(collect(keys(filter(x->(x.second["Run_Flag"] == true), ALEAF_setting["Simulation Configuration"]))))

    # simulation set
    np = nprocs()
    num_workers= length(workers())
    myhost = fetch(@spawnat myid() gethostname())
    i = 1
    nextidx() = (idx=i; i+=1; idx)

    # update total_num_required_simulation
    last_simulation_id = length(case_set)
    @aleaf_info "Total number of simulation batches: $last_simulation_id, Total number of workders: $num_workers"

    # Workers read the shared data folder directly: it is read-only during a run (after the
    # scenario-reduction/rep-day in-memory refactors), so no per-worker copy is needed.

    sim_status_list = Dict()
    
    @sync begin
        for p in workers()
            @async begin
                while true
                    
                    simulation_id = nextidx() 
                    if simulation_id > last_simulation_id
                        break
                    end
                    
                    case_id = case_set[simulation_id]

                    # each worker mutates its own copy of the settings
                    ALEAF_setting_for_worker = deepcopy(ALEAF_setting)

                    # run simulation (reads the shared data folder directly)
                    @aleaf_info "Worker $p will execute case ID: $case_id"
                    sim_status = remotecall_fetch(execute_LC_GTEP_case, p, ALEAF_setting_for_worker, case_id)
                    # free the worker's memory between cases (finalizers reclaim solver/GPU memory)
                    remotecall_fetch(GC.gc, p)
                    if sim_status == "completed"
                        @aleaf_info "Worker $p finished executing case ID: $case_id"
                    else
                        @aleaf_info "Worker $p failed to execute case ID: $case_id"
                    end
                    sim_status_list[string(p, "_", case_id)] = sim_status

                end
            end
        end
    end

    # check overall simulation status
    execution_status = "success"
    for (sim_id, sim_status) in sim_status_list
        if sim_status != "completed"
            execution_status = "failed"
        end
    end

    cleanup_multi_round_info_files(ALEAF_setting, case_set)

    return execution_status

end


# Seconds spent in RA calls made inside the multi-round expansion loop (reported as RA time, not expansion time).
const RA_TIME_IN_MULTI_ROUND = Ref(0.0)

function execute_LC_GTEP_case(ALEAF_setting::Dict{String,<:Any}, case_id)

    try
    
        case_id = parse(Int64, case_id)
        @aleaf_info "Current case = $case_id"

        # Run time 
        LC_GTEP_expansion_solve_time = 0.0
        LC_GTEP_operation_solve_time = 0.0 
        LC_GTEP_RA_solve_time = 0.0

        # Clean up output path
        clean_up_output_path(ALEAF_setting, case_id)

        # Execute and Report EXP
        start_time = time() 
        RA_TIME_IN_MULTI_ROUND[] = 0.0
        ALEAF_model_instance_LC_GTEP_expansion = nothing
        if ALEAF_setting["Simulation Configuration"][string(case_id)]["Run_expansion_flag"] == true
            ALEAF_model_instance_LC_GTEP_expansion = execute_LC_GTEP_expansion_model(ALEAF_setting, case_id)
        end
        # RA calls made inside the multi-round loop count as RA time, not expansion time
        LC_GTEP_RA_solve_time = round(RA_TIME_IN_MULTI_ROUND[], digits=2)
        LC_GTEP_expansion_solve_time = round(time() - start_time - RA_TIME_IN_MULTI_ROUND[], digits=2)

        # Execute and Report OP
        start_time = time()
        if ALEAF_setting["Simulation Configuration"][string(case_id)]["Run_operation_flag"] == true 

            @aleaf_info "============== Start ALEAF LC_GTEP Operational Model Runs"

            if ALEAF_setting["Simulation Configuration"][string(case_id)]["Run_expansion_flag"] == true
                # run operation after solving expansion
                execute_LC_GTEP_operation_model_after_expansion_run(ALEAF_model_instance_LC_GTEP_expansion, ALEAF_setting, case_id)
            else
                # no expansion model was run first
                if ALEAF_setting["Simulation Configuration"][string(case_id)]["Use_predefined_expansion_data_for_OP_flag"] == true
                    # run OP using the specified expansion results
                    execute_LC_GTEP_operation_model_using_predefined_expansion_data(ALEAF_setting, case_id)
                elseif ALEAF_setting["Simulation Configuration"][string(case_id)]["Use_predefined_expansion_data_for_OP_flag"] == false
                    # run OP without considering expansion
                    execute_LC_GTEP_operation_model(ALEAF_setting, case_id)
                end
            end
            @aleaf_info "[ALEAF LC_GTEP Operational Model]: Done."
        end
        LC_GTEP_operation_solve_time = round(time() - start_time, digits=2)

        # Solve RA
        start_time = time()
        if (ALEAF_setting["Simulation Configuration"][string(case_id)]["Run_RA_flag"] == true) 
        
            if (ALEAF_setting["Simulation Configuration"][string(case_id)]["Run_expansion_flag"] == true) && (ALEAF_setting["Planning Design"]["multi_round_solution_process_flag"] == true)
                # no need to re-run RA in multi-round solution process
            else
                if ALEAF_setting["Simulation Configuration"][string(case_id)]["Use_predefined_expansion_data_for_RA_flag"] == false
                    # run RA using the existing network defined in the database
                    RA_Info = execute_ALEAF_RA_model(ALEAF_setting, case_id)
                else
                    # read pre-defined network data
                    pre_defined_network_data_filename = ALEAF_setting["RA Setting"]["predefined_expansion_data_file_name_for_RA"]
                    external_GTEP_multi_round_info = JSON.parse(open(pre_defined_network_data_filename))
                    RA_Info = execute_ALEAF_RA_model_using_external_data(ALEAF_setting, case_id; external_GTEP_multi_round_info)
                end
            end
            LC_GTEP_RA_solve_time += round(time() - start_time, digits=2)
        end

        # report run time
        report_GTEP_run_time(ALEAF_setting, case_id, LC_GTEP_expansion_solve_time, LC_GTEP_operation_solve_time, LC_GTEP_RA_solve_time)
        @aleaf_info "[ALEAF LC_GTEP]: Run time (sec): Expansion: $LC_GTEP_expansion_solve_time, Operation: $LC_GTEP_operation_solve_time, RA: $LC_GTEP_RA_solve_time"
    
    catch
        rethrow()
        return "failed"
    end

    return "completed"
    
end


function execute_LC_GTEP_expansion_model(ALEAF_setting, case_id)

    # Enhanced-hybrid transmission expansion requires B-theta; warn and fall back otherwise.
    if get(ALEAF_setting["Simulation Setting"], "transmission_expansion_hybrid_flag", false) == true &&
        ALEAF_setting["Simulation Setting"]["power_flow_mode_flag"] != "B-theta"
        @aleaf_info "[ALEAF LC_GTEP]: WARNING transmission_expansion_hybrid_flag is set but power_flow_mode_flag = $(ALEAF_setting["Simulation Setting"]["power_flow_mode_flag"]); the hybrid formulation is B-theta only and will be ignored (standard formulation used)."
    end

    let pd = ALEAF_setting["Planning Design"]
        first_y = pd["first_stage_year_value"]; n_st = pd["num_stages_value"]; len_st = pd["num_years_per_stage_value"]
        @aleaf_info "[ALEAF LC_GTEP]: Horizon $(first_y)-$(first_y + n_st * len_st - 1): $n_st stage(s) of $len_st year(s), stage years = $([first_y + (k-1) * len_st for k in 1:n_st]); load-growth base year = $(pd["base_year_value"]), costs discounted to $(pd["dollar_year_value"])"
        stated_end = get(pd, "horizon_end_year_value", nothing)
        if stated_end isa Number && stated_end != first_y + n_st * len_st - 1
            @aleaf_info "[ALEAF LC_GTEP]: WARNING horizon_end_year_value = $stated_end is information only and differs from the derived horizon end $(first_y + n_st * len_st - 1); the derived value is used."
        end
    end

    start_time = time()
    ALEAF_model_instance_LC_GTEP_expansion = Abstract_ALEAF_Model

    if ALEAF_setting["Planning Design"]["multi_round_solution_process_flag"] == false
        ALEAF_model_instance_LC_GTEP_expansion = execute_LC_GTEP_direct(ALEAF_setting, case_id)
    else
        ALEAF_model_instance_LC_GTEP_expansion = execute_LC_GTEP_direct_multi_round(ALEAF_setting, case_id)
    end
    
    # Export expansion results
    LC_GTEP_expansion_solve_time = round(time() - start_time, digits=2)
    @aleaf_info "[ALEAF LC_GTEP]: Export expansion solution"
    export_LC_GTEP_result_EXP(ALEAF_setting, ALEAF_model_instance_LC_GTEP_expansion, LC_GTEP_expansion_solve_time, case_id)

    return ALEAF_model_instance_LC_GTEP_expansion

end


function execute_LC_GTEP_operation_model(ALEAF_setting, case_id)

    @aleaf_info "[ALEAF LC_GTEP Operational Model]: Running operation model without expansion results"
    
    OP_solutions = Dict()
    result_LC_GTEP_expansion = Dict{String,Any}()

    if ALEAF_setting["Simulation Setting"]["run_operation_in_parallel_flag"] == true
        OP_solutions = build_solve_LC_GTEP_operational_model_distributed(ALEAF_setting, result_LC_GTEP_expansion, case_id)
    else
        OP_solutions = build_solve_LC_GTEP_operational_model(ALEAF_setting, result_LC_GTEP_expansion, case_id)
    end

    # report ALEAF_model_instance results ------------------------------------------------
    @aleaf_info "[ALEAF LC_GTEP]: Export operation solution"
    export_LC_GTEP_result_OP(ALEAF_setting, OP_solutions, case_id, result_LC_GTEP_expansion)

    return nothing

end


function execute_LC_GTEP_operation_model_after_expansion_run(ALEAF_model_instance_LC_GTEP_expansion, ALEAF_setting, case_id)

    # ALEAF_model_instance build-------------------------------------------------
    @aleaf_info "[ALEAF LC_GTEP Operational Model]: Running operation model using expansion results from expansion model run"
    
    OP_solutions = Dict()
    result_LC_GTEP_expansion = Dict{String,Any}()
    recorded_investment_decisions = Dict{String,Any}()

    # record expansion decision
    result_LC_GTEP_expansion = ALEAF_model_instance_LC_GTEP_expansion.solution["1"][:nw]["1"]

    # recorded_investment_decisions for scenario generation in OP
    ids_i = [(i) for (i) in get_index(ALEAF_model_instance_LC_GTEP_expansion, :gen_index, 0)]
    ids_y = [(y) for (y) in get_index(ALEAF_model_instance_LC_GTEP_expansion, :planning_stages, 0)]
    ids_n = [(k) for (k) in get_index(ALEAF_model_instance_LC_GTEP_expansion, :bus, 0)]

    for year in ids_y
        recorded_investment_decisions[string(year)] = Dict{String, Any}()
        recorded_investment_decisions[string(year)]["Wind_TotalMW"] = Dict(string(i) => 0.0 for i in ids_n)   
        recorded_investment_decisions[string(year)]["PV_TotalMW"] = Dict(string(i) => 0.0 for i in ids_n) 
        recorded_investment_decisions[string(year)]["RTPV_TotalMW"] = Dict(string(i) => 0.0 for i in ids_n) 

        for i in ids_i
            bus_idx = parameter(ALEAF_model_instance_LC_GTEP_expansion, 0, :gen_index, "bus_idx", i)
            tech_idx = parameter(ALEAF_model_instance_LC_GTEP_expansion, 0, :gen_index, "genco_tech_id", i)    
            EXUNITS = parameter(ALEAF_model_instance_LC_GTEP_expansion, bus_idx, :gen_bus, tech_idx, "EXUNITS")
            profile_type = parameter(ALEAF_model_instance_LC_GTEP_expansion, 0, :gen_index, "Profile_Type", i)
            CAP = parameter(ALEAF_model_instance_LC_GTEP_expansion, bus_idx, :gen_bus, tech_idx, "CAP")

            U_G_iy = result_LC_GTEP_expansion["solution"]["expansion"][string("(", i, ", ", year, ")")]["u_G_iy"]

            # WIND <- onshore wind only
            if profile_type == "wind_ons"
                recorded_investment_decisions[string(year)]["Wind_TotalMW"][string(bus_idx)] += CAP * U_G_iy
            end

            # PV + RTPV, No CSP
            if profile_type == "pv"
                recorded_investment_decisions[string(year)]["PV_TotalMW"][string(bus_idx)] += CAP * U_G_iy
            end

            if profile_type == "rtpv"
                recorded_investment_decisions[string(year)]["RTPV_TotalMW"][string(bus_idx)] += CAP * U_G_iy
            end

        end

        # Precompute per-data-column VRE weights from the expansion reference (which carries the run-bus
        # finest-region map), so the coarse repday selection can re-key exactly without the run network.
        exp_region_map = get(ALEAF_model_instance_LC_GTEP_expansion.ref[:nw][0], :profile_data_region_map, nothing)
        recorded_investment_decisions[string(year)]["__col_weights"] = build_repday_selection_col_weights(
            ALEAF_model_instance_LC_GTEP_expansion.ref[:nw][0][:bus], exp_region_map,
            recorded_investment_decisions[string(year)]["Wind_TotalMW"],
            recorded_investment_decisions[string(year)]["PV_TotalMW"],
            recorded_investment_decisions[string(year)]["RTPV_TotalMW"])
    end

    if ALEAF_setting["Simulation Setting"]["run_operation_in_parallel_flag"] == true
        OP_solutions = build_solve_LC_GTEP_operational_model_distributed(ALEAF_setting, result_LC_GTEP_expansion, case_id; recorded_investment_decisions)
    else
        OP_solutions = build_solve_LC_GTEP_operational_model(ALEAF_setting, result_LC_GTEP_expansion, case_id; recorded_investment_decisions)
    end

    # report ALEAF_model_instance results ------------------------------------------------
    @aleaf_info "[ALEAF LC_GTEP]: Export operation solution"
    export_LC_GTEP_result_OP(ALEAF_setting, OP_solutions, case_id, result_LC_GTEP_expansion)

    return nothing

end


function execute_LC_GTEP_operation_model_using_predefined_expansion_data(ALEAF_setting, case_id)

    # define pre_defined_expansion_result_file
    data_path = joinpath(pwd(), "data", ALEAF_setting["Simulation Setting"]["test_system_name"], "expansion_results")
    pre_defined_expansion_result_filename = ALEAF_setting["Simulation Configuration"][string(case_id)]["predefined_expansion_data_file_name_for_OP"]
    pre_defined_expansion_result_path = joinpath(data_path, pre_defined_expansion_result_filename)

    @aleaf_info "[ALEAF LC_GTEP Operational Model]: Running operation model using predefined expansion data from $pre_defined_expansion_result_path"

    # parameters
    round_id = 1

    # get pre_defined_expansion_result_data
    pre_defined_expansion_result_data = JSON.parse(open(pre_defined_expansion_result_path))
    result_LC_GTEP_expansion = Dict{String,Any}()
    result_LC_GTEP_expansion = pre_defined_expansion_result_data[string(round_id)]["solution"]["1"]["nw"]["1"]

    # recorded_investment_decisions for scenario generation in OP
    recorded_investment_decisions = Dict{String,Any}()
    recorded_investment_decisions["1"] = pre_defined_expansion_result_data[string(round_id)]["recorded_investment_decisions"]   # only for year 1 because we only record the decision of round_id

    # Solve OPERATIONAL Model
    OP_solutions = Dict()
    if ALEAF_setting["Simulation Setting"]["run_operation_in_parallel_flag"] == true
        OP_solutions = build_solve_LC_GTEP_operational_model_distributed(ALEAF_setting, result_LC_GTEP_expansion, case_id; recorded_investment_decisions)
    else
        OP_solutions = build_solve_LC_GTEP_operational_model(ALEAF_setting, result_LC_GTEP_expansion, case_id; recorded_investment_decisions)
    end

    # report ALEAF_model_instance results ------------------------------------------------
    @aleaf_info "[ALEAF LC_GTEP]: Export operation solution"
    export_LC_GTEP_result_OP(ALEAF_setting, OP_solutions, case_id, result_LC_GTEP_expansion)
end


function get_dual_LC_GTEP_expansion_model!(am::Abstract_ALEAF_Model, ALEAF_setting::Dict{String,<:Any}, network_data::Dict{String,<:Any}; current_round_ids_y = [])

    # collect dual values
    JuMP_model = am.model[:nw][am.cnw][1]
    decomp_id = 1
    am.solution["1"][:nw][string(decomp_id)]["solution"]["dual"] = collect_dual_expansion_GTEP!(JuMP_model, am, decomp_id; round_ids_y = current_round_ids_y)

end


function report_GTEP_run_time(ALEAF_setting, case_id, LC_GTEP_expansion_solve_time, LC_GTEP_operation_solve_time, LC_GTEP_RA_solve_time)
    
    output_path = define_output_path(ALEAF_setting, case_id)
    case_name = ALEAF_setting["Simulation Configuration"][string(case_id)]["Case_ID"]
    
    run_time_file_name = joinpath(output_path, string(case_name, "__simulation_run_time.csv"))

    run_time_df = DataFrame(
        "Case_ID" => [case_name],
        "Expansion_Run_Time_s" => [LC_GTEP_expansion_solve_time],
        "Operation_Run_Time_s" => [LC_GTEP_operation_solve_time],
        "RA_Run_Time_s" => [LC_GTEP_RA_solve_time]
    )

    CSV.write(run_time_file_name, run_time_df)
end



function export_LC_GTEP_result_OP(ALEAF_setting::Dict{String,<:Any}, OP_solutions, case_id, result_LC_GTEP_expansion)
    
    # filename: LC_GEP_TestCase_Scenario_Time 
    case = ALEAF_setting["Simulation Setting"]["test_system_name"]
    case_name = ALEAF_setting["Simulation Configuration"][string(case_id)]["Case_ID"]
    time = Dates.format(now(), "yyyy_mm_dd_HH_MM") 

    output_path = define_output_path(ALEAF_setting, case_id)
    filename = string("ALEAF_LC_GTEP_OP_",case,"_",case_name,"_",time,".json")

    # output dict
    ALEAF_LC_GTEP_solution = Dict{String, Any}()

    if ALEAF_setting["Simulation Configuration"][string(case_id)]["Run_operation_flag"] == true
        # distributed path returns reduced aggregates (no per-hour solutions / network_data);
        # serial path still returns full daily_solutions + network_data.
        pu_applied = get(OP_solutions, "pu_applied", false)
        ALEAF_LC_GTEP_solution["operational model result"] = get(OP_solutions, "daily_solutions", Dict{String,Any}())
        ALEAF_LC_GTEP_solution["operational model system reference"] = OP_solutions["system_reference"]
        haskey(OP_solutions, "network_data") && (ALEAF_LC_GTEP_solution["network data"] = OP_solutions["network_data"])
        ALEAF_LC_GTEP_solution["setting"] = OP_solutions["setting"]["1"]
        if pu_applied
            ALEAF_LC_GTEP_solution["pu_applied"] = true
            ALEAF_LC_GTEP_solution["annual_gen_info_OP"] = OP_solutions["annual_gen_info_OP"]
            ALEAF_LC_GTEP_solution["scarcity_totals"] = OP_solutions["scarcity_totals"]
            ALEAF_LC_GTEP_solution["objectives"] = OP_solutions["objectives"]
            ALEAF_LC_GTEP_solution["build_decisions"] = OP_solutions["build_decisions"]
        end
    end

    if ALEAF_setting["Simulation Setting"]["export_model_reference_json_operation_flag"] == true
        stringdata = JSON.json(ALEAF_LC_GTEP_solution)
        open(joinpath(output_path, filename), "w") do f
            write(f, stringdata)
        end
        @aleaf_info string("", joinpath(output_path, filename), ": file saved")
    end

    # report output
    report_LCO_GTEP_result_OP!(ALEAF_setting, ALEAF_LC_GTEP_solution, result_LC_GTEP_expansion)
end


function clean_up_output_path(ALEAF_setting, case_id)

    output_path = define_output_path(ALEAF_setting, case_id)

    # a multi-round continue run either resumes (needs the report scratch + prior outputs) or warns and
    # returns on a completed checkpoint (needs the prior outputs); in both cases keep the dir intact
    if ALEAF_setting["Planning Design"]["multi_round_solution_process_flag"] == true &&
       ALEAF_setting["Planning Design"]["continue_from_previous_run_flag"] == true
        ckpt = joinpath(output_path, "GTEP_multi_round_info.json")
        if isfile(ckpt)
            try
                if JSON.parse(open(ckpt))["status"] in ("in progress", "completed")
                    return
                end
            catch
            end
        end
    end

    # remove existing files
    original_path = pwd()
    cd(output_path)
    remove_path = []
    for sce_id in 1:length(readdir())
        file_path = readdir(join=true)[sce_id]
        if !(endswith(file_path, ".lp") || endswith(file_path, ".json"))
            push!(remove_path, file_path)
        end
    end
    for remove_path_id in remove_path
        rm(remove_path_id; force=true, recursive=true)
    end
    cd(original_path)

end

function export_LC_GTEP_result_EXP(ALEAF_setting::Dict{String,<:Any}, am, LC_GTEP_expansion_solve_time, case_id)
    
    # filename: LC_GEP_TestCase_Scenario_Time 
    case = ALEAF_setting["Simulation Setting"]["test_system_name"]
    case_name = ALEAF_setting["Simulation Configuration"][string(case_id)]["Case_ID"]
    time = Dates.format(now(), "yyyy_mm_dd_HH_MM") 

    if ALEAF_setting["Simulation Configuration"][string(case_id)]["Run_expansion_flag"] == true
        output_path = am.data["output_path"]
    end

    filename = string("ALEAF_LC_GTEP_EXP_",case,"_",case_name,"_",time,".json")

    # output dict
    ALEAF_LC_GTEP_solution = Dict{String, Any}()

    if ALEAF_setting["Simulation Configuration"][string(case_id)]["Run_expansion_flag"] == true
        
        ALEAF_LC_GTEP_solution["expansion model result"] = am.solution["1"][:nw]
        ALEAF_LC_GTEP_solution["expansion model system reference"] = am.ref[:nw]
        ALEAF_LC_GTEP_solution["network data"] = am.data
        ALEAF_LC_GTEP_solution["setting"] = am.setting

        ALEAF_LC_GTEP_solution["setting"]["RA Setting"] = ALEAF_setting["RA Setting"]

    end

    # empty aleaf model instance to save memory
    am = Dict()
    
    if ALEAF_setting["Simulation Setting"]["export_model_reference_json_expansion_flag"] == true
        stringdata = JSON.json(ALEAF_LC_GTEP_solution)
        open(joinpath(output_path, filename), "w") do f
            write(f, stringdata)
        end
        @aleaf_info string("", joinpath(output_path, filename), ": file saved")
    end

    # report output; multi-round does its reporting per round inside the driver (converts once)
    if ALEAF_setting["Planning Design"]["multi_round_solution_process_flag"] == false
        report_LCO_GTEP_result_EXP!(ALEAF_setting, ALEAF_LC_GTEP_solution)
    end

end


# Report toggles live in the `Simulation Setting` sheet; a missing or blank cell means
# "write it", so existing workbooks keep their current all-reports behavior.
function report_enabled(ALEAF_setting, flag_name)
    v = get(ALEAF_setting["Simulation Setting"], flag_name, true)
    (v === missing || v === nothing) && return true
    # tolerate native booleans and text cells; only an explicit false disables the report
    v isa AbstractString && return !(lowercase(strip(v)) in ("false", "0", "no", "n", "f"))
    return v != false
end


function report_LCO_GTEP_result_EXP!(ALEAF_setting::Dict{String,<:Any}, result::Dict{String, Any})
    
    # apply per unit
    apply_pu_expansion_result_GTEP!(ALEAF_setting, result)

    if report_enabled(ALEAF_setting, "report_expansion_flag")
        report_result_expansion_gen_GTEP(result["expansion model result"], result["expansion model system reference"]; setting=result["setting"])
        report_result_expansion_line_GTEP(result["expansion model result"], result["expansion model system reference"], result["setting"])
    end

    if report_enabled(ALEAF_setting, "report_scarcity_EXP_flag")
        report_result_scarcity_ens_GTEP(result["expansion model result"], result["expansion model system reference"], result["setting"])
        report_result_scarcity_reserve_GTEP(result["expansion model result"], result["expansion model system reference"], result["setting"])
    end

    if report_enabled(ALEAF_setting, "report_power_flow_EXP_flag")
        report_result_power_flow_expansion_GTEP(result["expansion model result"], result["expansion model system reference"], result["setting"])
    end

    dispatch_EXP_enabled = report_enabled(ALEAF_setting, "report_dispatch_EXP_flag")
    summary_EXP_enabled = report_enabled(ALEAF_setting, "report_summary_EXP_flag")

    if dispatch_EXP_enabled
        report_result_dispatch_expansion_GTEP(result["expansion model result"], result["expansion model system reference"], result["setting"])
    elseif summary_EXP_enabled
        # summaries need annual_gen_info + annual_scarcity_info, both in-memory side effects of the
        # dispatch pass; write_csv=false keeps __market_EXP.csv off when dispatch reporting is off.
        report_result_dispatch_expansion_GTEP(result["expansion model result"], result["expansion model system reference"], result["setting"]; write_csv=false)
    end

    if summary_EXP_enabled
        report_result_tech_summary_expansion_GTEP(result["expansion model result"], result["expansion model system reference"], result["setting"])
        report_result_system_summary_expansion_GTEP(result["expansion model result"], result["expansion model system reference"], result["setting"])
    end

    report_result_repdays_GTEP(result["expansion model system reference"], result["setting"])
        
end


function report_LCO_GTEP_result_OP!(ALEAF_setting::Dict{String,<:Any}, result::Dict{String, Any}, result_LC_GTEP_expansion::Dict{String, Any})

    dispatch_OP_enabled = report_enabled(ALEAF_setting, "report_dispatch_OP_flag")
    summary_OP_enabled  = report_enabled(ALEAF_setting, "report_summary_OP_flag")
    power_flow_OP_enabled = report_enabled(ALEAF_setting, "report_power_flow_OP_flag")

    # Distributed path: workers already applied pu and streamed their own per-day-group dispatch /
    # market / policy / demand / power-flow files. Skip the master re-doing that work; instead feed the
    # summaries the reduced aggregates and concatenate the per-dg part files into the canonical outputs.
    if get(result, "pu_applied", false) == true

        # summaries read annual_gen_info_OP off the system reference
        result["operational model system reference"]["annual_gen_info_OP"] = result["annual_gen_info_OP"]

        if summary_OP_enabled
            report_result_tech_summary_operation_GTEP(result, result_LC_GTEP_expansion)
            report_result_system_summary_operation_GTEP(result, result_LC_GTEP_expansion)
        end

        # merge the per-dg worker part files into the single canonical files downstream expects
        concatenate_distributed_OP_files!(result, ALEAF_setting; dispatch=dispatch_OP_enabled, market=(dispatch_OP_enabled || summary_OP_enabled), power_flow=power_flow_OP_enabled)

        report_result_repdays_OP_GTEP(result["operational model system reference"], result["setting"])
        return nothing
    end

    # -----------------------------------------------------------------------------
    # Serial path (unchanged): full solutions in memory
    # -----------------------------------------------------------------------------
    # apply per unit
    apply_pu_operation_result_GTEP!(ALEAF_setting, result)

    if dispatch_OP_enabled
        report_result_dispatch_operation_GTEP(result)
    elseif summary_OP_enabled
        # summaries read annual_gen_info_OP and __market_OP.csv (dispatch-report side effects)
        report_result_dispatch_operation_GTEP(result; generate_file=false, market_only=true)
    end

    if power_flow_OP_enabled
        report_result_power_flow_operation_GTEP(result, result_LC_GTEP_expansion)
    end

    if summary_OP_enabled
        report_result_tech_summary_operation_GTEP(result, result_LC_GTEP_expansion)
        report_result_system_summary_operation_GTEP(result, result_LC_GTEP_expansion)
    end

    report_result_repdays_OP_GTEP(result["operational model system reference"], result["setting"])
end


# Merge each worker's per-day-group part files (suffix "__dg<id>") into the single canonical file
# downstream expects. Streams line-by-line (never loads a whole file), writes the header once from the
# first available part, appends the remaining parts skipping their header, then deletes the parts.
function _concat_parts!(final_path::String, part_paths::Vector{String})
    existing = filter(isfile, part_paths)
    isempty(existing) && return
    open(final_path, "w") do out
        first = true
        for p in existing
            open(p, "r") do inp
                if first
                    write(out, read(inp))   # header + rows from the first part
                    first = false
                else
                    readline(inp)           # drop header
                    while !eof(inp)
                        write(out, readline(inp)); write(out, '\n')
                    end
                end
            end
            rm(p; force=true)
        end
    end
end

# Collect the per-(year,dg) part files for one canonical base: files named "<case><base>__y<..>_dg<..>.csv".
# Matches on the "<base>__y" prefix so the canonical "<base>.csv" is never picked up as its own part.
function _find_op_parts(output_path::String, case_name::String, base::String)
    stem = string(case_name, base, "__y")
    parts = String[]
    for f in readdir(output_path)
        (startswith(f, stem) && endswith(f, ".csv")) && push!(parts, joinpath(output_path, f))
    end
    sort!(parts)
    return parts
end

function concatenate_distributed_OP_files!(result::Dict{String,Any}, ALEAF_setting::Dict{String,<:Any}; dispatch::Bool, market::Bool, power_flow::Bool)
    reference = result["operational model system reference"]
    setting   = result["setting"]
    output_path = parameter(reference["1"], 0, :output_path)
    case_name   = parameter(reference["1"], 0, :case_name)

    ids_y = try
        sort!(parse.(Int, collect(get_index(reference["1"], :planning_stages, 0))))
    catch
        sort!(collect(get_index(reference["1"], :planning_stages, 0)))
    end

    final(base) = joinpath(output_path, string(case_name, base, ".csv"))
    do_base(base) = _concat_parts!(final(base), _find_op_parts(output_path, case_name, base))

    if dispatch
        for y in ids_y
            do_base(string("__dispatch_OP_year_", y))
        end
        do_base("__policy_slack_OP")
        do_base("__demand_response_OP")
    end
    if market
        do_base("__market_OP")
    end
    if power_flow
        do_base("__power_flow_OP")
    end
end


# Scale one day-group's solution dict in place; shared by the serial reporter and each
# distributed worker so pu conversion is identical regardless of where reporting happens.
function apply_pu_to_solution!(sol::Dict, pu_power_base, pu_econ_base)

    # power flow
    if haskey(sol, "powerflow")
        for idx in keys(sol["powerflow"])
            if haskey(sol["powerflow"][idx], "f_kdhty")
                sol["powerflow"][idx]["f_kdhty"] *= pu_power_base
            end
        end
    end

    # storage
    for idx in keys(sol["expansion"])
        if haskey(sol["expansion"][idx], "u_ESE_iy") # this unit is storage
            sol["expansion"][idx]["u_ESE_iy"] *= pu_power_base
        end
    end

    # scarcity
    for idx in keys(sol["scarcity"])
        sol["scarcity"][idx]["ens_ndhty"] *= pu_power_base
    end

    # storage
    for idx in keys(sol["storage"])
        sol["storage"][idx]["chg_idhty"] *= pu_power_base
        sol["storage"][idx]["soc_idhty"] *= pu_power_base
    end

    # scarcity reserve
    for idx in keys(sol["scarcity_reserve"])
        sr = sol["scarcity_reserve"][idx]
        if haskey(sr, "rns_cont_zdhty");      sr["rns_cont_zdhty"]      *= pu_power_base end
        if haskey(sr, "demand_reserve_zdhty");sr["demand_reserve_zdhty"]*= pu_power_base end
        if haskey(sr, "rns_flex_up_zdhty"); sr["rns_flex_up_zdhty"] *= pu_power_base end
        if haskey(sr, "rns_flex_dn_zdhty"); sr["rns_flex_dn_zdhty"] *= pu_power_base end
    end

    # dispatch
    for idx in keys(sol["dispatch"])
        sol["dispatch"][idx]["g_idhty"] *= pu_power_base
    end

    # reserve
    if haskey(sol, "reserve")
        for idx in keys(sol["reserve"])
            if haskey(sol["reserve"][idx], "spin_idhty")
                sol["reserve"][idx]["spin_idhty"] *= pu_power_base
            end
            if haskey(sol["reserve"][idx], "flex_up_idhty")
                sol["reserve"][idx]["flex_up_idhty"] *= pu_power_base
            end
            if haskey(sol["reserve"][idx], "flex_dn_idhty")
                sol["reserve"][idx]["flex_dn_idhty"] *= pu_power_base
            end
            if haskey(sol["reserve"][idx], "reg_dn_idhty")
                sol["reserve"][idx]["reg_dn_idhty"] *= pu_power_base
            end
            if haskey(sol["reserve"][idx], "reg_up_idhty")
                sol["reserve"][idx]["reg_up_idhty"] *= pu_power_base
            end
            if haskey(sol["reserve"][idx], "nonspin_idhty")
                sol["reserve"][idx]["nonspin_idhty"] *= pu_power_base
            end
        end
    end

    # demand response
    # --> will be applied when reporting due to index issues

    # hybrid
    if haskey(sol, "hybrid")
        for idx in keys(sol["hybrid"])
            sol["hybrid"][idx]["g_G_ES_idhty"] *= pu_power_base
        end
    end

    # duals
    for idx in keys(sol["dual"])
        d = sol["dual"][idx]
        haskey(d, "LMP")            && (d["LMP"]            *= pu_econ_base)
        haskey(d, "Reg_Up_Price")   && (d["Reg_Up_Price"]   *= pu_econ_base)
        haskey(d, "Reg_Dn_Price")   && (d["Reg_Dn_Price"]   *= pu_econ_base)
        haskey(d, "Cont_Res_Price") && (d["Cont_Res_Price"] *= pu_econ_base)
        haskey(d, "Flex_Up_Price")  && (d["Flex_Up_Price"]  *= pu_econ_base)
        haskey(d, "Flex_Dn_Price")  && (d["Flex_Dn_Price"]  *= pu_econ_base)
        haskey(d, "CERT_MC")        && (d["CERT_MC"]        *= pu_econ_base)
        haskey(d, "REC_MC")         && (d["REC_MC"]         *= pu_econ_base)
        haskey(d, "CES_MC")         && (d["CES_MC"]         *= pu_econ_base)
    end
end


function apply_pu_operation_result_GTEP!(ALEAF_setting::Dict{String,<:Any}, result::Dict{String, Any})

    pu_power_base = ALEAF_setting["Simulation Setting"]["per_unit_base_value"]
    pu_econ_base = ALEAF_setting["Simulation Setting"]["per_unit_econ_base_value"] / pu_power_base    # per unit economic base

    for y in keys(result["operational model result"])
        for d in keys(result["operational model result"][y])
            apply_pu_to_solution!(result["operational model result"][y][d]["solution"], pu_power_base, pu_econ_base)
        end
    end
end


function apply_pu_expansion_result_GTEP!(ALEAF_setting::Dict{String,<:Any}, result::Dict{String, Any})

    pu_power_base = ALEAF_setting["Simulation Setting"]["per_unit_base_value"]
    pu_econ_base = ALEAF_setting["Simulation Setting"]["per_unit_econ_base_value"] / pu_power_base    # per unit economic base

    ids_decomp_group = [1]

    for decomp_group in ids_decomp_group

        # power flow
        if haskey(result["expansion model result"][string(decomp_group)]["solution"], "powerflow")
            for idx in keys(result["expansion model result"][string(decomp_group)]["solution"]["powerflow"])
                pf = result["expansion model result"][string(decomp_group)]["solution"]["powerflow"][idx]
                haskey(pf, "f_kdhty") && (pf["f_kdhty"] *= pu_power_base)
                # enhanced-hybrid increment: keep f_exp in the same MW units as f_kdhty
                haskey(pf, "f_exp_kdhty") && (pf["f_exp_kdhty"] *= pu_power_base)
            end
        end

        # storage
        for idx in keys(result["expansion model result"][string(decomp_group)]["solution"]["expansion"])
            if haskey(result["expansion model result"][string(decomp_group)]["solution"]["expansion"][idx], "u_ESE_iy") # this unit is storage
                result["expansion model result"][string(decomp_group)]["solution"]["expansion"][idx]["u_ESE_iy"] *= pu_power_base 
            end
        end

        # scarcity
        for idx in keys(result["expansion model result"][string(decomp_group)]["solution"]["scarcity"])
            result["expansion model result"][string(decomp_group)]["solution"]["scarcity"][idx]["ens_ndhty"] *= pu_power_base
        end

        # storage
        for idx in keys(result["expansion model result"][string(decomp_group)]["solution"]["storage"])
            result["expansion model result"][string(decomp_group)]["solution"]["storage"][idx]["chg_idhty"] *= pu_power_base
            result["expansion model result"][string(decomp_group)]["solution"]["storage"][idx]["soc_idhty"] *= pu_power_base
        end

        # scarcity reserve
        for idx in keys(result["expansion model result"][string(decomp_group)]["solution"]["scarcity_reserve"])
            sr = result["expansion model result"][string(decomp_group)]["solution"]["scarcity_reserve"][idx]
            if haskey(sr, "rns_cont_zdhty");      sr["rns_cont_zdhty"]      *= pu_power_base end
            if haskey(sr, "demand_reserve_zdhty");sr["demand_reserve_zdhty"]*= pu_power_base end
            if haskey(sr, "rns_flex_up_zdhty"); sr["rns_flex_up_zdhty"] *= pu_power_base end
            if haskey(sr, "rns_flex_dn_zdhty"); sr["rns_flex_dn_zdhty"] *= pu_power_base end
        end

        # dispatch
        for idx in keys(result["expansion model result"][string(decomp_group)]["solution"]["dispatch"])
            result["expansion model result"][string(decomp_group)]["solution"]["dispatch"][idx]["g_idhty"] *= pu_power_base
        end

        # reserve
        if haskey(result["expansion model result"][string(decomp_group)]["solution"], "reserve")
            for idx in keys(result["expansion model result"][string(decomp_group)]["solution"]["reserve"])
                if haskey(result["expansion model result"][string(decomp_group)]["solution"]["reserve"][idx], "spin_idhty")
                    result["expansion model result"][string(decomp_group)]["solution"]["reserve"][idx]["spin_idhty"] *= pu_power_base
                end
                if haskey(result["expansion model result"][string(decomp_group)]["solution"]["reserve"][idx], "flex_up_idhty")
                    result["expansion model result"][string(decomp_group)]["solution"]["reserve"][idx]["flex_up_idhty"] *= pu_power_base
                end
                if haskey(result["expansion model result"][string(decomp_group)]["solution"]["reserve"][idx], "flex_dn_idhty")
                    result["expansion model result"][string(decomp_group)]["solution"]["reserve"][idx]["flex_dn_idhty"] *= pu_power_base
                end
                if haskey(result["expansion model result"][string(decomp_group)]["solution"]["reserve"][idx], "reg_dn_idhty")
                    result["expansion model result"][string(decomp_group)]["solution"]["reserve"][idx]["reg_dn_idhty"] *= pu_power_base
                end
                if haskey(result["expansion model result"][string(decomp_group)]["solution"]["reserve"][idx], "reg_up_idhty")
                    result["expansion model result"][string(decomp_group)]["solution"]["reserve"][idx]["reg_up_idhty"] *= pu_power_base
                end
                if haskey(result["expansion model result"][string(decomp_group)]["solution"]["reserve"][idx], "nonspin_idhty")
                    result["expansion model result"][string(decomp_group)]["solution"]["reserve"][idx]["nonspin_idhty"] *= pu_power_base
                end
            end
        end

        # demand response
        # --> will be applied when reporting due to index issues

        # hybrid 
        if haskey(result["expansion model result"][string(decomp_group)]["solution"], "hybrid")
            for idx in keys(result["expansion model result"][string(decomp_group)]["solution"]["hybrid"])
                result["expansion model result"][string(decomp_group)]["solution"]["hybrid"][idx]["g_G_ES_idhty"] *= pu_power_base
            end
        end

        # duals
        for idx in keys(result["expansion model result"][string(decomp_group)]["solution"]["dual"])
            if haskey(result["expansion model result"][string(decomp_group)]["solution"]["dual"][idx], "LMP")
                result["expansion model result"][string(decomp_group)]["solution"]["dual"][idx]["LMP"] *= pu_econ_base
            end

            if haskey(result["expansion model result"][string(decomp_group)]["solution"]["dual"][idx], "Reg_Up_Price")
                result["expansion model result"][string(decomp_group)]["solution"]["dual"][idx]["Reg_Up_Price"] *= pu_econ_base                 
            end

            if haskey(result["expansion model result"][string(decomp_group)]["solution"]["dual"][idx], "Reg_Dn_Price")
                result["expansion model result"][string(decomp_group)]["solution"]["dual"][idx]["Reg_Dn_Price"] *= pu_econ_base 
            end

            if haskey(result["expansion model result"][string(decomp_group)]["solution"]["dual"][idx], "Cont_Res_Price")
                result["expansion model result"][string(decomp_group)]["solution"]["dual"][idx]["Cont_Res_Price"] *= pu_econ_base 
            end

            if haskey(result["expansion model result"][string(decomp_group)]["solution"]["dual"][idx], "Flex_Up_Price")
                result["expansion model result"][string(decomp_group)]["solution"]["dual"][idx]["Flex_Up_Price"] *= pu_econ_base 
            end

            if haskey(result["expansion model result"][string(decomp_group)]["solution"]["dual"][idx], "Flex_Dn_Price")
                result["expansion model result"][string(decomp_group)]["solution"]["dual"][idx]["Flex_Dn_Price"] *= pu_econ_base 
            end
            
            if haskey(result["expansion model result"][string(decomp_group)]["solution"]["dual"][idx], "CERT_MC")
                result["expansion model result"][string(decomp_group)]["solution"]["dual"][idx]["CERT_MC"] *= pu_econ_base 
            end

            if haskey(result["expansion model result"][string(decomp_group)]["solution"]["dual"][idx], "REC_MC")
                result["expansion model result"][string(decomp_group)]["solution"]["dual"][idx]["REC_MC"] *= pu_econ_base 
            end

            if haskey(result["expansion model result"][string(decomp_group)]["solution"]["dual"][idx], "CES_MC")
                result["expansion model result"][string(decomp_group)]["solution"]["dual"][idx]["CES_MC"] *= pu_econ_base 
            end
        
        end
    end


end



function get_annual_gen_info(reference)
    try
        return reference["0"]["annual_gen_info"]
    catch
        return reference[0][:annual_gen_info]
    end
end


function get_annual_scarcity_info(reference)
    try
        return reference["0"]["annual_scarcity_info"]
    catch
        return reference[0][:annual_scarcity_info]
    end
end


# ---------------------------------------------------------------------------------------------
# Column names of the hourly / stage report files (one place, so EXP and OP files stay identical)
# ---------------------------------------------------------------------------------------------
const REPORT_LABEL_RENAMES = Dict{String,String}(
    # identifiers and time
    "Scenario" => "Case_ID", "year" => "Stage", "day" => "Rep_Day", "hour" => "Hour", "time" => "Sub_Period",
    "Stochastic_scenario_ID" => "Stochastic_Scenario_ID", "NumDays" => "Days_Represented",
    "unit_id" => "Unit_ID", "PLANT_NAME" => "Plant_Name", "UnitGroup" => "Unit_Group",
    "bus_id" => "Bus_ID", "bus_idx" => "Bus_ID", "bus_name" => "Bus_Name", "node" => "Bus_ID", "zone" => "Reserve_Zone_ID",
    # unit dispatch
    "u_G_iy" => "Units_In_Service", "u_ESE_iy" => "Storage_Energy_MWh", "ICAP" => "Installed_Capacity_MW",
    "sto_c_idhty" => "Storage_Charging_Flag", "c_idhty" => "Units_Committed", "su_idhty" => "Units_Started",
    "g_idhty" => "Generation_MW", "reg_up_idhty" => "Reserve_RegUp_MW", "reg_dn_idhty" => "Reserve_RegDn_MW",
    "spin_idhty" => "Reserve_Spin_MW", "nonspin_idhty" => "Reserve_NSpin_MW",
    "flex_up_idhty" => "Reserve_FlexUp_MW", "flex_dn_idhty" => "Reserve_FlexDn_MW",
    "curt_idhty" => "Curtailment_MW", "chg_idhty" => "Charge_MW", "soc_idhty" => "SOC_MWh",
    "hybrid_chg_idhty" => "Hybrid_Charge_MW", "hybrid_type" => "Hybrid_Type", "Fuel_Type" => "Fuel",
    "Fuel_Cost" => "Fuel_Cost_USD", "MC" => "Marginal_Cost_USD_per_MWh", "vre_shape" => "VRE_Capacity_Factor",
    # market
    "load" => "Load_MW", "scarcity_E" => "Unserved_Energy_MW", "scarcity_SPIN" => "Spin_Shortfall_MW",
    "Demand_Reserve" => "Demand_Reserve_MW", "scarcity_NSPIN" => "NSpin_Shortfall_MW",
    "scarcity_FU" => "FlexUp_Shortfall_MW", "scarcity_FD" => "FlexDown_Shortfall_MW",
    "LMP" => "LMP_USD_per_MWh", "RCP_RU" => "Price_RegUp_USD_per_MW", "RCP_RD" => "Price_RegDn_USD_per_MW",
    "RCP_Spin" => "Price_Spin_USD_per_MW", "RCP_NSpin" => "Price_NSpin_USD_per_MW",
    "RCP_FU" => "Price_FlexUp_USD_per_MW", "RCP_FD" => "Price_FlexDn_USD_per_MW", "ens_indicator_ndhty" => "ENS_Flag",
    # policy and shortfall files
    "slack_CEG_ndy" => "Clean_Energy_Target_Slack", "slack_RPS_ny" => "RPS_Target_Slack",
    "ens_ndhty" => "Unserved_Energy_MW", "rns_spin_zdhty" => "Spin_Shortfall_MW", "demand_reserve_zdhty" => "Demand_Reserve_MW",
    "rns_nonspin_zdhty" => "NSpin_Shortfall_MW", "rns_flex_up_zdhty" => "FlexUp_Shortfall_MW", "rns_flex_dn_zdhty" => "FlexDown_Shortfall_MW",
    # expansion, lines and power flow
    "u_new_G_iy" => "New_Units", "u_ret_G_iy" => "Retired_Units", "u_new_ESH_iy" => "New_Storage_Unit_Hours",
    "line_id" => "Line_ID", "f_bus" => "From_Bus_ID", "t_bus" => "To_Bus_ID", "f_bus_name" => "From_Bus_Name", "t_bus_name" => "To_Bus_Name",
    "f_region" => "From_Region", "t_region" => "To_Region", "merged_line_UIDs" => "Merged_Line_UIDs",
    "orinal_rate" => "Original_Rating_MW", "original_rate" => "Original_Rating_MW", "final_rate" => "Final_Rating_MW",
    "expansion" => "Cumulative_Expansion_Fraction", "u_new_T_ky" => "New_Expansion_Fraction", "u_T_ky" => "Cumulative_Expansion_Fraction",
    "flow" => "Flow_MW", "LMP_from_bus" => "LMP_From_Bus_USD_per_MWh", "LMP_to_bus" => "LMP_To_Bus_USD_per_MWh",
    "congestion_flag" => "Congested_Flag", "wheeling_cost" => "Wheeling_Cost", "line_length" => "Line_Length", "length" => "Length",
    "f_exp" => "Expansion_Flow_MW", "flow_total" => "Total_Flow_MW", "kvl_residual" => "KVL_Residual_MW",
    # representative days
    "RepDay_id" => "Rep_Day", "Daygroup_ID" => "Day_Group_ID", "Numdays_Daygroup" => "Days_in_Group", "Reference_Day" => "Day_of_Year",
    # demand response / large flexible load
    "UNITGROUP" => "Unit_Group", "UNIT_CATEGORY" => "Unit_Category", "UNIT_REPORT_LABEL_1" => "Unit_Report_Label_1",
    "UNIT_REPORT_LABEL_2" => "Unit_Report_Label_2", "CAP" => "Unit_Capacity_MW", "INTERCON_LIM" => "Interconnection_Limit_MW",
    "lfl_lt" => "LFL_Load_MW", "lfl_DR_lt" => "LFL_DR_MW",
    "lfl_g_G_LFL_lt" => "LFL_Gen_to_Load_MW", "lfl_g_G_Grid_lt" => "LFL_Gen_to_Grid_MW", "lfl_g_G_ES_lt" => "LFL_Gen_to_Storage_MW",
    "lfl_g_ES_LFL_lt" => "LFL_Storage_to_Load_MW", "lfl_g_ES_Grid_lt" => "LFL_Storage_to_Grid_MW",
    "lfl_chg_Grid_ES_lt" => "LFL_Grid_to_Storage_MW", "lfl_soc_lt" => "LFL_Storage_SOC_MWh",
    # RA per-scenario files and water management
    "status" => "Unit_Available_Flag", "risk_status" => "Outage_Event_Flag", "g_ref_idht" => "Reference_Generation_MW",
    "chg_ref_idhty" => "Reference_Charge_MW", "soc_ref_idhty" => "Reference_SOC_MWh", "ens" => "Unserved_Energy_MW",
    "demand" => "Load_MW", "total_ICAP_loss" => "Outage_Capacity_MW", "segment" => "Segment", "water_use_ijdhty" => "Water_Use",
    # run time
    "Expansion Run Time" => "Expansion_Run_Time_s", "Operation Run Time" => "Operation_Run_Time_s", "RA Run Time" => "RA_Run_Time_s",
)
for k in 1:5
    REPORT_LABEL_RENAMES["lfl_seg_lt_$k"] = "LFL_DR_Segment_$(k)_MW"
    REPORT_LABEL_RENAMES["lfl_ind_lt_$k"] = "LFL_DR_Segment_$(k)_Active"
end

# Renamed labels; a calendar "Year" column is inserted right after "Stage" (rows must add it: see stage_calendar_year).
function report_labels(labels; overrides = Dict{String,String}(), add_year::Bool = true)
    out = String[get(overrides, l, get(REPORT_LABEL_RENAMES, l, l)) for l in labels]
    if add_year
        i = findfirst(==("Stage"), out)
        i === nothing || insert!(out, i + 1, "Year")
    end
    return out
end

ra_report_names(names) = Symbol.(report_labels(String.(names); add_year = false))

stage_calendar_year(setting, stage) = setting["Planning Design"]["first_stage_year_value"] + (stage - 1) * setting["Planning Design"]["num_years_per_stage_value"]

# ---------------------------------------------------------------------------------------------
# Summary-table presentation (shared by the EXP and OP writers)
# ---------------------------------------------------------------------------------------------

# Values below this magnitude are solver noise (e.g. 1e-8 "new units", 6e-5 MWh unserved energy) and are written as 0.
const REPORT_NOISE_TOL = 1e-3
const _NOISE_EXEMPT_COLUMNS = ("PRM", "CAPCRED")

_clean_noise(v) = (v isa AbstractFloat && abs(v) < REPORT_NOISE_TOL) ? 0.0 : v

# Dollar columns use a tolerance of $1 (e.g. 5e-5 MWh unserved energy x $9000/MWh is $0.5 of noise).
_is_cost_label(l) = l isa AbstractString && (occursin(r"(Cost|cost|Penalty|_PV|_real|_Committed|^Gen_ITC|^Generation_PTC)", l))

function zero_report_noise!(M::AbstractMatrix; exempt = _NOISE_EXEMPT_COLUMNS)
    for c in 1:size(M, 2)
        M[1, c] in exempt && continue
        tol = _is_cost_label(M[1, c]) ? 1.0 : REPORT_NOISE_TOL
        for r in 2:size(M, 1)
            v = M[r, c]
            v isa AbstractFloat && abs(v) < tol && (M[r, c] = 0.0)
        end
    end
    return M
end

function _rename_labels!(M::AbstractMatrix, renames)
    for c in 1:size(M, 2)
        haskey(renames, M[1, c]) && (M[1, c] = renames[M[1, c]])
    end
    return M
end

function _insert_column(M::AbstractMatrix, pos::Int, label, values)
    out = Array{Any}(undef, size(M, 1), size(M, 2) + 1)
    out[:, 1:pos-1] = M[:, 1:pos-1]
    out[1, pos] = label
    out[2:end, pos] = values
    out[:, pos+1:end] = M[:, pos:end]
    return out
end

# Costs whose payment stream spans several years: the by-stage table used to charge a stage with the whole
# stream of the units it builds. The by-stage columns now hold what is paid within the stage (they sum to the
# by-year rows); the whole-stream value stays available as <name>_Committed.
const _COMMITTED_COST_COLUMNS = ("Gen_Investment_Cost", "Gen_ITC", "Trans_Investment_Cost", "Generation_PTC")

# Returns (by_stage, by_year) tables ready to write. Inputs carry a header row.
function finalize_system_summaries(system_list::AbstractMatrix, extensive_list::AbstractMatrix, dollar_year)
    S = Array{Any}(copy(system_list))
    E = Array{Any}(copy(extensive_list))
    col(M, name) = findfirst(==(name), M[1, :])

    stage_col_E = col(E, "Stage")
    stage_col_S = col(S, "Stage")
    E_stage = E[2:end, stage_col_E]
    S_stage = S[2:end, stage_col_S]
    paid_in_stage(name_E) = [sum(Float64(E[1 + i, col(E, name_E)]) for i in eachindex(E_stage) if E_stage[i] == st) for st in S_stage]

    # committed (whole-stream) columns are appended at the end; the original columns become paid-in-stage
    committed_labels = String[]
    committed_values = Vector{Any}[]
    for name in _COMMITTED_COST_COLUMNS
        c = col(S, name)
        (c === nothing || col(E, name * "_CF") === nothing) && continue
        push!(committed_labels, name * "_Committed")
        push!(committed_values, Any[S[1 + i, c] for i in eachindex(S_stage)])
        S[2:end, c] = paid_in_stage(name * "_CF")
    end
    ct = col(S, "total_system_cost")
    if ct !== nothing && col(E, "total_system_cost_CF") !== nothing
        S[2:end, ct] = paid_in_stage("total_system_cost_CF")
    end
    for (lab, vals) in zip(committed_labels, committed_values)
        S = hcat(S, Any[lab; vals])
    end

    # by-year: no per-year copies of the multi-year objective or the constant year length
    keep = [c for c in 1:size(E, 2) if !(E[1, c] in ("ObjectiveValue", "Number of Years"))]
    E = E[:, keep]

    # names
    _rename_labels!(S, Dict("Scenario" => "Case_ID", "Start Year" => "Start_Year", "Number of Years" => "Years_in_Stage",
                            "Total_system_cost_OP" => "Operating_Cost"))
    E_ren = Dict{Any,Any}("Scenario" => "Case_ID", "total_system_cost_OP_CF" => "Operating_Cost_PV")
    for c in 1:size(E, 2)
        l = E[1, c]
        (l isa AbstractString) && endswith(l, "_CF") && !haskey(E_ren, l) && (E_ren[l] = l[1:end-3] * "_PV")
    end
    _rename_labels!(E, E_ren)

    # discount year
    S = _insert_column(S, col(S, "Years_in_Stage") + 1, "Discount_Year", fill(dollar_year, size(S, 1) - 1))
    E = _insert_column(E, col(E, "Year") + 1, "Discount_Year", fill(dollar_year, size(E, 1) - 1))

    zero_report_noise!(S)
    zero_report_noise!(E)
    return S, E
end

# Technology summary: same renames and noise cleaning; units that exist only as solver noise are blanked out.
function finalize_tech_summary(tech_list::AbstractMatrix)
    T = Array{Any}(copy(tech_list))
    col(name) = findfirst(==(name), T[1, :])
    cap_cols = [c for c in (col("ICAP"), col("ICap_New"), col("ICap_Ret")) if c !== nothing]
    first_metric = col("TotalUnits")
    if first_metric !== nothing && !isempty(cap_cols)
        keep_cols = Set(c for c in (col("CAPCRED"),) if c !== nothing)
        for r in 2:size(T, 1)
            all(c -> !(T[r, c] isa Number) || abs(T[r, c]) < REPORT_NOISE_TOL, cap_cols) || continue
            for c in first_metric:size(T, 2)
                c in keep_cols && continue
                T[r, c] isa Number && (T[r, c] = 0.0)
            end
        end
    end
    _rename_labels!(T, Dict("Scenario" => "Case_ID", "Start Year" => "Start_Year", "Number of Years" => "Years_in_Stage",
                            "PLANT_NAME" => "Plant_Name", "UnitGroup" => "Unit_Group",
                            "Gen_Investment_Cost" => "Gen_Investment_Cost_Committed", "Gen_ITC" => "Gen_ITC_Committed"))
    zero_report_noise!(T)
    return T
end

# Emit a constant-dollar (undiscounted) sibling of an annual summary. The by-year *_PV columns are per-year
# present values (discounted to the dollar year); reversing the single discount factor restores the real cost
# incurred each year in constant dollar-year dollars. Non-_PV columns (physical quantities) are copied verbatim.
# Adds System_Cost_per_MWh_real = total system cost / Demand_MWh for each year.
function write_constant_dollar_annual_summary(annual_file_name::AbstractString, annual_output_list, discount_rate, current_year)
    labels = annual_output_list[1, :]
    year_col = findfirst(==("Year"), labels)
    cost_cols = findall(l -> (l isa AbstractString) && endswith(l, "_PV"), labels)
    (year_col === nothing || isempty(cost_cols)) && return

    out = copy(annual_output_list)
    for r in 2:size(out, 1)
        yr = tryparse(Float64, string(out[r, year_col]))
        yr === nothing && continue
        undiscount = (1 + discount_rate)^(yr - current_year)
        for c in cost_cols
            v = out[r, c]
            v isa Number && (out[r, c] = v * undiscount)
        end
    end

    tc = findfirst(==("total_system_cost_PV"), out[1, :])
    dm = findfirst(==("Dispatched_Load_MWh"), out[1, :])
    if tc !== nothing && dm !== nothing
        out = _insert_column(out, size(out, 2) + 1, "System_Cost_per_MWh_real",
                             [(out[r, dm] isa Number && out[r, dm] > 0) ? out[r, tc] / out[r, dm] : 0.0 for r in 2:size(out, 1)])
    end

    # in this file the costs are undiscounted: relabel and drop the discount-year column
    for c in 1:size(out, 2)
        (out[1, c] isa AbstractString) && endswith(out[1, c], "_PV") && (out[1, c] = out[1, c][1:end-3] * "_real")
    end
    dy = findfirst(==("Discount_Year"), out[1, :])
    dy !== nothing && (out = out[:, [c for c in 1:size(out, 2) if c != dy]])

    constant_file_name = replace(annual_file_name, ".csv" => string("_real_", Int(round(current_year)), "usd.csv"))
    writedlm(constant_file_name, out, ",")
end


function report_result_tech_summary_expansion_GTEP(result::Dict{String, Any}, reference, setting, generate_file::Bool=true)

    pu_power_base = setting["Simulation Setting"]["per_unit_base_value"]
    pu_econ_base = setting["Simulation Setting"]["per_unit_econ_base_value"] / pu_power_base    # per unit economic base

    output_path = parameter(reference, 0, :output_path)
    tech_file_name = joinpath(output_path, string(parameter(reference, 0, :case_name), "__tech_summary_by_stage_EXP.csv"))
    extensive_tech_file_name = joinpath(output_path, string(parameter(reference, 0, :case_name), "__tech_summary_by_year_EXP.csv"))
            
    # Write Label 
    tech_label_list = [
        "Scenario",
        "Stage",
        "Start Year",
        "Number of Years",
        "PLANT_NAME",
        "Bus_ID",
        "Bus_Name",
        "Parent_Bus_Name",
        "Region_Name",
        "Tech_ID",
        "UnitGroup",
        "Unit_Category",
        "Unit_Report_Label_1",
        "Unit_Report_Label_2",
        "Fuel",
        "TotalUnits",
        "NewUnits",
        "RetUnits",
        "ICAP",
        "UCAP",
        "ICap_New",
        "ICap_Ret",
        "UCap_New",
        "UCap_Ret",
        "Storage_MWh",
        "Storage_Hr",
        "Generation",
        "Curtail",
        "Storage_Charge_MWh",
        "Reserve_RegUp",
        "Reserve_RegDn",
        "Reserve_Spin",
        "Reserve_NSpin",
        "Reserve_FlexUp",
        "Reserve_FlexDn",
        "Generation_Cost",
        "Charge_Cost",
        "Regulation_Cost",
        "Spin_Cost",
        "Nspin_Cost",
        "Flex_Cost",
        "UnitRevenue_E",
        "UnitRevenue_AS",
        "UnitRevenue_CRED",
        "UnitRevenue",
        "UnitProfit",
        "FuelConsumption",
        "FuelCost",
        "FOM",
        "CAPCRED",
        "Reference_Annual_Gen_Investment_Cost",
        "Gen_Investment_Cost",
        "Gen_ITC",
    ]

    annual_tech_label_list = [
        "Scenario",
        "Stage",
        "Start Year",
        "Number of Years",
        "Bus_ID",
        "Bus_Name",
        "Parent_Bus_Name",
        "Region_Name",
        "Tech_ID",
        "UnitGroup",
        "Unit_Category",
        "Unit_Report_Label_1",
        "Unit_Report_Label_2",
        "Fuel",
        "TotalUnits",
        "NewUnits",
        "RetUnits",
        "ICAP",
        "UCAP",
        "ICap_New",
        "ICap_Ret",
        "UCap_New",
        "UCap_Ret",
        "Storage_MWh",
        "Storage_Hr",
        "Generation",
        "Curtail",
        "Storage_Charge_MWh",
        "Reserve_RegUp",
        "Reserve_RegDn",
        "Reserve_Spin",
        "Reserve_NSpin",
        "Reserve_FlexUp",
        "Reserve_FlexDn",
        "Generation_Cost",
        "Charge_Cost",
        "Regulation_Cost",
        "Spin_Cost",
        "Nspin_Cost",
        "Flex_Cost",
        "UnitRevenue_E",
        "UnitRevenue_AS",
        "UnitRevenue_CRED",
        "UnitRevenue",
        "UnitProfit",
        "FuelConsumption",
        "FuelCost",
        "FOM",
        "CAPCRED",
    ]

    ids_y = []
    ids_i = []
    ids_n = []
    ids_z = []
    ids_p = []
    ids_k = []
    ids_d = []
    annual_gen_info = []
    try
        ids_y = [(y) for (y) in sort!(parse.(Int, (collect(get_index(reference, :planning_stages, 0)))))]
        ids_i = [(k) for (k) in sort!(parse.(Int, (collect(get_index(reference, :gen_index, 0)))))]
        ids_n = [(k) for (k) in sort!(parse.(Int, (collect(get_index(reference, :bus, 0)))))]
        ids_z = [(i) for (i) in sort!(parse.(Int, (collect(keys(reference[0][:zone]["reserve"])))))]
        ids_p = [(i) for (i) in sort!(parse.(Int, (collect(keys(reference[0][:zone]["policy"])))))]
        ids_d = [(k) for (k) in sort!(parse.(Int, (collect(get_index(reference, :repdays, 0)))))]
        annual_gen_info = get_annual_gen_info(reference)
        ids_k = [(k) for (k) in sort!(parse.(Int, (collect(get_index(reference, :branch, 0))))) if parameter(reference, 0, :branch, "model_flag", k) == true]
    catch
        ids_y = [(y) for (y) in sort!(collect(get_index(reference, :planning_stages, 0)))]
        ids_i = [(k) for (k) in sort!(collect(get_index(reference, :gen_index, 0)))]
        ids_n = [(k) for (k) in sort!(collect(get_index(reference, :bus, 0)))]
        ids_z = [(i) for (i) in sort!(collect(keys(reference[0][:zone]["reserve"])))]
        ids_p = [(i) for (i) in sort!(collect(keys(reference[0][:zone]["policy"])))]
        ids_d = [(k) for (k) in sort!(collect(get_index(reference, :repdays, 0)))]
        annual_gen_info = get_annual_gen_info(reference)
        ids_k = [(k) for (k) in sort!(collect(get_index(reference, :branch, 0))) if parameter(reference, 0, :branch, "model_flag", k) == true]
    end

    # total number of techs
    num_tech = 0
    num_bus_tech = 0
    num_tech_per_bus = [0]
    for bus_idx in ids_n
        try
            num_tech += length(keys(reference[string(bus_idx)]["gen_bus"]))
            num_bus_tech += length(keys(reference[string(bus_idx)]["gen_bus"]))
        catch
            num_tech += length(keys(reference[bus_idx][:gen_bus]))
            num_bus_tech += length(keys(reference[bus_idx][:gen_bus]))
        end
        push!(num_tech_per_bus, num_bus_tech)
    end

    tech_output_list = Array{Any}(undef, length(ids_y)*num_tech+1, length(tech_label_list))
    tech_output_list[1,:] = tech_label_list

    annual_tech_output_list = Array{Any}(undef, length(ids_y)*num_tech+1, length(annual_tech_label_list))
    annual_tech_output_list[1,:] = annual_tech_label_list

    # parameters
    Scenario = parameter(reference, 0, :case_name)
    current_year = setting["Planning Design"]["dollar_year_value"] # Discount factor relative to the current dollor year
    discount_rate = setting["Planning Design"]["discount_rate_value"]

    # Annualized investment cost for transmission lines
    for year_id in eachindex(ids_y)
        
        year_of_year_id = parameter(reference, 0, :planning_stages, "year", year_id)
        stage_length = parameter(reference, 0, :planning_stages, "stage_length", year_id) 
        
        num_stage = setting["Planning Design"]["num_stages_value"]
        remaining_years = (num_stage + 1 - year_id) * stage_length

        # Plant level data collection                
        @sync for bus_idx in ids_n

            @spawn begin

                Bus_Name = parameter(reference, 0, :bus, "bus_i", bus_idx)
                Parent_Bus_Name = reference[0][:bus][bus_idx]["region_config"]["parent_bus_name"]
                Region_Name = reference[0][:bus][bus_idx]["region_config"]["region_name"]


                ids_tech = []
                try
                    ids_tech = [(k) for (k) in sort!(parse.(Int, (collect( keys(reference[string(bus_idx)]["gen_bus"])))))]
                catch
                    ids_tech = [(k) for (k) in sort!(collect(keys(reference[bus_idx][:gen_bus])))]
                end

                ids_tech_of_this_bus = [i for i in 1:length(ids_tech)]

                for tech_id in ids_tech_of_this_bus

                    tech_idx = ids_tech[tech_id]
                    
                    gen_idx = parameter(reference, bus_idx, :gen_bus, tech_idx, "gen_idx")
                    unitgroup = parameter(reference, bus_idx, :gen_bus, tech_idx, "UNITGROUP")
                    Unit_Category = parameter(reference, bus_idx, :gen_bus, tech_idx, "UNIT_CATEGORY")
                    Unit_Report_Label_1 = parameter(reference, bus_idx, :gen_bus, tech_idx, "UNIT_REPORT_LABEL_1")
                    Unit_Report_Label_2 = parameter(reference, bus_idx, :gen_bus, tech_idx, "UNIT_REPORT_LABEL_2")

                    # plant name: based on NEW_Gen_UID or PLANT_NAME
                    plant_name = parameter(reference, bus_idx, :gen_bus, tech_idx, "NEW_Gen_UID")
                    if haskey(reference[bus_idx][:gen_bus][tech_idx], "PLANT_NAME")
                        plant_name = parameter(reference, bus_idx, :gen_bus, tech_idx, "PLANT_NAME")
                    end
                    
                    Fuel = parameter(reference, bus_idx, :gen_bus, tech_idx, "FUEL")
                    Emission_rate = parameter(reference, bus_idx, :gen_bus, tech_idx, "Emission_CO2") / pu_power_base
                    capacity = parameter(reference, bus_idx, :gen_bus, tech_idx, "CAP") * pu_power_base
                    FOM = parameter(reference, bus_idx, :gen_bus, tech_idx, "Annual_FOM")[string(year_id)] * pu_econ_base
                    CRP_i = parameter(reference, bus_idx, :gen_bus, tech_idx, "crpyears")

                    CAPCRED = parameter(reference, bus_idx, :gen_bus, tech_idx, "CAPCRED")
                    if CAPCRED isa String
                        data_identifier = reference[0][:zone]["capacity_credit"]["data_identifier"]
                        data_region = reference[0][:bus][bus_idx]["region_mapping_info"][data_identifier]
                        data_region = data_region isa AbstractString ? data_region : string(data_region)
                        CAPCRED = reference[0][:zone]["capacity_credit"][data_region][CAPCRED]
                    end

                    # get updated CAPCRED if available
                    if year_id > 1 # only check if this dict is not empty (i.e., this dict will be empty in year 1)
                        if (setting["Planning Design"]["multi_round_solution_process_flag"] == true) && (setting["Simulation Configuration"]["update_CAPCRED_in_each_round_of_Expansion_Flag"] == true)
                            if parameter(reference, bus_idx, :gen_bus, tech_idx, "ELCC_Flag") == true
                                try
                                    CAPCRED = parameter(reference, 0, :multi_round_info, "RA_Info", year_id)["ELCC_result_applied"][unitgroup]
                                catch
                                    @aleaf_info "[ALEAF LC_GTEP]: Failed to find a new capacity credit of $unitgroup in the report_result_tech_summary_expansion_GTEP function"
                                end
                            end
                        end
                    end
                        
                    # Collect data
                    idx_string = string("(", gen_idx, ", ", year_id, ")")

                    tag = string("Ret_", year_of_year_id)
                    planned_retirement = parameter(reference, bus_idx, :gen_bus, tech_idx, "Planned_Retirement")[tag]

                    TotalUnits = result["1"]["solution"]["expansion"][idx_string]["u_G_iy"]
                    NewUnits = result["1"]["solution"]["expansion"][idx_string]["u_new_G_iy"]  
                    RetUnits = result["1"]["solution"]["expansion"][idx_string]["u_ret_G_iy"]
                    
                    ICAP = capacity * TotalUnits 
                    UCAP = CAPCRED * ICAP

                    ICap_New = capacity * NewUnits
                    ICap_Ret = capacity * RetUnits
                    UCap_New = CAPCRED * ICap_New
                    UCap_Ret = CAPCRED * ICap_Ret

                    Storage_Hr = 0.0
                    u_ESE = 0.0
                    u_new_ESH = 0.0
                    if Unit_Category == "STORAGE"
                        if TotalUnits > 0
                            u_ESE = result["1"]["solution"]["expansion"][idx_string]["u_ESE_iy"]
                            u_new_ESH = result["1"]["solution"]["expansion"][idx_string]["u_new_ESH_iy"]
                            # MW floor avoids divide-by-near-zero on round-off-tiny storage units
                            if ICAP > 1e-3
                                Storage_Hr = u_ESE / ICAP
                            end
                        end
                    end
                    
                    ######### Investment and Retirement costs 
                    
                    # Determine the cost duration (investment payments)
                    total_gen_invc = 0.0
                    total_ITC = 0.0
                    total_retirement_cost = 0.0

                    payment_duration = min(CRP_i, remaining_years)

                    # Get the investment costs
                    investment_cost = parameter(reference, bus_idx, :gen_bus, tech_idx, "Annual_INVC")[string(year_id)] * pu_econ_base
                    storage_investment_cost = parameter(reference, bus_idx, :gen_bus, tech_idx, "Annual_STO_INV")[string(year_id)] * pu_econ_base
                    DECC = parameter(reference, bus_idx, :gen_bus, tech_idx, "DECC") * pu_econ_base

                    # Get ITC 
                    itc = 0.0
                    if (setting["Simulation Configuration"]["ITC_Flag"] == true) && (parameter(reference, bus_idx, :gen_bus, tech_idx, "ITC Flag") == true)
                        itc_year = min(2050, year_of_year_id)
                        itc = parameter(reference, bus_idx, :gen_bus, tech_idx, "ITC")[string(itc_year)] # %
                    end

                    # Collect data
                    u_new_G_iy = NewUnits
                    u_new_ESH_iy = u_new_ESH

                    # base investment cost
                    base_INVC = 1000 * investment_cost * capacity * u_new_G_iy / (1.0 - itc)
                    if Unit_Category == "STORAGE"
                        base_INVC += 1000 * storage_investment_cost * capacity * u_new_ESH_iy / (1.0 - itc)
                    end

                    for future_y in 1:payment_duration

                        future_year = year_of_year_id + future_y - 1    # Actual future year 
                        discount_factor = (1 + discount_rate)^-(future_year - current_year) # Discount factor relative to the current dollor year
        
                        gen_invc = 0.0
                        
                        # Expansion Cost (new units only) 
                        gen_invc += 1000 * investment_cost * capacity * u_new_G_iy / (1.0 - itc)
        
                        # Storage Duration Cost ($/kWh)
                        if Unit_Category == "STORAGE"
                            gen_invc += 1000 * storage_investment_cost * capacity * u_new_ESH_iy / (1.0 - itc)
                        end
        
                        # update the total investment cost
                        total_gen_invc += discount_factor * gen_invc
        
                        # Investment Tax Credits
                        total_ITC += itc * discount_factor * gen_invc

                    end
                    
                    # FOM
                    total_FOM = 0.0
                    annual_FOM = 0.0
                    if (FOM > 0) && (capacity > 0)
        
                        for future_y in 1:stage_length
                        
                            future_year = year_of_year_id + future_y - 1 
                            discount_factor_year_id = (1 + discount_rate)^-(future_year - current_year)
        
                            # Fixed OM Cost ($/kW-Year)
                            total_FOM += discount_factor_year_id * 1000 * FOM * ICAP

                            if future_year == 1
                                annual_FOM = discount_factor_year_id * 1000 * FOM * ICAP
                            end
                        end
                    end
                                        
                    # Retirement Cost
                    if DECC > 0
                        discount_factor_year_id = (1 + discount_rate)^-(year_of_year_id - current_year)
                        retirement_cost = 1000 * DECC * capacity * RetUnits
                        total_retirement_cost = discount_factor_year_id * retirement_cost
                    end
                    
                    # Get pre-calculated annual values from the dispatch output (already multiplied by stage_length)
                    Generation = annual_gen_info[gen_idx][year_id]["Generation"]
                    Curtail = annual_gen_info[gen_idx][year_id]["Curtail"]

                    Storage_Charge = annual_gen_info[gen_idx][year_id]["Storage_Charge"]
                    Reserve_RU = annual_gen_info[gen_idx][year_id]["Reserve_RU"]
                    Reserve_RD = annual_gen_info[gen_idx][year_id]["Reserve_RD"]
                    Reserve_Spin = annual_gen_info[gen_idx][year_id]["Reserve_Spin"]
                    Reserve_NSpin = annual_gen_info[gen_idx][year_id]["Reserve_NSpin"]
                    Reserve_FU = annual_gen_info[gen_idx][year_id]["Reserve_FU"]
                    Reserve_FD = annual_gen_info[gen_idx][year_id]["Reserve_FD"]
                    
                    UnitRevenue_E = annual_gen_info[gen_idx][year_id]["UnitRevenue_E"]
                    UnitRevenue_AS = annual_gen_info[gen_idx][year_id]["UnitRevenue_AS"]
                    UnitRevenue_CRED = annual_gen_info[gen_idx][year_id]["UnitRevenue_CRED"]
                    UnitRevenue = annual_gen_info[gen_idx][year_id]["UnitRevenue"]

                    Noload_Cost = annual_gen_info[gen_idx][year_id]["Noload_Cost"]
                    StartUp_Cost = annual_gen_info[gen_idx][year_id]["StartUp_Cost"]

                    Generation_Cost = annual_gen_info[gen_idx][year_id]["Generation_Cost"]
                    Charge_Cost = annual_gen_info[gen_idx][year_id]["Charge_Cost"]
                    
                    Reserve_Reg_Cost = annual_gen_info[gen_idx][year_id]["Reserve_Reg_Cost"]
                    Reserve_Spin_Cost = annual_gen_info[gen_idx][year_id]["Reserve_Spin_Cost"]
                    Reserve_NSpin_Cost = annual_gen_info[gen_idx][year_id]["Reserve_NSpin_Cost"]
                    Reserve_Flex_Cost = annual_gen_info[gen_idx][year_id]["Reserve_Flex_Cost"]

                    Fuel_Cost = annual_gen_info[gen_idx][year_id]["Fuel_Cost"]
                    Fuel_Consumption = annual_gen_info[gen_idx][year_id]["Fuel_Consumption"]

                    # Calculate unit revenue
                    total_revenue_E = 0.0
                    annual_revenue_E = 0.0

                    total_revenue_AS = 0.0
                    annual_revenue_AS = 0.0

                    total_revenue_CRED = 0.0
                    annual_revenue_CRED = 0.0

                    total_revenue = 0.0
                    annual_revenue = 0.0


                    for future_y in 1:stage_length
                        future_year = year_of_year_id + future_y - 1 
                        discount_factor = (1 + discount_rate)^-(future_year - current_year)

                        total_revenue_E += discount_factor * (UnitRevenue_E / stage_length)
                        total_revenue_AS += discount_factor * (UnitRevenue_AS / stage_length)
                        total_revenue_CRED += discount_factor * (UnitRevenue_CRED / stage_length)
                        total_revenue += discount_factor * (UnitRevenue / stage_length)

                        if future_year == 1
                            annual_revenue_E = discount_factor * (UnitRevenue_E / stage_length)
                            annual_revenue_AS = discount_factor * (UnitRevenue_AS / stage_length)
                            annual_revenue_CRED = discount_factor * (UnitRevenue_CRED / stage_length)
                            annual_revenue = discount_factor * (UnitRevenue / stage_length)
                        end
                    end
                    
                    # Calculate Unit Profit
                    total_profit = 0.0
                    annual_profit = 0.0
                    UnitProfit = annual_gen_info[gen_idx][year_id]["UnitProfit"] - Noload_Cost - StartUp_Cost - FOM # (already multiplied by stage_length)
                    for future_y in 1:stage_length
                        future_year = year_of_year_id + future_y - 1 
                        discount_factor = (1 + discount_rate)^-(future_year - current_year)

                        total_profit += discount_factor * (UnitProfit / stage_length)

                        if future_y == 1
                            annual_profit = discount_factor * (UnitProfit / stage_length)
                        end
                    end
                    

                    # Gen and Charge Cost 
                    total_gen_cost = 0.0
                    annual_gen_cost = 0.0

                    total_charge_cost = 0.0
                    annual_charge_cost = 0.0
                    for future_y in 1:stage_length
                        future_year = year_of_year_id + future_y - 1 
                        discount_factor = (1 + discount_rate)^-(future_year - current_year)

                        total_gen_cost += discount_factor * (Generation_Cost / stage_length)
                        total_charge_cost += discount_factor * (Charge_Cost / stage_length)

                        if future_y == 1
                            annual_gen_cost = discount_factor * (Generation_Cost / stage_length)
                            annual_charge_cost = discount_factor * (Charge_Cost / stage_length)
                        end
                    end

                    # AS Cost 
                    total_reg_cost = 0.0
                    annual_reg_cost = 0.0

                    total_spin_cost = 0.0
                    annual_spin_cost = 0.0

                    total_nspin_cost = 0.0                    
                    annual_nspin_cost = 0.0

                    total_flex_cost = 0.0
                    annual_flex_cost = 0.0

                    for future_y in 1:stage_length
                        future_year = year_of_year_id + future_y - 1 
                        discount_factor = (1 + discount_rate)^-(future_year - current_year)

                        total_reg_cost += discount_factor * (Reserve_Reg_Cost / stage_length)
                        total_spin_cost += discount_factor * (Reserve_Spin_Cost / stage_length)
                        total_nspin_cost += discount_factor * (Reserve_NSpin_Cost / stage_length)
                        total_flex_cost += discount_factor * (Reserve_Flex_Cost / stage_length)

                        if future_y == 1
                            annual_reg_cost = discount_factor * (Reserve_Reg_Cost / stage_length)
                            annual_spin_cost = discount_factor * (Reserve_Spin_Cost / stage_length)
                            annual_nspin_cost = discount_factor * (Reserve_NSpin_Cost / stage_length)
                            annual_flex_cost = discount_factor * (Reserve_Flex_Cost / stage_length)
                        end
                    end

                    # Fuel Cost
                    total_fuel_cost = 0.0
                    annual_fuel_cost = 0.0
                    for future_y in 1:stage_length
                        future_year = year_of_year_id + future_y - 1 
                        discount_factor = (1 + discount_rate)^-(future_year - current_year)

                        total_fuel_cost += discount_factor * (Fuel_Cost / stage_length)

                        if future_year == 1
                            annual_fuel_cost = discount_factor * (Fuel_Cost / stage_length)
                        end
                    end

                    if typeof(ids_tech[1]) == String
                        tech_idx = string(tech_idx)
                    end

                    row_id = (year_id-1)*num_tech + num_tech_per_bus[bus_idx] + tech_id + 1

                    # Update the list
                    tech_output_list[row_id, :] = [
                            Scenario, 
                            year_id, 
                            year_of_year_id, 
                            stage_length, 
                            plant_name,
                            bus_idx, 
                            Bus_Name, 
                            Parent_Bus_Name,
                            Region_Name, 
                            tech_idx, 
                            unitgroup,
                            Unit_Category,
                            Unit_Report_Label_1,
                            Unit_Report_Label_2,
                            Fuel,
                            TotalUnits, 
                            NewUnits, 
                            RetUnits, 
                            ICAP, 
                            UCAP, 
                            ICap_New, 
                            ICap_Ret, 
                            UCap_New, 
                            UCap_Ret, 
                            u_ESE, 
                            Storage_Hr, 
                            Generation, 
                            Curtail,
                            Storage_Charge, 
                            Reserve_RU, 
                            Reserve_RD, 
                            Reserve_Spin, 
                            Reserve_NSpin, 
                            Reserve_FU, 
                            Reserve_FD, 
                            total_gen_cost,
                            total_charge_cost,
                            total_reg_cost,
                            total_spin_cost,                            
                            total_nspin_cost,
                            total_flex_cost,
                            total_revenue_E, 
                            total_revenue_AS, 
                            total_revenue_CRED, 
                            total_revenue, 
                            total_profit, 
                            Fuel_Consumption, 
                            total_fuel_cost, 
                            total_FOM, 
                            CAPCRED, 
                            base_INVC, 
                            total_gen_invc, 
                            total_ITC
                        ]

                    annual_tech_output_list[row_id, :] = [
                            Scenario, 
                            year_id, 
                            year_of_year_id, 
                            stage_length, 
                            bus_idx, 
                            Bus_Name, 
                            Parent_Bus_Name,
                            Region_Name, 
                            tech_idx, 
                            unitgroup,
                            Unit_Category,
                            Unit_Report_Label_1,
                            Unit_Report_Label_2,
                            Fuel,
                            TotalUnits, 
                            NewUnits, 
                            RetUnits, 
                            ICAP, 
                            UCAP, 
                            ICap_New, 
                            ICap_Ret, 
                            UCap_New, 
                            UCap_Ret, 
                            u_ESE, 
                            Storage_Hr, 
                            Generation / stage_length, 
                            Curtail / stage_length,
                            Storage_Charge / stage_length, 
                            Reserve_RU / stage_length, 
                            Reserve_RD / stage_length, 
                            Reserve_Spin / stage_length, 
                            Reserve_NSpin / stage_length, 
                            Reserve_FU / stage_length, 
                            Reserve_FD / stage_length, 
                            annual_gen_cost,
                            annual_charge_cost,
                            annual_reg_cost,
                            annual_spin_cost,
                            annual_nspin_cost,
                            annual_flex_cost,
                            annual_revenue_E, 
                            annual_revenue_AS,
                            annual_revenue_CRED,
                            annual_revenue,
                            annual_profit, 
                            Fuel_Consumption / stage_length, 
                            annual_fuel_cost, 
                            annual_FOM, 
                            CAPCRED
                        ]
                    


                end
            end
        end
        
    end
    writedlm(tech_file_name, finalize_tech_summary(tech_output_list), ",")

    @aleaf_info "[ALEAF LC_GTEP]: - Tech summary reporting,"
end


# System-wide annual nominal demand (MWh, per-unit), mirroring
# get_annual_coincident_peak_demand_with_growth but summing hourly load instead of taking the max.
function get_annual_total_demand_with_growth(network_data::Dict{Symbol,<:Any}, year)
    demand_df = network_data[:time_series_data]["load"]
    region_map = get(network_data, :profile_data_region_map, nothing)
    hourly = nothing
    for bus_data in values(network_data[:bus])
        original_load = bus_data["aggregation_info"]["original_load_(bus_i, MW)"]
        for region_id in bus_data["aggregation_info"]["aggregated_regions_bus_i"]
            w = original_load[region_id] * get_annual_load_growth_factor(network_data, year, region_id)
            col = demand_df[!, Symbol(profile_data_region(region_map, "load", region_id))]
            hourly === nothing && (hourly = zeros(Float64, length(col)))
            hourly .+= w .* col
        end
    end
    return hourly === nothing ? 0.0 : sum(hourly)
end


# Policy-slack value (clean-energy target / RPS) of one index from an expansion solution. If the solution carries
# no "slack" group at all (for example a resume from a checkpoint written before slack was stored), the penalty
# cost is reported as 0 with a one-time warning instead of aborting the run at the end of the reporting stage.
const _WARNED_MISSING_POLICY_SLACK = Ref(false)
function policy_slack_value(solution::Dict, idx::AbstractString, name::AbstractString)
    slack = get(solution, "slack", nothing)
    if slack === nothing
        if !_WARNED_MISSING_POLICY_SLACK[]
            _WARNED_MISSING_POLICY_SLACK[] = true
            @aleaf_warn "Expansion solution has no policy slack results (older checkpoint?); clean-energy target / RPS slack costs are reported as 0."
        end
        return 0.0
    end
    return slack[idx][name]
end

function report_result_system_summary_expansion_GTEP(result::Dict{String, Any}, reference, setting, generate_file::Bool=true)

    pu_power_base = setting["Simulation Setting"]["per_unit_base_value"]
    pu_econ_ref_base = setting["Simulation Setting"]["per_unit_econ_base_value"]
    pu_econ_base = setting["Simulation Setting"]["per_unit_econ_base_value"] / pu_power_base    # per unit economic base

    output_path = parameter(reference, 0, :output_path)
    system_file_name = joinpath(output_path, string(parameter(reference, 0, :case_name), "__system_summary_by_stage_EXP.csv"))
    extensive_system_file_name = joinpath(output_path, string(parameter(reference, 0, :case_name), "__system_summary_by_year_EXP.csv"))
    
    system_label_list = [
        "Scenario",
        "Stage",
        "Start Year",
        "Number of Years",
        "Gen_Investment_Cost",
        "Gen_ITC",
        "Trans_Investment_Cost",
        "Trans_FOM_Cost",
        "Gen_Retirement_Cost",
        "FOM_Cost",
        "Generation_PTC",
        "Fuel_Cost",
        "VOM_Cost",
        "Commitment_Cost",
        "Regulation_Cost",
        "Spin_Cost",
        "Nspin_Cost",
        "Flex_Cost",
        "ENS_Cost",
        "RNS_Spin_Cost",
        "RNS_NSpin_Cost",
        "RNS_Flex_Cost",
        "CTAX_cost",
        "CEGT_Penalty",
        "RPS_Penalty",
        "total_system_cost",
        "ObjectiveValue",
        "Generation",
        "TnD_Loss",
        "Storage_Charge_MWh",
        "Reserve_RegUp",
        "Reserve_RegDn",
        "Reserve_Spin",
        "Reserve_NSpin",
        "Reserve_FlexUp",
        "Reserve_FlexDn",
        "ENS",
        "RNS_Spin",
        "RNS_NSpin",
        "RNS_Flex",
        "Emission",
        "PRM",
        "Peak_Demand_MW",
        "UCAP_MW",
        "Installed_Capacity_MW",
        "Annual_Input_Demand_MWh",
        "Dispatched_Load_MWh",
    ]

    extensive_system_label_list = [
        "Scenario",
        "Stage",
        "Year",
        "Number of Years",
        "Gen_Investment_Cost_CF",
        "Gen_ITC_CF",
        "Trans_Investment_Cost_CF",
        "Trans_FOM_Cost_CF",
        "Gen_Retirement_Cost_CF",
        "FOM_Cost_CF",
        "Generation_PTC_CF",
        "Fuel_Cost_CF",
        "VOM_Cost_CF",
        "Commitment_Cost_CF",
        "Regulation_Cost_CF",
        "Spin_Cost_CF",
        "Nspin_Cost_CF",
        "Flex_Cost_CF",
        "ENS_Cost_CF",
        "RNS_Spin_Cost_CF",
        "RNS_NSpin_Cost_CF",
        "RNS_Flex_Cost_CF",
        "CTAX_cost_CF",
        "CEGT_Penalty_CF",
        "RPS_Penalty_CF",
        "total_system_cost_CF",
        "ObjectiveValue",
        "Generation",
        "TnD_Loss",
        "Storage_Charge_MWh",
        "Reserve_RegUp",
        "Reserve_RegDn",
        "Reserve_Spin",
        "Reserve_NSpin",
        "Reserve_FlexUp",
        "Reserve_FlexDn",
        "ENS",
        "RNS_Spin",
        "RNS_NSpin",
        "RNS_Flex",
        "Emission",
        "PRM",
        "Peak_Demand_MW",
        "UCAP_MW",
        "Installed_Capacity_MW",
        "Annual_Input_Demand_MWh",
        "Dispatched_Load_MWh",
    ]

    ids_y = []
    ids_i = []
    ids_n = []
    ids_z = []
    ids_p = []
    ids_k = []
    ids_d = []
    annual_gen_info = []
    try
        ids_y = [(y) for (y) in sort!(parse.(Int, (collect(get_index(reference, :planning_stages, 0)))))]
        ids_i = [(k) for (k) in sort!(parse.(Int, (collect(get_index(reference, :gen_index, 0)))))]
        ids_n = [(k) for (k) in sort!(parse.(Int, (collect(get_index(reference, :bus, 0)))))]
        ids_z = [(i) for (i) in sort!(parse.(Int, (collect(keys(reference[0][:zone]["reserve"])))))]
        ids_p = [(i) for (i) in sort!(parse.(Int, (collect(keys(reference[0][:zone]["policy"])))))]
        ids_d = [(k) for (k) in sort!(parse.(Int, (collect(get_index(reference, :repdays, 0)))))]
        annual_gen_info = get_annual_gen_info(reference)
        ids_k = [(k) for (k) in sort!(parse.(Int, (collect(get_index(reference, :branch, 0))))) if parameter(reference, 0, :branch, "model_flag", k) == true]
    catch
        ids_y = [(y) for (y) in sort!(collect(get_index(reference, :planning_stages, 0)))]
        ids_i = [(k) for (k) in sort!(collect(get_index(reference, :gen_index, 0)))]
        ids_n = [(k) for (k) in sort!(collect(get_index(reference, :bus, 0)))]
        ids_z = [(i) for (i) in sort!(collect(keys(reference[0][:zone]["reserve"])))]
        ids_p = [(i) for (i) in sort!(collect(keys(reference[0][:zone]["policy"])))]
        ids_d = [(k) for (k) in sort!(collect(get_index(reference, :repdays, 0)))]
        annual_gen_info = get_annual_gen_info(reference)
        ids_k = [(k) for (k) in sort!(collect(get_index(reference, :branch, 0))) if parameter(reference, 0, :branch, "model_flag", k) == true]
    end
    annual_scarcity_info = get_annual_scarcity_info(reference)

    system_output_list = Array{Any}(undef, length(ids_y)+1, length(system_label_list))
    system_output_list[1,:] = system_label_list

    extensive_system_output_list = Array{Any}(undef, length(ids_y) * parameter(reference, 0, :planning_stages, "stage_length", 1)+1, length(extensive_system_label_list))
    extensive_system_output_list[1,:] = extensive_system_label_list

    # parameters
    Scenario = parameter(reference, 0, :case_name)
    current_year = setting["Planning Design"]["dollar_year_value"] # Discount factor relative to the current dollor year
    discount_rate = setting["Planning Design"]["discount_rate_value"]

    # total cost arrays
    annual_total_invc = zeros(Float64, length(ids_y))

    total_gen_invc_vector = zeros(Float64, length(ids_y) * parameter(reference, 0, :planning_stages, "stage_length", 1))
    total_gen_ITC_vector = zeros(Float64, length(ids_y) * parameter(reference, 0, :planning_stages, "stage_length", 1))
    total_trans_invc_vector = zeros(Float64, length(ids_y) * parameter(reference, 0, :planning_stages, "stage_length", 1))
    total_trans_fom_vector = zeros(Float64, length(ids_y) * parameter(reference, 0, :planning_stages, "stage_length", 1))
    total_gen_retc_vector = zeros(Float64, length(ids_y) * parameter(reference, 0, :planning_stages, "stage_length", 1))
    total_gen_FOM_vector = zeros(Float64, length(ids_y) * parameter(reference, 0, :planning_stages, "stage_length", 1))

    total_gen_cost_vector = zeros(Float64, length(ids_y) * parameter(reference, 0, :planning_stages, "stage_length", 1))
    total_fuel_cost_vector= zeros(Float64, length(ids_y) * parameter(reference, 0, :planning_stages, "stage_length", 1))
    total_VOM_cost_vector = zeros(Float64, length(ids_y) * parameter(reference, 0, :planning_stages, "stage_length", 1))
    total_commitment_cost_vector = zeros(Float64, length(ids_y) * parameter(reference, 0, :planning_stages, "stage_length", 1))

    total_regulation_cost_vector = zeros(Float64, length(ids_y) * parameter(reference, 0, :planning_stages, "stage_length", 1))
    total_spin_cost_vector = zeros(Float64, length(ids_y) * parameter(reference, 0, :planning_stages, "stage_length", 1))
    total_non_spin_cost_vector = zeros(Float64, length(ids_y) * parameter(reference, 0, :planning_stages, "stage_length", 1))
    total_flex_cost_vector = zeros(Float64, length(ids_y) * parameter(reference, 0, :planning_stages, "stage_length", 1))

    total_carbon_tax_vector = zeros(Float64, length(ids_y) * parameter(reference, 0, :planning_stages, "stage_length", 1))
    total_emissions_vector = zeros(Float64, length(ids_y) * parameter(reference, 0, :planning_stages, "stage_length", 1))
    total_ptc_vector = zeros(Float64, length(ids_y) * parameter(reference, 0, :planning_stages, "stage_length", 1))

    total_scarcity_cost_ENS_vector = zeros(Float64, length(ids_y) * parameter(reference, 0, :planning_stages, "stage_length", 1))
    total_scarcity_cost_SPIN_vector = zeros(Float64, length(ids_y) * parameter(reference, 0, :planning_stages, "stage_length", 1))
    total_scarcity_cost_NSPIN_vector = zeros(Float64, length(ids_y) * parameter(reference, 0, :planning_stages, "stage_length", 1))
    total_scarcity_cost_Flex_vector = zeros(Float64, length(ids_y) * parameter(reference, 0, :planning_stages, "stage_length", 1))

    total_scarcity_ENS_vector = zeros(Float64, length(ids_y) * parameter(reference, 0, :planning_stages, "stage_length", 1))
    total_scarcity_SPIN_vector = zeros(Float64, length(ids_y) * parameter(reference, 0, :planning_stages, "stage_length", 1))
    total_scarcity_NSPIN_vector = zeros(Float64, length(ids_y) * parameter(reference, 0, :planning_stages, "stage_length", 1))
    total_scarcity_FLEX_vector = zeros(Float64, length(ids_y) * parameter(reference, 0, :planning_stages, "stage_length", 1))

    total_CEGT_slack_vector = zeros(Float64, length(ids_y) * parameter(reference, 0, :planning_stages, "stage_length", 1))
    total_RPS_slack_vector = zeros(Float64, length(ids_y) * parameter(reference, 0, :planning_stages, "stage_length", 1))

    total_generation_vector = zeros(Float64, length(ids_y) * parameter(reference, 0, :planning_stages, "stage_length", 1))
    total_tnd_loss_vector = zeros(Float64, length(ids_y) * parameter(reference, 0, :planning_stages, "stage_length", 1))
    total_storage_charge_mwh_vector = zeros(Float64, length(ids_y) * parameter(reference, 0, :planning_stages, "stage_length", 1))
    total_reg_up_vector = zeros(Float64, length(ids_y) * parameter(reference, 0, :planning_stages, "stage_length", 1))
    total_reg_dn_vector = zeros(Float64, length(ids_y) * parameter(reference, 0, :planning_stages, "stage_length", 1))
    total_spin_vector = zeros(Float64, length(ids_y) * parameter(reference, 0, :planning_stages, "stage_length", 1))
    total_nspin_vector = zeros(Float64, length(ids_y) * parameter(reference, 0, :planning_stages, "stage_length", 1))
    total_flex_up_vector = zeros(Float64, length(ids_y) * parameter(reference, 0, :planning_stages, "stage_length", 1))
    total_flex_dn_vector = zeros(Float64, length(ids_y) * parameter(reference, 0, :planning_stages, "stage_length", 1))

    total_prm_vector = zeros(Float64, length(ids_y) * parameter(reference, 0, :planning_stages, "stage_length", 1))
    total_peak_demand_vector = zeros(Float64, length(ids_y) * parameter(reference, 0, :planning_stages, "stage_length", 1))
    total_ucap_vector = zeros(Float64, length(ids_y) * parameter(reference, 0, :planning_stages, "stage_length", 1))
    total_icap_vector = zeros(Float64, length(ids_y) * parameter(reference, 0, :planning_stages, "stage_length", 1))
    total_demand_vector = zeros(Float64, length(ids_y) * parameter(reference, 0, :planning_stages, "stage_length", 1))
    total_dispatched_load_vector = zeros(Float64, length(ids_y) * parameter(reference, 0, :planning_stages, "stage_length", 1))

    # start reporting for each year 
    for year_id in eachindex(ids_y)
        
        year_of_year_id = parameter(reference, 0, :planning_stages, "year", year_id)
        stage_length = parameter(reference, 0, :planning_stages, "stage_length", year_id) 
        
        num_stage = setting["Planning Design"]["num_stages_value"]
        remaining_years = (num_stage + 1 - year_id) * stage_length
        
        # Investment costs for generators + tax credits 
        total_gen_invc = 0.0    # 1) Upfront total investment costs 
        total_ITC = 0.0         # 2) Investment Tax Credits

        for i in ids_i

            bus_idx = parameter(reference, 0, :gen_index, "bus_idx", i)
            tech_idx = parameter(reference, 0, :gen_index, "genco_tech_id", i)
            capacity = parameter(reference, bus_idx, :gen_bus, tech_idx, "CAP") * pu_power_base
            CRP_i = parameter(reference, bus_idx, :gen_bus, tech_idx, "crpyears")

            # Determine the cost duration (investment payments)
            payment_duration = min(CRP_i, remaining_years)

            # Get the investment costs
            investment_cost = parameter(reference, bus_idx, :gen_bus, tech_idx, "Annual_INVC")[string(year_id)] * pu_econ_base
            storage_investment_cost = parameter(reference, bus_idx, :gen_bus, tech_idx, "Annual_STO_INV")[string(year_id)] * pu_econ_base

            # Get ITC 
            itc = 0.0
            if (setting["Simulation Configuration"]["ITC_Flag"] == true) && (parameter(reference, bus_idx, :gen_bus, tech_idx, "ITC Flag") == true)
                itc_year = min(2050, year_of_year_id)
                itc = parameter(reference, bus_idx, :gen_bus, tech_idx, "ITC")[string(itc_year)] # %
            end

            # Collect data
            idx_string = string("(", i, ", ", year_id, ")")
            u_new_G_iy = result["1"]["solution"]["expansion"][idx_string]["u_new_G_iy"]  
            u_new_ESH_iy = 0.0
            if parameter(reference, 0, :gen_index, "UNIT_CATEGORY", i) == "STORAGE"
                u_new_ESH_iy = result["1"]["solution"]["expansion"][idx_string]["u_new_ESH_iy"]
            end
            
            for future_y in 1:payment_duration

                future_year = year_of_year_id + future_y - 1    # Actual future year 
                discount_factor = (1 + discount_rate)^-(future_year - current_year) # Discount factor relative to the current dollor year

                gen_invc = 0.0
                
                # Expansion Cost (new units only) 
                gen_invc += 1000 * investment_cost * capacity * u_new_G_iy / (1.0 - itc)

                # Storage Duration Cost ($/kWh)
                if parameter(reference, 0, :gen_index, "UNIT_CATEGORY", i) == "STORAGE"
                    gen_invc += 1000 * storage_investment_cost * capacity * u_new_ESH_iy / (1.0 - itc)
                end

                # update the total investment cost
                total_gen_invc += discount_factor * gen_invc

                # Investment Tax Credits
                total_ITC += itc * discount_factor * gen_invc

                # update annual total_gen_invc_vector
                total_gen_invc_vector[(year_id-1)*stage_length + future_y] += discount_factor * gen_invc
                total_gen_ITC_vector[(year_id-1)*stage_length + future_y] += itc * discount_factor * gen_invc
            end
        end


        # Transmission investment costs
        total_upfront_invc_line = 0.0 # transmission expansion cost in each stage (thus multiplied by stage_length)
        total_upfront_fom_line = 0.0  # per-stage transmission FOM (existing baseline + expansion)
        
        trans_crpyears = setting["Planning Design"]["transmission_investment_CRP_value"]
        payment_duration = min(trans_crpyears, remaining_years)

        for k in ids_k
            
            rate_a = parameter(reference, 0, :branch, "rate_a", k) * pu_power_base
            length = get(reference[0][:branch][k], "length", get(reference[0][:branch][k], "Length", 0.0))
            transmission_expansion_cost = parameter(reference, 0, :branch, "transmission_expansion_cost", k) * pu_econ_base
            # DC-tie expansion cost is per-MW (no mile factor); AC lines are per-MW-mile
            len_factor = get(reference[0][:branch][k], "dc_line", false) == true ? 1.0 : length

            transmission_investment_cost_basis = transmission_expansion_cost * rate_a * len_factor
            
            for future_y in 1:payment_duration

                future_year = year_of_year_id + future_y - 1    # Actual future year 
                discount_factor = (1 + discount_rate)^-(future_year - current_year) # Discount factor relative to the current dollor year

                # Transmission expansion cost ($/mile)
                if parameter(reference, 0, :branch, "expansion_flag", k) == true
                    invc_line = transmission_investment_cost_basis * result["1"]["solution"]["expansion_line"][string("(", k, ", ", year_id, ")")]["u_new_T_ky"]
                    total_upfront_invc_line += discount_factor * invc_line
                    total_trans_invc_vector[(year_id-1)*stage_length + future_y] += discount_factor * invc_line
                end

            end

            # Transmission FOM on the full in-service grid: existing baseline (constant) + expansion (u_T_ky).
            # Charged every operating year of the stage, mirroring the rate_a*(1+u_T_ky) flow limit.
            transmission_fom_cost = parameter(reference, 0, :branch, "transmission_fom_cost", k) * pu_econ_base
            fom_basis = transmission_fom_cost * rate_a * len_factor
            u_T_ky_value = parameter(reference, 0, :branch, "expansion_flag", k) == true ? result["1"]["solution"]["expansion"][string("(", k, ", ", year_id, ")")]["u_T_ky"] : 0.0
            for future_y in 1:stage_length
                future_year = year_of_year_id + future_y - 1
                discount_factor = (1 + discount_rate)^-(future_year - current_year)
                fom_annual = discount_factor * fom_basis * (1 + u_T_ky_value)
                total_upfront_fom_line += fom_annual
                total_trans_fom_vector[(year_id-1)*stage_length + future_y] += fom_annual
            end
        end
       

        # Retirement costs for generators 
        total_retirement_cost = 0.0
        discount_factor_year_id = (1 + discount_rate)^-(year_of_year_id - current_year)

        for i in ids_i
            bus_idx = parameter(reference, 0, :gen_index, "bus_idx", i)
            tech_idx = parameter(reference, 0, :gen_index, "genco_tech_id", i)

            # Decommissioning Cost: k$/MW
            DECC = parameter(reference, bus_idx, :gen_bus, tech_idx, "DECC") * pu_econ_base

            if DECC > 0
                capacity = parameter(reference, bus_idx, :gen_bus, tech_idx, "CAP") * pu_power_base
                retirement_cost = 1000 * DECC * capacity * result["1"]["solution"]["expansion"][string("(", i, ", ", year_id, ")")]["u_ret_G_iy"]     
                total_retirement_cost += discount_factor_year_id * retirement_cost
                total_gen_retc_vector[(year_id-1)*stage_length + 1] += discount_factor_year_id * retirement_cost
            end
        end

        # FOM cost
        total_FOM = 0.0
        for i in ids_i
            bus_idx = parameter(reference, 0, :gen_index, "bus_idx", i)
            tech_idx = parameter(reference, 0, :gen_index, "genco_tech_id", i)

            capacity = parameter(reference, bus_idx, :gen_bus, tech_idx, "CAP") * pu_power_base
            FOM = parameter(reference, bus_idx, :gen_bus, tech_idx, "Annual_FOM")[string(year_id)] * pu_econ_base

            # Collect data
            u_G_iy = result["1"]["solution"]["expansion"][string("(", i, ", ", year_id, ")")]["u_G_iy"]  

            if (FOM > 0) && (capacity > 0)
        
                for future_y in 1:stage_length
                
                    future_year = year_of_year_id + future_y - 1 
                    discount_factor_year_id = (1 + discount_rate)^-(future_year - current_year)

                    # Fixed OM Cost ($/kW-Year)
                    total_FOM += discount_factor_year_id * 1000 * FOM * capacity * u_G_iy
                    total_gen_FOM_vector[(year_id-1)*stage_length + future_y] += discount_factor_year_id * 1000 * FOM * capacity * u_G_iy
                end
        
            end
        end

        # Generation, Ancillary Services Costs, carbon tax, and emissions
        ctax = setting["Simulation Configuration"]["CTAX"] * pu_econ_ref_base

        total_gen_cost = 0.0
        total_fuel_cost = 0.0
        total_VOM_cost = 0.0
        total_commitment_cost = 0.0

        total_regulation_cost = 0.0
        total_spin_cost = 0.0
        total_non_spin_cost = 0.0
        total_flex_cost = 0.0

        total_carbon_tax = 0.0
        total_emissions = 0.0

        total_generation = 0.0
        total_storage_charge_mwh = 0.0
        total_reserve_reg_up = 0.0
        total_reserve_reg_dn = 0.0
        total_reserve_spin = 0.0
        total_reserve_nspin = 0.0
        total_reserve_flex_up = 0.0
        total_reserve_flex_dn = 0.0

        for i in ids_i
    
            bus_idx = parameter(reference, 0, :gen_index, "bus_idx", i)
            tech_idx = parameter(reference, 0, :gen_index, "genco_tech_id", i)
            unit_group = parameter(reference, 0, :gen_index, "UNIT_GROUP", i)
            Annual_VOM = parameter(reference, bus_idx, :gen_bus, tech_idx, "Annual_VOM")[string(year_id)] * pu_econ_base 
            emission_rate = parameter(reference, bus_idx, :gen_bus, tech_idx, "Emission_CO2") / pu_power_base
            
            for future_y in 1:stage_length

                future_year = year_of_year_id + future_y - 1 
                discount_factor = (1 + discount_rate)^-(future_year - current_year)
                    
                # Generation 
                generation = annual_gen_info[i][year_id]["Generation"] / stage_length # annual_gen_info is pre-multiplied by stage_length*num_days; divide by stage_length for per-sub-year values
                geneneration_cost = annual_gen_info[i][year_id]["Generation_Cost"] / stage_length # annual_gen_info is pre-multiplied by stage_length*num_days; divide by stage_length for per-sub-year values
                fuel_cost = annual_gen_info[i][year_id]["Fuel_Cost"] / stage_length # annual_gen_info is pre-multiplied by stage_length*num_days; divide by stage_length for per-sub-year values
                VOM_cost = generation * Annual_VOM

                total_gen_cost += discount_factor * geneneration_cost
                total_fuel_cost += discount_factor * fuel_cost
                total_VOM_cost += discount_factor * VOM_cost
                total_generation += generation

                total_gen_cost_vector[(year_id-1)*stage_length + future_y] += discount_factor * geneneration_cost
                total_fuel_cost_vector[(year_id-1)*stage_length + future_y] += discount_factor * fuel_cost
                total_VOM_cost_vector[(year_id-1)*stage_length + future_y] += discount_factor * VOM_cost
                total_generation_vector[(year_id-1)*stage_length + future_y] += generation

                # Storage charge
                storage_charge_mwh = annual_gen_info[i][year_id]["Storage_Charge"] / stage_length # annual_gen_info is pre-multiplied by stage_length*num_days; divide by stage_length for per-sub-year values
                total_storage_charge_mwh += storage_charge_mwh
                total_storage_charge_mwh_vector[(year_id-1)*stage_length + future_y] += storage_charge_mwh

                # Ancillary Services Costs
                reserve_RU = annual_gen_info[i][year_id]["Reserve_RU"] / stage_length
                reserve_RD = annual_gen_info[i][year_id]["Reserve_RD"] / stage_length
                reserve_Spin = annual_gen_info[i][year_id]["Reserve_Spin"] / stage_length
                reserve_NSpin = annual_gen_info[i][year_id]["Reserve_NSpin"] / stage_length
                reserve_FU = annual_gen_info[i][year_id]["Reserve_FU"] / stage_length
                reserve_FD = annual_gen_info[i][year_id]["Reserve_FD"] / stage_length

                total_reserve_reg_up += reserve_RU
                total_reserve_reg_dn += reserve_RD
                total_reserve_spin += reserve_Spin
                total_reserve_nspin += reserve_NSpin
                total_reserve_flex_up += reserve_FU
                total_reserve_flex_dn += reserve_FD

                total_reg_up_vector[(year_id-1)*stage_length + future_y] += reserve_RU
                total_reg_dn_vector[(year_id-1)*stage_length + future_y] += reserve_RD
                total_spin_vector[(year_id-1)*stage_length + future_y] += reserve_Spin
                total_nspin_vector[(year_id-1)*stage_length + future_y] += reserve_NSpin
                total_flex_up_vector[(year_id-1)*stage_length + future_y] += reserve_FU
                total_flex_dn_vector[(year_id-1)*stage_length + future_y] += reserve_FD

                reg_cost = annual_gen_info[i][year_id]["Reserve_Reg_Cost"] / stage_length
                spin_cost = annual_gen_info[i][year_id]["Reserve_Spin_Cost"] / stage_length
                nspin_cost = annual_gen_info[i][year_id]["Reserve_NSpin_Cost"] / stage_length
                flex_cost = annual_gen_info[i][year_id]["Reserve_Flex_Cost"] / stage_length

                total_regulation_cost += discount_factor * reg_cost
                total_spin_cost += discount_factor * spin_cost
                total_non_spin_cost += discount_factor * nspin_cost
                total_flex_cost += discount_factor * flex_cost

                total_regulation_cost_vector[(year_id-1)*stage_length + future_y] += discount_factor * reg_cost
                total_spin_cost_vector[(year_id-1)*stage_length + future_y] += discount_factor * spin_cost
                total_non_spin_cost_vector[(year_id-1)*stage_length + future_y] += discount_factor * nspin_cost
                total_flex_cost_vector[(year_id-1)*stage_length + future_y] += discount_factor * flex_cost
                
                # Commitment Cost 
                noload_cost = annual_gen_info[i][year_id]["Noload_Cost"] / stage_length
                startup_cost = annual_gen_info[i][year_id]["StartUp_Cost"] / stage_length
                total_commitment_cost += discount_factor * (noload_cost + startup_cost)

                total_commitment_cost_vector[(year_id-1)*stage_length + future_y] += discount_factor * (noload_cost + startup_cost)
                
                # Carbon Tax and Emissions  
                total_carbon_tax += discount_factor * generation * emission_rate * ctax
                total_emissions += generation * emission_rate

                total_carbon_tax_vector[(year_id-1)*stage_length + future_y] += discount_factor * generation * emission_rate * ctax
                total_emissions_vector[(year_id-1)*stage_length + future_y] += generation * emission_rate

            end

        end


        # Generation Credit
        total_ptc = 0.0
        PTC_final_year_for_existing_assets = 5
        PTC_final_year_for_new_assets = 10

        if setting["Simulation Configuration"]["PTC_Flag"] == true
            for i in ids_i
                bus_idx = parameter(reference, 0, :gen_index, "bus_idx", i)
                tech_idx = parameter(reference, 0, :gen_index, "genco_tech_id", i)
                Tech_Type = parameter(reference, bus_idx, :gen_bus, tech_idx, "Tech_Type")
                unit_group = parameter(reference, 0, :gen_index, "UNIT_GROUP", i)
                profile_type = parameter(reference, 0, :gen_index, "Profile_Type", i)
                capacity = parameter(reference, bus_idx, :gen_bus, tech_idx, "CAP") * pu_power_base
                PMAX = parameter(reference, bus_idx, :gen_bus, tech_idx, "PMAX")
                
                if parameter(reference, bus_idx, :gen_bus, tech_idx, "PTC Flag") == true

                    if (Tech_Type == "Existing") && (unit_group != "nuclear")
                        
                        count = (year_id-1)*stage_length + 1
                    
                        for future_y in 1:stage_length

                            future_year = year_of_year_id + future_y - 1 
                            discount_factor = (1 + discount_rate)^-(future_year - current_year)
        
                            if count > PTC_final_year_for_existing_assets
                                break
                            end
        
                            ptc_year = min(2050, future_year)
                            ptc = parameter(reference, bus_idx, :gen_bus, tech_idx, "PTC")[string(ptc_year)] * 10 * pu_econ_base

                            # get generation
                            generation = annual_gen_info[i][year_id]["Generation"] / stage_length # annual_gen_info is pre-multiplied by stage_length*num_days; divide by stage_length for per-sub-year values

                            total_ptc += discount_factor * generation * ptc
                            total_ptc_vector[(year_id-1)*stage_length + future_y] += discount_factor * generation * ptc
        
                            count += 1
                        end

                    elseif (Tech_Type == "Existing") && (unit_group == "nuclear") 
                        
                        for future_y in 1:stage_length

                            future_year = year_of_year_id + future_y - 1 
                            discount_factor = (1 + discount_rate)^-(future_year - current_year)
        
                            ptc_year = min(2050, future_year)
                            ptc = parameter(reference, bus_idx, :gen_bus, tech_idx, "PTC")[string(ptc_year)] * 10 * pu_econ_base

                            # get generation
                            generation = annual_gen_info[i][year_id]["Generation"] / stage_length # annual_gen_info is pre-multiplied by stage_length*num_days; divide by stage_length for per-sub-year values

                            total_ptc += discount_factor * generation * ptc
                            total_ptc_vector[(year_id-1)*stage_length + future_y] += discount_factor * generation * ptc
                        end

                    elseif Tech_Type == "New" 
                        
                        ptc_year = min(2050, year_of_year_id)   # ptc is fixed to the investment year ptc
                        ptc = parameter(reference, bus_idx, :gen_bus, tech_idx, "PTC")[string(ptc_year)] * 10 * pu_econ_base 

                        ptc_duration = min(PTC_final_year_for_new_assets, remaining_years)

                        # Collect data
                        idx_string = string("(", i, ", ", year_id, ")")
                        u_new_G_iy = result["1"]["solution"]["expansion"][idx_string]["u_new_G_iy"]  
        
                        if u_new_G_iy > 0   # only if the unit is built
                            for future_y in 1:ptc_duration # for the (next 10 years or remaining years)

                                future_year = year_of_year_id + future_y - 1 
                                discount_factor = (1 + discount_rate)^-(future_year - current_year)
    
                                # Estimated annual generation in year y 
                                vre_data = reference[0][:bus][bus_idx]["vre_aggregated_data"][string(year_id)]
    
                                # Estimated generation per MW from the shape class (Profile_Type)
                                annual_generation_per_MW = 0.0
                                if profile_type != "NA"
                                    annual_generation_per_MW = vre_data[profile_type * "_shape"]
                                else    # nuclear
                                    annual_generation_per_MW = 8760 * PMAX # derated by FOR
                                end
    
                                total_ptc += discount_factor * annual_generation_per_MW * ptc * capacity * PMAX * u_new_G_iy
                                total_ptc_vector[(year_id-1)*stage_length + future_y] += discount_factor * annual_generation_per_MW * ptc * capacity * PMAX * u_new_G_iy
    
                            end
                        end
                        
                    
                    end
                end
            end
        end

        # Scarcity 

        total_scarcity_cost_ENS = 0.0
        total_scarcity_cost_SPIN = 0.0
        total_scarcity_cost_NSPIN = 0.0
        total_scarcity_cost_Flex = 0.0

        total_scarcity_ENS = 0.0
        total_scarcity_SPIN = 0.0
        total_scarcity_NSPIN = 0.0
        total_scarcity_FLEX = 0.0
        
        for future_y in 1:stage_length

            future_year = year_of_year_id + future_y - 1 
            discount_factor = (1 + discount_rate)^-(future_year - current_year)

            year_scarcity  = get(annual_scarcity_info, year_id, nothing)
            scarcity_E     = year_scarcity === nothing ? 0.0 : year_scarcity["scarcity_E"]
            scarcity_SPIN  = year_scarcity === nothing ? 0.0 : year_scarcity["scarcity_SPIN"]
            scarcity_NSPIN = year_scarcity === nothing ? 0.0 : year_scarcity["scarcity_NSPIN"]
            scarcity_FU    = year_scarcity === nothing ? 0.0 : year_scarcity["scarcity_FU"]
            scarcity_FD    = year_scarcity === nothing ? 0.0 : year_scarcity["scarcity_FD"]

            total_scarcity_cost_ENS += discount_factor * scarcity_E * setting["Simulation Configuration"]["VOLL"] * pu_econ_base
            total_scarcity_cost_SPIN += discount_factor * scarcity_SPIN * setting["Simulation Configuration"]["SRSP"] * pu_econ_base
            total_scarcity_cost_NSPIN += discount_factor * scarcity_NSPIN * setting["Simulation Configuration"]["NSRSP"] * pu_econ_base
            total_scarcity_cost_Flex += discount_factor * (scarcity_FU + scarcity_FD) * setting["Simulation Configuration"]["FLEXRSP"] * pu_econ_base

            total_scarcity_cost_ENS_vector[(year_id-1)*stage_length + future_y] = discount_factor * scarcity_E * setting["Simulation Configuration"]["VOLL"] * pu_econ_base
            total_scarcity_cost_SPIN_vector[(year_id-1)*stage_length + future_y] = discount_factor * scarcity_SPIN * setting["Simulation Configuration"]["SRSP"] * pu_econ_base
            total_scarcity_cost_NSPIN_vector[(year_id-1)*stage_length + future_y] = discount_factor * scarcity_NSPIN * setting["Simulation Configuration"]["NSRSP"] * pu_econ_base
            total_scarcity_cost_Flex_vector[(year_id-1)*stage_length + future_y] = discount_factor * (scarcity_FU + scarcity_FD) * setting["Simulation Configuration"]["FLEXRSP"] * pu_econ_base

            total_scarcity_ENS += scarcity_E
            total_scarcity_SPIN += scarcity_SPIN
            total_scarcity_NSPIN += scarcity_NSPIN
            total_scarcity_FLEX += (scarcity_FU + scarcity_FD)

            total_scarcity_ENS_vector[(year_id-1)*stage_length + future_y] = scarcity_E
            total_scarcity_SPIN_vector[(year_id-1)*stage_length + future_y] = scarcity_SPIN
            total_scarcity_NSPIN_vector[(year_id-1)*stage_length + future_y] = scarcity_NSPIN
            total_scarcity_FLEX_vector[(year_id-1)*stage_length + future_y] = (scarcity_FU + scarcity_FD)
            
        end

        # Planning Reserve Margin
        total_ucap = 0.0
        total_icap = 0.0
        for i in ids_i

            bus_idx = parameter(reference, 0, :gen_index, "bus_idx", i)
            tech_idx = parameter(reference, 0, :gen_index, "genco_tech_id", i)
            unit_group = parameter(reference, 0, :gen_index, "UNIT_GROUP", i)
            CAPCRED = parameter(reference, bus_idx, :gen_bus, tech_idx, "CAPCRED")
            if CAPCRED isa String
                data_identifier = reference[0][:zone]["capacity_credit"]["data_identifier"]
                data_region = reference[0][:bus][bus_idx]["region_mapping_info"][data_identifier]
                data_region = data_region isa AbstractString ? data_region : string(data_region)
                CAPCRED = reference[0][:zone]["capacity_credit"][data_region][CAPCRED]
            end

            capacity = parameter(reference, bus_idx, :gen_bus, tech_idx, "CAP") * pu_power_base

            # Collect data
            u_G_iy = result["1"]["solution"]["expansion"][string("(", i, ", ", year_id, ")")]["u_G_iy"]  

            total_ucap += capacity * u_G_iy * CAPCRED
            total_icap += capacity * u_G_iy
        end
        
        demand_mwh_year = get_annual_total_demand_with_growth(reference[0], year_of_year_id) * pu_power_base   # annual demand (MWh) of the stage's dispatched year
        # load actually dispatched: representative-day weighted (differs from the annual input demand)
        dispatched_load_year = let ys = get(annual_scarcity_info, year_id, nothing); ys === nothing ? 0.0 : get(ys, "Load_MWh", 0.0) end
        # PRM uses the coincident peak (per-unit); scale to MW to match total_ucap (MW).
        Peak_Demand = get_annual_coincident_peak_demand_with_growth(reference[0], year_of_year_id) * pu_power_base
        PRM = (total_ucap / (Peak_Demand)) - 1

        enforce_tnd_loss = setting["Planning Design"]["enforce_transmission_loss_flag"] == true
        tnd_loss_pct = setting["Planning Design"]["transmission_loss_percent_value"] * 0.01
        total_tnd_loss = 0.0

        for future_y in 1:stage_length

            future_year = year_of_year_id + future_y - 1
            peak_demand = get_annual_coincident_peak_demand_with_growth(reference[0], future_year) * pu_power_base

            total_prm_vector[(year_id-1)*stage_length + future_y] = (total_ucap / (peak_demand)) - 1
            total_peak_demand_vector[(year_id-1)*stage_length + future_y] = peak_demand
            total_ucap_vector[(year_id-1)*stage_length + future_y] = total_ucap
            total_icap_vector[(year_id-1)*stage_length + future_y] = total_icap
            total_demand_vector[(year_id-1)*stage_length + future_y] = demand_mwh_year
            total_dispatched_load_vector[(year_id-1)*stage_length + future_y] = dispatched_load_year

            # Annual nominal demand (MWh) at flat system-wide loss rate; zero unless the loss is enforced.
            if enforce_tnd_loss
                # loss on the stage-dispatched demand so the annual report energy-balances
                annual_demand = get_annual_total_demand_with_growth(reference[0], year_of_year_id) * pu_power_base
                tnd_loss = annual_demand * tnd_loss_pct
                total_tnd_loss_vector[(year_id-1)*stage_length + future_y] = tnd_loss
                total_tnd_loss += tnd_loss
            end
        end

        # Policy Slacks: CEGT 
        total_CEGT_slack = 0.0
        clean_energy_generation_penalty = setting["Simulation Configuration"]["Clean_Energy_Generation_Penalty"] * pu_econ_base
        
        if setting["Simulation Configuration"]["Clean_Energy_Generation_Target_EXP_Flag"] == true
            for future_y in 1:stage_length
                    
                future_year = year_of_year_id + future_y - 1 
                discount_factor_year_id = (1 + discount_rate)^-(future_year - current_year)

                CEGT_slack = 0.0

                if setting["Simulation Configuration"]["Clean_Energy_Generation_Target_EXP_Type"] == "Annual"
                    for n in ids_p
                        idx_string = string("(", n, ", ", year_id, ")")
                        CEGT_slack += discount_factor_year_id * clean_energy_generation_penalty * policy_slack_value(result["1"]["solution"], idx_string, "slack_CEG_ny") * pu_power_base                    
                    end

                elseif setting["Simulation Configuration"]["Clean_Energy_Generation_Target_EXP_Type"] == "Daygroup"
                    for n in ids_p
                        for d in ids_d
                            idx_string = string("(", n, ", ", d, ", ", year_id, ")")
                            CEGT_slack += discount_factor_year_id * clean_energy_generation_penalty * parameter(reference, 0, :repdays, "NumDays", d) * policy_slack_value(result["1"]["solution"], idx_string, "slack_CEG_ndy") * pu_power_base                                      
                        end
                    end
                end
                
                total_CEGT_slack += CEGT_slack
                total_CEGT_slack_vector[(year_id-1)*stage_length + future_y] += CEGT_slack
            end
        end

        # Policy Slacks: RPS 
        total_RPS_slack = 0.0
        if (setting["Simulation Configuration"]["RPS_Flag"] == true) && (setting["Simulation Configuration"]["Allow_Alternative_RPS_Compliance_Flag"] == true)

            rps_penalty = setting["Simulation Configuration"]["RPS_Penalty"] * pu_econ_base
        
            for future_y in 1:stage_length
                    
                future_year = year_of_year_id + future_y - 1 
                discount_factor_year_id = (1 + discount_rate)^-(future_year - current_year)

                RPS_slack = 0.0
                for n in ids_p
                    idx_string = string("(", n, ", ", year_id, ")")
                    RPS_slack += discount_factor_year_id * rps_penalty * policy_slack_value(result["1"]["solution"], idx_string, "slack_RPS_ny") * pu_power_base                        
                end
                
                total_RPS_slack += RPS_slack
                total_RPS_slack_vector[(year_id-1)*stage_length + future_y] += RPS_slack
            end
        end

        # Objective Value
        objective_value = 0.0
        try 
            objective_value = result["1"]["Objective Values"][string(year_id)] * setting["Simulation Setting"]["per_unit_econ_base_value"]
        catch
            objective_value = result["1"]["Objective Values"] * setting["Simulation Setting"]["per_unit_econ_base_value"]  
        end

        if isa(objective_value, Dict)
            objective_value = 0.0   
        end

        # calculate PTC of year_id 
        PTC_of_year_y = 0.0   # we need to calculate the PTC for the year_id because the total_PTC value includes PTCs from all the future years (10 years)
        for future_y in 1:stage_length
            PTC_of_year_y += total_ptc_vector[(year_id-1)*stage_length + future_y]
        end

        total_system_cost = total_gen_invc + total_upfront_invc_line + total_upfront_fom_line + total_retirement_cost + total_FOM
        total_system_cost += total_gen_cost + total_commitment_cost + total_regulation_cost + total_spin_cost + total_non_spin_cost + total_flex_cost 
        total_system_cost += total_scarcity_cost_ENS + total_scarcity_cost_SPIN + total_scarcity_cost_NSPIN + total_scarcity_cost_Flex 
        total_system_cost += total_carbon_tax + total_CEGT_slack + total_RPS_slack
        total_system_cost -= total_ITC
        total_system_cost -= PTC_of_year_y
        
        # update system output list
        system_output_list[year_id + 1, :] = [
            Scenario, 
            year_id, 
            year_of_year_id, 
            stage_length, 
            total_gen_invc, 
            -total_ITC, 
            total_upfront_invc_line,
            total_upfront_fom_line,
            total_retirement_cost, 
            total_FOM,
            -PTC_of_year_y, 
            total_fuel_cost, 
            total_VOM_cost, 
            total_commitment_cost, 
            total_regulation_cost, 
            total_spin_cost, 
            total_non_spin_cost, 
            total_flex_cost, 
            total_scarcity_cost_ENS, 
            total_scarcity_cost_SPIN, 
            total_scarcity_cost_NSPIN, 
            total_scarcity_cost_Flex, 
            total_carbon_tax, 
            total_CEGT_slack, 
            total_RPS_slack, 
            total_system_cost,
            objective_value,
            total_generation,
            total_tnd_loss,
            total_storage_charge_mwh,
            total_reserve_reg_up,
            total_reserve_reg_dn,
            total_reserve_spin,
            total_reserve_nspin,
            total_reserve_flex_up,
            total_reserve_flex_dn,
            total_scarcity_ENS,
            total_scarcity_SPIN,
            total_scarcity_NSPIN,
            total_scarcity_FLEX,
            total_emissions,
            PRM,
            Peak_Demand,
            total_ucap,
            total_icap,
            demand_mwh_year * stage_length,
            dispatched_load_year * stage_length
        ]

        # update extensive system output list
        for future_y in 1:stage_length

            year_pointer = (year_id-1)*stage_length + future_y

            total_system_cost = total_gen_invc_vector[year_pointer] +
                    total_trans_invc_vector[year_pointer] +
                    total_trans_fom_vector[year_pointer] +
                    total_gen_retc_vector[year_pointer] +
                    total_gen_FOM_vector[year_pointer] +
                    total_gen_cost_vector[year_pointer] +
                    total_commitment_cost_vector[year_pointer] +
                    total_regulation_cost_vector[year_pointer] +
                    total_spin_cost_vector[year_pointer] +
                    total_non_spin_cost_vector[year_pointer] +
                    total_flex_cost_vector[year_pointer] +
                    total_scarcity_cost_ENS_vector[year_pointer] +
                    total_scarcity_cost_SPIN_vector[year_pointer] +
                    total_scarcity_cost_NSPIN_vector[year_pointer] +
                    total_scarcity_cost_Flex_vector[year_pointer] +
                    total_carbon_tax_vector[year_pointer] +
                    total_CEGT_slack_vector[year_pointer] +
                    total_RPS_slack_vector[year_pointer] -
                    total_gen_ITC_vector[year_pointer] -
                    total_ptc_vector[year_pointer]
            
            extensive_objective_value = 0.0
            if (objective_value > 0) && (future_y == 1)
                extensive_objective_value = objective_value
            end

            extensive_system_output_list[year_pointer + 1, :] = [
                Scenario, 
                year_id, 
                year_of_year_id + future_y - 1,
                1, # number of years is always 1 for extensive output
                total_gen_invc_vector[year_pointer], 
                -total_gen_ITC_vector[year_pointer], 
                total_trans_invc_vector[year_pointer],
                total_trans_fom_vector[year_pointer],
                total_gen_retc_vector[year_pointer], 
                total_gen_FOM_vector[year_pointer],
                -total_ptc_vector[year_pointer], 
                total_fuel_cost_vector[year_pointer], 
                total_VOM_cost_vector[year_pointer], 
                total_commitment_cost_vector[year_pointer], 
                total_regulation_cost_vector[year_pointer], 
                total_spin_cost_vector[year_pointer], 
                total_non_spin_cost_vector[year_pointer], 
                total_flex_cost_vector[year_pointer],
                total_scarcity_cost_ENS_vector[year_pointer], 
                total_scarcity_cost_SPIN_vector[year_pointer], 
                total_scarcity_cost_NSPIN_vector[year_pointer], 
                total_scarcity_cost_Flex_vector[year_pointer],
                total_carbon_tax_vector[year_pointer], 
                total_CEGT_slack_vector[year_pointer], 
                total_RPS_slack_vector[year_pointer], 
                total_system_cost,
                extensive_objective_value,
                total_generation_vector[year_pointer],
                total_tnd_loss_vector[year_pointer],
                total_storage_charge_mwh_vector[year_pointer],
                total_reg_up_vector[year_pointer],
                total_reg_dn_vector[year_pointer],
                total_spin_vector[year_pointer],
                total_nspin_vector[year_pointer],
                total_flex_up_vector[year_pointer],
                total_flex_dn_vector[year_pointer],
                total_scarcity_ENS_vector[year_pointer],
                total_scarcity_SPIN_vector[year_pointer],
                total_scarcity_NSPIN_vector[year_pointer],
                total_scarcity_FLEX_vector[year_pointer],
                total_emissions_vector[year_pointer],
                total_prm_vector[year_pointer],
                total_peak_demand_vector[year_pointer],
                total_ucap_vector[year_pointer],
                total_icap_vector[year_pointer],
                total_demand_vector[year_pointer],
                total_dispatched_load_vector[year_pointer]
            ]
        end

    end

    # Write the system output list to a csv file
    system_output_list, extensive_system_output_list = finalize_system_summaries(system_output_list, extensive_system_output_list, current_year)
    writedlm(system_file_name, system_output_list, ",")
    writedlm(extensive_system_file_name, extensive_system_output_list, ",")
    write_constant_dollar_annual_summary(extensive_system_file_name, extensive_system_output_list, discount_rate, current_year)

    @aleaf_info "[ALEAF LC_GTEP]: - System summary reporting,"
end


# Accessors that read day-group-invariant build decisions from either the full in-memory solution
# (serial path, dg "1") or the reduced `build_decisions` payload (distributed path).
function _op_expansion_decisions(result, year_id)
    if get(result, "pu_applied", false) == true
        return result["build_decisions"][string(year_id)]["expansion"]
    end
    return result["operational model result"][string(year_id)]["1"]["solution"]["expansion"]
end

function _op_slack_decisions(result, year_id)
    if get(result, "pu_applied", false) == true
        return result["build_decisions"][string(year_id)]["slack"]
    end
    return result["operational model result"][string(year_id)]["1"]["solution"]["slack"]
end


function report_result_tech_summary_operation_GTEP(result, result_LC_GTEP_expansion)

    reference = result["operational model system reference"]
    setting = result["setting"]

    pu_power_base = setting["Simulation Setting"]["per_unit_base_value"]
    pu_econ_base = setting["Simulation Setting"]["per_unit_econ_base_value"] / pu_power_base    # per unit economic base

    output_path = parameter(reference["1"], 0, :output_path)
    tech_file_name = joinpath(output_path, string(parameter(reference["1"], 0, :case_name), "__tech_summary_by_stage_OP.csv"))
    
    # Write Label
    tech_label_list = [
        "Scenario",
        "Stage",
        "Start Year",
        "Number of Years",
        "PLANT_NAME",
        "Bus_ID",
        "Bus_Name",
        "Parent_Bus_Name",
        "Region_Name",
        "Tech_ID",
        "UnitGroup",
        "Unit_Category",
        "Unit_Report_Label_1",
        "Unit_Report_Label_2",
        "Fuel",
        "TotalUnits",
        "NewUnits",
        "RetUnits",
        "ICAP",
        "UCAP",
        "ICap_New", 
        "ICap_Ret",
        "UCap_New",
        "UCap_Ret", 
        "Storage_MWh", 
        "Storage_Hr", 
        "Generation", 
        "Curtail", 
        "Storage_Charge_MWh", 
        "Reserve_RegUp", 
        "Reserve_RegDn", 
        "Reserve_Spin", 
        "Reserve_NSpin",
        "Reserve_FlexUp", 
        "Reserve_FlexDn", 
        "Generation_Cost", 
        "Charge_Cost", 
        "Regulation_Cost", 
        "Spin_Cost", 
        "Nspin_Cost", 
        "Flex_Cost", 
        "UnitRevenue_E", 
        "UnitRevenue_AS", 
        "UnitRevenue_CRED", 
        "UnitRevenue", 
        "UnitProfit", 
        "FuelConsumption", 
        "FuelCost", 
        "FOM", 
        "CAPCRED", 
        "Reference_Annual_Gen_Investment_Cost",	
        "Gen_Investment_Cost",	
        "Gen_ITC"
        ]
     
    ids_y = []
    ids_i = []
    ids_n = []
    ids_z = []
    ids_p = []
    ids_k = []
    ids_d = []
    annual_gen_info = reference["annual_gen_info_OP"]
    try
        ids_y = [(y) for (y) in sort!(parse.(Int, (collect(get_index(reference["1"], :planning_stages, 0)))))]
        ids_i = [(k) for (k) in sort!(parse.(Int, (collect(get_index(reference["1"], :gen_index, 0)))))]
        ids_n = [(k) for (k) in sort!(parse.(Int, (collect(get_index(reference["1"], :bus, 0)))))]
        ids_z = [(i) for (i) in sort!(parse.(Int, (collect(keys(reference["1"][0][:zone]["reserve"])))))]
        ids_p = [(i) for (i) in sort!(parse.(Int, (collect(keys(reference["1"][0][:zone]["policy"])))))]
        ids_d = [(k) for (k) in sort!(parse.(Int, (collect(get_index(reference["1"], :repdays, 0)))))]
        ids_k = [(k) for (k) in sort!(parse.(Int, (collect(get_index(reference["1"], :branch, 0))))) if parameter(reference["1"], 0, :branch, "model_flag", k) == true]
    catch
        ids_y = [(y) for (y) in sort!(collect(get_index(reference["1"], :planning_stages, 0)))]
        ids_i = [(k) for (k) in sort!(collect(get_index(reference["1"], :gen_index, 0)))]
        ids_n = [(k) for (k) in sort!(collect(get_index(reference["1"], :bus, 0)))]
        ids_z = [(i) for (i) in sort!(parse.(Int, (collect(keys(reference["1"][0][:zone]["reserve"])))))]
        ids_p = [(i) for (i) in sort!(parse.(Int, (collect(keys(reference["1"][0][:zone]["policy"])))))]
        ids_d = [(k) for (k) in sort!(collect(get_index(reference["1"], :repdays, 0)))]
        ids_k = [(k) for (k) in sort!(collect(get_index(reference["1"], :branch, 0))) if parameter(reference["1"], 0, :branch, "model_flag", k) == true]
    end
    # __market_OP.csv is unused here; on the distributed path it does not exist yet (concatenated
    # after summaries), so skip the read entirely.
    if get(result, "pu_applied", false) != true
        market_outcome_df = DataFrame(CSV.File(joinpath(output_path, string(parameter(reference["1"], 0, :case_name), "__market_OP.csv"))))
    end

    # total number of techs
    num_tech = 0
    num_bus_tech = 0
    num_tech_per_bus = [0]
    for bus_idx in ids_n
        try
            num_tech += length(keys(reference["1"][string(bus_idx)]["gen_bus"]))
            num_bus_tech += length(keys(reference["1"][string(bus_idx)]["gen_bus"]))
        catch
            num_tech += length(keys(reference["1"][bus_idx][:gen_bus]))
            num_bus_tech += length(keys(reference["1"][bus_idx][:gen_bus]))
        end
        push!(num_tech_per_bus, num_bus_tech)
    end

    tech_output_list = Array{Any}(undef, length(ids_y)*num_tech+1, length(tech_label_list))
    tech_output_list[1,:] = tech_label_list

    # parameters
    Scenario = parameter(reference["1"], 0, :case_name)
    current_year = setting["Planning Design"]["dollar_year_value"] # Discount factor relative to the current dollor year
    discount_rate = setting["Planning Design"]["discount_rate_value"]   
    total_day_groups = setting["Simulation Configuration"]["NDAY_Groups_OP"]
    
    for year_id in eachindex(ids_y)
        
        reference_of_this_year = reference[string(year_id)]
        
        year_of_year_id = parameter(reference_of_this_year, 0, :planning_stages, "year", year_id)
        stage_length = parameter(reference_of_this_year, 0, :planning_stages, "stage_length", year_id) 
        
        num_stage = setting["Planning Design"]["num_stages_value"]
        remaining_years = (num_stage + 1 - year_id) * stage_length

        # Plant level data collection                
        @sync for bus_idx in ids_n

            @spawn begin

                Bus_Name = parameter(reference_of_this_year, 0, :bus, "bus_i", bus_idx)
                Region_Name = reference_of_this_year[0][:bus][bus_idx]["region_config"]["region_name"]
                Parent_Bus_Name = reference_of_this_year[0][:bus][bus_idx]["region_config"]["parent_bus_name"]

                ids_tech = []
                try
                    ids_tech = [(k) for (k) in sort!(parse.(Int, (collect( keys(reference_of_this_year[string(bus_idx)]["gen_bus"])))))]
                catch
                    ids_tech = [(k) for (k) in sort!(collect(keys(reference_of_this_year[bus_idx][:gen_bus])))]
                end

                ids_tech_of_this_bus = [i for i in 1:length(ids_tech)]

                for tech_id in ids_tech_of_this_bus

                    tech_idx = ids_tech[tech_id]
                    
                    gen_idx = parameter(reference_of_this_year, bus_idx, :gen_bus, tech_idx, "gen_idx")
                    unitgroup = parameter(reference_of_this_year, bus_idx, :gen_bus, tech_idx, "UNITGROUP")
                    Unit_Category = parameter(reference_of_this_year, bus_idx, :gen_bus, tech_idx, "UNIT_CATEGORY")
                    Unit_Report_Label_1 = parameter(reference_of_this_year, bus_idx, :gen_bus, tech_idx, "UNIT_REPORT_LABEL_1")
                    Unit_Report_Label_2 = parameter(reference_of_this_year, bus_idx, :gen_bus, tech_idx, "UNIT_REPORT_LABEL_2")

                    # plant name: based on NEW_Gen_UID or PLANT_NAME
                    plant_name = parameter(reference_of_this_year, bus_idx, :gen_bus, tech_idx, "NEW_Gen_UID")
                    if haskey(reference_of_this_year[bus_idx][:gen_bus][tech_idx], "PLANT_NAME")
                        plant_name = parameter(reference_of_this_year, bus_idx, :gen_bus, tech_idx, "PLANT_NAME")
                    end
                    
                    Fuel = parameter(reference_of_this_year, bus_idx, :gen_bus, tech_idx, "FUEL")
                    Emission_rate = parameter(reference_of_this_year, bus_idx, :gen_bus, tech_idx, "Emission_CO2") / pu_power_base
                    capacity = parameter(reference_of_this_year, bus_idx, :gen_bus, tech_idx, "CAP") * pu_power_base
                    FOM = parameter(reference_of_this_year, bus_idx, :gen_bus, tech_idx, "Annual_FOM")[string(year_id)] * pu_econ_base
                    CRP_i = parameter(reference_of_this_year, bus_idx, :gen_bus, tech_idx, "crpyears")

                    CAPCRED = parameter(reference_of_this_year, bus_idx, :gen_bus, tech_idx, "CAPCRED")
                    if CAPCRED isa String
                        data_identifier = reference_of_this_year[0][:zone]["capacity_credit"]["data_identifier"]
                        data_region = reference_of_this_year[0][:bus][bus_idx]["region_mapping_info"][data_identifier]
                        data_region = data_region isa AbstractString ? data_region : string(data_region)
                        CAPCRED = reference_of_this_year[0][:zone]["capacity_credit"][data_region][CAPCRED]
                    end

                    # get updated CAPCRED if available
                    if year_id > 1 # only check if this dict is not empty (i.e., this dict will be empty in year 1)
                        if (setting["Planning Design"]["multi_round_solution_process_flag"] == true) && (setting["Simulation Configuration"]["update_CAPCRED_in_each_round_of_Expansion_Flag"] == true)
                            if parameter(reference_of_this_year, bus_idx, :gen_bus, tech_idx, "ELCC_Flag") == true
                                try
                                    CAPCRED = parameter(reference_of_this_year, 0, :multi_round_info, "RA_Info", year_id)["ELCC_result_applied"][unitgroup]
                                catch
                                    @aleaf_info "[ALEAF LC_GTEP]: Failed to find a new capacity credit of $unitgroup in the report_result_tech_summary_expansion_GTEP function"
                                end
                            end
                        end
                    end

                    # Collect data
                    idx_string = string("(", gen_idx, ", ", year_id, ")")

                    tag = string("Ret_", year_of_year_id)
                    planned_retirement = parameter(reference_of_this_year, bus_idx, :gen_bus, tech_idx, "Planned_Retirement")[tag]

                    exp_dec = _op_expansion_decisions(result, year_id)
                    TotalUnits = exp_dec[idx_string]["u_G_iy"]
                    NewUnits = 0.0
                    RetUnits = 0.0
                    if length(result_LC_GTEP_expansion) > 0
                        NewUnits = result_LC_GTEP_expansion["solution"]["expansion"][idx_string]["u_new_G_iy"]  
                        RetUnits = result_LC_GTEP_expansion["solution"]["expansion"][idx_string]["u_ret_G_iy"]
                    end
                    
                    ICAP = capacity * TotalUnits 
                    UCAP = CAPCRED * ICAP

                    ICap_New = capacity * NewUnits
                    ICap_Ret = capacity * RetUnits
                    UCap_New = CAPCRED * ICap_New
                    UCap_Ret = CAPCRED * ICap_Ret

                    Storage_Hr = 0.0
                    u_ESE = 0.0
                    u_new_ESH = 0.0
                    if Unit_Category == "STORAGE"
                        if TotalUnits > 0
                            u_ESE = exp_dec[idx_string]["u_ESE_iy"]
                            # MW floor avoids divide-by-near-zero on round-off-tiny storage units
                            if ICAP > 1e-3
                                Storage_Hr = u_ESE / ICAP
                            end
                            if length(result_LC_GTEP_expansion) > 0
                                u_new_ESH = result_LC_GTEP_expansion["solution"]["expansion"][idx_string]["u_new_ESH_iy"]
                            end
                        end
                    end

                    ######### Investment and Retirement costs 
                    
                    # Determine the cost duration (investment payments)
                    total_gen_invc = 0.0
                    total_ITC = 0.0
                    total_retirement_cost = 0.0

                    payment_duration = min(CRP_i, remaining_years)

                    # Get the investment costs
                    investment_cost = parameter(reference_of_this_year, bus_idx, :gen_bus, tech_idx, "Annual_INVC")[string(year_id)] * pu_econ_base
                    storage_investment_cost = parameter(reference_of_this_year, bus_idx, :gen_bus, tech_idx, "Annual_STO_INV")[string(year_id)] * pu_econ_base
                    DECC = parameter(reference_of_this_year, bus_idx, :gen_bus, tech_idx, "DECC") * pu_econ_base

                    # Get ITC 
                    itc = 0.0
                    if (setting["Simulation Configuration"]["ITC_Flag"] == true) && (parameter(reference_of_this_year, bus_idx, :gen_bus, tech_idx, "ITC Flag") == true)
                        itc_year = min(2050, year_of_year_id)
                        itc = parameter(reference_of_this_year, bus_idx, :gen_bus, tech_idx, "ITC")[string(itc_year)] # %
                    end

                    # Collect data
                    u_new_G_iy = NewUnits
                    u_new_ESH_iy = u_new_ESH

                    # base investment cost
                    base_INVC = 1000 * investment_cost * capacity * u_new_G_iy / (1.0 - itc)
                    if Unit_Category == "STORAGE"
                        base_INVC += 1000 * storage_investment_cost * capacity * u_new_ESH_iy / (1.0 - itc)
                    end

                    for future_y in 1:payment_duration

                        future_year = year_of_year_id + future_y - 1    # Actual future year 
                        discount_factor = (1 + discount_rate)^-(future_year - current_year) # Discount factor relative to the current dollor year
        
                        gen_invc = 0.0
                        
                        # Expansion Cost (new units only) 
                        gen_invc += 1000 * investment_cost * capacity * u_new_G_iy / (1.0 - itc)
        
                        # Storage Duration Cost ($/kWh)
                        if Unit_Category == "STORAGE"
                            gen_invc += 1000 * storage_investment_cost * capacity * u_new_ESH_iy / (1.0 - itc)
                        end
        
                        # update the total investment cost
                        total_gen_invc += discount_factor * gen_invc
        
                        # Investment Tax Credits
                        total_ITC += itc * discount_factor * gen_invc

                    end

                    # FOM
                    total_FOM = 0.0
                    if (FOM > 0) && (capacity > 0)
        
                        for future_y in 1:stage_length
                        
                            future_year = year_of_year_id + future_y - 1 
                            discount_factor_year_id = (1 + discount_rate)^-(future_year - current_year)
        
                            # Fixed OM Cost ($/kW-Year)
                            total_FOM += discount_factor_year_id * 1000 * FOM * ICAP
                        end
                    end
                                        
                    # Retirement Cost
                    if DECC > 0
                        discount_factor_year_id = (1 + discount_rate)^-(year_of_year_id - current_year)
                        retirement_cost = 1000 * DECC * capacity * RetUnits
                        total_retirement_cost = discount_factor_year_id * retirement_cost
                    end

                    # Get pre-calculated annual values from the dispatch output (already multiplied by stage_length)
                    Generation = annual_gen_info[gen_idx][year_id]["Generation"]
                    Curtail = annual_gen_info[gen_idx][year_id]["Curtail"]

                    Storage_Charge = annual_gen_info[gen_idx][year_id]["Storage_Charge"]
                    Reserve_RU = annual_gen_info[gen_idx][year_id]["Reserve_RU"]
                    Reserve_RD = annual_gen_info[gen_idx][year_id]["Reserve_RD"]
                    Reserve_Spin = annual_gen_info[gen_idx][year_id]["Reserve_Spin"]
                    Reserve_NSpin = annual_gen_info[gen_idx][year_id]["Reserve_NSpin"]
                    Reserve_FU = annual_gen_info[gen_idx][year_id]["Reserve_FU"]
                    Reserve_FD = annual_gen_info[gen_idx][year_id]["Reserve_FD"]

                    UnitRevenue_E = annual_gen_info[gen_idx][year_id]["UnitRevenue_E"]
                    UnitRevenue_AS = annual_gen_info[gen_idx][year_id]["UnitRevenue_AS"]
                    UnitRevenue_CRED = annual_gen_info[gen_idx][year_id]["UnitRevenue_CRED"]
                    UnitRevenue = annual_gen_info[gen_idx][year_id]["UnitRevenue"]

                    Noload_Cost = annual_gen_info[gen_idx][year_id]["Noload_Cost"]
                    StartUp_Cost = annual_gen_info[gen_idx][year_id]["StartUp_Cost"]

                    Generation_Cost = annual_gen_info[gen_idx][year_id]["Generation_Cost"]
                    Charge_Cost = annual_gen_info[gen_idx][year_id]["Charge_Cost"]

                    Reserve_Reg_Cost = annual_gen_info[gen_idx][year_id]["Reserve_Reg_Cost"]
                    Reserve_Spin_Cost = annual_gen_info[gen_idx][year_id]["Reserve_Spin_Cost"]
                    Reserve_NSpin_Cost = annual_gen_info[gen_idx][year_id]["Reserve_NSpin_Cost"]
                    Reserve_Flex_Cost = annual_gen_info[gen_idx][year_id]["Reserve_Flex_Cost"]

                    Fuel_Cost = annual_gen_info[gen_idx][year_id]["Fuel_Cost"]
                    Fuel_Consumption = annual_gen_info[gen_idx][year_id]["Fuel_Consumption"]

                    # Calculate unit revenue
                    total_revenue_E = 0.0
                    total_revenue_AS = 0.0
                    total_revenue_CRED = 0.0
                    total_revenue = 0.0
                    for future_y in 1:stage_length
                        future_year = year_of_year_id + future_y - 1 
                        discount_factor = (1 + discount_rate)^-(future_year - current_year)

                        total_revenue_E += discount_factor * (UnitRevenue_E / stage_length)
                        total_revenue_AS += discount_factor * (UnitRevenue_AS / stage_length)
                        total_revenue_CRED += discount_factor * (UnitRevenue_CRED / stage_length)
                        total_revenue += discount_factor * (UnitRevenue / stage_length)
                    end

                    # Calculate Unit Profit
                    total_profit = 0.0
                    UnitProfit = annual_gen_info[gen_idx][year_id]["UnitProfit"] - Noload_Cost - StartUp_Cost - FOM # (already multiplied by stage_length)
                    for future_y in 1:stage_length
                        future_year = year_of_year_id + future_y - 1 
                        discount_factor = (1 + discount_rate)^-(future_year - current_year)

                        total_profit += discount_factor * (UnitProfit / stage_length)
                    end

                    # Gen and Charge Cost 
                    total_gen_cost = 0.0
                    total_charge_cost = 0.0
                    for future_y in 1:stage_length
                        future_year = year_of_year_id + future_y - 1 
                        discount_factor = (1 + discount_rate)^-(future_year - current_year)

                        total_gen_cost += discount_factor * (Generation_Cost / stage_length)
                        total_charge_cost += discount_factor * (Charge_Cost / stage_length)
                    end

                    # AS Cost 
                    total_reg_cost = 0.0
                    total_spin_cost = 0.0
                    total_nspin_cost = 0.0
                    total_flex_cost = 0.0

                    for future_y in 1:stage_length
                        future_year = year_of_year_id + future_y - 1 
                        discount_factor = (1 + discount_rate)^-(future_year - current_year)

                        total_reg_cost += discount_factor * (Reserve_Reg_Cost / stage_length)
                        total_spin_cost += discount_factor * (Reserve_Spin_Cost / stage_length)
                        total_nspin_cost += discount_factor * (Reserve_NSpin_Cost / stage_length)
                        total_flex_cost += discount_factor * (Reserve_Flex_Cost / stage_length)

                    end

                    # Fuel Cost
                    total_fuel_cost = 0.0
                    for future_y in 1:stage_length
                        future_year = year_of_year_id + future_y - 1 
                        discount_factor = (1 + discount_rate)^-(future_year - current_year)

                        total_fuel_cost += discount_factor * (Fuel_Cost / stage_length)
                    end

                    row_id = (year_id-1)*num_tech + num_tech_per_bus[bus_idx] + tech_id + 1


                    if typeof(ids_tech[1]) == String
                        tech_idx = string(tech_idx)
                    end

                    # Update the list
                    tech_output_list[row_id, :] = [
                        Scenario, 
                        year_id, 
                        year_of_year_id, 
                        stage_length, 
                        plant_name,
                        bus_idx, 
                        Bus_Name, 
                        Parent_Bus_Name, 
                        Region_Name, 
                        tech_idx, 
                        unitgroup,
                        Unit_Category,
                        Unit_Report_Label_1,
                        Unit_Report_Label_2,
                        Fuel,
                        TotalUnits, 
                        NewUnits, 
                        RetUnits, 
                        ICAP, 
                        UCAP, 
                        ICap_New, 
                        ICap_Ret, 
                        UCap_New, 
                        UCap_Ret, 
                        u_ESE, 
                        Storage_Hr, 
                        Generation, 
                        Curtail, 
                        Storage_Charge, 
                        Reserve_RU, 
                        Reserve_RD, 
                        Reserve_Spin, 
                        Reserve_NSpin, 
                        Reserve_FU, 
                        Reserve_FD, 
                        total_gen_cost, 
                        total_charge_cost, 
                        total_reg_cost, 
                        total_spin_cost, 
                        total_nspin_cost, 
                        total_flex_cost,
                        total_revenue_E, 
                        total_revenue_AS, 
                        total_revenue_CRED, 
                        total_revenue, 
                        total_profit, 
                        Fuel_Consumption, 
                        total_fuel_cost, 
                        total_FOM, 
                        CAPCRED, 
                        base_INVC, 
                        total_gen_invc, 
                        total_ITC
                        ]
                                                
                end
            end
        end
    end
    writedlm(tech_file_name, finalize_tech_summary(tech_output_list), ",")

    @aleaf_info "[ALEAF LC_GTEP]: - Tech Summary reporting,"
end



function report_result_system_summary_operation_GTEP(result, result_LC_GTEP_expansion)

    reference = result["operational model system reference"]
    setting = result["setting"]

    pu_power_base = setting["Simulation Setting"]["per_unit_base_value"]
    pu_econ_ref_base = setting["Simulation Setting"]["per_unit_econ_base_value"]
    pu_econ_base = setting["Simulation Setting"]["per_unit_econ_base_value"] / pu_power_base    # per unit economic base

    output_path = parameter(reference["1"], 0, :output_path)
    system_file_name = joinpath(output_path, string(parameter(reference["1"], 0, :case_name), "__system_summary_by_stage_OP.csv"))
    extensive_system_file_name = joinpath(output_path, string(parameter(reference["1"], 0, :case_name), "__system_summary_by_year_OP.csv"))
    
    # Write Label 
    system_label_list = [
        "Scenario",
        "Stage",
        "Start Year",
        "Number of Years",
        "Gen_Investment_Cost",
        "Gen_ITC",
        "Trans_Investment_Cost",
        "Trans_FOM_Cost",
        "Gen_Retirement_Cost",
        "FOM_Cost",
        "Generation_PTC",
        "Fuel_Cost",
        "VOM_Cost",
        "Commitment_Cost",
        "Regulation_Cost",
        "Spin_Cost",
        "Nspin_Cost",
        "Flex_Cost",
        "ENS_Cost",
        "RNS_Spin_Cost",
        "RNS_NSpin_Cost",
        "RNS_Flex_Cost",
        "CTAX_cost",
        "CEGT_Penalty",
        "RPS_Penalty",
        "total_system_cost",
        "Total_system_cost_OP",
        "ObjectiveValue",
        "Generation",
        "TnD_Loss",
        "Storage_Charge_MWh",
        "Reserve_RegUp",
        "Reserve_RegDn",
        "Reserve_Spin",
        "Reserve_NSpin",
        "Reserve_FlexUp",
        "Reserve_FlexDn",
        "ENS",
        "RNS_Spin",
        "RNS_NSpin",
        "RNS_Flex",
        "Emission",
        "PRM",
        "Peak_Demand_MW",
        "UCAP_MW",
        "Installed_Capacity_MW",
        "Annual_Input_Demand_MWh",
        "Dispatched_Load_MWh",
    ]

    extensive_system_label_list = [
        "Scenario",
        "Stage",
        "Year",
        "Number of Years",
        "Gen_Investment_Cost_CF",
        "Gen_ITC_CF",
        "Trans_Investment_Cost_CF",
        "Trans_FOM_Cost_CF",
        "Gen_Retirement_Cost_CF",
        "FOM_Cost_CF",
        "Generation_PTC_CF",
        "Fuel_Cost_CF",
        "VOM_Cost_CF",
        "Commitment_Cost_CF",
        "Regulation_Cost_CF",
        "Spin_Cost_CF",
        "Nspin_Cost_CF",
        "Flex_Cost_CF",
        "ENS_Cost_CF",
        "RNS_Spin_Cost_CF",
        "RNS_NSpin_Cost_CF",
        "RNS_Flex_Cost_CF",
        "CTAX_cost_CF",
        "CEGT_Penalty_CF",
        "RPS_Penalty_CF",
        "total_system_cost_CF",
        "total_system_cost_OP_CF",
        "ObjectiveValue",
        "Generation",
        "TnD_Loss",
        "Storage_Charge_MWh",
        "Reserve_RegUp",
        "Reserve_RegDn",
        "Reserve_Spin",
        "Reserve_NSpin",
        "Reserve_FlexUp",
        "Reserve_FlexDn",
        "ENS",
        "RNS_Spin",
        "RNS_NSpin",
        "RNS_Flex",
        "Emission",
        "PRM",
        "Peak_Demand_MW",
        "UCAP_MW",
        "Installed_Capacity_MW",
        "Annual_Input_Demand_MWh",
        "Dispatched_Load_MWh",
    ]


    ids_y = []
    ids_i = []
    ids_n = []
    ids_z = []
    ids_p = []
    ids_k = []
    ids_d = []
    annual_gen_info = reference["annual_gen_info_OP"]
    try
        ids_y = [(y) for (y) in sort!(parse.(Int, (collect(get_index(reference["1"], :planning_stages, 0)))))]
        ids_i = [(k) for (k) in sort!(parse.(Int, (collect(get_index(reference["1"], :gen_index, 0)))))]
        ids_n = [(k) for (k) in sort!(parse.(Int, (collect(get_index(reference["1"], :bus, 0)))))]
        ids_z = [(i) for (i) in sort!(parse.(Int, (collect(keys(reference["1"][0][:zone]["reserve"])))))]
        ids_p = [(i) for (i) in sort!(parse.(Int, (collect(keys(reference["1"][0][:zone]["policy"])))))]
        ids_d = [(k) for (k) in sort!(parse.(Int, (collect(get_index(reference["1"], :repdays, 0)))))]
        ids_k = [(k) for (k) in sort!(parse.(Int, (collect(get_index(reference["1"], :branch, 0))))) if parameter(reference["1"], 0, :branch, "model_flag", k) == true]
    catch
        ids_y = [(y) for (y) in sort!(collect(get_index(reference["1"], :planning_stages, 0)))]
        ids_i = [(k) for (k) in sort!(collect(get_index(reference["1"], :gen_index, 0)))]
        ids_n = [(k) for (k) in sort!(collect(get_index(reference["1"], :bus, 0)))]
        ids_z = [(i) for (i) in sort!(parse.(Int, (collect(keys(reference["1"][0][:zone]["reserve"])))))]
        ids_p = [(i) for (i) in sort!(parse.(Int, (collect(keys(reference["1"][0][:zone]["policy"])))))]
        ids_d = [(k) for (k) in sort!(collect(get_index(reference["1"], :repdays, 0)))]
        ids_k = [(k) for (k) in sort!(collect(get_index(reference["1"], :branch, 0))) if parameter(reference["1"], 0, :branch, "model_flag", k) == true]
    end
    # Distributed path reads per-year scarcity totals from the reduced payload (the concatenated
    # __market_OP.csv does not exist yet); serial path sums them from the CSV as before.
    op_pu_applied = get(result, "pu_applied", false) == true
    if !op_pu_applied
        market_outcome_df = DataFrame(CSV.File(joinpath(output_path, string(parameter(reference["1"], 0, :case_name), "__market_OP.csv"))))
    end

    system_output_list = Array{Any}(undef, length(ids_y)+1, length(system_label_list))
    system_output_list[1,:] = system_label_list

    extensive_system_output_list = Array{Any}(undef, length(ids_y) * parameter(reference["1"], 0, :planning_stages, "stage_length", 1)+1, length(extensive_system_label_list))
    extensive_system_output_list[1,:] = extensive_system_label_list

    # parameters
    Scenario = parameter(reference["1"], 0, :case_name)
    current_year = setting["Planning Design"]["dollar_year_value"] # Discount factor relative to the current dollor year
    discount_rate = setting["Planning Design"]["discount_rate_value"]
    total_day_groups = setting["Simulation Configuration"]["NDAY_Groups_OP"]

    # total cost arrays
    annual_total_invc = zeros(Float64, length(ids_y))

    total_gen_invc_vector = zeros(Float64, length(ids_y) * parameter(reference["1"], 0, :planning_stages, "stage_length", 1))
    total_gen_ITC_vector = zeros(Float64, length(ids_y) * parameter(reference["1"], 0, :planning_stages, "stage_length", 1))
    total_trans_invc_vector = zeros(Float64, length(ids_y) * parameter(reference["1"], 0, :planning_stages, "stage_length", 1))
    total_trans_fom_vector = zeros(Float64, length(ids_y) * parameter(reference["1"], 0, :planning_stages, "stage_length", 1))
    total_gen_retc_vector = zeros(Float64, length(ids_y) * parameter(reference["1"], 0, :planning_stages, "stage_length", 1))
    total_gen_FOM_vector = zeros(Float64, length(ids_y) * parameter(reference["1"], 0, :planning_stages, "stage_length", 1))

    total_gen_cost_vector = zeros(Float64, length(ids_y) * parameter(reference["1"], 0, :planning_stages, "stage_length", 1))
    total_fuel_cost_vector= zeros(Float64, length(ids_y) * parameter(reference["1"], 0, :planning_stages, "stage_length", 1))
    total_VOM_cost_vector = zeros(Float64, length(ids_y) * parameter(reference["1"], 0, :planning_stages, "stage_length", 1))
    total_commitment_cost_vector = zeros(Float64, length(ids_y) * parameter(reference["1"], 0, :planning_stages, "stage_length", 1))

    total_regulation_cost_vector = zeros(Float64, length(ids_y) * parameter(reference["1"], 0, :planning_stages, "stage_length", 1))
    total_spin_cost_vector = zeros(Float64, length(ids_y) * parameter(reference["1"], 0, :planning_stages, "stage_length", 1))
    total_non_spin_cost_vector = zeros(Float64, length(ids_y) * parameter(reference["1"], 0, :planning_stages, "stage_length", 1))
    total_flex_cost_vector = zeros(Float64, length(ids_y) * parameter(reference["1"], 0, :planning_stages, "stage_length", 1))

    total_carbon_tax_vector = zeros(Float64, length(ids_y) * parameter(reference["1"], 0, :planning_stages, "stage_length", 1))
    total_emissions_vector = zeros(Float64, length(ids_y) * parameter(reference["1"], 0, :planning_stages, "stage_length", 1))
    total_ptc_vector = zeros(Float64, length(ids_y) * parameter(reference["1"], 0, :planning_stages, "stage_length", 1))

    total_scarcity_cost_ENS_vector = zeros(Float64, length(ids_y) * parameter(reference["1"], 0, :planning_stages, "stage_length", 1))
    total_scarcity_cost_SPIN_vector = zeros(Float64, length(ids_y) * parameter(reference["1"], 0, :planning_stages, "stage_length", 1))
    total_scarcity_cost_NSPIN_vector = zeros(Float64, length(ids_y) * parameter(reference["1"], 0, :planning_stages, "stage_length", 1))
    total_scarcity_cost_Flex_vector = zeros(Float64, length(ids_y) * parameter(reference["1"], 0, :planning_stages, "stage_length", 1))

    total_scarcity_ENS_vector = zeros(Float64, length(ids_y) * parameter(reference["1"], 0, :planning_stages, "stage_length", 1))
    total_scarcity_SPIN_vector = zeros(Float64, length(ids_y) * parameter(reference["1"], 0, :planning_stages, "stage_length", 1))
    total_scarcity_NSPIN_vector = zeros(Float64, length(ids_y) * parameter(reference["1"], 0, :planning_stages, "stage_length", 1))
    total_scarcity_FLEX_vector = zeros(Float64, length(ids_y) * parameter(reference["1"], 0, :planning_stages, "stage_length", 1))

    total_CEGT_slack_vector = zeros(Float64, length(ids_y) * parameter(reference["1"], 0, :planning_stages, "stage_length", 1))
    total_RPS_slack_vector = zeros(Float64, length(ids_y) * parameter(reference["1"], 0, :planning_stages, "stage_length", 1))

    total_generation_vector = zeros(Float64, length(ids_y) * parameter(reference["1"], 0, :planning_stages, "stage_length", 1))
    total_tnd_loss_vector = zeros(Float64, length(ids_y) * parameter(reference["1"], 0, :planning_stages, "stage_length", 1))
    total_storage_charge_mwh_vector = zeros(Float64, length(ids_y) * parameter(reference["1"], 0, :planning_stages, "stage_length", 1))
    total_reg_up_vector = zeros(Float64, length(ids_y) * parameter(reference["1"], 0, :planning_stages, "stage_length", 1))
    total_reg_dn_vector = zeros(Float64, length(ids_y) * parameter(reference["1"], 0, :planning_stages, "stage_length", 1))
    total_spin_vector = zeros(Float64, length(ids_y) * parameter(reference["1"], 0, :planning_stages, "stage_length", 1))
    total_nspin_vector = zeros(Float64, length(ids_y) * parameter(reference["1"], 0, :planning_stages, "stage_length", 1))
    total_flex_up_vector = zeros(Float64, length(ids_y) * parameter(reference["1"], 0, :planning_stages, "stage_length", 1))
    total_flex_dn_vector = zeros(Float64, length(ids_y) * parameter(reference["1"], 0, :planning_stages, "stage_length", 1))

    total_prm_vector = zeros(Float64, length(ids_y) * parameter(reference["1"], 0, :planning_stages, "stage_length", 1))
    total_peak_demand_vector = zeros(Float64, length(ids_y) * parameter(reference["1"], 0, :planning_stages, "stage_length", 1))
    total_ucap_vector = zeros(Float64, length(ids_y) * parameter(reference["1"], 0, :planning_stages, "stage_length", 1))
    total_icap_vector = zeros(Float64, length(ids_y) * parameter(reference["1"], 0, :planning_stages, "stage_length", 1))
    total_demand_vector = zeros(Float64, length(ids_y) * parameter(reference["1"], 0, :planning_stages, "stage_length", 1))
    total_dispatched_load_vector = zeros(Float64, length(ids_y) * parameter(reference["1"], 0, :planning_stages, "stage_length", 1))
        
    # start reporting for each year 
    for year_id in eachindex(ids_y)

        reference_of_this_year = reference[string(year_id)]

        # day-group-invariant build decisions (dg "1"); source differs serial vs distributed
        exp_dec_y = _op_expansion_decisions(result, year_id)

        year_of_year_id = parameter(reference_of_this_year, 0, :planning_stages, "year", year_id)
        stage_length = parameter(reference_of_this_year, 0, :planning_stages, "stage_length", year_id)
        
        num_stage = setting["Planning Design"]["num_stages_value"]
        remaining_years = (num_stage + 1 - year_id) * stage_length
        
        # Investment costs for generators + tax credits 
        total_gen_invc = 0.0    # 1) Upfront total investment costs 
        total_ITC = 0.0         # 2) Investment Tax Credits

        for i in ids_i

            bus_idx = parameter(reference_of_this_year, 0, :gen_index, "bus_idx", i)
            tech_idx = parameter(reference_of_this_year, 0, :gen_index, "genco_tech_id", i)
            capacity = parameter(reference_of_this_year, bus_idx, :gen_bus, tech_idx, "CAP") * pu_power_base
            CRP_i = parameter(reference_of_this_year, bus_idx, :gen_bus, tech_idx, "crpyears")

            # Determine the cost duration (investment payments)
            payment_duration = min(CRP_i, remaining_years)

            # Get the investment costs
            investment_cost = parameter(reference_of_this_year, bus_idx, :gen_bus, tech_idx, "Annual_INVC")[string(year_id)] * pu_econ_base
            storage_investment_cost = parameter(reference_of_this_year, bus_idx, :gen_bus, tech_idx, "Annual_STO_INV")[string(year_id)] * pu_econ_base

            # Get ITC 
            itc = 0.0
            if (setting["Simulation Configuration"]["ITC_Flag"] == true) && (parameter(reference_of_this_year, bus_idx, :gen_bus, tech_idx, "ITC Flag") == true)
                itc_year = min(2050, year_of_year_id)
                itc = parameter(reference_of_this_year, bus_idx, :gen_bus, tech_idx, "ITC")[string(itc_year)] # %
            end

            # Collect data
            idx_string = string("(", i, ", ", year_id, ")")
            u_new_G_iy = 0.0
            u_new_ESH_iy = 0.0
            
            if haskey(result_LC_GTEP_expansion, "solution")
                u_new_G_iy = result_LC_GTEP_expansion["solution"]["expansion"][idx_string]["u_new_G_iy"]  

                if parameter(reference_of_this_year, 0, :gen_index, "UNIT_CATEGORY", i) == "STORAGE"
                    u_new_ESH_iy = result_LC_GTEP_expansion["solution"]["expansion"][idx_string]["u_new_ESH_iy"]
                end
            end
            
            for future_y in 1:payment_duration

                future_year = year_of_year_id + future_y - 1    # Actual future year 
                discount_factor = (1 + discount_rate)^-(future_year - current_year) # Discount factor relative to the current dollor year

                gen_invc = 0.0
                
                # Expansion Cost (new units only) 
                gen_invc += 1000 * investment_cost * capacity * u_new_G_iy / (1.0 - itc)

                # Storage Duration Cost ($/kWh)
                if parameter(reference_of_this_year, 0, :gen_index, "UNIT_CATEGORY", i) == "STORAGE"
                    gen_invc += 1000 * storage_investment_cost * capacity * u_new_ESH_iy / (1.0 - itc)
                end

                # update the total investment cost
                total_gen_invc += discount_factor * gen_invc

                # Investment Tax Credits
                total_ITC += itc * discount_factor * gen_invc

                # update annual total_gen_invc_vector
                total_gen_invc_vector[(year_id-1)*stage_length + future_y] += discount_factor * gen_invc
                total_gen_ITC_vector[(year_id-1)*stage_length + future_y] += itc * discount_factor * gen_invc
            end
        end


        # Transmission investment costs
        total_upfront_invc_line = 0.0 # transmission expansion cost in each stage (thus multiplied by stage_length)
        total_upfront_fom_line = 0.0  # per-stage transmission FOM (existing baseline + expansion)
        
        trans_crpyears = setting["Planning Design"]["transmission_investment_CRP_value"]
        payment_duration = min(trans_crpyears, remaining_years)

        for k in ids_k
            
            rate_a = parameter(reference_of_this_year, 0, :branch, "rate_a", k) * pu_power_base
            length = get(reference_of_this_year[0][:branch][k], "length", get(reference_of_this_year[0][:branch][k], "Length", 0.0))
            transmission_expansion_cost = parameter(reference_of_this_year, 0, :branch, "transmission_expansion_cost", k) * pu_econ_base
            # DC-tie expansion cost is per-MW (no mile factor); AC lines are per-MW-mile
            len_factor = get(reference_of_this_year[0][:branch][k], "dc_line", false) == true ? 1.0 : length

            transmission_investment_cost_basis = transmission_expansion_cost * rate_a * len_factor

            u_new_T_ky = 0.0
            if haskey(result_LC_GTEP_expansion, "solution")
                u_new_T_ky = result_LC_GTEP_expansion["solution"]["expansion_line"][string("(", k, ", ", year_id, ")")]["u_new_T_ky"]
            end
            
            for future_y in 1:payment_duration

                future_year = year_of_year_id + future_y - 1    # Actual future year 
                discount_factor = (1 + discount_rate)^-(future_year - current_year) # Discount factor relative to the current dollor year

                # Transmission expansion cost ($/mile)
                if parameter(reference_of_this_year, 0, :branch, "expansion_flag", k) == true
                    invc_line = transmission_investment_cost_basis * u_new_T_ky
                    total_upfront_invc_line += discount_factor * invc_line
                    total_trans_invc_vector[(year_id-1)*stage_length + future_y] += discount_factor * invc_line
                end

            end

            # Transmission FOM on the full in-service grid: existing baseline (constant) + expansion (u_T_ky).
            # Charged every operating year of the stage, mirroring the rate_a*(1+u_T_ky) flow limit.
            transmission_fom_cost = parameter(reference_of_this_year, 0, :branch, "transmission_fom_cost", k) * pu_econ_base
            fom_basis = transmission_fom_cost * rate_a * len_factor
            u_T_ky_value = 0.0
            if (parameter(reference_of_this_year, 0, :branch, "expansion_flag", k) == true) && haskey(result_LC_GTEP_expansion, "solution")
                u_T_ky_value = result_LC_GTEP_expansion["solution"]["expansion"][string("(", k, ", ", year_id, ")")]["u_T_ky"]
            end
            for future_y in 1:stage_length
                future_year = year_of_year_id + future_y - 1
                discount_factor = (1 + discount_rate)^-(future_year - current_year)
                fom_annual = discount_factor * fom_basis * (1 + u_T_ky_value)
                total_upfront_fom_line += fom_annual
                total_trans_fom_vector[(year_id-1)*stage_length + future_y] += fom_annual
            end
        end
       

        # Retirement costs for generators 
        total_retirement_cost = 0.0
        discount_factor_year_id = (1 + discount_rate)^-(year_of_year_id - current_year)

        for i in ids_i
            bus_idx = parameter(reference_of_this_year, 0, :gen_index, "bus_idx", i)
            tech_idx = parameter(reference_of_this_year, 0, :gen_index, "genco_tech_id", i)

            # Decommissioning Cost: k$/MW
            DECC = parameter(reference_of_this_year, bus_idx, :gen_bus, tech_idx, "DECC") * pu_econ_base

            u_ret_G_iy = 0.0
            if haskey(result_LC_GTEP_expansion, "solution")
                u_ret_G_iy = result_LC_GTEP_expansion["solution"]["expansion"][string("(", i, ", ", year_id, ")")]["u_ret_G_iy"]  
            end

            if DECC > 0
                capacity = parameter(reference_of_this_year, bus_idx, :gen_bus, tech_idx, "CAP") * pu_power_base
                retirement_cost = 1000 * DECC * capacity * u_ret_G_iy    
                total_retirement_cost += discount_factor_year_id * retirement_cost
                total_gen_retc_vector[(year_id-1)*stage_length + 1] += discount_factor_year_id * retirement_cost
            end
        end

        # FOM cost
        total_FOM = 0.0
        for i in ids_i
            bus_idx = parameter(reference_of_this_year, 0, :gen_index, "bus_idx", i)
            tech_idx = parameter(reference_of_this_year, 0, :gen_index, "genco_tech_id", i)

            capacity = parameter(reference_of_this_year, bus_idx, :gen_bus, tech_idx, "CAP") * pu_power_base
            FOM = parameter(reference_of_this_year, bus_idx, :gen_bus, tech_idx, "Annual_FOM")[string(year_id)] * pu_econ_base

            # Collect data
            u_G_iy = exp_dec_y[string("(", i, ", ", year_id, ")")]["u_G_iy"]

            if (FOM > 0) && (capacity > 0)
        
                for future_y in 1:stage_length
                
                    future_year = year_of_year_id + future_y - 1 
                    discount_factor_year_id = (1 + discount_rate)^-(future_year - current_year)

                    # Fixed OM Cost ($/kW-Year)
                    total_FOM += discount_factor_year_id * 1000 * FOM * capacity * u_G_iy
                    total_gen_FOM_vector[(year_id-1)*stage_length + future_y] += discount_factor_year_id * 1000 * FOM * capacity * u_G_iy
                end
        
            end
        end

        # Generation, Ancillary Services Costs, carbon tax, and emissions
        ctax = setting["Simulation Configuration"]["CTAX"] * pu_econ_ref_base

        total_gen_cost = 0.0
        total_fuel_cost = 0.0
        total_VOM_cost = 0.0
        total_commitment_cost = 0.0

        total_regulation_cost = 0.0
        total_spin_cost = 0.0
        total_non_spin_cost = 0.0
        total_flex_cost = 0.0

        total_carbon_tax = 0.0
        total_emissions = 0.0

        total_generation = 0.0
        total_storage_charge_mwh = 0.0
        total_reserve_reg_up = 0.0
        total_reserve_reg_dn = 0.0
        total_reserve_spin = 0.0
        total_reserve_nspin = 0.0
        total_reserve_flex_up = 0.0
        total_reserve_flex_dn = 0.0
        
        for i in ids_i
    
            bus_idx = parameter(reference_of_this_year, 0, :gen_index, "bus_idx", i)
            tech_idx = parameter(reference_of_this_year, 0, :gen_index, "genco_tech_id", i)
            unit_group = parameter(reference_of_this_year, 0, :gen_index, "UNIT_GROUP", i)
            Annual_VOM = parameter(reference_of_this_year, bus_idx, :gen_bus, tech_idx, "Annual_VOM")[string(year_id)] * pu_econ_base 
            emission_rate = parameter(reference_of_this_year, bus_idx, :gen_bus, tech_idx, "Emission_CO2") / pu_power_base

            for future_y in 1:stage_length

                future_year = year_of_year_id + future_y - 1 
                discount_factor = (1 + discount_rate)^-(future_year - current_year)
                    
                # Generation
                generation = annual_gen_info[i][year_id]["Generation"] / stage_length # annual_gen_info is pre-multiplied by stage_length*num_days; divide by stage_length for per-sub-year values
                geneneration_cost = annual_gen_info[i][year_id]["Generation_Cost"] / stage_length # annual_gen_info is pre-multiplied by stage_length*num_days; divide by stage_length for per-sub-year values
                fuel_cost = annual_gen_info[i][year_id]["Fuel_Cost"] / stage_length # annual_gen_info is pre-multiplied by stage_length*num_days; divide by stage_length for per-sub-year values
                VOM_cost = generation * Annual_VOM

                total_gen_cost += discount_factor * geneneration_cost
                total_fuel_cost += discount_factor * fuel_cost
                total_VOM_cost += discount_factor * VOM_cost
                total_generation += generation

                total_gen_cost_vector[(year_id-1)*stage_length + future_y] += discount_factor * geneneration_cost
                total_fuel_cost_vector[(year_id-1)*stage_length + future_y] += discount_factor * fuel_cost
                total_VOM_cost_vector[(year_id-1)*stage_length + future_y] += discount_factor * VOM_cost
                total_generation_vector[(year_id-1)*stage_length + future_y] += generation

                # Storage charge
                storage_charge_mwh = annual_gen_info[i][year_id]["Storage_Charge"] / stage_length # annual_gen_info is pre-multiplied by stage_length*num_days; divide by stage_length for per-sub-year values
                total_storage_charge_mwh += storage_charge_mwh
                total_storage_charge_mwh_vector[(year_id-1)*stage_length + future_y] += storage_charge_mwh

                # Ancillary Services Costs
                reserve_RU = annual_gen_info[i][year_id]["Reserve_RU"] / stage_length
                reserve_RD = annual_gen_info[i][year_id]["Reserve_RD"] / stage_length
                reserve_Spin = annual_gen_info[i][year_id]["Reserve_Spin"] / stage_length
                reserve_NSpin = annual_gen_info[i][year_id]["Reserve_NSpin"] / stage_length
                reserve_FU = annual_gen_info[i][year_id]["Reserve_FU"] / stage_length
                reserve_FD = annual_gen_info[i][year_id]["Reserve_FD"] / stage_length

                total_reserve_reg_up += reserve_RU
                total_reserve_reg_dn += reserve_RD
                total_reserve_spin += reserve_Spin
                total_reserve_nspin += reserve_NSpin
                total_reserve_flex_up += reserve_FU
                total_reserve_flex_dn += reserve_FD

                total_reg_up_vector[(year_id-1)*stage_length + future_y] += reserve_RU
                total_reg_dn_vector[(year_id-1)*stage_length + future_y] += reserve_RD
                total_spin_vector[(year_id-1)*stage_length + future_y] += reserve_Spin
                total_nspin_vector[(year_id-1)*stage_length + future_y] += reserve_NSpin
                total_flex_up_vector[(year_id-1)*stage_length + future_y] += reserve_FU
                total_flex_dn_vector[(year_id-1)*stage_length + future_y] += reserve_FD

                reg_cost = annual_gen_info[i][year_id]["Reserve_Reg_Cost"] / stage_length
                spin_cost = annual_gen_info[i][year_id]["Reserve_Spin_Cost"] / stage_length
                nspin_cost = annual_gen_info[i][year_id]["Reserve_NSpin_Cost"] / stage_length
                flex_cost = annual_gen_info[i][year_id]["Reserve_Flex_Cost"] / stage_length

                total_regulation_cost += discount_factor * reg_cost
                total_spin_cost += discount_factor * spin_cost
                total_non_spin_cost += discount_factor * nspin_cost
                total_flex_cost += discount_factor * flex_cost

                total_regulation_cost_vector[(year_id-1)*stage_length + future_y] += discount_factor * reg_cost
                total_spin_cost_vector[(year_id-1)*stage_length + future_y] += discount_factor * spin_cost
                total_non_spin_cost_vector[(year_id-1)*stage_length + future_y] += discount_factor * nspin_cost
                total_flex_cost_vector[(year_id-1)*stage_length + future_y] += discount_factor * flex_cost
                
                # Commitment Cost 
                noload_cost = annual_gen_info[i][year_id]["Noload_Cost"] / stage_length
                startup_cost = annual_gen_info[i][year_id]["StartUp_Cost"] / stage_length
                total_commitment_cost += discount_factor * (noload_cost + startup_cost)

                total_commitment_cost_vector[(year_id-1)*stage_length + future_y] += discount_factor * (noload_cost + startup_cost)
                
                # Carbon Tax and Emissions  
                total_carbon_tax += discount_factor * generation * emission_rate * ctax
                total_emissions += generation * emission_rate

                total_carbon_tax_vector[(year_id-1)*stage_length + future_y] += discount_factor * generation * emission_rate * ctax
                total_emissions_vector[(year_id-1)*stage_length + future_y] += generation * emission_rate

            end

        end


        # Generation Credit
        total_ptc = 0.0
        PTC_final_year_for_existing_assets = 5
        PTC_final_year_for_new_assets = 10

        if setting["Simulation Configuration"]["PTC_Flag"] == true
            for i in ids_i
                bus_idx = parameter(reference_of_this_year, 0, :gen_index, "bus_idx", i)
                tech_idx = parameter(reference_of_this_year, 0, :gen_index, "genco_tech_id", i)
                Tech_Type = parameter(reference_of_this_year, bus_idx, :gen_bus, tech_idx, "Tech_Type")
                unit_group = parameter(reference_of_this_year, 0, :gen_index, "UNIT_GROUP", i)
                profile_type = parameter(reference_of_this_year, 0, :gen_index, "Profile_Type", i)
                capacity = parameter(reference_of_this_year, bus_idx, :gen_bus, tech_idx, "CAP") * pu_power_base
                PMAX = parameter(reference_of_this_year, bus_idx, :gen_bus, tech_idx, "PMAX")

                if parameter(reference_of_this_year, bus_idx, :gen_bus, tech_idx, "PTC Flag") == true

                    if (Tech_Type == "Existing") && (unit_group != "nuclear")
                        
                        count = (year_id-1)*stage_length + 1
                    
                        for future_y in 1:stage_length

                            future_year = year_of_year_id + future_y - 1 
                            discount_factor = (1 + discount_rate)^-(future_year - current_year)
        
                            if count > PTC_final_year_for_existing_assets
                                break
                            end
        
                            ptc_year = min(2050, future_year)
                            ptc = parameter(reference_of_this_year, bus_idx, :gen_bus, tech_idx, "PTC")[string(ptc_year)] * 10 * pu_econ_base

                            # get generation
                            generation = annual_gen_info[i][year_id]["Generation"] / stage_length # annual_gen_info is pre-multiplied by stage_length*num_days; divide by stage_length for per-sub-year values

                            total_ptc += discount_factor * generation * ptc
                            total_ptc_vector[(year_id-1)*stage_length + future_y] += discount_factor * generation * ptc
        
                            count += 1
                        end

                    elseif (Tech_Type == "Existing") && (unit_group == "nuclear") 
                        
                        for future_y in 1:stage_length

                            future_year = year_of_year_id + future_y - 1 
                            discount_factor = (1 + discount_rate)^-(future_year - current_year)
        
                            ptc_year = min(2050, future_year)
                            ptc = parameter(reference_of_this_year, bus_idx, :gen_bus, tech_idx, "PTC")[string(ptc_year)] * 10 * pu_econ_base

                            # get generation
                            generation = annual_gen_info[i][year_id]["Generation"] / stage_length # annual_gen_info is pre-multiplied by stage_length*num_days; divide by stage_length for per-sub-year values

                            total_ptc += discount_factor * generation * ptc
                            total_ptc_vector[(year_id-1)*stage_length + future_y] += discount_factor * generation * ptc
                        end

                    elseif Tech_Type == "New" 
                        
                        ptc_year = min(2050, year_of_year_id)   # ptc is fixed to the investment year ptc
                        ptc = parameter(reference_of_this_year, bus_idx, :gen_bus, tech_idx, "PTC")[string(ptc_year)] * 10 * pu_econ_base 

                        ptc_duration = min(PTC_final_year_for_new_assets, remaining_years)

                        # Collect data
                        idx_string = string("(", i, ", ", year_id, ")")
                        u_new_G_iy = 0.0
                        if haskey(result_LC_GTEP_expansion, "solution")
                            u_new_G_iy = result_LC_GTEP_expansion["solution"]["expansion"][idx_string]["u_new_G_iy"]
                        end
        
                        if u_new_G_iy > 0   # only if the unit is built
                            for future_y in 1:ptc_duration # for the (next 10 years or remaining years)

                                future_year = year_of_year_id + future_y - 1 
                                discount_factor = (1 + discount_rate)^-(future_year - current_year)
    
                                # Estimated annual generation in year y 
                                vre_data = reference_of_this_year[0][:bus][bus_idx]["vre_aggregated_data"][string(year_id)]
    
                                # Estimated generation per MW from the shape class (Profile_Type)
                                annual_generation_per_MW = 0.0
                                if profile_type != "NA"
                                    annual_generation_per_MW = vre_data[profile_type * "_shape"]
                                else    # nuclear
                                    annual_generation_per_MW = 8760 * PMAX # derated by FOR
                                end
    
                                total_ptc += discount_factor * annual_generation_per_MW * ptc * capacity * PMAX * u_new_G_iy
                                total_ptc_vector[(year_id-1)*stage_length + future_y] += discount_factor * annual_generation_per_MW * ptc * capacity * PMAX * u_new_G_iy
    
                            end
                        end
                        
                    
                    end
                end
            end
        end

        # Scarcity 

        total_scarcity_cost_ENS = 0.0
        total_scarcity_cost_SPIN = 0.0
        total_scarcity_cost_NSPIN = 0.0
        total_scarcity_cost_Flex = 0.0

        total_scarcity_ENS = 0.0
        total_scarcity_SPIN = 0.0
        total_scarcity_NSPIN = 0.0
        total_scarcity_FLEX = 0.0
        
        for future_y in 1:stage_length

            future_year = year_of_year_id + future_y - 1 
            discount_factor = (1 + discount_rate)^-(future_year - current_year)

            if op_pu_applied
                st_y = get(result["scarcity_totals"], year_id, Dict{String,Float64}())
                scarcity_E     = get(st_y, "scarcity_E", 0.0)
                scarcity_SPIN  = get(st_y, "scarcity_SPIN", 0.0)
                scarcity_NSPIN = get(st_y, "scarcity_NSPIN", 0.0)
                scarcity_FU    = get(st_y, "scarcity_FU", 0.0)
                scarcity_FD    = get(st_y, "scarcity_FD", 0.0)
            else
                scarcity_E = sum(market_outcome_df[(market_outcome_df.Stage .== float(year_id)), :Unserved_Energy_MW] .* market_outcome_df[(market_outcome_df.Stage .== float(year_id)), :Days_Represented])
                scarcity_SPIN = sum(market_outcome_df[(market_outcome_df.Stage .== float(year_id)), :Spin_Shortfall_MW] .* market_outcome_df[(market_outcome_df.Stage .== float(year_id)), :Days_Represented])
                scarcity_NSPIN = sum(market_outcome_df[(market_outcome_df.Stage .== float(year_id)), :NSpin_Shortfall_MW] .* market_outcome_df[(market_outcome_df.Stage .== float(year_id)), :Days_Represented])
                scarcity_FU = sum(market_outcome_df[(market_outcome_df.Stage .== float(year_id)), :FlexUp_Shortfall_MW] .* market_outcome_df[(market_outcome_df.Stage .== float(year_id)), :Days_Represented])
                scarcity_FD = sum(market_outcome_df[(market_outcome_df.Stage .== float(year_id)), :FlexDown_Shortfall_MW] .* market_outcome_df[(market_outcome_df.Stage .== float(year_id)), :Days_Represented])
            end

            total_scarcity_cost_ENS += discount_factor * scarcity_E * setting["Simulation Configuration"]["VOLL"] * pu_econ_base
            total_scarcity_cost_SPIN += discount_factor * scarcity_SPIN * setting["Simulation Configuration"]["SRSP"] * pu_econ_base
            total_scarcity_cost_NSPIN += discount_factor * scarcity_NSPIN * setting["Simulation Configuration"]["NSRSP"] * pu_econ_base
            total_scarcity_cost_Flex += discount_factor * (scarcity_FU + scarcity_FD) * setting["Simulation Configuration"]["FLEXRSP"] * pu_econ_base

            total_scarcity_cost_ENS_vector[(year_id-1)*stage_length + future_y] = discount_factor * scarcity_E * setting["Simulation Configuration"]["VOLL"] * pu_econ_base
            total_scarcity_cost_SPIN_vector[(year_id-1)*stage_length + future_y] = discount_factor * scarcity_SPIN * setting["Simulation Configuration"]["SRSP"] * pu_econ_base
            total_scarcity_cost_NSPIN_vector[(year_id-1)*stage_length + future_y] = discount_factor * scarcity_NSPIN * setting["Simulation Configuration"]["NSRSP"] * pu_econ_base
            total_scarcity_cost_Flex_vector[(year_id-1)*stage_length + future_y] = discount_factor * (scarcity_FU + scarcity_FD) * setting["Simulation Configuration"]["FLEXRSP"] * pu_econ_base

            total_scarcity_ENS += scarcity_E
            total_scarcity_SPIN += scarcity_SPIN
            total_scarcity_NSPIN += scarcity_NSPIN
            total_scarcity_FLEX += (scarcity_FU + scarcity_FD)

            total_scarcity_ENS_vector[(year_id-1)*stage_length + future_y] = scarcity_E
            total_scarcity_SPIN_vector[(year_id-1)*stage_length + future_y] = scarcity_SPIN
            total_scarcity_NSPIN_vector[(year_id-1)*stage_length + future_y] = scarcity_NSPIN
            total_scarcity_FLEX_vector[(year_id-1)*stage_length + future_y] = (scarcity_FU + scarcity_FD)
            
        end

        # Planning Reserve Margin
        total_ucap = 0.0
        total_icap = 0.0
        for i in ids_i

            bus_idx = parameter(reference_of_this_year, 0, :gen_index, "bus_idx", i)
            tech_idx = parameter(reference_of_this_year, 0, :gen_index, "genco_tech_id", i)
            unit_group = parameter(reference_of_this_year, 0, :gen_index, "UNIT_GROUP", i)
            CAPCRED = parameter(reference_of_this_year, bus_idx, :gen_bus, tech_idx, "CAPCRED")
            if CAPCRED isa String
                data_identifier = reference_of_this_year[0][:zone]["capacity_credit"]["data_identifier"]
                data_region = reference_of_this_year[0][:bus][bus_idx]["region_mapping_info"][data_identifier]
                data_region = data_region isa AbstractString ? data_region : string(data_region)
                CAPCRED = reference_of_this_year[0][:zone]["capacity_credit"][data_region][CAPCRED]
            end

            capacity = parameter(reference_of_this_year, bus_idx, :gen_bus, tech_idx, "CAP") * pu_power_base

            # Collect data
            u_G_iy = exp_dec_y[string("(", i, ", ", year_id, ")")]["u_G_iy"]

            total_ucap += capacity * u_G_iy * CAPCRED
            total_icap += capacity * u_G_iy
        end
        
        demand_mwh_year = get_annual_total_demand_with_growth(reference_of_this_year[0], year_of_year_id) * pu_power_base   # annual demand (MWh) of the stage's dispatched year
        # load actually dispatched: representative-day weighted (differs from the annual input demand)
        dispatched_load_year = if op_pu_applied
            get(get(result["scarcity_totals"], year_id, Dict{String,Float64}()), "Load_MWh", 0.0)
        else
            yr_rows = market_outcome_df.Stage .== float(year_id)
            sum(market_outcome_df[yr_rows, :Load_MW] .* market_outcome_df[yr_rows, :Days_Represented])
        end
        # PRM uses the coincident peak (per-unit); scale to MW to match total_ucap (MW).
        Peak_Demand = get_annual_coincident_peak_demand_with_growth(reference_of_this_year[0], year_of_year_id) * pu_power_base
        PRM = (total_ucap / (Peak_Demand)) - 1

        enforce_tnd_loss = setting["Planning Design"]["enforce_transmission_loss_flag"] == true
        tnd_loss_pct = setting["Planning Design"]["transmission_loss_percent_value"] * 0.01
        total_tnd_loss = 0.0

        for future_y in 1:stage_length

            future_year = year_of_year_id + future_y - 1
            peak_demand = get_annual_coincident_peak_demand_with_growth(reference_of_this_year[0], future_year) * pu_power_base

            total_prm_vector[(year_id-1)*stage_length + future_y] = (total_ucap / (peak_demand)) - 1
            total_peak_demand_vector[(year_id-1)*stage_length + future_y] = peak_demand
            total_ucap_vector[(year_id-1)*stage_length + future_y] = total_ucap
            total_icap_vector[(year_id-1)*stage_length + future_y] = total_icap
            total_demand_vector[(year_id-1)*stage_length + future_y] = demand_mwh_year
            total_dispatched_load_vector[(year_id-1)*stage_length + future_y] = dispatched_load_year

            # Annual nominal demand (MWh) at flat system-wide loss rate; zero unless the loss is enforced.
            if enforce_tnd_loss
                # loss on the stage-dispatched demand so the annual report energy-balances
                annual_demand = get_annual_total_demand_with_growth(reference_of_this_year[0], year_of_year_id) * pu_power_base
                tnd_loss = annual_demand * tnd_loss_pct
                total_tnd_loss_vector[(year_id-1)*stage_length + future_y] = tnd_loss
                total_tnd_loss += tnd_loss
            end
        end

        # Policy Slacks: CEGT 
        total_CEGT_slack = 0.0
        clean_energy_generation_penalty = setting["Simulation Configuration"]["Clean_Energy_Generation_Penalty"] * pu_econ_base
        
        if setting["Simulation Configuration"]["Clean_Energy_Generation_Target_OP_Flag"] == true
            for future_y in 1:stage_length
                    
                future_year = year_of_year_id + future_y - 1 
                discount_factor_year_id = (1 + discount_rate)^-(future_year - current_year)

                CEGT_slack = 0.0

                slack_dec_y = _op_slack_decisions(result, year_id)
                for n in ids_p
                    for d in ids_d
                        idx_string = string("(", n, ", ", d, ", ", year_id, ")")
                        CEGT_slack += discount_factor_year_id * clean_energy_generation_penalty * parameter(reference_of_this_year, 0, :repdays, "NumDays", d) * slack_dec_y[idx_string]["slack_CEG_ndy"] * pu_power_base
                    end
                end
                
                total_CEGT_slack += CEGT_slack
                total_CEGT_slack_vector[(year_id-1)*stage_length + future_y] += CEGT_slack
            end
        end

        # Policy Slacks: RPS 
        total_RPS_slack = 0.0
        
        # Objective Value
        objective_value = 0.0
        for future_y in 1:stage_length 

            # apply discounting to the objective value because the operational model objective function does not apply discounting
            future_year = year_of_year_id + future_y - 1 
            discount_factor_year_id = (1 + discount_rate)^-(future_year - current_year)

            objectives_year = op_pu_applied ? result["objectives"][string(year_id)] : Dict(dg => result["operational model result"][string(year_id)][dg]["objective"] for dg in keys(result["operational model result"][string(year_id)]))
            for day_group_id in keys(objectives_year)
                objective_value += discount_factor_year_id * objectives_year[day_group_id] * setting["Simulation Setting"]["per_unit_econ_base_value"]
            end
        end

        
        # calculate PTC of year_id 
        PTC_of_year_y = 0.0   # we need to calculate the PTC for the year_id because the total_PTC value includes PTCs from all the future years (10 years)
        for future_y in 1:stage_length
            PTC_of_year_y += total_ptc_vector[(year_id-1)*stage_length + future_y]
        end

        total_system_cost = total_gen_invc + total_upfront_invc_line + total_upfront_fom_line + total_retirement_cost + total_FOM
        total_system_cost += total_gen_cost + total_commitment_cost + total_regulation_cost + total_spin_cost + total_non_spin_cost + total_flex_cost 
        total_system_cost += total_scarcity_cost_ENS + total_scarcity_cost_SPIN + total_scarcity_cost_NSPIN + total_scarcity_cost_Flex 
        total_system_cost += total_carbon_tax + total_CEGT_slack + total_RPS_slack
        total_system_cost -= total_ITC
        total_system_cost -= PTC_of_year_y

        total_system_cost_OP = total_system_cost - (total_gen_invc + total_upfront_invc_line + total_upfront_fom_line + total_retirement_cost + total_FOM + total_RPS_slack) + total_ITC
        
        # update system output list
        system_output_list[year_id + 1, :] = [
            Scenario, 
            year_id, 
            year_of_year_id, 
            stage_length, 
            total_gen_invc, 
            -total_ITC, 
            total_upfront_invc_line,
            total_upfront_fom_line,
            total_retirement_cost, 
            total_FOM,
            -PTC_of_year_y, 
            total_fuel_cost, 
            total_VOM_cost, 
            total_commitment_cost, 
            total_regulation_cost, 
            total_spin_cost, 
            total_non_spin_cost, 
            total_flex_cost, 
            total_scarcity_cost_ENS, 
            total_scarcity_cost_SPIN, 
            total_scarcity_cost_NSPIN, 
            total_scarcity_cost_Flex, 
            total_carbon_tax, 
            total_CEGT_slack, 
            total_RPS_slack, 
            total_system_cost,
            total_system_cost_OP,
            objective_value,
            total_generation,
            total_tnd_loss,
            total_storage_charge_mwh,
            total_reserve_reg_up,
            total_reserve_reg_dn,
            total_reserve_spin,
            total_reserve_nspin,
            total_reserve_flex_up,
            total_reserve_flex_dn,
            total_scarcity_ENS,
            total_scarcity_SPIN,
            total_scarcity_NSPIN,
            total_scarcity_FLEX,
            total_emissions,
            PRM,
            Peak_Demand,
            total_ucap,
            total_icap,
            demand_mwh_year * stage_length,
            dispatched_load_year * stage_length
        ]

        # update extensive system output list
        for future_y in 1:stage_length

            year_pointer = (year_id-1)*stage_length + future_y

            total_system_cost = total_gen_invc_vector[year_pointer] +
                    total_trans_invc_vector[year_pointer] +
                    total_trans_fom_vector[year_pointer] +
                    total_gen_retc_vector[year_pointer] +
                    total_gen_FOM_vector[year_pointer] +
                    total_gen_cost_vector[year_pointer] +
                    total_commitment_cost_vector[year_pointer] +
                    total_regulation_cost_vector[year_pointer] +
                    total_spin_cost_vector[year_pointer] +
                    total_non_spin_cost_vector[year_pointer] +
                    total_flex_cost_vector[year_pointer] +
                    total_scarcity_cost_ENS_vector[year_pointer] +
                    total_scarcity_cost_SPIN_vector[year_pointer] +
                    total_scarcity_cost_NSPIN_vector[year_pointer] +
                    total_scarcity_cost_Flex_vector[year_pointer] +
                    total_carbon_tax_vector[year_pointer] +
                    total_CEGT_slack_vector[year_pointer] +
                    total_RPS_slack_vector[year_pointer] -
                    total_gen_ITC_vector[year_pointer] -
                    total_ptc_vector[year_pointer]

            total_system_cost_OP = 
                    total_gen_cost_vector[year_pointer] +
                    total_commitment_cost_vector[year_pointer] +
                    total_regulation_cost_vector[year_pointer] +
                    total_spin_cost_vector[year_pointer] +
                    total_non_spin_cost_vector[year_pointer] +
                    total_flex_cost_vector[year_pointer] +
                    total_scarcity_cost_ENS_vector[year_pointer] +
                    total_scarcity_cost_SPIN_vector[year_pointer] +
                    total_scarcity_cost_NSPIN_vector[year_pointer] +
                    total_scarcity_cost_Flex_vector[year_pointer] +
                    total_carbon_tax_vector[year_pointer] +
                    total_CEGT_slack_vector[year_pointer] -                    
                    total_ptc_vector[year_pointer]
            
            extensive_objective_value = 0.0
            if (objective_value > 0) && (future_y == 1)
                extensive_objective_value = objective_value
            end

            extensive_system_output_list[year_pointer + 1, :] = [
                Scenario, 
                year_id, 
                year_of_year_id + future_y - 1,
                1, # number of years is always 1 for extensive output
                total_gen_invc_vector[year_pointer], 
                -total_gen_ITC_vector[year_pointer], 
                total_trans_invc_vector[year_pointer],
                total_trans_fom_vector[year_pointer],
                total_gen_retc_vector[year_pointer], 
                total_gen_FOM_vector[year_pointer],
                -total_ptc_vector[year_pointer], 
                total_fuel_cost_vector[year_pointer], 
                total_VOM_cost_vector[year_pointer], 
                total_commitment_cost_vector[year_pointer], 
                total_regulation_cost_vector[year_pointer], 
                total_spin_cost_vector[year_pointer], 
                total_non_spin_cost_vector[year_pointer], 
                total_flex_cost_vector[year_pointer],
                total_scarcity_cost_ENS_vector[year_pointer], 
                total_scarcity_cost_SPIN_vector[year_pointer], 
                total_scarcity_cost_NSPIN_vector[year_pointer], 
                total_scarcity_cost_Flex_vector[year_pointer],
                total_carbon_tax_vector[year_pointer], 
                total_CEGT_slack_vector[year_pointer], 
                total_RPS_slack_vector[year_pointer], 
                total_system_cost,
                total_system_cost_OP,
                extensive_objective_value,
                total_generation_vector[year_pointer],
                total_tnd_loss_vector[year_pointer],
                total_storage_charge_mwh_vector[year_pointer],
                total_reg_up_vector[year_pointer],
                total_reg_dn_vector[year_pointer],
                total_spin_vector[year_pointer],
                total_nspin_vector[year_pointer],
                total_flex_up_vector[year_pointer],
                total_flex_dn_vector[year_pointer],
                total_scarcity_ENS_vector[year_pointer],
                total_scarcity_SPIN_vector[year_pointer],
                total_scarcity_NSPIN_vector[year_pointer],
                total_scarcity_FLEX_vector[year_pointer],
                total_emissions_vector[year_pointer],
                total_prm_vector[year_pointer],
                total_peak_demand_vector[year_pointer],
                total_ucap_vector[year_pointer],
                total_icap_vector[year_pointer],
                total_demand_vector[year_pointer],
                total_dispatched_load_vector[year_pointer]
            ]
        end

    end

    # Write the system output list to a csv file
    system_output_list, extensive_system_output_list = finalize_system_summaries(system_output_list, extensive_system_output_list, current_year)
    writedlm(system_file_name, system_output_list, ",")
    writedlm(extensive_system_file_name, extensive_system_output_list, ",")
    write_constant_dollar_annual_summary(extensive_system_file_name, extensive_system_output_list, discount_rate, current_year)

    @aleaf_info "[ALEAF LC_GTEP]: - System summary reporting,"

    
end


function report_result_power_flow_expansion_GTEP(result::Dict{String, Any}, reference, setting; years=nothing, out_dir=nothing)

    # -----------------------------
    # Scalars / flags
    # -----------------------------
    pu_power_base = setting["Simulation Setting"]["per_unit_base_value"]
    ph_flag      = false
    multi_round  = (setting["Planning Design"]["multi_round_solution_process_flag"] == true)

    Scenario = parameter(reference, 0, :case_name)

    # -----------------------------
    # Paths
    # -----------------------------
    output_path = parameter(reference, 0, :output_path)
    dir = out_dir === nothing ? output_path : out_dir
    out_dir === nothing || mkpath(dir)
    file_name   = joinpath(dir, string(Scenario, "__power_flow_EXP.csv"))

    # -----------------------------
    # Labels
    # -----------------------------
    label_list = [
        "Scenario","year","day","hour","time","line_id","f_bus","t_bus","f_region","t_region",
        "orinal_rate","expansion","final_rate","flow","LMP_from_bus","LMP_to_bus","congestion_flag",
        "wheeling_cost","line_length","NumDays",
        # Enhanced-hybrid audit (0 for non-hybrid corridors): increment flow, true corridor flow,
        # and the KVL residual |f_exp − u·f_base| (per-unit) — 0 at u=0 and u=ū, positive only mid-build.
        "f_exp","flow_total","kvl_residual",
    ]

    function get_bus_label(reference_data, bus_ref)
        try
            return parameter(reference_data, 0, :bus, "bus_i", bus_ref)
        catch
            for (_, bus_data) in reference_data[0][:bus]
                if get(bus_data, "bus_i", nothing) == bus_ref
                    return bus_data["bus_i"]
                end
            end
            return bus_ref
        end
    end

    # -----------------------------
    # Indices
    # -----------------------------
    ids_h = collect(setting["run_H"])
    ids_t = collect(setting["run_T"])

    ids_y = try
        sort!(parse.(Int, collect(get_index(reference, :planning_stages, 0))))
    catch
        sort!(collect(get_index(reference, :planning_stages, 0)))
    end
    local_ids_y = years === nothing ? ids_y : years

    ids_d = try
        sort!(parse.(Int, collect(get_index(reference, :repdays, 0))))
    catch
        sort!(collect(get_index(reference, :repdays, 0)))
    end

    ids_k_all = try
        sort!(parse.(Int, collect(get_index(reference, :branch, 0))))
    catch
        sort!(collect(get_index(reference, :branch, 0)))
    end

    # keep only modeled lines
    ids_k = [k for k in ids_k_all if parameter(reference, 0, :branch, "model_flag", k) == true]

    ids_decomp_group = ph_flag ? collect(get_index(reference, :repday_groups, 0)) : [1]

    # -----------------------------
    # Cache branch + bus metadata (k-only / bus-only)
    # -----------------------------
    bus_name = Dict{Int,Any}()
    bus_index_lookup = Dict{Any,Int}()
    for b in try
        sort!(parse.(Int, collect(get_index(reference, :bus, 0))))
    catch
        sort!(collect(get_index(reference, :bus, 0)))
    end
        bus_name[b] = get_bus_label(reference, b)
        bus_index_lookup[get(reference[0][:bus][b], "bus_i", b)] = b
    end

    f_bus_of_k   = Dict{Int,Int}()
    t_bus_of_k   = Dict{Int,Int}()
    rate0_of_k   = Dict{Int,Float64}()
    length_of_k  = Dict{Int,Float64}()
    f_reg_of_k   = Dict{Int,Any}()
    t_reg_of_k   = Dict{Int,Any}()

    for k in ids_k
        fb = parameter(reference, 0, :branch, "f_bus", k)
        tb = parameter(reference, 0, :branch, "t_bus", k)
        f_bus_of_k[k]  = fb
        t_bus_of_k[k]  = tb
        rate0_of_k[k]  = parameter(reference, 0, :branch, "rate_a", k) * pu_power_base
        length_of_k[k] = get(reference[0][:branch][k], "length", get(reference[0][:branch][k], "Length", 0.0))

        # from/to endpoints keyed by bus_i label (falls back to the generated bus label)
        f_reg_of_k[k] = get(bus_name, fb, get_bus_label(reference, fb))
        t_reg_of_k[k] = get(bus_name, tb, get_bus_label(reference, tb))
    end

    # -----------------------------
    # Open file + write header once
    # -----------------------------
    io = open(file_name, "w")
    println(io, join(report_labels(label_list), ","))

    file_lock = ReentrantLock()

    # -----------------------------
    # Main loops (parallelize over day)
    # -----------------------------
    for year in local_ids_y

        # cumulative u_T_ky matches the flow limit rate_a*(1+u_T_ky); round-local u_new_T_ky would understate it in myopic multi-round runs
        exp_u_T_ky = Dict{Int,Float64}()
        for k in ids_k
            idx_ky = "($k, $year)"
            exp_u_T_ky[k] = result["1"]["solution"]["expansion"][idx_ky]["u_T_ky"]
        end

        for decomp_group in ids_decomp_group
            grpkey = string(decomp_group)
            sol = result[grpkey]["solution"]

            # repdays for group
            ids_d_new = if ph_flag
                collect(parameter(reference, 0, :repday_groups, "Day_Idx_List", decomp_group))
            else
                ids_d
            end

            @threads for dd in eachindex(ids_d_new)
                day_id = ids_d_new[dd]

                # thread-local buffer
                buf = IOBuffer()

                NumDays = if multi_round
                    parameter(reference, 0, :multi_round_info, "planning_stages", year)[string(day_id)]["NumDays"]
                else
                    parameter(reference, 0, :repdays, "NumDays", day_id)
                end

                for hour_id in ids_h, time_id in ids_t
                    for line_id in ids_k

                        fb = f_bus_of_k[line_id]
                        tb = t_bus_of_k[line_id]

                        idx_kdhty = "($line_id, $day_id, $hour_id, $time_id, $year)"

                        orinal_rate = rate0_of_k[line_id]
                        expansion   = exp_u_T_ky[line_id]
                        final_rate  = (1 + expansion) * orinal_rate
                        line_length = length_of_k[line_id]

                        pf_sol = sol["powerflow"][idx_kdhty]
                        flow = pf_sol["f_kdhty"]
                        # Hybrid split flow: f_exp exists only for eligible AC corridors. For
                        # non-hybrid corridors (DC ties, flag off) the audit columns are trivial:
                        # f_exp=0, flow_total=flow, and the KVL residual is undefined -> 0.
                        if haskey(pf_sol, "f_exp_kdhty")
                            f_exp_val = pf_sol["f_exp_kdhty"]
                            flow_total = flow + f_exp_val
                            kvl_residual = f_exp_val - expansion * flow
                        else
                            f_exp_val = 0.0
                            flow_total = flow
                            kvl_residual = 0.0
                        end

                        fb_idx = get(bus_index_lookup, fb, fb)
                        tb_idx = get(bus_index_lookup, tb, tb)
                        idx_from = "($fb_idx, $day_id, $hour_id, $time_id, $year)"
                        idx_to   = "($tb_idx, $day_id, $hour_id, $time_id, $year)"

                        LMP_from_bus = sol["dual"][idx_from]["LMP"]
                        LMP_to_bus   = sol["dual"][idx_to]["LMP"]

                        congestion_flag = abs(flow) >= final_rate

                        wheeling_cost = 0.0  # keep as before

                        row = Any[
                            Scenario, year, stage_calendar_year(setting, year), day_id, hour_id, time_id, line_id,
                            fb, tb, f_reg_of_k[line_id], t_reg_of_k[line_id],
                            orinal_rate, expansion, final_rate, flow,
                            LMP_from_bus, LMP_to_bus, congestion_flag,
                            wheeling_cost, line_length, NumDays,
                            f_exp_val, flow_total, kvl_residual
                        ]
                        println(buf, join(row, ","))
                    end
                end

                # flush buffer (lock)
                lock(file_lock)
                try
                    write(io, take!(buf))
                finally
                    unlock(file_lock)
                end
            end
        end
    end

    close(io)

    @aleaf_info "[ALEAF LC_GTEP]: - Power flow solution reporting,"
    return nothing
end


function report_result_power_flow_operation_GTEP(result::Dict{String, Any}, result_LC_GTEP_expansion; generate_file::Bool=true, day_group_subset::Union{Nothing,Vector{Int}}=nothing, file_suffix::String="")

    reference = result["operational model system reference"]
    setting = result["setting"]

    pu_power_base = setting["Simulation Setting"]["per_unit_base_value"]

    output_path = parameter(reference["1"], 0, :output_path)
    file_name = string(parameter(reference["1"], 0, :case_name), "__power_flow_OP", file_suffix, ".csv")
    file_name = joinpath(output_path, file_name)

    # Write Label 
    label_list = [
        "Scenario",
        "year",
        "day",
        "hour",
        "time",
        "line_id",
        "f_bus",
        "t_bus",
        "f_region",
        "t_region",
        "orinal_rate",
        "expansion",
        "final_rate",
        "flow",
        "LMP_from_bus",
        "LMP_to_bus",
        "congestion_flag",
        "wheeling_cost",
        "line_length",
        "NumDays",
    ]

    # Write Outputs
    ids_y = []
    ids_k = []
    ids_h = [(h) for (h) in setting["run_H"]]
    ids_t = [(h) for (h) in setting["run_T"]]
    try
        ids_y = [(y) for (y) in sort!(parse.(Int, (collect(get_index(reference["1"], :planning_stages, 0)))))]
        ids_k = [(k) for (k) in sort!(parse.(Int, (collect(get_index(reference["1"], :branch, 0))))) if parameter(reference["1"], 0, :branch, "model_flag", k) == true]
    catch
        ids_y = [(y) for (y) in sort!(collect(get_index(reference["1"], :planning_stages, 0)))]
        ids_k = [(k) for (k) in sort!(collect(get_index(reference["1"], :branch, 0))) if parameter(reference["1"], 0, :branch, "model_flag", k) == true]
    end

    # Stream rows straight to disk (buffered per day-group) instead of preallocating one giant
    # Array{Any} covering every (line, day, hour, time) row. At 365 day-groups that array is tens of
    # GB of boxed cells and OOMs the master. This mirrors the dispatch reporter's streaming pattern.
    @inline function csv_escape(x)
        x === nothing && return ""
        x isa Number && return string(x)
        s = String(x)
        return occursin(r"[,\n\"]", s) ? "\"" * replace(s, "\"" => "\"\"") * "\"" : s
    end

    io = open(file_name, "w")
    println(io, join(report_labels(label_list), ","))
    io_lock = ReentrantLock()

    # grouping days
    Scenario = parameter(reference["1"], 0, :case_name)
    total_day_groups = setting["Simulation Configuration"]["NDAY_Groups_OP"]
    day_group_ids = day_group_subset === nothing ? collect(1:total_day_groups) : day_group_subset

    op_result_pf = result["operational model result"]
    for year in ids_y
        # workers hold a single-year view; skip stages this view does not solve
        haskey(op_result_pf, string(year)) || continue

        reference_of_this_year = reference[string(year)]
        bus_label_lookup = Dict{Any, Any}()
        bus_index_lookup = Dict{Any, Any}()
        for (bus_idx, bus_data) in reference_of_this_year[0][:bus]
            bus_label = get(bus_data, "bus_i", bus_idx)
            bus_label_lookup[bus_idx] = bus_label
            bus_label_lookup[bus_label] = bus_label
            bus_index_lookup[bus_idx] = bus_idx
            bus_index_lookup[bus_label] = bus_idx
        end


        @sync for day_group_id in day_group_ids

            @spawn begin

                buf = IOBuffer()  # per-task buffer flushed to the shared file under io_lock

                start_day = parameter(reference_of_this_year, 0, :repday_groups, "Start_Day_Id", day_group_id)
                end_day = parameter(reference_of_this_year, 0, :repday_groups, "End_Day_Id", day_group_id)

                # define days group index
                ids_d =  parameter(reference_of_this_year, 0, :repday_groups, "Day_Idx_List", day_group_id)

                for day_id in ids_d

                    NumDays = parameter(reference_of_this_year, 0, :repdays, "NumDays", day_id)

                    for hour_id in ids_h
                        for time_id in ids_t
                            for line_id in ids_k

                                idx_kdhty = string("(", line_id, ", ", day_id, ", ", hour_id, ", ", time_id, ", ", year, ")")
                                idx_ky = string("(", line_id, ", ", year, ")")

                                f_bus = parameter(reference_of_this_year, 0, :branch, "f_bus", line_id)
                                t_bus = parameter(reference_of_this_year, 0, :branch, "t_bus", line_id)
                                f_bus_idx = get(bus_index_lookup, f_bus, f_bus)
                                t_bus_idx = get(bus_index_lookup, t_bus, t_bus)

                                f_region = get(bus_label_lookup, f_bus, f_bus)
                                t_region = get(bus_label_lookup, t_bus, t_bus)

                                branch_data = reference_of_this_year[0][:branch][line_id]
                                line_length = get(branch_data, "length", get(branch_data, "Length", 0.0))
                                
                                original_rate = parameter(reference_of_this_year, 0, :branch, "rate_a", line_id) * pu_power_base       

                                if length(result_LC_GTEP_expansion) > 0
                                    # cumulative u_T_ky matches the flow limit rate_a*(1+u_T_ky)
                                    expansion_result = result_LC_GTEP_expansion["solution"]["expansion"][idx_ky]["u_T_ky"]
                                    final_rate = (1 + expansion_result) * original_rate
                                else
                                    expansion_result = 0.0
                                    final_rate = original_rate
                                end

                                flow = result["operational model result"][string(year)][string(day_group_id)]["solution"]["powerflow"][idx_kdhty]["f_kdhty"]

                                idx_string_ndhty_from_bus = string("(", f_bus_idx, ", ", day_id, ", ", hour_id, ", ", time_id, ", ", year,")")
                                idx_string_ndhty_to_bus = string("(", t_bus_idx, ", ", day_id, ", ", hour_id, ", ", time_id, ", ", year,")")

                                LMP_from_bus = result["operational model result"][string(year)][string(day_group_id)]["solution"]["dual"][idx_string_ndhty_from_bus]["LMP"]
                                LMP_to_bus = result["operational model result"][string(year)][string(day_group_id)]["solution"]["dual"][idx_string_ndhty_to_bus]["LMP"]

                                # congestion flag
                                congestion_flag = false
                                if abs(flow) >= final_rate
                                    congestion_flag = true
                                end

                                # wheeling cost disabled (would be abs(flow)*line_length*wheeling_cost_value*0.01)
                                wheeling_cost = 0.0

                                println(buf, join(csv_escape.(Any[
                                    Scenario,
                                    year,
                                    stage_calendar_year(setting, year),
                                    day_id,
                                    hour_id,
                                    time_id,
                                    line_id,
                                    f_bus,
                                    t_bus,
                                    f_region,
                                    t_region,
                                    original_rate,
                                    expansion_result,
                                    final_rate,
                                    flow,
                                    LMP_from_bus,
                                    LMP_to_bus,
                                    congestion_flag,
                                    wheeling_cost,
                                    line_length,
                                    NumDays,
                                ]), ","))


                            end
                        end
                    end
                end

                # flush this day-group's rows to the shared file
                lock(io_lock); write(io, take!(buf)); unlock(io_lock)
            end
        end
    end
    close(io)

    # @aleaf_info "[ALEAF LC_GTEP]: - Power flow solution reporting,"
end


# market_only writes just __market_OP.csv (needed by summaries) without the bulky per-year files
function report_result_dispatch_operation_GTEP(result::Dict{String, Any}; generate_file::Bool=true, market_only::Bool=false, day_group_subset::Union{Nothing,Vector{Int}}=nothing, file_suffix::String="")

    reference = result["operational model system reference"]
    setting   = result["setting"]

    # -----------------------------
    # CSV escaping helper
    # -----------------------------
    @inline function csv_escape(x)
        x === nothing && return ""
        if x isa Number
            return string(x)
        end
        s = String(x)
        if occursin(r"[,\n\"]", s)
            return "\"" * replace(s, "\"" => "\"\"") * "\""
        else
            return s
        end
    end

    # -----------------------------
    # Paths
    # -----------------------------
    output_path = parameter(reference["1"], 0, :output_path)
    case_name   = parameter(reference["1"], 0, :case_name)

    dispatch_file_prefix = joinpath(output_path, string(case_name, "__dispatch_OP_year_"))
    market_file_name     = joinpath(output_path, string(case_name, "__market_OP", file_suffix, ".csv"))
    policy_file_name     = joinpath(output_path, string(case_name, "__policy_slack_OP", file_suffix, ".csv"))
    demand_file_name     = joinpath(output_path, string(case_name, "__demand_response_OP", file_suffix, ".csv"))

    # -----------------------------
    # Labels
    # -----------------------------
    dispatch_label_list = [
        "Scenario","year","day","hour","time","unit_id","PLANT_NAME","bus_id","Bus_Name","Parent_Bus_Name","Region_Name",
        "Tech_ID","UnitGroup","Unit_Category","Unit_Report_Label_1","Unit_Report_Label_2","u_G_iy","u_ESE_iy","ICAP","sto_c_idhty","c_idhty","su_idhty",
        "g_idhty","reg_up_idhty","reg_dn_idhty","spin_idhty","nonspin_idhty","flex_up_idhty","flex_dn_idhty","curt_idhty",
        "chg_idhty","soc_idhty","hybrid_chg_idhty","hybrid_type","Inertia","MC","NumDays"
    ]

    market_label_list = [
        "Scenario","year","day","hour","time","bus_id","Bus_Name","Parent_Bus_Name","Region_Name",
        "load","scarcity_E","scarcity_SPIN","Demand_Reserve","scarcity_NSPIN","scarcity_FU","scarcity_FD",
        "LMP","RCP_RU","RCP_RD","RCP_Spin","RCP_NSpin","RCP_FU","RCP_FD","NumDays"
    ]

    policy_label_list = ["Scenario","year","day","bus_id","slack_CEG_ndy"]

    demand_response_label_list = [
        "Scenario","year","day","hour","time","PLANT_NAME","bus_idx","bus_name","Region_Name","Parent_Bus_Name","Online_Year",
        "UNITGROUP","UNIT_CATEGORY","UNIT_REPORT_LABEL_1","UNIT_REPORT_LABEL_2","CAP","INTERCON_LIM","Integer_Flag","Daily_DR_Limit_MWh","Num_DR_Segments",
        "Pct_MW_1","Pct_MW_2","Pct_MW_3","Pct_MW_4","Pct_MW_5","Price_1","Price_2","Price_3","Price_4","Price_5",
        "Hybrid_Gen","Hybrid_Gen_CAP","Hybrid_ES","Hybrid_ES_CAP",
        "lfl_lt","lfl_DR_lt","lfl_seg_lt_1","lfl_seg_lt_2","lfl_seg_lt_3","lfl_seg_lt_4","lfl_seg_lt_5",
        "lfl_ind_lt_1","lfl_ind_lt_2","lfl_ind_lt_3","lfl_ind_lt_4","lfl_ind_lt_5",
        "lfl_g_G_LFL_lt","lfl_g_G_Grid_lt","lfl_g_G_ES_lt",
        "lfl_g_ES_LFL_lt","lfl_g_ES_Grid_lt","lfl_chg_Grid_ES_lt","lfl_soc_lt"
    ]

    # -----------------------------
    # Indices
    # -----------------------------
    ids_h = collect(setting["run_H"])
    ids_t = collect(setting["run_T"])

    ids_y = try
        sort!(parse.(Int, collect(get_index(reference["1"], :planning_stages, 0))))
    catch
        sort!(collect(get_index(reference["1"], :planning_stages, 0)))
    end

    ids_i = try
        sort!(parse.(Int, collect(get_index(reference["1"], :gen_index, 0))))
    catch
        sort!(collect(get_index(reference["1"], :gen_index, 0)))
    end

    ids_n = try
        sort!(parse.(Int, collect(get_index(reference["1"], :bus, 0))))
    catch
        sort!(collect(get_index(reference["1"], :bus, 0)))
    end

    ids_d = try
        sort!(parse.(Int, collect(get_index(reference["1"], :repdays, 0))))
    catch
        sort!(collect(get_index(reference["1"], :repdays, 0)))
    end

    ids_p = try
        sort!(parse.(Int, collect(keys(reference["1"][0][:zone]["policy"]))))
    catch
        sort!(parse.(Int, collect(keys(reference["1"][0][:zone]["policy"]))))
    end

    ids_lfl = try
        sort!(parse.(Int, collect(get_index(reference["1"], :demand, 0))))
    catch
        sort!(collect(get_index(reference["1"], :demand, 0)))
    end

    dispatch_label_list = report_labels(dispatch_label_list)
    market_label_list = report_labels(market_label_list)
    policy_label_list = report_labels(policy_label_list; overrides = Dict("bus_id" => "Policy_Zone_ID"))
    demand_response_label_list = report_labels(demand_response_label_list)

    total_day_groups = setting["Simulation Configuration"]["NDAY_Groups_OP"]
    day_group_ids = day_group_subset === nothing ? collect(1:total_day_groups) : day_group_subset

    # -----------------------------
    # Scalars / flags
    # -----------------------------
    pu_power_base = setting["Simulation Setting"]["per_unit_base_value"]
    pu_econ_base  = setting["Simulation Setting"]["per_unit_econ_base_value"] / pu_power_base
    Scenario      = case_name

    reserve_cost_type = setting["Planning Design"]["reserve_cost_type_flag"]
    # Minimum regulation cost floor ($/MWh). base_reg here is already real $/MWh
    # (= reg_cost_param * MC, with MC = Annual_MC * pu_econ_base), so use the raw $/MWh value.
    min_reg_cost_pu   = get(setting["Planning Design"], "min_regulation_cost_value", 0.0)
    dispatch_mode_uc  = (setting["Simulation Configuration"]["Dispatch_Mode_in_OP"] == "Unit Commitment")
    ptc_enabled       = (setting["Simulation Configuration"]["PTC_Flag"] == true)

    op_result = result["operational model result"]

    # -----------------------------
    # Annual totals (scalars)
    # -----------------------------
    metric_names = (
        "Generation","Curtail","Storage_Charge","Reserve_RU","Reserve_RD","Reserve_Spin","Reserve_NSpin","Reserve_FU","Reserve_FD",
        "UnitRevenue_E","UnitRevenue_AS","UnitRevenue_CRED","UnitRevenue","UnitProfit",
        "Noload_Cost","StartUp_Cost","Operating_Cost","Generation_Cost","Charge_Cost","Fuel_Cost","Fuel_Consumption",
        "Reserve_Reg_Cost","Reserve_Spin_Cost","Reserve_NSpin_Cost","Reserve_Flex_Cost"
    )

    Annual_Gen_Info = Dict{Int,Dict{Int,Dict{String,Float64}}}()
    for i in ids_i
        Annual_Gen_Info[i] = Dict{Int,Dict{String,Float64}}()
        for y in ids_y
            Annual_Gen_Info[i][y] = Dict(m => 0.0 for m in metric_names)
        end
    end
    annual_lock = ReentrantLock()

    # Per-year NumDays-weighted scarcity totals (mirror the __market_OP.csv sum used by the system
    # summary) so the distributed path can obtain them without re-reading the concatenated CSV.
    scarcity_metric_names = ("scarcity_E","scarcity_SPIN","scarcity_NSPIN","scarcity_FU","scarcity_FD","Load_MWh")
    Annual_Scarcity_Info = Dict(y => Dict(m => 0.0 for m in scarcity_metric_names) for y in ids_y)

    # -----------------------------
    # Open files + headers (streaming)
    # -----------------------------
    market_io  = nothing
    policy_io  = nothing
    demand_io  = nothing
    dispatch_ios = Dict{Int,IO}()
    dispatch_locks = Dict{Int,ReentrantLock}()

    market_lock = ReentrantLock()
    policy_lock = ReentrantLock()
    demand_lock = ReentrantLock()

    if generate_file || market_only
        market_io = open(market_file_name, "w");  println(market_io, join(market_label_list, ","))
    end
    if generate_file
        policy_io = open(policy_file_name, "w");  println(policy_io, join(policy_label_list, ","))
        demand_io = open(demand_file_name, "w");  println(demand_io, join(demand_response_label_list, ","))

        for y in ids_y
            # workers hold a single-year view; only open files for stages actually solved
            haskey(op_result, string(y)) || continue
            io = open(string(dispatch_file_prefix, y, file_suffix, ".csv"), "w")
            println(io, join(dispatch_label_list, ","))
            dispatch_ios[y] = io
            dispatch_locks[y] = ReentrantLock()
        end
    end

    # ensure file handles close even on error
    try
        # -----------------------------
        # Main loops
        # -----------------------------
        for year in ids_y
            # workers hold a single-year view; skip stages this view does not solve
            haskey(op_result, string(year)) || continue
            reference_of_this_year = reference[string(year)]
            stage_length = parameter(reference_of_this_year, 0, :planning_stages, "stage_length", year)

            # -----------------------------
            # Precompute reserve zone per bus (per year)
            # -----------------------------
            bus_reserve_zone_idx = Dict{Int,Any}()
            for bus_idx in ids_n
                for zone_id in keys(reference_of_this_year[0][:zone]["reserve"])
                    if bus_idx in reference_of_this_year[0][:zone]["reserve"][zone_id]["aggregation_info"]["zone_bus_idx"]
                        bus_reserve_zone_idx[bus_idx] = zone_id
                        break
                    end
                end
            end

            # -----------------------------
            # Cache bus metadata (per year)
            # -----------------------------
            bus_name        = Dict{Int,Any}()
            region_name     = Dict{Int,Any}()
            parent_bus_name = Dict{Int,Any}()
            agg_bus_ids     = Dict{Int,Any}()
            agg_bus_load    = Dict{Int,Any}()

            for b in ids_n
                bus_name[b]        = parameter(reference_of_this_year, 0, :bus, "bus_i", b)
                region_name[b]     = reference_of_this_year[0][:bus][b]["region_config"]["region_name"]
                parent_bus_name[b] = reference_of_this_year[0][:bus][b]["region_config"]["parent_bus_name"]
                agg_bus_ids[b]     = reference_of_this_year[0][:bus][b]["aggregation_info"]["aggregated_regions_bus_i"]
                agg_bus_load[b]    = reference_of_this_year[0][:bus][b]["aggregation_info"]["original_load_(bus_i, MW)"]
            end

            # -----------------------------
            # Cache generator metadata (per year; depends on i)
            # -----------------------------
            bus_of_i         = Dict{Int,Int}()
            tech_of_i        = Dict{Int,Int}()
            unitgroup_of_i   = Dict{Int,Any}()
            unitcat_of_i     = Dict{Int,Any}()
            unitlabel1_of_i  = Dict{Int,Any}()
            unitlabel2_of_i  = Dict{Int,Any}()
            plant_of_i       = Dict{Int,Any}()
            hybrid_type_of_i = Dict{Int,Any}()
            cap_of_i_mw      = Dict{Int,Float64}()
            pmax_of_i        = Dict{Int,Float64}()
            pf_of_i          = Dict{Int,Float64}()
            inertiaH_of_i    = Dict{Int,Float64}()
            hr_of_i          = Dict{Int,Float64}()
            timeseries_tag_of_i = Dict{Int,Any}()
            fixed_profile_type_of_i = Dict{Int,Any}()

            has_commitment_of_i        = Dict{Int,Bool}()
            storage_commitment_of_i    = Dict{Int,Bool}()
            ptc_flag_of_i              = Dict{Int,Bool}()

            reg_cost_param   = Dict{Int,Float64}()
            spin_cost_param  = Dict{Int,Float64}()
            nspin_cost_param = Dict{Int,Float64}()
            flex_cost_param  = Dict{Int,Float64}()

            nlc_param = Dict{Int,Float64}()
            suc_param = Dict{Int,Float64}()

            for i in ids_i
                b    = parameter(reference_of_this_year, 0, :gen_index, "bus_idx", i)
                tech = parameter(reference_of_this_year, 0, :gen_index, "genco_tech_id", i)

                bus_of_i[i]       = b
                tech_of_i[i]      = tech
                unitgroup_of_i[i] = parameter(reference_of_this_year, 0, :gen_index, "UNIT_GROUP", i)

                ucat  = parameter(reference_of_this_year, b, :gen_bus, tech, "UNIT_CATEGORY")
                unitcat_of_i[i]  = ucat
                unitlabel1_of_i[i] = parameter(reference_of_this_year, b, :gen_bus, tech, "UNIT_REPORT_LABEL_1")
                unitlabel2_of_i[i] = parameter(reference_of_this_year, b, :gen_bus, tech, "UNIT_REPORT_LABEL_2")

                pname = parameter(reference_of_this_year, b, :gen_bus, tech, "NEW_Gen_UID")
                if haskey(reference_of_this_year[b][:gen_bus][tech], "PLANT_NAME")
                    pname = parameter(reference_of_this_year, b, :gen_bus, tech, "PLANT_NAME")
                end
                plant_of_i[i] = pname

                hybrid_type_of_i[i] = parameter(reference_of_this_year, 0, :gen_index, "hybrid_type", i)

                cap_of_i_mw[i]   = parameter(reference_of_this_year, b, :gen_bus, tech, "CAP") * pu_power_base
                pmax_of_i[i]     = parameter(reference_of_this_year, b, :gen_bus, tech, "PMAX")
                pf_of_i[i]       = parameter(reference_of_this_year, b, :gen_bus, tech, "Power_Factor")
                inertiaH_of_i[i] = parameter(reference_of_this_year, b, :gen_bus, tech, "Inertia_Constant")
                hr_of_i[i]       = parameter(reference_of_this_year, b, :gen_bus, tech, "HR")
                timeseries_tag_of_i[i] = parameter(reference_of_this_year, b, :gen_bus, tech, "Timeseries_Tag")

                fixed_profile_type_of_i[i] = parameter(reference_of_this_year, b, :gen_bus, tech, "Profile_Type")

                has_commitment_of_i[i]     = (parameter(reference_of_this_year, b, :gen_bus, tech, "Commitment") == true)
                storage_commitment_of_i[i] = (parameter(reference_of_this_year, b, :gen_bus, tech, "Storage Commitment") == true)
                ptc_flag_of_i[i]           = ptc_enabled && (parameter(reference_of_this_year, b, :gen_bus, tech, "PTC Flag") == true)

                reg_cost_param[i]   = float(parameter(reference_of_this_year, b, :gen_bus, tech, "reg_cost"))
                spin_cost_param[i]  = float(parameter(reference_of_this_year, b, :gen_bus, tech, "spin_cost"))
                nspin_cost_param[i] = float(parameter(reference_of_this_year, b, :gen_bus, tech, "nspin_cost"))
                flex_cost_param[i]  = float(parameter(reference_of_this_year, b, :gen_bus, tech, "flex_cost"))

                nlc_param[i] = float(parameter(reference_of_this_year, b, :gen_bus, tech, "NLC"))
                suc_param[i] = float(parameter(reference_of_this_year, b, :gen_bus, tech, "SUC"))
            end

            function get_fixed_profile_shape(i::Int, d::Int, h::Int, t::Int, y::Int)
                type = fixed_profile_type_of_i[i]
                type == "NA" && return 0.0

                b = bus_of_i[i]
                Timeseries_Tag = timeseries_tag_of_i[i]

                if Timeseries_Tag == "LOCAL"
                    type_key = if type == "wind_ons"
                        "wind_ons_shape"
                    elseif type == "wind_ofs"
                        "wind_ofs_shape"
                    elseif type == "pv"
                        "pv_shape"
                    elseif type == "hydro"
                        "hydro_shape"
                    else
                        string(type, "_shape")
                    end
                    return get_vre_zdt_shape(reference_of_this_year[0], y, b, d, h, type_key)
                else
                    if type == "wind_ons"
                        return parameter(reference_of_this_year, 0, :planning_stages, "repdays", "data", "wind_ons", d, h, t, y)[Timeseries_Tag]
                    elseif type == "wind_ofs"
                        return parameter(reference_of_this_year, 0, :planning_stages, "repdays", "data", "wind_ofs", d, h, t, y)[Timeseries_Tag]
                    elseif type == "csp"
                        return parameter(reference_of_this_year, 0, :planning_stages, "repdays", "data", "csp", d, h, t, y)[Timeseries_Tag]
                    end
                end

                return 0.0
            end

            # -----------------------------
            # Cache large-load static metadata (per year)
            # -----------------------------
            lfl_static = Dict{Int,NamedTuple}()
            for lfl_idx in ids_lfl
                b   = parse(Int, parameter(reference_of_this_year, 0, :demand, "bus_idx", lfl_idx))
                CAP = parameter(reference_of_this_year, 0, :demand, "CAP", lfl_idx) * pu_power_base

                lfl_static[lfl_idx] = (
                    PLANT_NAME = parameter(reference_of_this_year, 0, :demand, "PLANT_NAME", lfl_idx),
                    bus_idx = b,
                    bus_name = bus_name[b],
                    Region_Name = region_name[b],
                    Parent_Bus_Name = parent_bus_name[b],
                    Online_Year = parameter(reference_of_this_year, 0, :demand, "Online_Year", lfl_idx),
                    UNITGROUP = parameter(reference_of_this_year, 0, :demand, "UNITGROUP", lfl_idx),
                    UNIT_CATEGORY = parameter(reference_of_this_year, 0, :demand, "UNIT_CATEGORY", lfl_idx),
                    UNIT_REPORT_LABEL_1 = parameter(reference_of_this_year, 0, :demand, "UNIT_REPORT_LABEL_1", lfl_idx),
                    UNIT_REPORT_LABEL_2 = parameter(reference_of_this_year, 0, :demand, "UNIT_REPORT_LABEL_2", lfl_idx),
                    CAP = CAP,
                    INTERCON_LIM = parameter(reference_of_this_year, 0, :demand, "INTERCON_LIM", lfl_idx) * pu_power_base,
                    Integer_Flag = parameter(reference_of_this_year, 0, :demand, "Integer_Flag", lfl_idx),
                    Daily_DR_Limit_MWh = parameter(reference_of_this_year, 0, :demand, "Daily_DR_Limit_MWh", lfl_idx) * pu_power_base,
                    Num_DR_Segments = parameter(reference_of_this_year, 0, :demand, "Num_DR_Segments", lfl_idx),
                    Pct_MW_1 = parameter(reference_of_this_year, 0, :demand, "Pct_MW_1", lfl_idx) * CAP,
                    Pct_MW_2 = parameter(reference_of_this_year, 0, :demand, "Pct_MW_2", lfl_idx) * CAP,
                    Pct_MW_3 = parameter(reference_of_this_year, 0, :demand, "Pct_MW_3", lfl_idx) * CAP,
                    Pct_MW_4 = parameter(reference_of_this_year, 0, :demand, "Pct_MW_4", lfl_idx) * CAP,
                    Pct_MW_5 = parameter(reference_of_this_year, 0, :demand, "Pct_MW_5", lfl_idx) * CAP,
                    Price_1 = parameter(reference_of_this_year, 0, :demand, "Price_1", lfl_idx) * pu_econ_base,
                    Price_2 = parameter(reference_of_this_year, 0, :demand, "Price_2", lfl_idx) * pu_econ_base,
                    Price_3 = parameter(reference_of_this_year, 0, :demand, "Price_3", lfl_idx) * pu_econ_base,
                    Price_4 = parameter(reference_of_this_year, 0, :demand, "Price_4", lfl_idx) * pu_econ_base,
                    Price_5 = parameter(reference_of_this_year, 0, :demand, "Price_5", lfl_idx) * pu_econ_base,
                    Hybrid_Gen = parameter(reference_of_this_year, 0, :demand, "Hybrid_Gen", lfl_idx),
                    Hybrid_Gen_CAP = parameter(reference_of_this_year, 0, :demand, "Hybrid_Gen_CAP", lfl_idx),
                    Hybrid_ES = parameter(reference_of_this_year, 0, :demand, "Hybrid_ES", lfl_idx),
                    Hybrid_ES_CAP = parameter(reference_of_this_year, 0, :demand, "Hybrid_ES_CAP", lfl_idx),
                )
            end

            # -----------------------------
            # Loop day groups (OP decomposition)
            # -----------------------------
            for day_group_id in day_group_ids
                sol = op_result[string(year)][string(day_group_id)]["solution"]

                ids_d_new = collect(parameter(reference_of_this_year, 0, :repday_groups, "Day_Idx_List", day_group_id))

                Threads.@threads for dd in eachindex(ids_d_new)
                    day_id = ids_d_new[dd]

                    # thread-local buffers
                    buf_dispatch = IOBuffer()
                    buf_market   = IOBuffer()
                    buf_policy   = IOBuffer()
                    buf_demand   = IOBuffer()

                    NumDays = parameter(reference_of_this_year, 0, :repdays, "NumDays", day_id)
                    scale   = NumDays * stage_length

                    # thread-local annual accumulation for this day (no locks)
                    localA = Dict{Int,Dict{String,Float64}}()
                    @inline function get_localA(i)
                        get!(localA, i) do
                            Dict(m => 0.0 for m in metric_names)
                        end
                    end

                    # thread-local scarcity accumulation (NumDays-weighted, matches market-CSV sum)
                    localS = Dict(m => 0.0 for m in scarcity_metric_names)

                    # -----------------------------
                    # hour/time loops
                    # -----------------------------
                    for hour_id in ids_h, time_id in ids_t

                        # -----------------------------
                        # Dispatch (generators)
                        # -----------------------------
                        for i in ids_i
                            b    = bus_of_i[i]
                            tech = tech_of_i[i]

                            idx_idhty = "($i, $day_id, $hour_id, $time_id, $year)"

                            u_G_iy = sol["expansion"]["($i, $year)"]["u_G_iy"]
                            CAP    = cap_of_i_mw[i]
                            ICAP   = CAP * u_G_iy

                            disp = sol["dispatch"][idx_idhty]
                            reserve_disp = get(sol["reserve"], idx_idhty, Dict{String,Any}())
                            g_idhty       = disp["g_idhty"]
                            reg_up_idhty  = get(reserve_disp, "reg_up_idhty", 0.0)
                            reg_dn_idhty  = get(reserve_disp, "reg_dn_idhty", 0.0)
                            spin_idhty    = get(reserve_disp, "spin_idhty", 0.0)
                            flex_up_idhty = get(reserve_disp, "flex_up_idhty", 0.0)
                            flex_dn_idhty = get(reserve_disp, "flex_dn_idhty", 0.0)
                            nonspin_idhty = get(reserve_disp, "nonspin_idhty", 0.0)

                            curt_idhty = 0.0
                            shape = get_fixed_profile_shape(i, day_id, hour_id, time_id, year)
                            if shape > 0.0
                                availability = shape * u_G_iy * CAP * pmax_of_i[i]
                                curt_idhty = max(availability - g_idhty, 0.0)
                            end

                            c_idhty  = 0.0
                            su_idhty = 0.0
                            if haskey(sol, "commitment") && haskey(sol["commitment"], idx_idhty)
                                com = sol["commitment"][idx_idhty]
                                c_idhty  = com["c_idhty"]
                                su_idhty = com["su_idhty"]
                            end

                            Unit_Category = unitcat_of_i[i]

                            u_ESE_iy      = 0.0
                            chg_idhty     = 0.0
                            soc_idhty     = 0.0
                            sto_c_idhty   = 0.0

                            if Unit_Category == "STORAGE"
                                u_ESE_iy  = sol["expansion"]["($i, $year)"]["u_ESE_iy"]
                                sto       = sol["storage"][idx_idhty]
                                chg_idhty = sto["chg_idhty"]
                                soc_idhty = sto["soc_idhty"]
                                if storage_commitment_of_i[i]
                                    sto_c_idhty = sol["storage_commitment"][idx_idhty]["sto_c_idhty"]
                                end
                            end

                            hybrid_type = hybrid_type_of_i[i]
                            hybrid_chg_idhty = 0.0
                            if hybrid_type == "GEN"
                                hybrid_chg_idhty = sol["hybrid"][idx_idhty]["g_G_ES_idhty"]
                            end

                            MC = parameter(reference_of_this_year, b, :gen_bus, tech, "Annual_MC")[string(year)][day_id] * pu_econ_base

                            inertia = 0.0
                            if haskey(sol, "commitment")
                                inertia = (CAP / pf_of_i[i] * inertiaH_of_i[i]) * c_idhty
                            end

                            if generate_file
                                row = Any[
                                    Scenario, year, stage_calendar_year(setting, year), day_id, hour_id, time_id, i, plant_of_i[i], b,
                                    bus_name[b], parent_bus_name[b], region_name[b],
                                    tech, unitgroup_of_i[i], Unit_Category, unitlabel1_of_i[i], unitlabel2_of_i[i],
                                    u_G_iy, u_ESE_iy, ICAP, sto_c_idhty, c_idhty, su_idhty,
                                    g_idhty, reg_up_idhty, reg_dn_idhty, spin_idhty, nonspin_idhty,
                                    flex_up_idhty, flex_dn_idhty, curt_idhty, chg_idhty, soc_idhty,
                                    hybrid_chg_idhty, hybrid_type, inertia, MC, NumDays
                                ]
                                println(buf_dispatch, join(csv_escape.(row), ","))
                            end

                            # -----------------------------
                            # Prices + revenues/costs for annuals
                            # -----------------------------
                            rz = bus_reserve_zone_idx[b]
                            idx_ndhty = "($b, $day_id, $hour_id, $time_id, $year)"
                            idx_zdhty = "($rz, $day_id, $hour_id, $time_id, $year)"

                            dualn = sol["dual"][idx_ndhty]
                            dualz = sol["dual"][idx_zdhty]

                            LMP = dualn["LMP"]
                            RCP_RU   = dualz["Reg_Up_Price"]
                            RCP_RD   = dualz["Reg_Dn_Price"]
                            RCP_Cont = dualz["Cont_Res_Price"]
                            RCP_NSPIN= dualz["NonSpin_Res_Price"]
                            RCP_FU   = dualz["Flex_Up_Price"]
                            RCP_FD   = dualz["Flex_Dn_Price"]

                            gen_cost = MC * g_idhty

                            # PTC (OP: only current stage year)
                            gen_credit = 0.0
                            if ptc_flag_of_i[i]
                                ptc_year = min(2050, parameter(reference_of_this_year, 0, :planning_stages, "year", year))
                                ptc_tbl  = parameter(reference_of_this_year, b, :gen_bus, tech, "PTC")
                                gen_credit = (ptc_tbl[string(ptc_year)] * 10 * pu_econ_base) * g_idhty
                            end

                            base_reg  = reg_cost_param[i]
                            base_spin = spin_cost_param[i]
                            base_nsp  = nspin_cost_param[i]
                            base_flex = flex_cost_param[i]

                            if reserve_cost_type == "percentage"
                                base_reg  *= MC
                                base_spin *= MC
                                base_nsp  *= MC
                                base_flex *= MC
                            else
                                base_reg  *= pu_econ_base
                                base_spin *= pu_econ_base
                                base_nsp  *= pu_econ_base
                                base_flex *= pu_econ_base
                            end

                            # Apply the minimum regulation cost floor so reported reg cost
                            # matches what the objective optimized (regulation never free).
                            base_reg = max(base_reg, min_reg_cost_pu)

                            reg_cost   = base_reg  * (reg_up_idhty + reg_dn_idhty)
                            cont_cost  = base_spin * spin_idhty
                            nspin_cost = base_nsp  * nonspin_idhty
                            flex_cost  = base_flex * (flex_up_idhty + flex_dn_idhty)

                            noloadc = 0.0
                            su_cost = 0.0
                            if dispatch_mode_uc && (Unit_Category != "STORAGE") && has_commitment_of_i[i]
                                noloadc = CAP * c_idhty * nlc_param[i] * pu_econ_base
                                su_cost = su_idhty * suc_param[i] * pu_econ_base
                            end

                            Revenue_E    = LMP * g_idhty
                            Revenue_AS   = RCP_RU*reg_up_idhty + RCP_RD*reg_dn_idhty + RCP_Cont*spin_idhty +
                                           RCP_NSPIN*nonspin_idhty + RCP_FU*flex_up_idhty + RCP_FD*flex_dn_idhty
                            Revenue_CRED = gen_credit
                            Revenue      = Revenue_E + Revenue_AS + Revenue_CRED

                            cost = gen_cost + reg_cost + cont_cost + flex_cost + nspin_cost + noloadc + su_cost

                            charge_cost = 0.0
                            if Unit_Category == "STORAGE"
                                cost += LMP * chg_idhty
                                charge_cost = LMP * chg_idhty
                            end

                            FuelConsumption = g_idhty * hr_of_i[i]
                            fuel_cost = parameter(reference_of_this_year, b, :gen_bus, tech, "Annual_FC")[string(year)][day_id] * pu_econ_base
                            FuelCost  = FuelConsumption * fuel_cost

                            A = get_localA(i)
                            A["Generation"]         += g_idhty * scale
                            A["Curtail"]            += curt_idhty * scale
                            A["Storage_Charge"]     += chg_idhty * scale
                            A["Reserve_RU"]         += reg_up_idhty * scale
                            A["Reserve_RD"]         += reg_dn_idhty * scale
                            A["Reserve_Spin"]       += spin_idhty * scale
                            A["Reserve_NSpin"]      += nonspin_idhty * scale
                            A["Reserve_FU"]         += flex_up_idhty * scale
                            A["Reserve_FD"]         += flex_dn_idhty * scale

                            A["UnitRevenue_E"]      += Revenue_E * scale
                            A["UnitRevenue_AS"]     += Revenue_AS * scale
                            A["UnitRevenue_CRED"]   += Revenue_CRED * scale
                            A["UnitRevenue"]        += Revenue * scale

                            A["UnitProfit"]         += (Revenue - cost) * scale
                            A["Noload_Cost"]        += noloadc * scale
                            A["StartUp_Cost"]       += su_cost * scale
                            A["Operating_Cost"]     += cost * scale

                            A["Generation_Cost"]    += gen_cost * scale
                            A["Charge_Cost"]        += charge_cost * scale
                            A["Fuel_Cost"]          += FuelCost * scale
                            A["Fuel_Consumption"]   += FuelConsumption * scale

                            A["Reserve_Reg_Cost"]   += base_reg  * (reg_up_idhty + reg_dn_idhty) * scale
                            A["Reserve_Spin_Cost"]  += base_spin * spin_idhty * scale
                            A["Reserve_NSpin_Cost"] += base_nsp  * nonspin_idhty * scale
                            A["Reserve_Flex_Cost"]  += base_flex * (flex_up_idhty + flex_dn_idhty) * scale
                        end # i loop

                        # -----------------------------
                        # Market (buses)
                        # -----------------------------
                        if generate_file || market_only
                            emitted_reserve_zones = Set{Any}()   # zonal reserve shortfall is a per-zone quantity: report it once, not on every bus in the zone
                            for b in ids_n
                                rz = bus_reserve_zone_idx[b]

                                Demand = 0.0
                                load_tbl = parameter(reference_of_this_year, 0, :planning_stages, "repdays", year)[string(day_id)]["data"][string(hour_id)]["1"]["load"]
                                region_map = get(reference_of_this_year[0], :profile_data_region_map, nothing)
                                for ba_id in agg_bus_ids[b]
                                    Demand += load_tbl[profile_data_region(region_map, "load", ba_id)] * agg_bus_load[b][ba_id] * get_load_growth_factor(reference_of_this_year[0], year, ba_id)
                                end
                                Demand *= pu_power_base

                                idx_ndhty = "($b, $day_id, $hour_id, $time_id, $year)"
                                idx_zdhty = "($rz, $day_id, $hour_id, $time_id, $year)"

                                dualn = sol["dual"][idx_ndhty]
                                dualz = sol["dual"][idx_zdhty]

                                LMP = dualn["LMP"]
                                RCP_RU   = dualz["Reg_Up_Price"]
                                RCP_RD   = dualz["Reg_Dn_Price"]
                                RCP_Cont = dualz["Cont_Res_Price"]
                                RCP_NSPIN= dualz["NonSpin_Res_Price"]
                                RCP_FU   = dualz["Flex_Up_Price"]
                                RCP_FD   = dualz["Flex_Dn_Price"]

                                scarcity_E     = get(sol["scarcity"][idx_ndhty], "ens_ndhty", 0.0)
                                if rz in emitted_reserve_zones
                                    scarcity_CONT = demand_reserve = scarcity_NSPIN = scarcity_FU = scarcity_FD = 0.0
                                else
                                    push!(emitted_reserve_zones, rz)
                                    scarcity_CONT  = get(sol["scarcity_reserve"][idx_zdhty], "rns_cont_zdhty", 0.0)
                                    demand_reserve = get(sol["scarcity_reserve"][idx_zdhty], "demand_reserve_zdhty", 0.0)
                                    scarcity_NSPIN = get(sol["scarcity_reserve"][idx_zdhty], "rns_nonspin_zdhty", 0.0)
                                    scarcity_FU    = get(sol["scarcity_reserve"][idx_zdhty], "rns_flex_up_zdhty", 0.0)
                                    scarcity_FD    = get(sol["scarcity_reserve"][idx_zdhty], "rns_flex_dn_zdhty", 0.0)
                                end

                                row = Any[
                                    Scenario, year, stage_calendar_year(setting, year), day_id, hour_id, time_id, b,
                                    bus_name[b], parent_bus_name[b], region_name[b],
                                    Demand, scarcity_E, scarcity_CONT, demand_reserve, scarcity_NSPIN, scarcity_FU, scarcity_FD,
                                    LMP, RCP_RU, RCP_RD, RCP_Cont, RCP_NSPIN, RCP_FU, RCP_FD,
                                    NumDays
                                ]
                                println(buf_market, join(csv_escape.(row), ","))

                                # NumDays-weighted, exactly the sum the system summary takes over the CSV
                                localS["scarcity_E"]     += scarcity_E     * NumDays
                                localS["scarcity_SPIN"]  += scarcity_CONT  * NumDays
                                localS["scarcity_NSPIN"] += scarcity_NSPIN * NumDays
                                localS["scarcity_FU"]    += scarcity_FU    * NumDays
                                localS["scarcity_FD"]    += scarcity_FD    * NumDays
                                localS["Load_MWh"]       += Demand * NumDays   # dispatched (representative-day weighted) load
                            end
                        end

                        # -----------------------------
                        # Demand response (large load)
                        # -----------------------------
                        if haskey(sol, "demand")
                            if generate_file
                                demand_solution = sol["demand"]
                                for lfl_idx in ids_lfl
                                    st = lfl_static[lfl_idx]
                                    idx_ldhty = "($lfl_idx, $day_id, $hour_id, $time_id, $year)"

                                    lfl_lt    = demand_solution[idx_ldhty]["lfl_ldhty"] * pu_power_base
                                    lfl_DR_lt = demand_solution[idx_ldhty]["lfl_DR_ldhty"] * pu_power_base

                                    # avoid allocations: fixed-size scalars
                                    seg1=0.0; seg2=0.0; seg3=0.0; seg4=0.0; seg5=0.0
                                    ind1=0;   ind2=0;   ind3=0;   ind4=0;   ind5=0

                                    for seg_id in 1:st.Num_DR_Segments
                                        idx = "($lfl_idx, $seg_id, $day_id, $hour_id, $time_id, $year)"
                                        rec = demand_solution[idx]
                                        v = rec["lfl_seg_lsdhty"] * pu_power_base
                                        w = (seg_id > 1) ? rec["lfl_ind_lsdhty"] : 1
                                        if seg_id == 1; seg1=v; ind1=w
                                        elseif seg_id == 2; seg2=v; ind2=w
                                        elseif seg_id == 3; seg3=v; ind3=w
                                        elseif seg_id == 4; seg4=v; ind4=w
                                        elseif seg_id == 5; seg5=v; ind5=w
                                        end
                                    end

                                    lfl_g_G_LFL = 0.0; lfl_g_G_Grid = 0.0; lfl_g_G_ES = 0.0
                                    if st.Hybrid_Gen != "NA"
                                        lfl_g_G_LFL  = demand_solution[idx_ldhty]["lfl_g_G_LFL_ldhty"] * pu_power_base
                                        lfl_g_G_Grid = demand_solution[idx_ldhty]["lfl_g_G_Grid_ldhty"] * pu_power_base
                                        if st.Hybrid_ES != "NA"
                                            lfl_g_G_ES = demand_solution[idx_ldhty]["lfl_g_G_ES_ldhty"] * pu_power_base
                                        end
                                    end

                                    lfl_g_ES_LFL = 0.0; lfl_g_ES_Grid = 0.0; lfl_chg_Grid_ES = 0.0; lfl_soc = 0.0
                                    if st.Hybrid_ES != "NA"
                                        lfl_g_ES_LFL    = demand_solution[idx_ldhty]["lfl_g_ES_LFL_ldhty"] * pu_power_base
                                        lfl_g_ES_Grid   = demand_solution[idx_ldhty]["lfl_g_ES_Grid_ldhty"] * pu_power_base
                                        lfl_chg_Grid_ES = demand_solution[idx_ldhty]["lfl_chg_Grid_ES_ldhty"] * pu_power_base
                                        lfl_soc         = demand_solution[idx_ldhty]["lfl_soc_ldhty"] * pu_power_base
                                    end

                                    row = Any[
                                        Scenario, year, stage_calendar_year(setting, year), day_id, hour_id, time_id,
                                        st.PLANT_NAME, st.bus_idx, st.bus_name, st.Region_Name, st.Parent_Bus_Name, st.Online_Year,
                                        st.UNITGROUP, st.UNIT_CATEGORY, st.UNIT_REPORT_LABEL_1, st.UNIT_REPORT_LABEL_2, st.CAP, st.INTERCON_LIM, st.Integer_Flag,
                                        st.Daily_DR_Limit_MWh, st.Num_DR_Segments,
                                        st.Pct_MW_1, st.Pct_MW_2, st.Pct_MW_3, st.Pct_MW_4, st.Pct_MW_5,
                                        st.Price_1, st.Price_2, st.Price_3, st.Price_4, st.Price_5,
                                        st.Hybrid_Gen, st.Hybrid_Gen_CAP, st.Hybrid_ES, st.Hybrid_ES_CAP,
                                        lfl_lt, lfl_DR_lt,
                                        seg1, seg2, seg3, seg4, seg5,
                                        ind1, ind2, ind3, ind4, ind5,
                                        lfl_g_G_LFL, lfl_g_G_Grid, lfl_g_G_ES,
                                        lfl_g_ES_LFL, lfl_g_ES_Grid, lfl_chg_Grid_ES, lfl_soc
                                    ]
                                    println(buf_demand, join(csv_escape.(row), ","))
                                end
                            end
                        end

                    end # hour/time loop

                    # -----------------------------
                    # Policy (once per day)
                    # -----------------------------
                    if generate_file
                        for z in ids_p
                            slack_CEG = 0.0
                            if setting["Simulation Configuration"]["Clean_Energy_Generation_Target_OP_Flag"] == true
                                idx_ndy = "($z, $day_id, $year)"
                                slack_CEG = sol["slack"][idx_ndy]["slack_CEG_ndy"] * pu_power_base
                            end
                            println(buf_policy, join(csv_escape.(Any[Scenario, year, stage_calendar_year(setting, year), day_id, z, slack_CEG]), ","))
                        end
                    end

                    # -----------------------------
                    # Merge annual totals (one lock per day)
                    # -----------------------------
                    lock(annual_lock)
                    for (i, Ai) in localA
                        A = Annual_Gen_Info[i][year]
                        for m in metric_names
                            A[m] += Ai[m]
                        end
                    end
                    S = Annual_Scarcity_Info[year]
                    for m in scarcity_metric_names
                        S[m] += localS[m]
                    end
                    unlock(annual_lock)

                    # -----------------------------
                    # Flush buffers (locks)
                    # -----------------------------
                    if generate_file || market_only
                        lock(market_lock);          write(market_io,          take!(buf_market));  unlock(market_lock)
                    end
                    if generate_file
                        lock(dispatch_locks[year]); write(dispatch_ios[year], take!(buf_dispatch)); unlock(dispatch_locks[year])
                        lock(demand_lock);          write(demand_io,          take!(buf_demand));  unlock(demand_lock)
                        lock(policy_lock);          write(policy_io,          take!(buf_policy));  unlock(policy_lock)
                    end

                end # threads day
            end # day_group
        end # year

    finally
        if generate_file || market_only
            try close(market_io) catch end
        end
        if generate_file
            try close(policy_io) catch end
            try close(demand_io) catch end
            for y in keys(dispatch_ios)
                try close(dispatch_ios[y]) catch end
            end
        end
    end

    # save annual_gen_info_OP (scalar totals)
    result["operational model system reference"]["annual_gen_info_OP"] = Annual_Gen_Info
    # per-year scarcity totals (distributed path reads these back instead of the concatenated CSV)
    result["operational model system reference"]["annual_scarcity_info_OP"] = Annual_Scarcity_Info

    # @aleaf_info "[ALEAF LC_GTEP]: - OP Dispatch solution reporting,"
    return nothing
end


function report_result_water_management_operation_GTEP(result::Dict{String, Any}; generate_file::Bool=true)
    
    pu_base = result["setting"]["Simulation Setting"]["per_unit_base_value"]
    
    dispatch = DataFrame(
                        year = Int[],
                        day = Int[], 
                        hour = Int[], 
                        time = Int[], 
                        segment = Int[], 

                        unit_id = Int[],
                        bus_id = Int[],
                        Unit_Report_Label_1 = String[],
                        Unit_Report_Label_2 = String[],

                        water_use_ijdhty = Float64[]
                        )


    # grouping days
    total_day_groups = result["setting"]["Simulation Configuration"]["NDAY_Groups_OP"]
    ids_i = result["operational model system reference"][0][:water_management]["gen_index"]
    ids_j = result["operational model system reference"][0][:water_management]["segment_index"]

    ids_y = [(y) for (y) in get_index(result["operational model system reference"], :planning_stages, 0)]
    ids_h = [(h) for (h) in result["setting"]["run_H"]]
    ids_t = [(t) for (t) in result["setting"]["run_T"]]

    for y in ids_y

        for day_group_id in 1:total_day_groups

            start_day = result["operational model system reference"][0][:repday_groups][day_group_id]["Start_Day_Id"]
            end_day = result["operational model system reference"][0][:repday_groups][day_group_id]["End_Day_Id"]
            
            # define days group index
            ids_d = result["operational model system reference"][0][:repday_groups][day_group_id]["Day_Idx_List"]

            for d in ids_d
                for h in ids_h
                    for t in ids_t
                        for i in ids_i

                            year = y
                            day = d
                            hour = h
                            time = t

                            unit_id = i
                            bus_idx = result["operational model system reference"][0][:gen_index][i]["bus_idx"]
                            tech_idx = result["operational model system reference"][0][:gen_index][i]["genco_tech_id"]
                            Unit_Report_Label_1 = parameter(result["operational model system reference"], bus_idx, :gen_bus, tech_idx, "UNIT_REPORT_LABEL_1")
                            Unit_Report_Label_2 = parameter(result["operational model system reference"], bus_idx, :gen_bus, tech_idx, "UNIT_REPORT_LABEL_2")

                            for j in ids_j

                                segment = j
                                idx_idhty = string("(", i, ", ", j, ", ", d, ", ", h, ", ", t, ", ", y, ")")
                                water_use_ijdhty = result["operational model result"][string(y)][string(day_group_id)]["solution"]["water_management"][idx_idhty]["water_use_ijdhty"]
                            

                                push!(dispatch, [
                                            year,
                                            day,
                                            hour,
                                            time,
                                            segment,

                                            unit_id,
                                            bus_idx,
                                            Unit_Report_Label_1,
                                            Unit_Report_Label_2,

                                            water_use_ijdhty
                                ])
                            end
                        end
                    end
                end
            end
        end
    end

    if generate_file == true

        output_path = result["operational model system reference"][0][:output_path]
        file_name = string(result["operational model system reference"][0][:case_name], "__water_management_OP.csv")
        full_path = joinpath(output_path, file_name)
        out = copy(dispatch)
        rename!(out, Symbol.(names(out)) .=> ra_report_names(names(out)))
        insertcols!(out, 2, :Year => [stage_calendar_year(result["setting"], y) for y in out.Stage])
        CSV.write(full_path, out)
    end

    return dispatch
end


# write_csv=false computes annual gen/scarcity metrics with no files; market_only keeps just the
# market piece while suppressing the bulky per-year files.
function report_result_dispatch_expansion_GTEP(result::Dict{String, Any}, reference, setting; write_csv::Bool = true, years=nothing, out_dir=nothing, market_only::Bool = false)

    # -----------------------------
    # CSV escaping helper
    # -----------------------------
    @inline function csv_escape(x)
        x === nothing && return ""
        if x isa Number
            return string(x)
        end
        s = String(x)
        if occursin(r"[,\n\"]", s)
            return "\"" * replace(s, "\"" => "\"\"") * "\""
        else
            return s
        end
    end

    # -----------------------------
    # Paths
    # -----------------------------
    output_path = parameter(reference, 0, :output_path)
    case_name   = parameter(reference, 0, :case_name)
    dir = out_dir === nothing ? output_path : out_dir
    out_dir === nothing || mkpath(dir)

    dispatch_file_prefix = joinpath(dir, string(case_name, "__dispatch_EXP_year_"))
    market_file_name     = joinpath(dir, string(case_name, "__market_EXP.csv"))
    policy_file_name     = joinpath(dir, string(case_name, "__policy_slack_EXP.csv"))
    demand_file_name     = joinpath(dir, string(case_name, "__demand_response_EXP.csv"))

    # -----------------------------
    # Labels
    # -----------------------------
    dispatch_label_list = [
        "Scenario","year","Stochastic_scenario_ID","day","hour","time","unit_id","PLANT_NAME","bus_id",
        "Bus_Name","Parent_Bus_Name","Region_Name","Tech_ID","UnitGroup","Unit_Category","Unit_Report_Label_1","Unit_Report_Label_2",
        "u_G_iy","u_ESE_iy","ICAP","sto_c_idhty","c_idhty","su_idhty","g_idhty","reg_up_idhty","reg_dn_idhty",
        "spin_idhty","nonspin_idhty","flex_up_idhty","flex_dn_idhty","curt_idhty","chg_idhty","soc_idhty",
        "hybrid_chg_idhty","hybrid_type", "Fuel_Type", "Fuel_Consumption_MMBtu", "Fuel_Cost", "Inertia","MC","NumDays",
        # the VRE availability shape behind curt_idhty, so reported curtailment can be checked against the input profile
        "vre_shape",
    ]

    market_label_list = [
        "Scenario","year","Stochastic_scenario_ID","day","hour","time","bus_id","Bus_Name","Parent_Bus_Name","Region_Name",
        "load","scarcity_E","scarcity_SPIN","Demand_Reserve","scarcity_NSPIN","scarcity_FU","scarcity_FD",
        "LMP","RCP_RU","RCP_RD","RCP_Spin","RCP_NSpin","RCP_FU","RCP_FD","NumDays","ens_indicator_ndhty",
    ]

    policy_label_list = ["Scenario","year","day","bus_id","slack_CEG_ndy","slack_RPS_ny"]

    demand_response_label_list = [
        "Scenario","year","day","hour","time","PLANT_NAME","bus_idx","bus_name","Region_Name","Parent_Bus_Name","Online_Year",
        "UNITGROUP","UNIT_CATEGORY","UNIT_REPORT_LABEL_1","UNIT_REPORT_LABEL_2","CAP","INTERCON_LIM","Integer_Flag","Daily_DR_Limit_MWh","Num_DR_Segments",
        "Pct_MW_1","Pct_MW_2","Pct_MW_3","Pct_MW_4","Pct_MW_5","Price_1","Price_2","Price_3","Price_4","Price_5",
        "Hybrid_Gen","Hybrid_Gen_CAP","Hybrid_ES","Hybrid_ES_CAP",
        "lfl_lt","lfl_DR_lt","lfl_seg_lt_1","lfl_seg_lt_2","lfl_seg_lt_3","lfl_seg_lt_4","lfl_seg_lt_5",
        "lfl_ind_lt_1","lfl_ind_lt_2","lfl_ind_lt_3","lfl_ind_lt_4","lfl_ind_lt_5",
        "lfl_g_G_LFL_lt","lfl_g_G_Grid_lt","lfl_g_G_ES_lt",
        "lfl_g_ES_LFL_lt","lfl_g_ES_Grid_lt","lfl_chg_Grid_ES_lt","lfl_soc_lt",
    ]

    dispatch_label_list = report_labels(dispatch_label_list)
    market_label_list = report_labels(market_label_list)
    policy_label_list = report_labels(policy_label_list; overrides = Dict("bus_id" => "Policy_Zone_ID"))
    demand_response_label_list = report_labels(demand_response_label_list)

    # -----------------------------
    # Indices
    # -----------------------------
    ids_h = collect(setting["run_H"])
    ids_t = collect(setting["run_T"])

    ids_y = try
        sort!(parse.(Int, collect(get_index(reference, :planning_stages, 0))))
    catch
        sort!(collect(get_index(reference, :planning_stages, 0)))
    end
    local_ids_y = years === nothing ? ids_y : years

    ids_i = try
        sort!(parse.(Int, collect(get_index(reference, :gen_index, 0))))
    catch
        sort!(collect(get_index(reference, :gen_index, 0)))
    end

    ids_n = try
        sort!(parse.(Int, collect(get_index(reference, :bus, 0))))
    catch
        sort!(collect(get_index(reference, :bus, 0)))
    end

    ids_d = try
        sort!(parse.(Int, collect(get_index(reference, :repdays, 0))))
    catch
        sort!(collect(get_index(reference, :repdays, 0)))
    end

    ids_p = try
        sort!(parse.(Int, collect(keys(reference[0][:zone]["policy"]))))
    catch
        sort!(parse.(Int, collect(keys(reference[0][:zone]["policy"]))))
    end

    ids_lfl = try
        sort!(parse.(Int, collect(get_index(reference, :demand, 0))))
    catch
        sort!(collect(get_index(reference, :demand, 0)))
    end

    ids_decomp_group = [1]

    # -----------------------------
    # Scalars
    # -----------------------------
    pu_power_base = setting["Simulation Setting"]["per_unit_base_value"]
    pu_econ_base  = setting["Simulation Setting"]["per_unit_econ_base_value"] / pu_power_base
    current_year  = setting["Planning Design"]["dollar_year_value"]
    discount_rate = setting["Planning Design"]["discount_rate_value"]
    Scenario      = case_name
    reserve_modeling_option = setting["Simulation Configuration"]["operating_reserve_modeling_option"]

    # -----------------------------
    # Precompute reserve zone per bus
    # -----------------------------
    bus_reserve_zone_idx = Dict{Int,Any}()
    for bus_idx in ids_n
        for zone_id in keys(reference[0][:zone]["reserve"])
            if bus_idx in reference[0][:zone]["reserve"][zone_id]["aggregation_info"]["zone_bus_idx"]
                bus_reserve_zone_idx[bus_idx] = zone_id
                break
            end
        end
    end

    # -----------------------------
    # Cache bus metadata
    # -----------------------------
    bus_name = Dict{Int,Any}()
    region_name = Dict{Int,Any}()
    parent_bus_name = Dict{Int,Any}()
    agg_bus_ids = Dict{Int,Any}()          # aggregated BA ids per bus
    agg_bus_load = Dict{Int,Any}()         # original_load vector/dict per bus

    for b in ids_n
        bus_name[b] = parameter(reference, 0, :bus, "bus_i", b)
        region_name[b] = reference[0][:bus][b]["region_config"]["region_name"]
        parent_bus_name[b] = reference[0][:bus][b]["region_config"]["parent_bus_name"]
        agg_bus_ids[b] = reference[0][:bus][b]["aggregation_info"]["aggregated_regions_bus_i"]
        agg_bus_load[b] = reference[0][:bus][b]["aggregation_info"]["original_load_(bus_i, MW)"]
    end

    # -----------------------------
    # Cache generator metadata (depends only on i)
    # -----------------------------
    bus_of_i    = Dict{Int,Int}()
    tech_of_i   = Dict{Int,Int}()
    unitgroup_of_i = Dict{Int,Any}()
    unitcat_of_i   = Dict{Int,Any}()
    unitlabel1_of_i = Dict{Int,Any}()
    unitlabel2_of_i = Dict{Int,Any}()
    plant_of_i     = Dict{Int,Any}()
    hybrid_type_of_i = Dict{Int,Any}()
    cap_of_i_mw    = Dict{Int,Float64}()
    pmax_of_i      = Dict{Int,Float64}()
    pf_of_i        = Dict{Int,Float64}()
    inertiaH_of_i  = Dict{Int,Float64}()
    hr_of_i        = Dict{Int,Float64}()
    fueltype_of_i   = Dict{Int,Any}()
    timeseries_tag_of_i = Dict{Int,Any}()
    fixed_profile_type_of_i = Dict{Int,Any}()

    # cache some flags/cost params used often
    has_commitment_of_i = Dict{Int,Bool}()
    storage_commitment_of_i = Dict{Int,Bool}()
    ptc_flag_of_i = Dict{Int,Bool}()

    # reserve cost bases (either absolute in $/MW or percent*MC later)
    reg_cost_param = Dict{Int,Float64}()
    spin_cost_param = Dict{Int,Float64}()
    nspin_cost_param = Dict{Int,Float64}()
    flex_cost_param = Dict{Int,Float64}()

    # UC costs
    nlc_param = Dict{Int,Float64}()
    suc_param = Dict{Int,Float64}()

    reserve_cost_type = setting["Planning Design"]["reserve_cost_type_flag"]
    # Minimum regulation cost floor ($/MWh). base_reg here is already real $/MWh
    # (= reg_cost_param * MC, with MC = Annual_MC * pu_econ_base), so use the raw $/MWh value.
    min_reg_cost_pu = get(setting["Planning Design"], "min_regulation_cost_value", 0.0)
    dispatch_mode_uc = (setting["Simulation Configuration"]["Dispatch_Mode_in_EXP"] == "Unit Commitment")
    ptc_enabled = (setting["Simulation Configuration"]["PTC_Flag"] == true)

    for i in ids_i
        b = parameter(reference, 0, :gen_index, "bus_idx", i)
        tech = parameter(reference, 0, :gen_index, "genco_tech_id", i)

        bus_of_i[i] = b
        tech_of_i[i] = tech
        unitgroup_of_i[i] = parameter(reference, 0, :gen_index, "UNIT_GROUP", i)

        ucat = parameter(reference, b, :gen_bus, tech, "UNIT_CATEGORY")
        unitcat_of_i[i] = ucat
        unitlabel1_of_i[i] = parameter(reference, b, :gen_bus, tech, "UNIT_REPORT_LABEL_1")
        unitlabel2_of_i[i] = parameter(reference, b, :gen_bus, tech, "UNIT_REPORT_LABEL_2")

        # plant name resolution once
        pname = parameter(reference, b, :gen_bus, tech, "NEW_Gen_UID")
        if haskey(reference[b][:gen_bus][tech], "PLANT_NAME")
            pname = parameter(reference, b, :gen_bus, tech, "PLANT_NAME")
        end
        plant_of_i[i] = pname

        hybrid_type_of_i[i] = parameter(reference, 0, :gen_index, "hybrid_type", i)

        cap_of_i_mw[i] = parameter(reference, b, :gen_bus, tech, "CAP") * pu_power_base
        pmax_of_i[i] = parameter(reference, b, :gen_bus, tech, "PMAX")
        pf_of_i[i] = parameter(reference, b, :gen_bus, tech, "Power_Factor")
        inertiaH_of_i[i] = parameter(reference, b, :gen_bus, tech, "Inertia_Constant")
        hr_of_i[i] = parameter(reference, b, :gen_bus, tech, "HR")
        timeseries_tag_of_i[i] = parameter(reference, b, :gen_bus, tech, "Timeseries_Tag")

        fixed_profile_type_of_i[i] = parameter(reference, b, :gen_bus, tech, "Profile_Type")

        has_commitment_of_i[i] = (parameter(reference, b, :gen_bus, tech, "Commitment") == true)
        storage_commitment_of_i[i] = (parameter(reference, b, :gen_bus, tech, "Storage Commitment") == true)
        ptc_flag_of_i[i] = ptc_enabled && (parameter(reference, b, :gen_bus, tech, "PTC Flag") == true)

        reg_cost_param[i]   = float(parameter(reference, b, :gen_bus, tech, "reg_cost"))
        spin_cost_param[i]  = float(parameter(reference, b, :gen_bus, tech, "spin_cost"))
        nspin_cost_param[i] = float(parameter(reference, b, :gen_bus, tech, "nspin_cost"))
        flex_cost_param[i]  = float(parameter(reference, b, :gen_bus, tech, "flex_cost"))

        nlc_param[i] = float(parameter(reference, b, :gen_bus, tech, "NLC"))
        suc_param[i] = float(parameter(reference, b, :gen_bus, tech, "SUC"))

        fueltype_of_i[i] = parameter(reference, b, :gen_bus, tech, "FUEL")
    end

    reserve_group_of_i = Dict{Int, Int}()
    if reserve_modeling_option == "aggregated"
        for i in ids_i
            z = parse(Int, bus_reserve_zone_idx[bus_of_i[i]])
            unit_group = unitgroup_of_i[i]
            group_idx = reference[0][:reserve_group_reverse_lookup][(unit_group, z)]
            reserve_group_of_i[i] = group_idx
        end
    end

    function get_fixed_profile_shape(i::Int, d::Int, h::Int, t::Int, y::Int)
        type = fixed_profile_type_of_i[i]
        type == "NA" && return 0.0

        b = bus_of_i[i]
        Timeseries_Tag = timeseries_tag_of_i[i]

        if Timeseries_Tag == "LOCAL"
            type_key = if type == "wind_ons"
                "wind_ons_shape"
            elseif type == "wind_ofs"
                "wind_ofs_shape"
            elseif type == "pv"
                "pv_shape"
            elseif type == "hydro"
                "hydro_shape"
            else
                string(type, "_shape")
            end
            return get_vre_zdt_shape(reference[0], y, b, d, h, type_key)
        else
            if type == "wind_ons"
                return parameter(reference, 0, :planning_stages, "repdays", "data", "wind_ons", d, h, t, y)[Timeseries_Tag]
            elseif type == "wind_ofs"
                return parameter(reference, 0, :planning_stages, "repdays", "data", "wind_ofs", d, h, t, y)[Timeseries_Tag]
            elseif type == "csp"
                return parameter(reference, 0, :planning_stages, "repdays", "data", "csp", d, h, t, y)[Timeseries_Tag]
            end
        end

        return 0.0
    end

    # -----------------------------
    # Cache large-load static metadata
    # -----------------------------
    lfl_static = Dict{Int,NamedTuple}()
    for lfl_idx in ids_lfl
        b = parse(Int, parameter(reference, 0, :demand, "bus_idx", lfl_idx))
        CAP = parameter(reference, 0, :demand, "CAP", lfl_idx) * pu_power_base
        lfl_static[lfl_idx] = (
            PLANT_NAME = parameter(reference, 0, :demand, "PLANT_NAME", lfl_idx),
            bus_idx = b,
            bus_name = bus_name[b],
            Region_Name = region_name[b],
            Parent_Bus_Name = parent_bus_name[b],
            Online_Year = parameter(reference, 0, :demand, "Online_Year", lfl_idx),
            UNITGROUP = parameter(reference, 0, :demand, "UNITGROUP", lfl_idx),
            UNIT_CATEGORY = parameter(reference, 0, :demand, "UNIT_CATEGORY", lfl_idx),
            UNIT_REPORT_LABEL_1 = parameter(reference, 0, :demand, "UNIT_REPORT_LABEL_1", lfl_idx),
            UNIT_REPORT_LABEL_2 = parameter(reference, 0, :demand, "UNIT_REPORT_LABEL_2", lfl_idx),
            CAP = CAP,
            INTERCON_LIM = parameter(reference, 0, :demand, "INTERCON_LIM", lfl_idx) * pu_power_base,
            Integer_Flag = parameter(reference, 0, :demand, "Integer_Flag", lfl_idx),
            Daily_DR_Limit_MWh = parameter(reference, 0, :demand, "Daily_DR_Limit_MWh", lfl_idx) * pu_power_base,
            Num_DR_Segments = parameter(reference, 0, :demand, "Num_DR_Segments", lfl_idx),
            Pct_MW_1 = parameter(reference, 0, :demand, "Pct_MW_1", lfl_idx) * CAP,
            Pct_MW_2 = parameter(reference, 0, :demand, "Pct_MW_2", lfl_idx) * CAP,
            Pct_MW_3 = parameter(reference, 0, :demand, "Pct_MW_3", lfl_idx) * CAP,
            Pct_MW_4 = parameter(reference, 0, :demand, "Pct_MW_4", lfl_idx) * CAP,
            Pct_MW_5 = parameter(reference, 0, :demand, "Pct_MW_5", lfl_idx) * CAP,
            Price_1 = parameter(reference, 0, :demand, "Price_1", lfl_idx) * pu_econ_base,
            Price_2 = parameter(reference, 0, :demand, "Price_2", lfl_idx) * pu_econ_base,
            Price_3 = parameter(reference, 0, :demand, "Price_3", lfl_idx) * pu_econ_base,
            Price_4 = parameter(reference, 0, :demand, "Price_4", lfl_idx) * pu_econ_base,
            Price_5 = parameter(reference, 0, :demand, "Price_5", lfl_idx) * pu_econ_base,
            Hybrid_Gen = parameter(reference, 0, :demand, "Hybrid_Gen", lfl_idx),
            Hybrid_Gen_CAP = parameter(reference, 0, :demand, "Hybrid_Gen_CAP", lfl_idx),
            Hybrid_ES = parameter(reference, 0, :demand, "Hybrid_ES", lfl_idx),
            Hybrid_ES_CAP = parameter(reference, 0, :demand, "Hybrid_ES_CAP", lfl_idx),
        )
    end

    # -----------------------------
    # Annual totals as scalars
    # -----------------------------
    metric_names = (
        "Generation","Curtail","Storage_Charge","Reserve_RU","Reserve_RD","Reserve_Spin","Reserve_NSpin","Reserve_FU","Reserve_FD",
        "UnitRevenue_E","UnitRevenue_AS","UnitRevenue_CRED","UnitRevenue","UnitProfit",
        "Noload_Cost","StartUp_Cost","Generation_Cost","Charge_Cost","Fuel_Cost","Fuel_Consumption",
        "Reserve_Reg_Cost","Reserve_Spin_Cost","Reserve_NSpin_Cost","Reserve_Flex_Cost"
    )
    Annual_Gen_Info = Dict{Int,Dict{Int,Dict{String,Float64}}}()
    for i in ids_i
        Annual_Gen_Info[i] = Dict{Int,Dict{String,Float64}}()
        for y in ids_y
            Annual_Gen_Info[i][y] = Dict(m => 0.0 for m in metric_names)
        end
    end
    annual_lock = ReentrantLock()

    # NumDays-weighted scarcity totals per year: the only values the system summary took from
    # __market_EXP.csv, so carrying them here removes a multi-million-row CSV round-trip.
    scarcity_metric_names = ("scarcity_E", "scarcity_SPIN", "scarcity_NSPIN", "scarcity_FU", "scarcity_FD", "Load_MWh")
    Annual_Scarcity_Info = Dict{Int,Dict{String,Float64}}(
        y => Dict(m => 0.0 for m in scarcity_metric_names) for y in ids_y)

    # -----------------------------
    # Open files + write headers once
    # -----------------------------
    # write_csv==false keeps the metric-accumulation path identical by routing all IO to devnull.
    # market_only writes only the market piece (bulky files -> devnull) while still computing annual info.
    write_market = write_csv || market_only
    write_bulky  = write_csv && !market_only
    market_io = write_market ? open(market_file_name, "w") : devnull;  println(market_io, join(market_label_list, ","))
    policy_io = write_bulky ? open(policy_file_name, "w") : devnull;  println(policy_io, join(policy_label_list, ","))
    demand_io = write_bulky ? open(demand_file_name, "w") : devnull;  println(demand_io, join(demand_response_label_list, ","))

    dispatch_ios = Dict{Int,IO}()
    for y in local_ids_y
        io = write_bulky ? open(string(dispatch_file_prefix, y, ".csv"), "w") : devnull
        println(io, join(dispatch_label_list, ","))
        dispatch_ios[y] = io
    end

    market_lock = ReentrantLock()
    policy_lock = ReentrantLock()
    demand_lock = ReentrantLock()
    dispatch_locks = Dict(y => ReentrantLock() for y in local_ids_y)

    try
        # -----------------------------
        # Main loops
        # -----------------------------
        for year in local_ids_y
            stage_length = parameter(reference, 0, :planning_stages, "stage_length", year)
            year_of_year_id = parameter(reference, 0, :planning_stages, "year", year)
            discount_factor = (1 + discount_rate)^-(year_of_year_id - current_year)

            for decomp_group in ids_decomp_group
                grpkey = string(decomp_group)
                sol = result[grpkey]["solution"]

                # repdays for group
                ids_d_new = ids_d

                @threads for dd in eachindex(ids_d_new)
                    day_id = ids_d_new[dd]

                    # thread-local buffers
                    buf_dispatch = IOBuffer()
                    buf_market   = IOBuffer()
                    buf_policy   = IOBuffer()
                    buf_demand   = IOBuffer()

                    # NumDays + Scenario_ID once per day
                    local NumDays::Float64
                    local Scenario_ID::Int
                    if setting["Planning Design"]["multi_round_solution_process_flag"] == true
                        info = parameter(reference, 0, :multi_round_info, "planning_stages", year)[string(day_id)]
                        NumDays = info["NumDays"]
                        Scenario_ID = (setting["Simulation Configuration"]["Stochastic_Expansion_Flag"] == true) ? info["Scenario_ID"] : 0
                    else
                        NumDays = parameter(reference, 0, :repdays, "NumDays", day_id)
                        Scenario_ID = (setting["Simulation Configuration"]["Stochastic_Expansion_Flag"] == true) ? parameter(reference, 0, :repdays, "Scenario_ID", day_id) : 0
                    end
                    scale = NumDays * stage_length

                    # thread-local annual accumulation for this day (no locks)
                    localA = Dict{Int,Dict{String,Float64}}()
                    localS = Dict{Int,Dict{String,Float64}}()
                    @inline function get_localA(i)
                        get!(localA, i) do
                            Dict(m => 0.0 for m in metric_names)
                        end
                    end

                    # -----------------------------
                    # Dispatch (generators)
                    # -----------------------------
                    for hour_id in ids_h, time_id in ids_t
                        for i in ids_i
                            b    = bus_of_i[i]
                            tech = tech_of_i[i]

                            idx_idhty = "($i, $day_id, $hour_id, $time_id, $year)"

                            # expansion
                            u_G_iy = sol["expansion"]["($i, $year)"]["u_G_iy"]
                            CAP    = cap_of_i_mw[i]
                            ICAP   = CAP * u_G_iy

                            disp = sol["dispatch"][idx_idhty]
                            g_idhty       = disp["g_idhty"]
                            reg_up_idhty  = 0.0
                            reg_dn_idhty  = 0.0
                            spin_idhty    = 0.0
                            flex_up_idhty = 0.0
                            flex_dn_idhty = 0.0
                            nonspin_idhty = 0.0

                            if reserve_modeling_option == "individual"
                                reserve_disp = get(sol["reserve"], idx_idhty, Dict{String,Any}())
                                reg_up_idhty  = get(reserve_disp, "reg_up_idhty", 0.0)
                                reg_dn_idhty  = get(reserve_disp, "reg_dn_idhty", 0.0)
                                spin_idhty    = get(reserve_disp, "spin_idhty", 0.0)
                                flex_up_idhty = get(reserve_disp, "flex_up_idhty", 0.0)
                                flex_dn_idhty = get(reserve_disp, "flex_dn_idhty", 0.0)
                                nonspin_idhty = get(reserve_disp, "nonspin_idhty", 0.0)
                            end

                            if reserve_modeling_option == "aggregated"
                                r = reserve_group_of_i[i]
                                group_info = reference[0][:reserve_group_lookup][r]
                                unit_group = group_info["unit_group"]
                                z = group_info["zone_idx"]
                                group_data = reference[0][:zone]["reserve"][string(z)]["aggregation_info"]["gen_technology_groups"][unit_group]
                                group_active_capacity = 0.0
                                for grouped_i in group_data["gen_idx"]
                                    grouped_CAP = cap_of_i_mw[grouped_i]
                                    grouped_PMAX = parameter(reference, bus_of_i[grouped_i], :gen_bus, tech_of_i[grouped_i], "PMAX")
                                    grouped_u_G_iy = sol["expansion"]["($grouped_i, $year)"]["u_G_iy"]
                                    group_active_capacity += grouped_CAP * grouped_PMAX * grouped_u_G_iy
                                end

                                local_PMAX = parameter(reference, b, :gen_bus, tech, "PMAX")
                                reserve_share = 0.0
                                if group_active_capacity > 0.0
                                    reserve_share = (CAP * local_PMAX * u_G_iy) / group_active_capacity
                                end
                                idx_rdhty = "($r, $day_id, $hour_id, $time_id, $year)"
                                reserve_disp = get(sol["reserve"], idx_rdhty, Dict{String,Any}())

                                reg_up_idhty  = get(reserve_disp, "reg_up_idhty", 0.0) * reserve_share
                                reg_dn_idhty  = get(reserve_disp, "reg_dn_idhty", 0.0) * reserve_share
                                spin_idhty    = get(reserve_disp, "spin_idhty", 0.0) * reserve_share
                                flex_up_idhty = get(reserve_disp, "flex_up_idhty", 0.0) * reserve_share
                                flex_dn_idhty = get(reserve_disp, "flex_dn_idhty", 0.0) * reserve_share
                                nonspin_idhty = get(reserve_disp, "nonspin_idhty", 0.0) * reserve_share
                            end

                            curt_idhty = 0.0
                            shape = get_fixed_profile_shape(i, day_id, hour_id, time_id, year)
                            if shape > 0.0
                                availability = shape * u_G_iy * CAP * pmax_of_i[i]
                                curt_idhty = max(availability - g_idhty, 0.0)
                            end

                            c_idhty  = 0.0
                            su_idhty = 0.0
                            if haskey(sol, "commitment") && haskey(sol["commitment"], idx_idhty)
                                com = sol["commitment"][idx_idhty]
                                c_idhty  = com["c_idhty"]
                                su_idhty = com["su_idhty"]
                            end

                            Unit_Category = unitcat_of_i[i]
                            u_ESE_iy = 0.0
                            chg_idhty = 0.0
                            soc_idhty = 0.0
                            sto_c_idhty = 0.0
                            if Unit_Category == "STORAGE"
                                u_ESE_iy = sol["expansion"]["($i, $year)"]["u_ESE_iy"]
                                sto = sol["storage"][idx_idhty]
                                chg_idhty = sto["chg_idhty"]
                                soc_idhty = sto["soc_idhty"]
                                if storage_commitment_of_i[i]
                                    sto_c_idhty = sol["storage_commitment"][idx_idhty]["sto_c_idhty"]
                                end
                            end

                            hybrid_type = hybrid_type_of_i[i]
                            hybrid_chg_idhty = 0.0
                            if hybrid_type == "GEN"
                                hybrid_chg_idhty = sol["hybrid"][idx_idhty]["g_G_ES_idhty"]
                            end

                            MC = parameter(reference, b, :gen_bus, tech, "Annual_MC")[string(year)][day_id] * pu_econ_base

                            inertia = 0.0
                            if haskey(sol, "commitment")
                                inertia = (CAP / pf_of_i[i] * inertiaH_of_i[i]) * c_idhty
                            end
                            
                            # fuel consumption and cost for this hour 
                            FuelType = fueltype_of_i[i]
                            FuelConsumption = g_idhty * hr_of_i[i]
                            fuel_cost = parameter(reference, b, :gen_bus, tech, "Annual_FC")[string(year)][day_id] * pu_econ_base
                            FuelCost = FuelConsumption * fuel_cost

                            # write dispatch row 
                            row = Any[
                                        Scenario, year, stage_calendar_year(setting, year), Scenario_ID, day_id, hour_id, time_id, i, plant_of_i[i], b,
                                        bus_name[b], parent_bus_name[b], region_name[b], tech, unitgroup_of_i[i],
                                        Unit_Category, unitlabel1_of_i[i], unitlabel2_of_i[i], u_G_iy, u_ESE_iy, ICAP, sto_c_idhty,
                                        c_idhty, su_idhty, g_idhty, reg_up_idhty, reg_dn_idhty, spin_idhty,
                                        nonspin_idhty, flex_up_idhty, flex_dn_idhty, curt_idhty, chg_idhty,
                                        soc_idhty, hybrid_chg_idhty, hybrid_type, FuelType, FuelConsumption, FuelCost,
                                        inertia, MC, NumDays, shape
                                    ]
                            println(buf_dispatch, join(csv_escape.(row), ","))
                            # -----------------------------

                            # market prices (for revenue/cost calcs)
                            rz = bus_reserve_zone_idx[b]
                            idx_ndhty = "($b, $day_id, $hour_id, $time_id, $year)"
                            idx_zdhty = "($rz, $day_id, $hour_id, $time_id, $year)"

                            dualn = sol["dual"][idx_ndhty]
                            dualz = sol["dual"][idx_zdhty]

                            LMP = dualn["LMP"] / discount_factor
                            RCP_RU   = dualz["Reg_Up_Price"] / discount_factor
                            RCP_RD   = dualz["Reg_Dn_Price"] / discount_factor
                            RCP_Cont = dualz["Cont_Res_Price"] / discount_factor
                            RCP_NSPIN= dualz["NonSpin_Res_Price"] / discount_factor
                            RCP_FU   = dualz["Flex_Up_Price"] / discount_factor
                            RCP_FD   = dualz["Flex_Dn_Price"] / discount_factor

                            # ---- revenue/cost ----
                            gen_cost = MC * g_idhty

                            gen_credit = 0.0
                            if ptc_flag_of_i[i]
                                base_year = parameter(reference, 0, :planning_stages, "year", year)
                                ptc_tbl = parameter(reference, b, :gen_bus, tech, "PTC")
                                for fy in 1:stage_length
                                    ptc_year = min(2050, base_year + fy - 1)
                                    ptc = ptc_tbl[string(ptc_year)] * 10 * pu_econ_base
                                    gen_credit += ptc * g_idhty
                                end
                            end

                            base_reg  = reg_cost_param[i]
                            base_spin = spin_cost_param[i]
                            base_nsp  = nspin_cost_param[i]
                            base_flex = flex_cost_param[i]

                            if reserve_cost_type == "percentage"
                                base_reg  *= MC
                                base_spin *= MC
                                base_nsp  *= MC
                                base_flex *= MC
                            else
                                base_reg  *= pu_econ_base
                                base_spin *= pu_econ_base
                                base_nsp  *= pu_econ_base
                                base_flex *= pu_econ_base
                            end

                            # Apply the minimum regulation cost floor so reported reg cost
                            # matches what the objective optimized (regulation never free).
                            base_reg = max(base_reg, min_reg_cost_pu)

                            reg_cost   = base_reg  * (reg_up_idhty + reg_dn_idhty)
                            cont_cost  = base_spin * spin_idhty
                            nspin_cost = base_nsp  * nonspin_idhty
                            flex_cost  = base_flex * (flex_up_idhty + flex_dn_idhty)

                            noloadc = 0.0
                            su_cost = 0.0
                            if dispatch_mode_uc && (Unit_Category != "STORAGE") && has_commitment_of_i[i]
                                noloadc = CAP * c_idhty * nlc_param[i] * pu_econ_base
                                su_cost = su_idhty * suc_param[i] * pu_econ_base
                            end

                            Revenue_E = LMP * g_idhty
                            Revenue_AS = RCP_RU*reg_up_idhty + RCP_RD*reg_dn_idhty + RCP_Cont*spin_idhty +
                                        RCP_NSPIN*nonspin_idhty + RCP_FU*flex_up_idhty + RCP_FD*flex_dn_idhty
                            Revenue_CRED = gen_credit
                            Revenue = Revenue_E + Revenue_AS + Revenue_CRED

                            cost = gen_cost + reg_cost + cont_cost + flex_cost + nspin_cost + noloadc + su_cost

                            charge_cost = 0.0
                            if Unit_Category == "STORAGE"
                                cost += LMP * chg_idhty
                                charge_cost += LMP * chg_idhty
                            end
                            
                            

                            # ---- accumulate annual totals 
                            A = get_localA(i)
                            A["Generation"]       += g_idhty * scale
                            A["Curtail"]          += curt_idhty * scale
                            A["Storage_Charge"]   += chg_idhty * scale
                            A["Reserve_RU"]       += reg_up_idhty * scale
                            A["Reserve_RD"]       += reg_dn_idhty * scale
                            A["Reserve_Spin"]     += spin_idhty * scale
                            A["Reserve_NSpin"]    += nonspin_idhty * scale
                            A["Reserve_FU"]       += flex_up_idhty * scale
                            A["Reserve_FD"]       += flex_dn_idhty * scale

                            A["UnitRevenue_E"]    += Revenue_E * scale
                            A["UnitRevenue_AS"]   += Revenue_AS * scale
                            A["UnitRevenue_CRED"] += Revenue_CRED * scale
                            A["UnitRevenue"]      += Revenue * scale

                            A["UnitProfit"]       += (Revenue - cost) * scale
                            A["Noload_Cost"]      += noloadc * scale
                            A["StartUp_Cost"]     += su_cost * scale

                            A["Generation_Cost"]  += gen_cost * scale
                            A["Charge_Cost"]      += charge_cost * scale
                            A["Fuel_Cost"]        += FuelCost * scale
                            A["Fuel_Consumption"] += FuelConsumption * scale

                            A["Reserve_Reg_Cost"]   += base_reg  * (reg_up_idhty + reg_dn_idhty) * scale
                            A["Reserve_Spin_Cost"]  += base_spin * spin_idhty * scale
                            A["Reserve_NSpin_Cost"] += base_nsp  * nonspin_idhty * scale
                            A["Reserve_Flex_Cost"]  += base_flex * (flex_up_idhty + flex_dn_idhty) * scale
                        end

                        # -----------------------------
                        # Market (buses)
                        # -----------------------------
                        emitted_reserve_zones = Set{Any}()   # zonal reserve shortfall is a per-zone quantity: report it once, not on every bus in the zone
                        for b in ids_n
                            rz = bus_reserve_zone_idx[b]

                            # Demand: dot product (cached BA list + load factors)
                            Demand = 0.0
                            region_map = get(reference[0], :profile_data_region_map, nothing)
                            if setting["Planning Design"]["multi_round_solution_process_flag"] == true
                                load_tbl = parameter(reference, 0, :multi_round_info, "planning_stages", year)[string(day_id)]["data"][string(hour_id)]["1"]["load"]
                                for ba_id in agg_bus_ids[b]
                                    Demand += load_tbl[profile_data_region(region_map, "load", ba_id)] * agg_bus_load[b][ba_id] * get_load_growth_factor(reference[0], year, ba_id)
                                end
                            else
                                load_tbl = parameter(reference, 0, :planning_stages, "repdays", year)[string(day_id)]["data"][string(hour_id)]["1"]["load"]
                                for ba_id in agg_bus_ids[b]
                                    Demand += load_tbl[profile_data_region(region_map, "load", ba_id)] * agg_bus_load[b][ba_id]
                                end
                            end
                            Demand *= pu_power_base

                            idx_ndhty = "($b, $day_id, $hour_id, $time_id, $year)"
                            idx_zdhty = "($rz, $day_id, $hour_id, $time_id, $year)"

                            dualn = sol["dual"][idx_ndhty]
                            dualz = sol["dual"][idx_zdhty]

                            LMP = dualn["LMP"] / discount_factor
                            RCP_RU   = dualz["Reg_Up_Price"] / discount_factor
                            RCP_RD   = dualz["Reg_Dn_Price"] / discount_factor
                            RCP_Cont = dualz["Cont_Res_Price"] / discount_factor
                            RCP_NSPIN= dualz["NonSpin_Res_Price"] / discount_factor
                            RCP_FU   = dualz["Flex_Up_Price"] / discount_factor
                            RCP_FD   = dualz["Flex_Dn_Price"] / discount_factor

                            scarcity_E     = get(sol["scarcity"][idx_ndhty], "ens_ndhty", 0.0)
                            if rz in emitted_reserve_zones
                                scarcity_CONT = demand_reserve = scarcity_NSPIN = scarcity_FU = scarcity_FD = 0.0
                            else
                                push!(emitted_reserve_zones, rz)
                                scarcity_CONT  = get(sol["scarcity_reserve"][idx_zdhty], "rns_cont_zdhty", 0.0)
                                demand_reserve = get(sol["scarcity_reserve"][idx_zdhty], "demand_reserve_zdhty", 0.0)
                                scarcity_NSPIN = get(sol["scarcity_reserve"][idx_zdhty], "rns_nonspin_zdhty", 0.0)
                                scarcity_FU    = get(sol["scarcity_reserve"][idx_zdhty], "rns_flex_up_zdhty", 0.0)
                                scarcity_FD    = get(sol["scarcity_reserve"][idx_zdhty], "rns_flex_dn_zdhty", 0.0)
                            end

                            ens_indicator = get(sol["scarcity"][idx_ndhty], "ens_indicator_ndhty", 0)

                            row = Any[
                                Scenario, year, stage_calendar_year(setting, year), Scenario_ID, day_id, hour_id, time_id, b,
                                bus_name[b], parent_bus_name[b], region_name[b],
                                Demand, scarcity_E, scarcity_CONT, demand_reserve, scarcity_NSPIN, scarcity_FU, scarcity_FD,
                                LMP, RCP_RU, RCP_RD, RCP_Cont, RCP_NSPIN, RCP_FU, RCP_FD,
                                NumDays, ens_indicator
                            ]
                            println(buf_market, join(csv_escape.(row), ","))

                            # Same values, same reserve-zone de-duplication, same weighting as the
                            # summary's former CSV pass; accumulated here so the file is not needed.
                            localS_y = get!(localS, year) do
                                Dict(m => 0.0 for m in scarcity_metric_names)
                            end
                            localS_y["scarcity_E"]     += scarcity_E * NumDays
                            localS_y["scarcity_SPIN"]  += scarcity_CONT * NumDays
                            localS_y["scarcity_NSPIN"] += scarcity_NSPIN * NumDays
                            localS_y["scarcity_FU"]    += scarcity_FU * NumDays
                            localS_y["scarcity_FD"]    += scarcity_FD * NumDays
                            localS_y["Load_MWh"]       += Demand * NumDays   # dispatched (representative-day weighted) load
                        end

                        # -----------------------------
                        # Demand response (large load)
                        # -----------------------------
                        if haskey(sol, "demand")
                            demand_solution = sol["demand"]
                            for lfl_idx in ids_lfl
                                st = lfl_static[lfl_idx]
                                idx_ldhty = "($lfl_idx, $day_id, $hour_id, $time_id, $year)"
                                lfl_lt    = demand_solution[idx_ldhty]["lfl_ldhty"] * pu_power_base
                                lfl_DR_lt = demand_solution[idx_ldhty]["lfl_DR_ldhty"] * pu_power_base

                                segs = ntuple(_->0.0, 5)
                                inds = ntuple(_->0, 5)

                                # build vectors without allocating
                                lfl_seg = Vector{Float64}(undef, 5); fill!(lfl_seg, 0.0)
                                lfl_ind = Vector{Float64}(undef, 5); fill!(lfl_ind, 0.0)

                                for seg_id in 1:st.Num_DR_Segments
                                    idx = "($lfl_idx, $seg_id, $day_id, $hour_id, $time_id, $year)"
                                    rec = demand_solution[idx]
                                    lfl_seg[seg_id] = rec["lfl_seg_lsdhty"] * pu_power_base
                                    if seg_id == 1
                                        lfl_ind[seg_id] = 1.0
                                    else
                                        val = rec["lfl_ind_lsdhty"]
                                        lfl_ind[seg_id] = st.Integer_Flag ? round(val) : val
                                    end
                                end

                                lfl_g_G_LFL = 0.0; lfl_g_G_Grid = 0.0; lfl_g_G_ES = 0.0
                                if st.Hybrid_Gen != "NA"
                                    lfl_g_G_LFL  = demand_solution[idx_ldhty]["lfl_g_G_LFL_ldhty"] * pu_power_base
                                    lfl_g_G_Grid = demand_solution[idx_ldhty]["lfl_g_G_Grid_ldhty"] * pu_power_base
                                    if st.Hybrid_ES != "NA"
                                        lfl_g_G_ES = demand_solution[idx_ldhty]["lfl_g_G_ES_ldhty"] * pu_power_base
                                    end
                                end

                                lfl_g_ES_LFL = 0.0; lfl_g_ES_Grid = 0.0; lfl_chg_Grid_ES = 0.0; lfl_soc = 0.0
                                if st.Hybrid_ES != "NA"
                                    lfl_g_ES_LFL    = demand_solution[idx_ldhty]["lfl_g_ES_LFL_ldhty"] * pu_power_base
                                    lfl_g_ES_Grid   = demand_solution[idx_ldhty]["lfl_g_ES_Grid_ldhty"] * pu_power_base
                                    lfl_chg_Grid_ES = demand_solution[idx_ldhty]["lfl_chg_Grid_ES_ldhty"] * pu_power_base
                                    lfl_soc         = demand_solution[idx_ldhty]["lfl_soc_ldhty"] * pu_power_base
                                end

                                row = Any[
                                    Scenario, year, stage_calendar_year(setting, year), day_id, hour_id, time_id,
                                    st.PLANT_NAME, st.bus_idx, st.bus_name, st.Region_Name, st.Parent_Bus_Name, st.Online_Year,
                                    st.UNITGROUP, st.UNIT_CATEGORY, st.UNIT_REPORT_LABEL_1, st.UNIT_REPORT_LABEL_2, st.CAP, st.INTERCON_LIM, st.Integer_Flag,
                                    st.Daily_DR_Limit_MWh, st.Num_DR_Segments,
                                    st.Pct_MW_1, st.Pct_MW_2, st.Pct_MW_3, st.Pct_MW_4, st.Pct_MW_5,
                                    st.Price_1, st.Price_2, st.Price_3, st.Price_4, st.Price_5,
                                    st.Hybrid_Gen, st.Hybrid_Gen_CAP, st.Hybrid_ES, st.Hybrid_ES_CAP,
                                    lfl_lt, lfl_DR_lt,
                                    lfl_seg[1], lfl_seg[2], lfl_seg[3], lfl_seg[4], lfl_seg[5],
                                    lfl_ind[1], lfl_ind[2], lfl_ind[3], lfl_ind[4], lfl_ind[5],
                                    lfl_g_G_LFL, lfl_g_G_Grid, lfl_g_G_ES,
                                    lfl_g_ES_LFL, lfl_g_ES_Grid, lfl_chg_Grid_ES, lfl_soc
                                ]
                                println(buf_demand, join(csv_escape.(row), ","))
                            end
                        end
                    end # hour/time loop

                    # -----------------------------
                    # Policy (once per day)
                    # -----------------------------
                    for z in ids_p
                        slack_CEG = 0.0
                        slack_RPS = 0.0

                        if setting["Simulation Configuration"]["Clean_Energy_Generation_Target_EXP_Flag"] == true
                            if setting["Simulation Configuration"]["Clean_Energy_Generation_Target_EXP_Type"] == "Daygroup"
                                idx_ndy = "($z, $day_id, $year)"
                                slack_CEG = sol["slack"][idx_ndy]["slack_CEG_ndy"] * pu_power_base
                            else
                                idx_ny = "($z, $year)"
                                slack_CEG = sol["slack"][idx_ny]["slack_CEG_ny"] * pu_power_base
                            end
                        end

                        if (setting["Simulation Configuration"]["RPS_Flag"] == true) &&
                        (setting["Simulation Configuration"]["Allow_Alternative_RPS_Compliance_Flag"] == true)
                            idx_ny = "($z, $year)"
                            slack_RPS = sol["slack"][idx_ny]["slack_RPS_ny"] * pu_power_base
                        end

                        println(buf_policy, join(Any[Scenario, year, stage_calendar_year(setting, year), day_id, z, slack_CEG, slack_RPS], ","))
                    end

                    # -----------------------------
                    # Merge annual totals (one lock per day)
                    # -----------------------------
                    lock(annual_lock)
                    for (i, Ai) in localA
                        A = Annual_Gen_Info[i][year]
                        for m in metric_names
                            A[m] += Ai[m]
                        end
                    end
                    for (y_s, Si) in localS
                        S = Annual_Scarcity_Info[y_s]
                        for m in scarcity_metric_names
                            S[m] += Si[m]
                        end
                    end
                    unlock(annual_lock)

                    # -----------------------------
                    # Flush buffers (locks)
                    # -----------------------------
                    lock(dispatch_locks[year]); write(dispatch_ios[year], take!(buf_dispatch)); unlock(dispatch_locks[year])
                    lock(market_lock);          write(market_io, take!(buf_market));          unlock(market_lock)
                    lock(demand_lock);          write(demand_io, take!(buf_demand));          unlock(demand_lock)
                    lock(policy_lock);          write(policy_io, take!(buf_policy));          unlock(policy_lock)
                end # threads day
            end # decomp_group
        end # year

    finally

        if write_market
            close(market_io)
        end
        if write_bulky
            close(policy_io)
            close(demand_io)
            for y in local_ids_y
                close(dispatch_ios[y])
            end
        end
    end

    # store annual info
    try
        reference["0"]["annual_gen_info"] = Annual_Gen_Info
        reference["0"]["annual_scarcity_info"] = Annual_Scarcity_Info
    catch
        reference[0][:annual_gen_info] = Annual_Gen_Info
        reference[0][:annual_scarcity_info] = Annual_Scarcity_Info
    end

    @aleaf_info "[ALEAF LC_GTEP]: - Dispatch solution reporting,"
end


function report_result_scarcity_reserve_GTEP(result::Dict{String, Any}, reference, setting; years=nothing, out_dir=nothing)

    output_path = parameter(reference, 0, :output_path)
    dir = out_dir === nothing ? output_path : out_dir
    out_dir === nothing || mkpath(dir)
    file_name = string(parameter(reference, 0, :case_name), "__reserve_shortfall_EXP.csv")
    file_name = joinpath(dir, file_name)

    # Write Label 
    label_list = ["Scenario", "year", "day", "hour", "time", "zone", "rns_spin_zdhty", "demand_reserve_zdhty", "rns_nonspin_zdhty", "rns_flex_up_zdhty", "rns_flex_dn_zdhty"]

    # Write Outputs
    ids_y = []
    ids_z = []
    ids_d = []
    ids_h = [(h) for (h) in setting["run_H"]]
    ids_t = [(h) for (h) in setting["run_T"]]
    try
        ids_y = [(y) for (y) in sort!(parse.(Int, (collect(get_index(reference, :planning_stages, 0)))))]
        ids_z = [(y) for (y) in sort!(parse.(Int, (collect(keys(reference[0][:zone]["reserve"])))))]
        ids_d = [(y) for (y) in sort!(parse.(Int, (collect(get_index(reference, :repdays, 0)))))]
    catch
        ids_y = [(y) for (y) in sort!(collect(get_index(reference, :planning_stages, 0)))]
        ids_z = [(i) for (i) in sort!(parse.(Int, (collect(keys(reference[0][:zone]["reserve"])))))]
        ids_d = [(i) for (i) in sort!(collect(get_index(reference, :repdays, 0)))]
    end
    local_ids_y = years === nothing ? ids_y : years

    ids_decomp_group = [1]

    Scenario = parameter(reference, 0, :case_name)

    # -----------------------------
    # Write file
    # -----------------------------
    io = open(file_name, "w")
    println(io, join(report_labels(label_list), ","))

    try
        for year in local_ids_y

            for decomp_group in ids_decomp_group

                ids_d_new = []
                try
                    ids_d_new = [(y) for (y) in sort!(parse.(Int, (collect(get_index(reference, :repdays, 0)))))]
                catch
                    ids_d_new = [(y) for (y) in sort!(collect(get_index(reference, :repdays, 0)))]
                end

                for day_id in ids_d_new

                    for zode_id in ids_z
                        for hour_id in ids_h
                            for time_id in ids_t

                                idx = string("(", zode_id, ", ", day_id, ", ", hour_id, ", ", time_id, ", ", year,")")
                                srz = get(result[string(decomp_group)]["solution"]["scarcity_reserve"], string(idx), Dict{String,Any}())
                                rns_cont_zdhty = get(srz, "rns_cont_zdhty", 0.0)
                                demand_reserve_zdhty = get(srz, "demand_reserve_zdhty", 0.0)
                                rns_nonspin_zdhty = get(srz, "rns_nonspin_zdhty", 0.0)
                                rns_flex_up_zdhty = get(srz, "rns_flex_up_zdhty", 0.0)
                                rns_flex_dn_zdhty = get(srz, "rns_flex_dn_zdhty", 0.0)

                                row_id = (year-1)*length(ids_d)*length(ids_z)*length(ids_h)*length(ids_t) + (day_id-1)*length(ids_z)*length(ids_h)*length(ids_t) + (hour_id-1)*length(ids_z)*length(ids_t) + (time_id-1)*length(ids_z) + (zode_id) + 1

                                row = Any[Scenario, year, stage_calendar_year(setting, year), day_id, hour_id, time_id, zode_id, rns_cont_zdhty, demand_reserve_zdhty, rns_nonspin_zdhty, rns_flex_up_zdhty, rns_flex_dn_zdhty]
                                println(io, join(map(_clean_noise, row), ","))
                            end
                        end
                    end
                end
            end
        end
    finally
        close(io)
    end

    @aleaf_info "[ALEAF LC_GTEP]: - Reserve scarcity solution reporting,"
    return nothing
end


function report_result_scarcity_ens_GTEP(result::Dict{String, Any}, reference, setting; years=nothing, out_dir=nothing)

    output_path = parameter(reference, 0, :output_path)
    dir = out_dir === nothing ? output_path : out_dir
    out_dir === nothing || mkpath(dir)
    file_name = string(parameter(reference, 0, :case_name), "__unserved_energy_EXP.csv")
    file_name = joinpath(dir, file_name)

    # Write Label 
    label_list = ["Scenario", "year", "day", "hour", "time", "node", "ens_ndhty"]

    # Write Outputs
    ids_y = []
    ids_n = []
    ids_d = []
    ids_h = [(h) for (h) in setting["run_H"]]
    ids_t = [(h) for (h) in setting["run_T"]]
    try
        ids_y = [(y) for (y) in sort!(parse.(Int, (collect(get_index(reference, :planning_stages, 0)))))]
        ids_n = [(y) for (y) in sort!(parse.(Int, (collect(get_index(reference, :bus, 0)))))]
        ids_d = [(y) for (y) in sort!(parse.(Int, (collect(get_index(reference, :repdays, 0)))))]
    catch
        ids_y = [(y) for (y) in sort!(collect(get_index(reference, :planning_stages, 0)))]
        ids_n = [(i) for (i) in sort!(collect(get_index(reference, :bus, 0)))]
        ids_d = [(i) for (i) in sort!(collect(get_index(reference, :repdays, 0)))]
    end
    local_ids_y = years === nothing ? ids_y : years

    ids_decomp_group = [1]

    Scenario = parameter(reference, 0, :case_name)

    # -----------------------------
    # Write file
    # -----------------------------
    io = open(file_name, "w")
    println(io, join(report_labels(label_list), ","))
    try
        for year in local_ids_y

            for decomp_group in ids_decomp_group

                ids_d_new = []
                try
                    ids_d_new = [(y) for (y) in sort!(parse.(Int, (collect(get_index(reference, :repdays, 0)))))]                        
                catch
                    ids_d_new = [(y) for (y) in sort!(collect(get_index(reference, :repdays, 0)))]
                end

                for day_id in ids_d_new

                    for node_id in ids_n
                        for hour_id in ids_h
                            for time_id in ids_t

                                idx = string("(", node_id, ", ", day_id, ", ", hour_id, ", ", time_id, ", ", year,")")
                                ens_ndhty = result[string(decomp_group)]["solution"]["scarcity"][string(idx)]["ens_ndhty"]
                                                
                                row = Any[Scenario, year, stage_calendar_year(setting, year), day_id, hour_id, time_id, node_id, ens_ndhty]
                                println(io, join(map(_clean_noise, row), ","))
                            end
                        end
                    end
                end
            end
        end
    finally
        close(io)
    end

    @aleaf_info "[ALEAF LC_GTEP]: - Scarcity solution reporting,"
    return nothing
end


function report_result_expansion_line_GTEP(result::Dict{String, Any}, reference, setting; years=nothing, out_dir=nothing)

    pu_power_base = setting["Simulation Setting"]["per_unit_base_value"]
    pu_econ_base = setting["Simulation Setting"]["per_unit_econ_base_value"] / pu_power_base    # per unit economic base

    output_path = parameter(reference, 0, :output_path)
    dir = out_dir === nothing ? output_path : out_dir
    out_dir === nothing || mkpath(dir)
    file_name = string(parameter(reference, 0, :case_name), "__line_expansion_EXP.csv")
    file_name = joinpath(dir, file_name)

    # Write Label
    label_list = ["Scenario", "line_id", "year", "f_bus", "t_bus", "f_bus_name", "t_bus_name", "merged_line_UIDs", "original_rate", "u_new_T_ky", "u_T_ky", "final_rate", "length"]

    # Write Outputs
    ids_y = []
    ids_k = []
    try
        ids_y = [(y) for (y) in sort!(parse.(Int, (collect(get_index(reference, :planning_stages, 0)))))]
        ids_k = [(k) for (k) in sort!(parse.(Int, (collect(get_index(reference, :branch, 0))))) if parameter(reference, 0, :branch, "model_flag", k) == true]
    catch
        ids_y = [(y) for (y) in sort!(collect(get_index(reference, :planning_stages, 0)))]
        ids_k = [(k) for (k) in sort!(collect(get_index(reference, :branch, 0))) if parameter(reference, 0, :branch, "model_flag", k) == true]
    end
    local_ids_y = years === nothing ? ids_y : years

    Scenario = parameter(reference, 0, :case_name)

    # Branch endpoints may still reference original bus_i values in the no-config path.
    function get_bus_label(reference_data, bus_ref)
        try
            return reference_data[0][:bus][bus_ref]["bus_i"]
        catch
        end

        for (_, bus_data) in reference_data[0][:bus]
            if get(bus_data, "bus_i", nothing) == bus_ref
                return bus_data["bus_i"]
            end
        end
        return bus_ref
    end

    # -----------------------------
    # Write file
    # -----------------------------
    io = open(file_name, "w")
    println(io, join(report_labels(label_list), ","))
    try
        for year in local_ids_y

            for line_id in ids_k

                idx = "($line_id, $year)"

                orinal_rate = parameter(reference, 0, :branch, "rate_a", line_id) * pu_power_base
                line_length = get(reference[0][:branch][line_id], "length", get(reference[0][:branch][line_id], "Length", 0.0))
                f_bus = parameter(reference, 0, :branch, "f_bus", line_id)
                t_bus = parameter(reference, 0, :branch, "t_bus", line_id)

                f_bus_name = get_bus_label(reference, f_bus)
                t_bus_name = get_bus_label(reference, t_bus)
                merged_line_uids = ""
                try
                    merged_line_uids = join(string.(parameter(reference, 0, :branch, "merged_line_uids", line_id)), "|")
                catch
                end

                expansion_result = result["1"]["solution"]["expansion_line"][string(idx)]["u_new_T_ky"]
                # cumulative u_T_ky matches the flow limit rate_a*(1+u_T_ky); round-local u_new_T_ky understates final in-service capacity in myopic multi-round runs
                cumulative_result = try
                    result["1"]["solution"]["expansion"][string(idx)]["u_T_ky"]
                catch
                    expansion_result
                end
                final_rate = (1 + cumulative_result) * orinal_rate

                row = Any[Scenario, line_id, year, stage_calendar_year(setting, year), f_bus, t_bus, f_bus_name, t_bus_name, merged_line_uids, orinal_rate, expansion_result, cumulative_result, final_rate, line_length]
                println(io, join(map(_clean_noise, row), ","))

            end
        end
    finally
        close(io)
    end

    @aleaf_info "[ALEAF LC_GTEP]: - Line expansion solution reporting,"
    return nothing
end


function report_result_expansion_gen_GTEP(result::Dict{String, Any}, reference; years=nothing, out_dir=nothing, setting=nothing)

    # -----------------------------
    # Paths
    # -----------------------------
    output_path = parameter(reference, 0, :output_path)
    case_name   = parameter(reference, 0, :case_name)
    dir = out_dir === nothing ? output_path : out_dir
    out_dir === nothing || mkpath(dir)
    file_name   = joinpath(dir, string(case_name, "__gen_expansion_EXP.csv"))

    # -----------------------------
    # Labels
    # -----------------------------
    label_list = report_labels(["Scenario","year","unit_id","PLANT_NAME","UnitGroup","Unit_Category","Bus_Name","Unit_Capacity_MW","u_new_G_iy","u_ret_G_iy","u_G_iy","New_Capacity_MW","Retired_Capacity_MW","Total_Capacity_MW","u_new_ESH_iy","u_ESE_iy"])

    # -----------------------------
    # Indices
    # -----------------------------
    ids_y = try
        sort!(parse.(Int, collect(get_index(reference, :planning_stages, 0))))
    catch
        sort!(collect(get_index(reference, :planning_stages, 0)))
    end

    ids_i = try
        sort!(parse.(Int, collect(get_index(reference, :gen_index, 0))))
    catch
        sort!(collect(get_index(reference, :gen_index, 0)))
    end
    local_ids_y = years === nothing ? ids_y : years

    Scenario = case_name

    # -----------------------------
    # Write file (streaming)
    # -----------------------------
    # per-unit identification (plant, group, category, bus, unit size in MW) so the file can be read on its own
    pu_power_base = setting === nothing ? 1.0 : setting["Simulation Setting"]["per_unit_base_value"]
    unit_meta = Dict{Any,Any}()
    for i in ids_i
        b = parameter(reference, 0, :gen_index, "bus_idx", i)
        tech = parameter(reference, 0, :gen_index, "genco_tech_id", i)
        plant = haskey(reference[b][:gen_bus][tech], "PLANT_NAME") ? parameter(reference, b, :gen_bus, tech, "PLANT_NAME") : parameter(reference, b, :gen_bus, tech, "NEW_Gen_UID")
        unit_meta[i] = (plant = plant,
                        unit_group = parameter(reference, 0, :gen_index, "UNIT_GROUP", i),
                        category = parameter(reference, b, :gen_bus, tech, "UNIT_CATEGORY"),
                        bus_name = parameter(reference, 0, :bus, "bus_i", b),
                        cap_mw = parameter(reference, b, :gen_bus, tech, "CAP") * pu_power_base)
    end

    io = open(file_name, "w")
    println(io, join(label_list, ","))
    try
        exp = result["1"]["solution"]["expansion"]

        for year in local_ids_y
            for unit_id in ids_i
                idx = "($unit_id, $year)"

                rec = exp[idx]
                u_new_G_iy = rec["u_new_G_iy"]
                u_ret_G_iy = rec["u_ret_G_iy"]
                u_G_iy     = rec["u_G_iy"]

                u_new_ESH_iy = get(rec, "u_new_ESH_iy", 0.0)
                u_ESE_iy     = get(rec, "u_ESE_iy",     0.0)

                meta = unit_meta[unit_id]
                cap_mw = meta.cap_mw
                row = Any[Scenario, year, parameter(reference, 0, :planning_stages, "year", year), unit_id,
                          meta.plant, meta.unit_group, meta.category, meta.bus_name, cap_mw,
                          u_new_G_iy, u_ret_G_iy, u_G_iy, u_new_G_iy * cap_mw, u_ret_G_iy * cap_mw, u_G_iy * cap_mw,
                          u_new_ESH_iy, u_ESE_iy]
                println(io, join(map(_clean_noise, row), ","))
            end
        end
    finally
        close(io)
    end

    @aleaf_info "[ALEAF LC_GTEP]: - Generation expansion solution reporting,"
    return nothing
end


function report_result_repdays_OP_GTEP(reference, setting)
    
    output_path = parameter(reference["1"], 0, :output_path)
    file_name = string(parameter(reference["1"], 0, :case_name), "__representative_days_OP.csv")
    file_name = joinpath(output_path, file_name)

    # Write Label 
    label_list = report_labels(["Scenario", "year", "RepDay_id", "NumDays", "Daygroup_ID", "Numdays_Daygroup", "Reference_Day"])

    # Write Outputs
    ids_y = []
    ids_d = []
    try
        ids_y = [(y) for (y) in sort!(parse.(Int, (collect(get_index(reference["1"], :planning_stages, 0)))))]
        ids_d = [(y) for (y) in sort!(parse.(Int, (collect(get_index(reference["1"], :repdays, 0)))))]
    catch
        ids_y = [(y) for (y) in sort!(collect(get_index(reference["1"], :planning_stages, 0)))]
        ids_d = [(i) for (i) in sort!(collect(get_index(reference["1"], :repdays, 0)))]
    end

    Scenario = parameter(reference["1"], 0, :case_name)

    output_list = Array{Any}(undef, length(ids_y)*length(ids_d)+1, length(label_list))
    output_list[1,:] = label_list
    for year in ids_y

        reference_of_this_year = reference[string(year)]

        ids_d_new = []
        try
            ids_d_new = [(y) for (y) in sort!(parse.(Int, (collect(get_index(reference_of_this_year, :repdays, 0)))))]                        
        catch
            ids_d_new = [(y) for (y) in sort!(collect(get_index(reference_of_this_year, :repdays, 0)))]
        end

        for day_id in ids_d_new

            NumDays = parameter(reference_of_this_year, 0, :repdays, "NumDays", day_id)
            Daygroup_ID = parameter(reference_of_this_year, 0, :repdays, "Day_Group_ID", day_id)
            Numdays_Daygroup = parameter(reference_of_this_year, 0, :repdays, "NumDays_Group", day_id)
            Reference_Day = parameter(reference_of_this_year, 0, :repdays, "Day", day_id)
        
            row_id = (year-1) * length(ids_d) + day_id + 1
            
            output_list[row_id, :] = [Scenario, year, stage_calendar_year(setting, year), day_id, NumDays, Daygroup_ID, Numdays_Daygroup, Reference_Day]

        end
    end
    writedlm(file_name, output_list, ",")

    @aleaf_info "[ALEAF LC_GTEP]: - Repday selection reporting,"
end


function report_result_repdays_GTEP(reference, setting; years=nothing, out_dir=nothing)

    output_path = parameter(reference, 0, :output_path)
    dir = out_dir === nothing ? output_path : out_dir
    out_dir === nothing || mkpath(dir)
    file_name = string(parameter(reference, 0, :case_name), "__representative_days_EXP.csv")
    file_name = joinpath(dir, file_name)

    # Write Label 
    label_list = report_labels(["Scenario", "year", "RepDay_id", "Stochastic_scenario_ID", "NumDays", "Daygroup_ID", "Numdays_Daygroup", "Reference_Day"])

    Scenario = parameter(reference, 0, :case_name)

    # Write Outputs
    ids_y = []
    ids_d = []
    try
        ids_y = [(y) for (y) in sort!(parse.(Int, (collect(get_index(reference, :planning_stages, 0)))))]
        ids_d = [(y) for (y) in sort!(parse.(Int, (collect(get_index(reference, :repdays, 0)))))]
    catch
        ids_y = [(y) for (y) in sort!(collect(get_index(reference, :planning_stages, 0)))]
        ids_d = [(i) for (i) in sort!(collect(get_index(reference, :repdays, 0)))]
    end
    local_ids_y = years === nothing ? ids_y : years

    ids_decomp_group = [1]

    output_list = Array{Any}(undef, length(ids_y)*length(ids_d)+1, length(label_list))
    output_list[1,:] = label_list
    for year in local_ids_y

        for decomp_group in ids_decomp_group
        
            ids_d_new = []
            try
                ids_d_new = [(y) for (y) in sort!(parse.(Int, (collect(get_index(reference, :repdays, 0)))))]                        
            catch
                ids_d_new = [(y) for (y) in sort!(collect(get_index(reference, :repdays, 0)))]
            end

            for day_id in ids_d_new

                if setting["Planning Design"]["multi_round_solution_process_flag"] == true
                    if setting["Simulation Configuration"]["Stochastic_Expansion_Flag"] == true
                        Scenario_ID = parameter(reference, 0, :multi_round_info, "planning_stages", year)[string(day_id)]["Scenario_ID"]
                    else
                        Scenario_ID = 0
                    end
                    NumDays = parameter(reference, 0, :multi_round_info, "planning_stages", year)[string(day_id)]["NumDays"]
                    Daygroup_ID = parameter(reference, 0, :multi_round_info, "planning_stages", year)[string(day_id)]["Day_Group_ID"]
                    Numdays_Daygroup = parameter(reference, 0, :multi_round_info, "planning_stages", year)[string(day_id)]["NumDays_Group"]
                    Reference_Day = parameter(reference, 0, :multi_round_info, "planning_stages", year)[string(day_id)]["Day"]
                
                    row_id = (year-1) * length(ids_d) + day_id + 1
                    
                    output_list[row_id, :] = [Scenario, year, stage_calendar_year(setting, year), day_id, Scenario_ID, NumDays, Daygroup_ID, Numdays_Daygroup, Reference_Day]
                
                else
                    if setting["Simulation Configuration"]["Stochastic_Expansion_Flag"] == true
                        Scenario_ID = parameter(reference, 0, :repdays, "Scenario_ID", day_id)
                    else
                        Scenario_ID = 0
                    end
                    NumDays = parameter(reference, 0, :repdays, "NumDays", day_id)
                    Daygroup_ID = parameter(reference, 0, :repdays, "Day_Group_ID", day_id)
                    Numdays_Daygroup = parameter(reference, 0, :repdays, "NumDays_Group", day_id)
                    Reference_Day = parameter(reference, 0, :repdays, "Day", day_id)
                
                    row_id = (year-1) * length(ids_d) + day_id + 1
                    
                    output_list[row_id, :] = [Scenario, year, stage_calendar_year(setting, year), day_id, Scenario_ID, NumDays, Daygroup_ID, Numdays_Daygroup, Reference_Day]
                end

            end
        end
    end
    writedlm(file_name, output_list, ",")

    @aleaf_info "[ALEAF LC_GTEP]: - Repday selection reporting,"
end


function build_solve_LC_GTEP_operational_model(ALEAF_setting::Dict{String,<:Any}, result_LC_GTEP_expansion::Dict, case_id::Int; recorded_investment_decisions::Dict{String,Any}=Dict{String,Any}())
    
    # Check expansion result
    expansion_fix = false
    if length(result_LC_GTEP_expansion) > 0
        expansion_fix = true
    end

    # initialize
    OP_data_record_for_reporting = Dict{String, Any}()
    OP_data_record_for_reporting["daily_solutions"] = Dict{String, Any}()
    OP_data_record_for_reporting["network_data"] = Dict{String, Any}()
    OP_data_record_for_reporting["system_reference"] = Dict{String, Any}()
    OP_data_record_for_reporting["setting"] = Dict{String, Any}()
    
    ids_y = [y for y in 1:ALEAF_setting["Planning Design"]["num_stages_value"]] 

    for y in ids_y

        # generate network data
        rid_y = (length(recorded_investment_decisions) > 0 && haskey(recorded_investment_decisions, string(y))) ? recorded_investment_decisions[string(y)] : Dict{String,Any}()
        # Pick repdays on a coarse selection network first (bounded pre-processing), then build the run
        # network with that assignment. Falls back to in-line selection if coarse yields nothing.
        pre_groups, pre_days = (Dict{String,Any}(), Dict{String,Any}())
        coarse = select_repdays_coarse(ALEAF_setting, case_id, "operation"; year_id=y, recorded_investment_decisions=rid_y)
        coarse !== nothing && ((pre_groups, pre_days) = coarse)
        if length(rid_y) > 0
            network_data = generate_networkdata_LC_GTEP(ALEAF_setting, case_id, "operation"; recorded_investment_decisions=rid_y, year_id=y, precomputed_repday_groups=pre_groups, precomputed_repdays=pre_days)
        else
            network_data = generate_networkdata_LC_GTEP(ALEAF_setting, case_id, "operation"; year_id=y, precomputed_repday_groups=pre_groups, precomputed_repdays=pre_days)
        end

        @aleaf_info "[ALEAF LC_GTEP Operational Model]: Generate network data."

        # grouping days
        total_day_groups = ALEAF_setting["Simulation Configuration"][string(case_id)]["NDAY_Groups_OP"]

        # build common ALEAF model instance structure
        ALEAF_model_instance = build_model_structure_GTEP(ALEAF_setting, network_data, Abstract_LC_GTEP_Model)
        @aleaf_info "[ALEAF LC_GTEP Operational Model]: Build common ALEAF model instance structure."

        # Add reference data. add_ref transfers the heavy hourly payload into ref[:nw][0] (default) instead
        # of deep-copying it, so the serial nodal master never doubles the largest structure (OOM fix).
        add_ref_LC_GTEP_model!(ALEAF_model_instance, ALEAF_setting; recorded_investment_decisions)
        @aleaf_info "[ALEAF LC_GTEP Operational Model]: Add reference data."

        # Disable unnecessary flags
        ALEAF_model_instance.setting["Planning Design"]["multi_round_solution_process_flag"] = false
        ALEAF_model_instance.setting["Simulation Configuration"]["transmission_expansion_flag"] = false

        # Free the redundant heavy hourly payload retained in network_data. `add_common_ref_GTEP!`
        # already deepcopied everything into the independent `ref[:nw][0]` (used by every solve + the
        # OP reporter, which reads the system reference, not this dict). The per-day-group rebuild below
        # only needs the light skeleton (bus / gen_technology / repday_groups), and the only other
        # reader of the stored network_data just counts `repdays`. So drop the per-repday
        # ["data"] tensors and the raw full-year time_series_data DataFrames — at nodal scale keeping
        # them alongside ref[:nw][0] doubled the largest structure and caused the master OOM. Guard with
        # the same JSON-export flag add_ref uses so full-reference JSON export is unaffected.
        keep_full_data = get(ALEAF_setting["Simulation Setting"], "export_model_reference_json_operation_flag", false) == true ||
                         get(ALEAF_setting["Simulation Setting"], "export_model_reference_json_expansion_flag", false) == true
        if !keep_full_data
            if haskey(network_data, "repdays")
                for (_, rd) in network_data["repdays"]
                    isa(rd, Dict) && haskey(rd, "data") && delete!(rd, "data")
                end
            end
            if haskey(network_data, "planning_stages")
                for (_, stage) in network_data["planning_stages"]
                    (isa(stage, Dict) && haskey(stage, "repdays")) || continue
                    for (_, rd) in stage["repdays"]
                        isa(rd, Dict) && haskey(rd, "data") && delete!(rd, "data")
                    end
                end
            end
            haskey(network_data, "time_series_data") && delete!(network_data, "time_series_data")
        end

        OP_data_record_for_reporting["daily_solutions"][string(y)] = Dict{String, Any}()
        # store the (now light) network_data skeleton; add_ref may have rebound am.data to a scalar-only
        # dict, so this remains the retained copy. Heavy hourly data lives only in ref[:nw][0] now.
        OP_data_record_for_reporting["network_data"][string(y)] = network_data
        OP_data_record_for_reporting["system_reference"][string(y)] = ALEAF_model_instance.ref[:nw] 
        OP_data_record_for_reporting["setting"][string(y)] = ALEAF_model_instance.setting

        # Solve Model
        num_solved_group = 0
        for day_group_id in 1:total_day_groups
            
            OP_data_record_for_reporting["daily_solutions"][string(y)][string(day_group_id)] = run_operation_model_for_day_group(y, day_group_id, ALEAF_setting, network_data, ALEAF_model_instance.ref, ALEAF_model_instance.setting, expansion_fix, result_LC_GTEP_expansion)
        
            num_solved_group += 1
            if mod(num_solved_group,5) == 0
                progress = round(num_solved_group/total_day_groups, digits=2) * 100 
                @aleaf_info "[ALEAF LC_GTEP Operational Model]: === Solved $num_solved_group day groups out of $total_day_groups day groups in year $y. ($progress % done - year $y) === "
            end

        end
    end

    return OP_data_record_for_reporting
end


function run_operation_model_for_day_group(y, day_group_id, ALEAF_setting, network_data, aleaf_model_ref, aleaf_model_setting, expansion_fix, result_LC_GTEP_expansion)
            
    start_day = network_data["repday_groups"][string(day_group_id)]["Start_Day_Id"]
    end_day = network_data["repday_groups"][string(day_group_id)]["End_Day_Id"]
    
    # build common ALEAF model instance structure
    start_time = time()
    _detailed_log = lowercase(string(get(ALEAF_setting["Simulation Setting"], "logging_level_value", "simple"))) == "detailed"
    _detailed_log && @aleaf_info "[OP build] === day_group=$day_group_id (days $start_day..$end_day): building model structure ==="   # detailed build timing
    ALEAF_model_instance = build_model_structure_GTEP(ALEAF_setting, network_data, Abstract_LC_GTEP_Model)
    _detailed_log && @aleaf_info "[OP build] day_group=$day_group_id: model structure built ($(round(time()-start_time, digits=1))s)"   # detailed build timing

    # Add reference data and setting
    ALEAF_model_instance.ref = aleaf_model_ref
    ALEAF_model_instance.setting= aleaf_model_setting

    # build LCO GEP model instance
    build_LCO_GTEP_operational_mode_instance!(ALEAF_model_instance, result_LC_GTEP_expansion, day_group_id, y; expansion_fix)
    _detailed_log && @aleaf_info "[OP build] day_group=$day_group_id: full model instance ready ($(round(time()-start_time, digits=1))s); handing to solver"   # detailed build timing

    # Define JuMP_model, solution_list, and solver setting
    JuMP_model = ALEAF_model_instance.model[:nw][ALEAF_model_instance.cnw][day_group_id]
    solution_list = ALEAF_model_instance.sol[:nw][ALEAF_model_instance.cnw][day_group_id]
    solver_setting = ALEAF_model_instance.setting["Solver Setting"]
    iteration = 1

    # Export ALEAF_model_instance model
    start_time_2 = time()
    if ALEAF_setting["Simulation Setting"]["export_model_lp_operation_flag"] == true
        output_path = parameter(ALEAF_model_instance, 0, :output_path)
        file_name = ALEAF_setting["Simulation Setting"]["model_lp_file_name_value"]
        file_name = string("Operation_", file_name)
        JuMP.write_to_file(JuMP_model, joinpath(output_path, file_name))
        @aleaf_info "[ALEAF LC_GTEP Operational Model]: Export ALEAF_model_instance model (lp format)."
    end
    preparation_time = round(time() - start_time, digits=2)

    # Solve ALEAF LC_GTEP Operational Model
    start_time = time()
    result_LC_GTEP_operation = solve_model_GTEP!(JuMP_model, day_group_id, solution_list, solver_setting; iteration, PH_flag=false)
    solution_time = round(time() - start_time, digits=2)

    # Get Dual Values
    start_time = time()
    dual_status = JuMP.dual_status(JuMP_model)
    if dual_status == MOI.NO_SOLUTION 
        # fix discrete variables and resolve
        JuMP.fix_discrete_variables(JuMP_model)    # fix discrete variables
        JuMP.optimize!(JuMP_model)
    end
    # collect dual values
    collect_dual_LC_GTEP_operation!(JuMP_model, ALEAF_model_instance, result_LC_GTEP_operation, day_group_id, y)
    dual_collection_time = round(time() - start_time, digits=2)

    # Calculate total time
    total_time = preparation_time + solution_time + dual_collection_time

    @aleaf_info "[ALEAF LC_GTEP Operational Model]: Solved Day Group id: $day_group_id, Start day: $start_day, End day: $end_day, Prep (sec): $preparation_time, Solution (sec): $solution_time, Dual (sec): $dual_collection_time, Total (sec): $total_time"

    # Clean up to save memory
    ALEAF_model_instance = 0.0 
    # ---------------------------------------------------------------------------------

    return result_LC_GTEP_operation
end


# Build the small "day assignment" payload shipped to distributed OP workers from the master's full
# network_data: the full repday_groups (Start/End day, Day_Idx_List, Probability, ... — a few KB) plus
# per-repday METADATA with the heavy hourly ["data"] stripped out (workers rebuild that for their own
# day from the raw timeseries). Index keys are preserved verbatim so worker solutions align with the
# master reporting reference.
function extract_repday_assignment_LC_GTEP(network_data::Dict{String,<:Any})
    repday_groups = deepcopy(network_data["repday_groups"])
    repdays_meta = Dict{String,Any}()
    for (k, v) in network_data["repdays"]
        m = Dict{String,Any}()
        for (kk, vv) in v
            kk == "data" && continue
            m[kk] = deepcopy(vv)
        end
        repdays_meta[k] = m
    end
    return repday_groups, repdays_meta
end


function build_solve_LC_GTEP_operational_model_distributed(ALEAF_setting::Dict{String,<:Any}, result_LC_GTEP_expansion::Dict, case_id::Int; recorded_investment_decisions::Dict{String,Any}=Dict{String,Any}())

    # Check expansion result
    expansion_fix = false
    if length(result_LC_GTEP_expansion) > 0
        expansion_fix = true
    end

    # initialize. Workers do their own per-day-group reporting and return only small aggregate
    # payloads; the master keeps the reduced aggregates (no per-hour solutions, no network_data).
    OP_data_record_for_reporting = Dict{String, Any}()
    OP_data_record_for_reporting["system_reference"] = Dict{String, Any}()
    OP_data_record_for_reporting["setting"] = Dict{String, Any}()
    OP_data_record_for_reporting["annual_gen_info_OP"] = Dict{Any, Any}()   # [gen_idx][year_id] => metrics
    OP_data_record_for_reporting["scarcity_totals"] = Dict{Any, Any}()      # [year_id] => scarcity metrics
    OP_data_record_for_reporting["objectives"] = Dict{String, Any}()        # [year] => Dict(dg => objective)
    OP_data_record_for_reporting["build_decisions"] = Dict{String, Any}()   # [year] => Dict("expansion","slack")
    OP_data_record_for_reporting["pu_applied"] = true

    # initialize
    ids_y = [y for y in 1:ALEAF_setting["Planning Design"]["num_stages_value"]]

    function solve_multi_years_in_parallel(y, ALEAF_setting, case_id, result_LC_GTEP_expansion, expansion_fix, precomputed_repday_groups, precomputed_repdays)

        # collect the small per-day-group payloads
        payloads = Dict{String, Any}()

        # grouping days
        total_day_groups = [i for i in 1:ALEAF_setting["Simulation Configuration"][string(case_id)]["NDAY_Groups_OP"]]
        np = nprocs()

        # The day-group-invariant build decisions are carried back by exactly one worker: the one
        # solving the smallest day-group id in this run (always present, = 1 for a full run).
        capture_build_decisions_dg = minimum(total_day_groups)

        # Fault-tolerant work queue: a worker death (ProcessExitedException/IOError — e.g. node fault or
        # OOM) requeues its day-group onto another live worker instead of aborting the whole run. A
        # thread-safe pending stack + counters are guarded by qlock; `remaining` only drops on success
        # or permanent failure, so tasks keep draining until every day-group is done or capped.
        pending = collect(Int, total_day_groups)
        attempts = Dict{Int,Int}(dg => 0 for dg in total_day_groups)
        failed = Int[]
        remaining = length(total_day_groups)
        max_attempts = 3
        qlock = ReentrantLock()

        @sync begin
            for p in workers()
                if p != myid() || np == 1
                    @async begin
                        while true
                            dg = 0
                            done = false
                            lock(qlock) do
                                if remaining <= 0
                                    done = true
                                elseif !isempty(pending)
                                    dg = pop!(pending)
                                end
                            end
                            done && break
                            if dg == 0
                                sleep(0.2); continue    # nothing free right now; items in flight elsewhere
                            end
                            try
                                # Ship only the small day assignment; the worker builds its own single-day
                                # network data, writes its per-dg output files, and returns a small payload.
                                res = remotecall_fetch(build_solve_GTEP_OP_subproblem, p, ALEAF_setting, case_id, y, dg, result_LC_GTEP_expansion, expansion_fix, precomputed_repday_groups, precomputed_repdays; recorded_investment_decisions=recorded_investment_decisions, capture_build_decisions_dg=capture_build_decisions_dg)
                                lock(qlock) do
                                    payloads[string(dg)] = res
                                    remaining -= 1
                                end
                            catch e
                                if e isa Distributed.ProcessExitedException || e isa Base.IOError
                                    resched = false
                                    lock(qlock) do
                                        attempts[dg] += 1
                                        if attempts[dg] < max_attempts
                                            push!(pending, dg); resched = true
                                        else
                                            push!(failed, dg); remaining -= 1
                                        end
                                    end
                                    @aleaf_warn "[ALEAF LC_GTEP Operational Model]: worker $p died on day group $dg ($(typeof(e))); $(resched ? "requeued (attempt $(attempts[dg]))" : "gave up after $max_attempts attempts"). Retiring this worker."
                                    break   # this worker is gone; end its task and let others drain the queue
                                else
                                    rethrow(e)   # real code error on a live worker — surface it immediately
                                end
                            end
                        end
                    end
                end
            end
        end

        # Completeness guard: catches both capped failures AND the case where every worker died leaving
        # requeued items with no one to run them. Never return partial results silently.
        undone = [dg for dg in total_day_groups if !haskey(payloads, string(dg))]
        if !isempty(undone)
            terminate_with_error(; msg="[ALEAF LC_GTEP Operational Model]: day group(s) $(sort(undone)) not completed (failed after $max_attempts attempts or all workers died); see per-worker logs")
        end

        return payloads
    end

    # Light reporting-only master reference: the master never solves the OP instance and its reporters
    # (tech/system summary, repday selection) read only repday/planning_stages metadata + per-bus/per-policy
    # annual vre_aggregated_data + the injected annual_gen_info_OP. So skip the two heavy per-repday builds on
    # the master: the hourly repday ["data"] (build_repday_hourly_data=false) and vre_aggregated_data_zdt +
    # :hydro_budget (add_ref light=true). Workers still build the FULL reference (they solve + report per day).
    # Keep the full reference only when a model-reference JSON export is requested.
    keep_full_data = get(ALEAF_setting["Simulation Setting"], "export_model_reference_json_operation_flag", false) == true ||
                     get(ALEAF_setting["Simulation Setting"], "export_model_reference_json_expansion_flag", false) == true
    light_ref = !keep_full_data

    for y in ids_y

        network_data = Dict{String, Any}()
        rid_y = (length(recorded_investment_decisions) > 0 && haskey(recorded_investment_decisions, string(y))) ? recorded_investment_decisions[string(y)] : Dict{String,Any}()

        # Pick repdays on a coarse selection network first, then feed the assignment into BOTH the master
        # reference build and the workers (identical to the master's own selection, just far cheaper).
        pre_groups, pre_days = (Dict{String,Any}(), Dict{String,Any}())
        coarse = select_repdays_coarse(ALEAF_setting, case_id, "operation"; year_id=y, recorded_investment_decisions=rid_y)
        coarse !== nothing && ((pre_groups, pre_days) = coarse)

        # Light master reference (see above): skip building the per-repday hourly data when light_ref is on.
        if length(rid_y) > 0
            network_data = generate_networkdata_LC_GTEP(ALEAF_setting, case_id, "operation"; recorded_investment_decisions=rid_y, year_id=y, build_repday_hourly_data=!light_ref, precomputed_repday_groups=pre_groups, precomputed_repdays=pre_days)
        else
            network_data = generate_networkdata_LC_GTEP(ALEAF_setting, case_id, "operation"; year_id=y, build_repday_hourly_data=!light_ref, precomputed_repday_groups=pre_groups, precomputed_repdays=pre_days)
        end
        @aleaf_info "[ALEAF LC_GTEP Operational Model]: Generate network data."

        # Dispatch the day-groups FIRST so workers start solving immediately. The small day assignment
        # is extracted from the master's network_data; each worker rebuilds ONLY its own day-group from
        # it, reusing the exact repday index keys so worker output files line up with the reference.
        precomputed_repday_groups, precomputed_repdays = extract_repday_assignment_LC_GTEP(network_data)
        payloads_y = solve_multi_years_in_parallel(y, ALEAF_setting, case_id, result_LC_GTEP_expansion, expansion_fix, precomputed_repday_groups, precomputed_repdays)

        # -----------------------------------------------------------------------------
        # Reduce the per-day-group payloads for this year
        # -----------------------------------------------------------------------------
        objectives_y = Dict{String, Any}()
        for (dg, pl) in payloads_y
            objectives_y[dg] = pl["objective"]

            # element-wise SUM annual_gen_info across day-groups (each is per gen_idx, per year)
            agi = pl["annual_gen_info_contrib"]
            if agi !== nothing
                for (gen_idx, per_year) in agi
                    dst_g = get!(OP_data_record_for_reporting["annual_gen_info_OP"], gen_idx, Dict{Any,Any}())
                    for (yr, metrics) in per_year
                        dst_m = get!(dst_g, yr, Dict{String,Float64}())
                        for (m, v) in metrics
                            dst_m[m] = get(dst_m, m, 0.0) + v
                        end
                    end
                end
            end

            # sum per-year scarcity totals across day-groups
            asi = pl["scarcity_totals"]
            if asi !== nothing
                for (yr, metrics) in asi
                    dst_s = get!(OP_data_record_for_reporting["scarcity_totals"], yr, Dict{String,Float64}())
                    for (m, v) in metrics
                        dst_s[m] = get(dst_s, m, 0.0) + v
                    end
                end
            end

            # day-group-invariant build decisions captured on dg 1 only
            if haskey(pl, "build_decisions")
                OP_data_record_for_reporting["build_decisions"][string(y)] = pl["build_decisions"]
            end
        end
        OP_data_record_for_reporting["objectives"][string(y)] = objectives_y

        # Build the reporting-only reference AFTER dispatch (off the workers' critical path).
        # build_jump_models=false skips the unused per-repday JuMP models: this instance is never
        # solved; reporting reads only .ref[:nw] / .setting (network_data no longer stored).
        # light=light_ref skips the heavy zdt aggregation + water-budget (see the note at the top of the loop).
        ALEAF_model_instance = build_model_structure_GTEP(ALEAF_setting, network_data, Abstract_LC_GTEP_Model; build_jump_models=false)
        add_ref_LC_GTEP_model!(ALEAF_model_instance, ALEAF_setting; recorded_investment_decisions=rid_y, light=light_ref)

        OP_data_record_for_reporting["system_reference"][string(y)] = ALEAF_model_instance.ref[:nw]
        OP_data_record_for_reporting["setting"][string(y)] = ALEAF_model_instance.setting

    end

    return OP_data_record_for_reporting
end


function build_solve_GTEP_OP_subproblem(ALEAF_setting, case_id, y, day_group_id, result_LC_GTEP_expansion, expansion_fix, precomputed_repday_groups, precomputed_repdays; recorded_investment_decisions::Dict{String,Any}=Dict{String,Any}(), capture_build_decisions_dg::Int=1)

    # Each worker builds the network data for ONLY its assigned day-group. Topology is cached
    # in-process (get_network_data_cached!), scenario reduction is skipped (the master's assignment is
    # reused), and only this day's repday + hourly data are materialized — so the worker footprint is
    # ~1/NDAY_Groups of the full-horizon build. Guard against an absent year key for multi-round runs.
    rid_y = (length(recorded_investment_decisions) > 0 && haskey(recorded_investment_decisions, string(y))) ? recorded_investment_decisions[string(y)] : Dict{String,Any}()
    network_data = generate_networkdata_LC_GTEP(ALEAF_setting, case_id, "operation"; year_id=y, recorded_investment_decisions=rid_y, target_day_group_id=day_group_id, precomputed_repday_groups=precomputed_repday_groups, precomputed_repdays=precomputed_repdays)

    start_day = network_data["repday_groups"][string(day_group_id)]["Start_Day_Id"]
    end_day = network_data["repday_groups"][string(day_group_id)]["End_Day_Id"]

    # build common ALEAF model instance structure
    ALEAF_model_instance = build_model_structure_GTEP(ALEAF_setting, network_data, Abstract_LC_GTEP_Model)
    ALEAF_model_instance_fixed = build_model_structure_GTEP(ALEAF_setting, network_data, Abstract_LC_GTEP_Model)

    # Add reference data (rid_y computed above); the per-day network was built from the master's
    # precomputed assignment, and add_ref applies VRE from recorded expansion decisions for this year.
    # transfer_heavy_data=false: two instances add_ref the SAME shared network_data, and the second
    # (fixed) instance also needs the heavy hourly payload in its nw[0]; keep the full deepcopy for both.
    # Worker footprint is ~1/NDAY_Groups so the deepcopy is affordable here.
    add_ref_LC_GTEP_model!(ALEAF_model_instance, ALEAF_setting; recorded_investment_decisions=rid_y, transfer_heavy_data=false)
    if ALEAF_model_instance.data["model_type"]["model_type_OP"] != "LP"  # only update the model instance if the model is not LP
        add_ref_LC_GTEP_model!(ALEAF_model_instance_fixed, ALEAF_setting; recorded_investment_decisions=rid_y, transfer_heavy_data=false)
    end

    # Disable transmission expansion
    ALEAF_model_instance.setting["Planning Design"]["multi_round_solution_process_flag"] = false
    ALEAF_model_instance.setting["Simulation Configuration"]["transmission_expansion_flag"] = false
    if ALEAF_model_instance.data["model_type"]["model_type_OP"] != "LP"  # only update the model instance if the model is not LP
        ALEAF_model_instance_fixed.setting["Simulation Configuration"]["transmission_expansion_flag"] = false
        ALEAF_model_instance_fixed.setting["Planning Design"]["multi_round_solution_process_flag"] = false
    end

    # build LCO GEP model instance 
    build_LCO_GTEP_operational_mode_instance!(ALEAF_model_instance, result_LC_GTEP_expansion, day_group_id, y; expansion_fix)

    # Define JuMP_model, solution_list, and solver setting
    JuMP_model = ALEAF_model_instance.model[:nw][ALEAF_model_instance.cnw][day_group_id]
    solution_list = ALEAF_model_instance.sol[:nw][ALEAF_model_instance.cnw][day_group_id]
    solver_setting = ALEAF_model_instance.setting["Solver Setting"]
    iteration = 1

    # Solve ALEAF LC_GTEP Operational Model
    result_LC_GTEP_operation = solve_model_GTEP!(JuMP_model, day_group_id, solution_list, solver_setting; iteration, PH_flag=false)
    @aleaf_info "[ALEAF LC_GTEP Operational Model]: Solved Year id: $y, Day Group id: $day_group_id, Start day: $start_day, End day: $end_day,"

    # Solve model again with fixed integer variables to get dual values
    if ALEAF_model_instance.data["model_type"]["model_type_OP"] == "LP"

        # collect dual values
        collect_dual_LC_GTEP_operation!(JuMP_model, ALEAF_model_instance, result_LC_GTEP_operation, day_group_id, y)

    else

        @aleaf_info "[ALEAF LC_GTEP Operational Model]: Relax MILP model to obtain duals"

        build_LCO_GTEP_operational_mode_instance!(ALEAF_model_instance_fixed, result_LC_GTEP_expansion, day_group_id, y; expansion_fix, operation_fix=true, result_operation = result_LC_GTEP_operation)
        JuMP_model_fixed = ALEAF_model_instance_fixed.model[:nw][ALEAF_model_instance_fixed.cnw][day_group_id]
        solution_list = ALEAF_model_instance_fixed.sol[:nw][ALEAF_model_instance_fixed.cnw][day_group_id]
        solver_setting = ALEAF_model_instance_fixed.setting["Solver Setting"]
        iteration = 1

        result_LC_GTEP_operation_fixed = solve_model_GTEP!(JuMP_model_fixed, day_group_id, solution_list, solver_setting; iteration, PH_flag=false)

        # collect dual values
        collect_dual_LC_GTEP_operation!(JuMP_model_fixed, ALEAF_model_instance, result_LC_GTEP_operation, day_group_id, y)
    end

    # ---------------------------------------------------------------------------------
    # Worker-side reporting: each worker streams ONLY its own day-group's dispatch/market/policy/
    # demand/power-flow rows into per-dg files (suffix "__dg<id>") on the shared filesystem, then
    # returns a small aggregate payload. The master never holds the per-hour solution.
    # ---------------------------------------------------------------------------------
    pu_power_base = ALEAF_setting["Simulation Setting"]["per_unit_base_value"]
    pu_econ_base  = ALEAF_setting["Simulation Setting"]["per_unit_econ_base_value"] / pu_power_base
    apply_pu_to_solution!(result_LC_GTEP_operation["solution"], pu_power_base, pu_econ_base)

    # minimal result-shaped view: one day-group, single solved year. The reporters derive indices from
    # reference["1"] and read per-year metadata from reference[string(year)]; both are topology/config
    # identical across stages, so alias every stage id to this worker's single-year ref[:nw]. op_result
    # holds only year y, and the reporters skip stages absent from op_result.
    single_ref = ALEAF_model_instance.ref[:nw]
    reference = Dict{String,Any}(string(s) => single_ref for s in 1:ALEAF_setting["Planning Design"]["num_stages_value"])
    opres     = Dict{String,Any}(string(y) => Dict{String,Any}(string(day_group_id) => result_LC_GTEP_operation))
    wresult   = Dict{String,Any}(
        "operational model system reference" => reference,
        "operational model result"           => opres,
        "setting"                            => ALEAF_model_instance.setting,
    )

    dispatch_OP_enabled = report_enabled(ALEAF_setting, "report_dispatch_OP_flag")
    summary_OP_enabled  = report_enabled(ALEAF_setting, "report_summary_OP_flag")
    # year in the suffix keeps year-less files (market/policy/demand/power_flow) from colliding across
    # stages; dispatch files already carry _year_<y> so the dg part alone would also be unique there.
    dg_suffix = "__y$(y)_dg$(day_group_id)"

    if dispatch_OP_enabled
        report_result_dispatch_operation_GTEP(wresult; day_group_subset=[day_group_id], file_suffix=dg_suffix)
    elseif summary_OP_enabled
        # summaries need annual_gen_info_OP + annual_scarcity_info_OP + __market_OP rows only
        report_result_dispatch_operation_GTEP(wresult; generate_file=false, market_only=true, day_group_subset=[day_group_id], file_suffix=dg_suffix)
    end

    if report_enabled(ALEAF_setting, "report_power_flow_OP_flag")
        report_result_power_flow_operation_GTEP(wresult, result_LC_GTEP_expansion; day_group_subset=[day_group_id], file_suffix=dg_suffix)
    end

    # annual_gen_info_OP / annual_scarcity_info_OP are stashed as top-level keys on `reference`
    # (result["operational model system reference"]) by the dispatch reporter
    agi = get(reference, "annual_gen_info_OP", nothing)
    asi = get(reference, "annual_scarcity_info_OP", nothing)

    payload = Dict{String,Any}(
        "objective"             => result_LC_GTEP_operation["objective"],
        "annual_gen_info_contrib" => agi,
        "scarcity_totals"       => asi,
    )
    # build decisions are day-group-invariant; capture them once (from dg 1) for the summaries
    # Build decisions are day-group-invariant; only the designated day-group carries them back to the
    # master (avoids shipping them from every worker). The driver designates the smallest day-group id
    # in the run, so this works for a full run (dg 1) and for any future day-group subset.
    if day_group_id == capture_build_decisions_dg
        payload["build_decisions"] = Dict{String,Any}(
            "expansion" => result_LC_GTEP_operation["solution"]["expansion"],
            "slack"     => get(result_LC_GTEP_operation["solution"], "slack", Dict{String,Any}()),
        )
    end
    # Worker memory hygiene: drop this day-group's heavy structures and force reclamation so native
    # CPLEX memory (freed only by finalizers on GC) and the Julia heap are returned BETWEEN the
    # day-groups this worker processes — otherwise they accumulate across successive remotecalls.
    # `payload` already holds the small pieces it needs (objective, aggregates, and references to the
    # build-decision sub-dicts). Mirrors the multi-round GC pattern.
    JuMP_model = nothing
    ALEAF_model_instance = nothing
    ALEAF_model_instance_fixed = nothing
    network_data = nothing
    wresult = nothing; reference = nothing; opres = nothing; single_ref = nothing
    result_LC_GTEP_operation = nothing
    GC.gc(); GC.gc()

    return payload
    # ---------------------------------------------------------------------------------


end



function collect_dual_expansion_GTEP!(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, decomp_id::Int; round_ids_y = [])

    dual_list = Dict{String, Any}()

    ids_d = [(d) for (d) in get_index(am, :repdays, 0)]
    ids_day_groups = [(c) for (c) in get_index(am, :repday_groups, 0)]

    if am.setting["Planning Design"]["multi_round_solution_process_flag"] == true
        ids_y = round_ids_y
    else
        ids_y = [(y) for (y) in get_index(am, :planning_stages, 0)]
    end

    ids_pr = [(z) for (z) in get_index(am, :zone, 0, "planning_reserve")]
    ids_z = [(z) for (z) in get_index(am, :zone, 0, "reserve")]
    ids_p = [(z) for (z) in get_index(am, :zone, 0, "policy")]
    ids_zdht = [(z,d,h,t) for (z) in ids_z for d in ids_d for h in am.setting["run_H"] for t in am.setting["run_T"]]
    ids_ndht = [(n,d,h,t) for (n) in get_index(am, :bus, 0) for (d) in ids_d for h in am.setting["run_H"] for t in am.setting["run_T"]]
    ids_ndy = [(n,d,y) for (n) in get_index(am, :bus, 0) for (d) in ids_d for y in ids_y]
    ids_pdy = [(p,d,y) for (p) in ids_p for (d) in ids_d for y in ids_y]
    ids_n = [(n) for (n) in get_index(am, :bus, 0)]

    for y in ids_y

        stage_length = am.ref[:nw][0][:planning_stages][y]["stage_length"]
        current_year = am.setting["Planning Design"]["base_year_value"]

        discount_factor = 0.0
        for future_y in 1:stage_length
            future_year = am.ref[:nw][0][:planning_stages][y]["year"] + future_y - 1 
            discount_factor += 1.0
        end

        for (n,d,h,t) in ids_ndht
            
            idx = (n,d,h,t,y)
            
            if !haskey(dual_list, string(idx))
                dual_list[string(idx)] = Dict{String, Any}()
            end

            num_days = parameter(am, 0, :repdays, "NumDays", d)

            # LMP
            const_name = "constraint_LoadBalance_ndhty_real"
            if am.setting["Simulation Setting"]["power_flow_mode_flag"] == "PTDF"
                const_name = "constraint_PTDF_LoadBalance_ndhty_real"
            end
            dual_list[string(idx)]["LMP"] = collect_dual_value_ndhty(JuMP_model, const_name, n, d, h, t, y) / (discount_factor * num_days)

        end
    end

    for y in ids_y

        stage_length = am.ref[:nw][0][:planning_stages][y]["stage_length"]
        current_year = am.setting["Planning Design"]["base_year_value"]

        discount_factor = 0.0
        for future_y in 1:stage_length
            future_year = am.ref[:nw][0][:planning_stages][y]["year"] + future_y - 1 
            discount_factor += 1.0
        end
    
        for (z,d,h,t) in ids_zdht
            
            idx = (z,d,h,t,y) 

            if !haskey(dual_list, string(idx))
                dual_list[string(idx)] = Dict{String, Any}()
            end

            num_days = parameter(am, 0, :repdays, "NumDays", d)

            # Regulation
            const_name = "constraint_reg_Up_zdhty"
            dual_list[string(idx)]["Reg_Up_Price"] = collect_dual_value_zdhty(JuMP_model, const_name, z, d, h, t, y) / (discount_factor * num_days)

            # Regulation
            const_name = "constraint_reg_Dn_zdhty"
            dual_list[string(idx)]["Reg_Dn_Price"] = collect_dual_value_zdhty(JuMP_model, const_name, z, d, h, t, y) / (discount_factor * num_days)

            # Spin and Non-spin
            if am.setting["Simulation Configuration"]["Dispatch_Mode_in_EXP"] == "Unit Commitment"
                # cont_R
                const_name = "constraint_R_Spin_zdhty"
                dual_list[string(idx)]["Cont_Res_Price"] = collect_dual_value_zdhty(JuMP_model, const_name, z, d, h, t, y) / (discount_factor * num_days)

                # cont_R
                const_name = "constraint_R_NonSpin_zdhty"
                dual_list[string(idx)]["NonSpin_Res_Price"] = collect_dual_value_zdhty(JuMP_model, const_name, z, d, h, t, y) / (discount_factor * num_days)
            else
                # cont_R
                const_name = "constraint_R_Cont_zdhty"
                dual_list[string(idx)]["Cont_Res_Price"] = collect_dual_value_zdhty(JuMP_model, const_name, z, d, h, t, y) / (discount_factor * num_days)

                dual_list[string(idx)]["NonSpin_Res_Price"] = 0.0
            end

            # Flex_Up
            const_name = "constraint_flex_Up_zdhty"
            dual_list[string(idx)]["Flex_Up_Price"] = collect_dual_value_zdhty(JuMP_model, const_name, z, d, h, t, y) / (discount_factor * num_days)

            # Flex_Down
            const_name = "constraint_flex_Dn_zdhty"
            dual_list[string(idx)]["Flex_Dn_Price"] = collect_dual_value_zdhty(JuMP_model, const_name, z, d, h, t, y) / (discount_factor * num_days)

        end
    end

    # Policy
    # Carbon emission reduction target constraint
    
    if am.setting["Simulation Configuration"]["Carbon_Emission_Reduction_Target_Flag"] == true
        
        for y in ids_y

            stage_length = am.ref[:nw][0][:planning_stages][y]["stage_length"]

            if am.ref[:nw][0][:planning_stages][y]["year"] >= am.setting["Simulation Configuration"]["Carbon_Emission_Reduction_Target_Start_Year"]
            
                for n in ids_p
                    idx = (n, y)

                    if !haskey(dual_list, string(idx))
                        dual_list[string(idx)] = Dict{String, Any}()
                    end

                    const_name = string("constraint_carbon_emission_reduction_target_ny", "_($n,$y)")
                    dual_list[string(idx)]["CERT_MC"] = safe_constraint_dual(JuMP_model, const_name) / stage_length

                end

            end
        end
    end

    # Clean Energy Generation Constraint
    if am.setting["Simulation Configuration"]["Clean_Energy_Generation_Target_EXP_Flag"] == true

        if am.setting["Simulation Configuration"]["Clean_Energy_Generation_Target_EXP_Type"] == "Daygroup"
            for idx in ids_pdy
                (n,d,y) = idx

                stage_length = am.ref[:nw][0][:planning_stages][y]["stage_length"]
                num_days = parameter(am, 0, :repdays, "NumDays", d)

                if am.ref[:nw][0][:planning_stages][y]["year"] >= am.setting["Simulation Configuration"]["Clean_Energy_Generation_Target_Start_Year"]
                    if !haskey(dual_list, string(idx))
                        dual_list[string(idx)] = Dict{String, Any}()
                    end

                    const_name = string("constraint_clean_energy_generation_ndy", "_($n,$d,$y)")
                    dual_list[string(idx)]["CES_MC"] = safe_constraint_dual(JuMP_model, const_name) / (stage_length * num_days)
                end
            end
        elseif am.setting["Simulation Configuration"]["Clean_Energy_Generation_Target_EXP_Type"] == "Annual"
            for y in ids_y
                for p in ids_p
                    idx = (p,y)
                    stage_length = am.ref[:nw][0][:planning_stages][y]["stage_length"]
                    
                    if am.ref[:nw][0][:planning_stages][y]["year"] >= am.setting["Simulation Configuration"]["Clean_Energy_Generation_Target_Start_Year"]
                        if !haskey(dual_list, string(idx))
                            dual_list[string(idx)] = Dict{String, Any}()
                        end

                        const_name = string("constraint_clean_energy_generation_ny", "_($p,$y)")
                        dual_list[string(idx)]["CES_MC"] = safe_constraint_dual(JuMP_model, const_name) / (stage_length)
                    end
                end
            end
        end
    end

    # RPS Constraint
    if am.setting["Simulation Configuration"]["RPS_Flag"] == true
        for y in ids_y

            stage_length = am.ref[:nw][0][:planning_stages][y]["stage_length"]
            
            for n in ids_p
                idx = (n, y)

                if !haskey(dual_list, string(idx))
                    dual_list[string(idx)] = Dict{String, Any}()
                end

                const_name = string("constraint_RPS_regional_ny", "_($n,$y)")
                dual_list[string(idx)]["REC_MC"] = safe_constraint_dual(JuMP_model, const_name) / stage_length

            end
        end
    end

    # RPM
    if am.setting["Simulation Configuration"]["enforce_min_reserve_margin_flag"] == true
        for y in ids_y

            stage_length = am.ref[:nw][0][:planning_stages][y]["stage_length"]
            
            for n in ids_pr
                idx = (n, y)
                
                if !haskey(dual_list, string(idx))
                    dual_list[string(idx)] = Dict{String, Any}()
                end

                const_name = string("constraint_PRM_y_p_real", "_($n,$y)")
                dual_list[string(idx)]["PRM_MC"] = safe_constraint_dual(JuMP_model, const_name) / stage_length
            end
        end
    end

    return dual_list

end


function safe_constraint_dual(JuMP_model::JuMP.AbstractModel, const_name::String)

    const_ref = JuMP.constraint_by_name(JuMP_model, const_name)
    if const_ref === nothing
        return 0.0
    end

    return (try JuMP.dual(const_ref) catch; 0.0 end)

end


function collect_dual_value_zdhty(JuMP_model::JuMP.AbstractModel, const_name::String, z::Int, d::Int, h::Int, t::Int, y::Int)

    const_name = string(const_name, "_($z,$d,$h,$t,$y)")
    dual_value = safe_constraint_dual(JuMP_model, const_name)

    return dual_value
end


function collect_dual_value_ndhty(JuMP_model::JuMP.AbstractModel, const_name::String, n::Int, d::Int, h::Int, t::Int, y::Int)

    const_name = string(const_name, "_($n,$d,$h,$t,$y)")
    dual_value = safe_constraint_dual(JuMP_model, const_name)

    return dual_value
end


function collect_dual_LC_GTEP_operation!(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, result_LC_GEP_operation::Dict{String,<:Any}, day_group_id::Int, y::Int)
    
    ids_d = [(d) for (d) in am.ref[:nw][0][:repday_groups][day_group_id]["Day_Idx_List"]]    
    result_LC_GEP_operation["solution"]["dual"] = Dict{String, Any}()

    ids_z = [(z) for (z) in get_index(am, :zone, 0, "reserve")]
    ids_p = [(z) for (z) in get_index(am, :zone, 0, "policy")]
    ids_zdhty = [(z,d,h,t,y) for (z) in ids_z for d in ids_d for h in am.setting["run_H"] for t in am.setting["run_T"]]
    ids_ndhty = [(n,d,h,t,y) for (n) in get_index(am, :bus, 0) for (d) in ids_d for h in am.setting["run_H"] for t in am.setting["run_T"]]
    ids_pdy = [(p,d,y) for (p) in ids_p for (d) in ids_d]

    for idx in ids_ndhty
        (n,d,h,t,y) = idx
        if !haskey(result_LC_GEP_operation["solution"]["dual"], string(idx))
            result_LC_GEP_operation["solution"]["dual"][string(idx)] = Dict{String, Any}()
        end

        num_days = parameter(am, 0, :repdays, "NumDays", d)
        # num_days = 1.0

        # LMP
        const_name = "constraint_LoadBalance_ndhty_real"
        if am.setting["Simulation Setting"]["power_flow_mode_flag"] == "PTDF"
            const_name = "constraint_PTDF_LoadBalance_ndhty_real"
        end
        result_LC_GEP_operation["solution"]["dual"][string(idx)]["LMP"] = collect_dual_value_ndhty(JuMP_model, const_name, n, d, h, t, y) / num_days
        
    end

    for idx in ids_zdhty
        (z,d,h,t,y) = idx
        if !haskey(result_LC_GEP_operation["solution"]["dual"], string(idx))
            result_LC_GEP_operation["solution"]["dual"][string(idx)] = Dict{String, Any}()
        end

        num_days = parameter(am, 0, :repdays, "NumDays", d)
        # num_days = 1.0

        # Regulation
        const_name = "constraint_reg_Up_zdhty"
        result_LC_GEP_operation["solution"]["dual"][string(idx)]["Reg_Up_Price"] = collect_dual_value_zdhty(JuMP_model, const_name, z, d, h, t, y) / num_days

        # Regulation
        const_name = "constraint_reg_Dn_zdhty"
        result_LC_GEP_operation["solution"]["dual"][string(idx)]["Reg_Dn_Price"] = collect_dual_value_zdhty(JuMP_model, const_name, z, d, h, t, y) / num_days

        # Spin and Non-spin
        if am.setting["Simulation Configuration"]["Dispatch_Mode_in_OP"] == "Unit Commitment"
            const_name = "constraint_R_Spin_zdhty"
            result_LC_GEP_operation["solution"]["dual"][string(idx)]["Cont_Res_Price"] = collect_dual_value_zdhty(JuMP_model, const_name, z, d, h, t, y) / num_days

            const_name = "constraint_R_NonSpin_zdhty"
            result_LC_GEP_operation["solution"]["dual"][string(idx)]["NonSpin_Res_Price"] = collect_dual_value_zdhty(JuMP_model, const_name, z, d, h, t, y) / num_days
        else
            const_name = "constraint_R_Cont_zdhty"
            result_LC_GEP_operation["solution"]["dual"][string(idx)]["Cont_Res_Price"] = collect_dual_value_zdhty(JuMP_model, const_name, z, d, h, t, y) / num_days

            result_LC_GEP_operation["solution"]["dual"][string(idx)]["NonSpin_Res_Price"] = 0
        end

        # Flex_Up
        const_name = "constraint_flex_Up_zdhty"
        result_LC_GEP_operation["solution"]["dual"][string(idx)]["Flex_Up_Price"] = collect_dual_value_zdhty(JuMP_model, const_name, z, d, h, t, y) / num_days

        # Flex_Down
        const_name = "constraint_flex_Dn_zdhty"
        result_LC_GEP_operation["solution"]["dual"][string(idx)]["Flex_Dn_Price"] = collect_dual_value_zdhty(JuMP_model, const_name, z, d, h, t, y) / num_days


    end

    # Clean Energy Generation Constraint
    if am.setting["Simulation Configuration"]["Clean_Energy_Generation_Target_OP_Flag"] == true
        for idx in ids_pdy
            (n,d,y) = idx

            num_days = parameter(am, 0, :repdays, "NumDays", d)
            # num_days = 1.0

            if !haskey(result_LC_GEP_operation["solution"]["dual"], string(idx))
                result_LC_GEP_operation["solution"]["dual"][string(idx)] = Dict{String, Any}()
            end
            
            if am.ref[:nw][0][:planning_stages][y]["year"] >= am.setting["Simulation Configuration"]["Clean_Energy_Generation_Target_Start_Year"]

                const_name = string("constraint_clean_energy_generation_ndy", "_($n,$d,$y)")
                result_LC_GEP_operation["solution"]["dual"][string(idx)]["CES_MC"] = (try JuMP.dual(JuMP.constraint_by_name(JuMP_model, const_name)) catch; 0.0 end) / (num_days)
            end
        end
    end
end


function execute_LC_GTEP_direct_multi_round(ALEAF_setting::Dict{String,<:Any}, case_id::Int)

    # ALEAF_model_instance build-------------------------------------------------
    iteration = 1   # no iteration in the direct mode
    decomp_group = 1    # Only one group when PH is not being used.
    
    # export round info file (JSON format) via atomic temp-then-rename
    function export_multi_round_info(output_path, GTEP_multi_round_info)
        output_path = parameter(ALEAF_model_instance, 0, :output_path)
        file_name = joinpath(output_path, "GTEP_multi_round_info.json")
        tmp_name = string(file_name, ".tmp")

        stringdata = JSON.json(GTEP_multi_round_info)
        open(tmp_name, "w") do f
            write(f, stringdata)
        end
        mv(tmp_name, file_name; force=true)
    end

    # keep only expansion, expansion_line and policy slack (raw) so RA and OP-from-predefined read the checkpoint unchanged
    function prune_checkpoint_solution(round_solution)
        pruned = Dict{String, Any}(string(iteration) => Dict{String, Any}("nw" => Dict{String, Any}(string(decomp_group) => Dict{String, Any}())))
        src = round_solution[string(iteration)][:nw][string(decomp_group)]
        dst = pruned[string(iteration)]["nw"][string(decomp_group)]
        dst["objective"] = src["objective"]
        dst["solution"] = Dict{String, Any}()
        # deepcopy so the later in-place per-unit conversion of the round solution keeps the checkpoint raw
        for cat in ("expansion", "expansion_line", "slack")
            haskey(src["solution"], cat) && (dst["solution"][cat] = deepcopy(src["solution"][cat]))
        end
        return pruned
    end

    # merge only converted expansion, expansion_line and policy slack so end-of-run summaries and OP-after-expansion
    # see the same converted expansion the standard report path produced (dispatch/dual are not carried)
    function merge_expansion_solution!(ALEAF_model_instance, round_solution, round_id)
        dst = ALEAF_model_instance.solution[string(iteration)][:nw][string(decomp_group)]
        src = round_solution[string(iteration)][:nw][string(decomp_group)]
        for cat in ("expansion", "expansion_line", "slack")
            haskey(src["solution"], cat) || continue
            haskey(dst["solution"], cat) || (dst["solution"][cat] = Dict{String, Any}())
            merge!(dst["solution"][cat], src["solution"][cat])
        end
    end

    # convert a pruned checkpoint round (raw expansion+expansion_line) and merge into am.solution.
    # only u_ESE_iy is per-unit-scaled in expansion, so replicate just that conversion here.
    function merge_pruned_round_expansion!(ALEAF_model_instance, pruned_round)
        pu_power_base = ALEAF_setting["Simulation Setting"]["per_unit_base_value"]
        src = deepcopy(pruned_round[string(iteration)]["nw"][string(decomp_group)]["solution"])
        for idx in keys(get(src, "expansion", Dict{String, Any}()))
            haskey(src["expansion"][idx], "u_ESE_iy") && (src["expansion"][idx]["u_ESE_iy"] *= pu_power_base)
        end
        dst = ALEAF_model_instance.solution[string(iteration)][:nw][string(decomp_group)]
        for cat in ("expansion", "expansion_line", "slack")
            haskey(src, cat) || continue
            haskey(dst["solution"], cat) || (dst["solution"][cat] = Dict{String, Any}())
            merge!(dst["solution"][cat], src[cat])
        end
    end

    # ALEAF_model_instance build-------------------------------------------------
    # generate network data (coarse repday selection first; this init build carries no decisions yet)
    pre_groups, pre_days = (Dict{String,Any}(), Dict{String,Any}())
    coarse = select_repdays_coarse(ALEAF_setting, case_id, "expansion")
    coarse !== nothing && ((pre_groups, pre_days) = coarse)
    network_data = generate_networkdata_LC_GTEP(ALEAF_setting, case_id, "expansion"; print_output_flag=false, precomputed_repday_groups=pre_groups, precomputed_repdays=pre_days)

    # build common ALEAF model instance structure. transfer_heavy_data=false: preserve exact current
    # deepcopy behavior for the multi-round EXP driver (OOM fix targets serial nodal OP only).
    ALEAF_model_instance = build_model_structure_GTEP(ALEAF_setting, network_data, Abstract_LC_GTEP_Model)
    add_ref_LC_GTEP_model!(ALEAF_model_instance, ALEAF_setting; transfer_heavy_data=false)
    ALEAF_model_instance.solution = Dict{String, Any}()

    # Initialize solution space    
    if !haskey(ALEAF_model_instance.solution, string(iteration))
        ALEAF_model_instance.solution[string(iteration)] = Dict{Symbol,Any}(:nw => Dict{String,Any}())             
        ALEAF_model_instance.solution[string(iteration)][:nw][string(decomp_group)] = Dict{String,Any}()             
        ALEAF_model_instance.solution[string(iteration)][:nw][string(decomp_group)]["solution"] = Dict{String,Any}()             
    end

    # Initialize Recorded (Accumulated) Investment Decisions
    recorded_investment_decisions = Dict{String, Any}()
    ids_y = [(y) for (y) in get_index(ALEAF_model_instance, :planning_stages, 0)]
    initialize_recorded_investment_decisions!(ALEAF_model_instance, recorded_investment_decisions)
    
    # Build multi-round info data
    GTEP_multi_round_info = Dict()

    # Define simulation Rounds
    num_decision_stages_per_round = ALEAF_setting["Planning Design"]["num_decision_stages_per_round_value"]
    num_lookahead_stages_per_round = ALEAF_setting["Planning Design"]["num_lookahead_stages_per_round_value"]
    if num_decision_stages_per_round < 1 || num_lookahead_stages_per_round < 0
        error("[ALEAF LC_GTEP Expansion Model]: invalid rolling-horizon setting: num_decision_stages_per_round_value = $num_decision_stages_per_round " *
              "(must be >= 1) and num_lookahead_stages_per_round_value = $num_lookahead_stages_per_round (must be >= 0).")
    end
    # each round solves decision + look-ahead stages; only the decision stages are committed
    num_stages_per_simulation_round = num_decision_stages_per_round + num_lookahead_stages_per_round
    num_overlaps_between_simulation_rounds = num_lookahead_stages_per_round

    # Per-round reporting state (report_parts_dir wipe deferred until resume state is known)
    output_path = parameter(ALEAF_model_instance, 0, :output_path)
    report_parts_dir = joinpath(output_path, "_report_parts")
    report_expansion_enabled  = report_enabled(ALEAF_setting, "report_expansion_flag")
    report_scarcity_enabled   = report_enabled(ALEAF_setting, "report_scarcity_EXP_flag")
    report_power_flow_enabled = report_enabled(ALEAF_setting, "report_power_flow_EXP_flag")
    report_dispatch_enabled   = report_enabled(ALEAF_setting, "report_dispatch_EXP_flag")
    report_summary_enabled    = report_enabled(ALEAF_setting, "report_summary_EXP_flag")
    accumulated_annual_gen_info = Dict{Int, Dict{Int, Dict{String, Float64}}}()  # gen -> stage -> metrics
    accumulated_scarcity_info = Dict{Int, Dict{String, Float64}}()                # stage -> scarcity totals
    stage_objective_values = Dict{String, Any}()                                 # stage -> owning round objective
    round_order = Int[]

    # container for per-decision-year network data (populated per round, before that round's reports)
    ALEAF_model_instance.ref[:nw][0][:multi_round_info] = Dict{Int64, Any}()
    for year in ids_y
        ALEAF_model_instance.ref[:nw][0][:multi_round_info][year] = Dict{String, Any}()
    end

    # Beginning of loop
    round_id = 1
    resuming = false

    # Resume an unfinished multi-round run from the schema_version=3 checkpoint (Phase 4)
    checkpoint_path = joinpath(parameter(ALEAF_model_instance, 0, :output_path), "GTEP_multi_round_info.json")
    if ALEAF_setting["Planning Design"]["continue_from_previous_run_flag"] == true && ispath(checkpoint_path)

        prior = nothing
        try
            prior = JSON.parse(open(checkpoint_path))
        catch e
            @aleaf_error "[ALEAF LC_GTEP Expansion Model]: Failed to parse GTEP_multi_round_info.json (corrupt checkpoint)."
            rethrow(e)
        end

        # guard: refuse to resume an incompatible or old-format (bloated, no schema_version) checkpoint.
        # v3 = repday-derived report caches follow each round's own days; v2 pieces on disk would mix with v3 ones
        if get(prior, "schema_version", nothing) != 3 ||
           get(prior, "case_id", nothing) != case_id ||
           get(prior, "num_stages_per_simulation_round", nothing) != num_stages_per_simulation_round ||
           get(prior, "num_overlaps_between_simulation_rounds", nothing) != num_overlaps_between_simulation_rounds ||
           get(prior, "max_planning_year", nothing) != maximum(ids_y)
            error("[ALEAF LC_GTEP Expansion Model]: GTEP_multi_round_info.json is incompatible with the current run "*
                  "(schema/case_id/round config/horizon mismatch, or an old-format file). "*
                  "Delete $checkpoint_path to start fresh.")
        end

        if prior["status"] == "completed"
            @aleaf_warn "[ALEAF LC_GTEP Expansion Model]: prior run already completed; nothing to resume (delete the checkpoint to start fresh)."
            return ALEAF_model_instance

        elseif prior["status"] == "in progress"

            @aleaf_info "[ALEAF LC_GTEP Expansion Model]: Multi-Round simulation [Continue from the last simulation]"

            resuming = true
            last = prior["last round"]
            round_id = last + 1
            GTEP_multi_round_info = deepcopy(prior)

            # restore accumulated investment decisions from the last committed round
            recorded_investment_decisions = deepcopy(prior[string(last)]["updated_investment_decisions"])

            # rebuild reporting accumulators for the pre-resume rounds (JSON keys are strings)
            round_order = collect(1:last)
            for r in 1:last
                for (i_s, per_year) in prior[string(r)]["annual_gen_info"]
                    i = parse(Int, i_s)
                    haskey(accumulated_annual_gen_info, i) || (accumulated_annual_gen_info[i] = Dict{Int, Dict{String, Float64}}())
                    for (y_s, metrics) in per_year
                        accumulated_annual_gen_info[i][parse(Int, y_s)] = Dict{String, Float64}(k => Float64(v) for (k, v) in metrics)
                    end
                end
                for (y_s, totals) in get(prior[string(r)], "scarcity_info", Dict{String, Any}())
                    accumulated_scarcity_info[parse(Int, y_s)] = Dict{String, Float64}(k => Float64(v) for (k, v) in totals)
                end
                round_objective = prior[string(r)]["objective"]
                for decision_year_id in prior[string(r)]["round_ids_y_decision"]
                    stage_objective_values[string(decision_year_id)] = round_objective
                end
                # restore converted expansion into am.solution and per-year RA_Info + repday meta for reports
                merge_pruned_round_expansion!(ALEAF_model_instance, prior[string(r)]["solution"])
                for decision_year_id in prior[string(r)]["round_ids_y_decision"]
                    ALEAF_model_instance.ref[:nw][0][:multi_round_info][decision_year_id]["RA_Info"] = prior[string(r)]["RA_Info"]
                    ALEAF_model_instance.ref[:nw][0][:multi_round_info][decision_year_id]["planning_stages"] = prior[string(r)]["repday_meta"][string(decision_year_id)]
                end
            end
        end
    end

    # wipe the scratch piece dir only on a fresh start; on resume the prior pieces must be preserved
    if !resuming
        ispath(report_parts_dir) && rm(report_parts_dir; force=true, recursive=true)
    end

    # Start Simulation

    while true
        
        @aleaf_info "[ALEAF LC_GTEP Expansion Model]: Multi-Round simulation; current round = $round_id"
        
        GTEP_multi_round_info[string(round_id)] = Dict()
        GTEP_multi_round_info[string(round_id)]["round_idx"] = round_id
        GTEP_multi_round_info[string(round_id)]["round_first_y"] = num_stages_per_simulation_round * (round_id - 1) - num_overlaps_between_simulation_rounds * (round_id - 1) + 1
        GTEP_multi_round_info[string(round_id)]["round_last_y"] = GTEP_multi_round_info[string(round_id)]["round_first_y"] + num_stages_per_simulation_round - 1
        GTEP_multi_round_info[string(round_id)]["RA_Info"] = Dict{String, Any}()

        # check last year
        if GTEP_multi_round_info[string(round_id)]["round_last_y"] > maximum(ids_y)
            GTEP_multi_round_info[string(round_id)]["round_last_y"] = maximum(ids_y)
        end

        # determine current round ids_y and ids_y_decision
        GTEP_multi_round_info[string(round_id)]["round_ids_y"] = [GTEP_multi_round_info[string(round_id)]["round_first_y"]:GTEP_multi_round_info[string(round_id)]["round_last_y"];]

        if GTEP_multi_round_info[string(round_id)]["round_last_y"] == maximum(ids_y)
            # The decision stages in the last round includes all the simulation years in the round
            GTEP_multi_round_info[string(round_id)]["round_ids_y_decision"] = [GTEP_multi_round_info[string(round_id)]["round_first_y"]:GTEP_multi_round_info[string(round_id)]["round_last_y"];]
        else
            GTEP_multi_round_info[string(round_id)]["round_ids_y_decision"] = [GTEP_multi_round_info[string(round_id)]["round_first_y"]:(GTEP_multi_round_info[string(round_id)]["round_first_y"]+num_decision_stages_per_round-1);]
        end
        
        # update recorded_investment_decisions and ELCC values
        GTEP_multi_round_info[string(round_id)]["recorded_investment_decisions"] = deepcopy(recorded_investment_decisions)    
        if round_id > 1 # we only have ELCC results from round_id 2
            if haskey(GTEP_multi_round_info[string(round_id-1)]["RA_Info"], "ELCC_result")
                GTEP_multi_round_info[string(round_id)]["RA_Info"]["ELCC_result_applied"] = GTEP_multi_round_info[string(round_id-1)]["RA_Info"]["ELCC_result"]
            end
        end

        # build and solve LCO GTEP model instance (expansion mode)
        years_to_report = string(GTEP_multi_round_info[string(round_id)]["round_ids_y"])
        decision_years_to_report = string(GTEP_multi_round_info[string(round_id)]["round_ids_y_decision"])
        @aleaf_info "[ALEAF LC_GTEP Expansion Model]: Modeled Years = $years_to_report"
        @aleaf_info "[ALEAF LC_GTEP Expansion Model]: Decision Years = $decision_years_to_report"

        GTEP_multi_round_info[string(round_id)]["solution"], GTEP_multi_round_info[string(round_id)]["network_data"] = build_solve_LCO_GTEP_expansion_multi_round_direct(ALEAF_setting, case_id, decomp_group, GTEP_multi_round_info[string(round_id)], round_id)
        @aleaf_info "[ALEAF LC_GTEP Expansion Model]: Build and solve GTEP model instance (multi-round expansion mode)."
        
        # keep the full solve result and decision years local; only a slim view is kept long-term
        round_solution = GTEP_multi_round_info[string(round_id)]["solution"]
        round_network_data = GTEP_multi_round_info[string(round_id)]["network_data"]
        round_ids_y_decision = GTEP_multi_round_info[string(round_id)]["round_ids_y_decision"]
        push!(round_order, round_id)

        # 1) record decision-year solutions from RAW values
        record_decision_year_solutions!(ALEAF_model_instance, round_solution, recorded_investment_decisions, round_ids_y_decision, iteration, decomp_group)
        GTEP_multi_round_info[string(round_id)]["updated_investment_decisions"] = deepcopy(recorded_investment_decisions)

        # 2) build the small checkpoint entry (raw expansion+expansion_line) before conversion
        GTEP_multi_round_info[string(round_id)]["solution"] = prune_checkpoint_solution(round_solution)
        # stage-keyed objective: each decision stage takes its owning round's objective
        round_objective = round_solution[string(iteration)][:nw][string(decomp_group)]["objective"]
        for decision_year_id in round_ids_y_decision
            stage_objective_values[string(decision_year_id)] = round_objective
        end

        # 3) convert this round's solution once for reporting (same conversion as the standard path)
        apply_pu_expansion_result_GTEP!(ALEAF_setting, Dict{String, Any}("expansion model result" => round_solution[string(iteration)][:nw]))

        # 4) merge converted expansion into am.solution (slim update; no dispatch/dual)
        merge_expansion_solution!(ALEAF_model_instance, round_solution, round_id)

        # perform RA analysis (before reports so multi_round_info RA_Info is populated for this round)
        if ALEAF_setting["Simulation Configuration"][string(case_id)]["Run_RA_flag"] == true

            _ra_t0 = time()
            RA_Info = execute_ALEAF_RA_model(ALEAF_setting, case_id; current_year=GTEP_multi_round_info[string(round_id)]["round_ids_y_decision"][1], recorded_investment_decisions = GTEP_multi_round_info[string(round_id)]["updated_investment_decisions"])
            RA_TIME_IN_MULTI_ROUND[] += time() - _ra_t0
            GTEP_multi_round_info[string(round_id)]["RA_Info"]["RA_metrics"] = RA_Info["RA_metrics"]
            if haskey(RA_Info, "ELCC_result")
                GTEP_multi_round_info[string(round_id)]["RA_Info"]["ELCC_result"] = RA_Info["ELCC_result"]
            end
            if haskey(RA_Info, "capacity_credit_result")
                GTEP_multi_round_info[string(round_id)]["RA_Info"]["capacity_credit_result"] = deepcopy(RA_Info["capacity_credit_result"])
            end

        end

        # 5) populate multi_round_info for this round's decision years, before its reports
        for decision_year_id in round_ids_y_decision
            ALEAF_model_instance.ref[:nw][0][:multi_round_info][decision_year_id] = round_network_data[string(decision_year_id)]
            ALEAF_model_instance.ref[:nw][0][:multi_round_info][decision_year_id]["RA_Info"] = GTEP_multi_round_info[string(round_id)]["RA_Info"]
        end

        # 5b) align the reporting instance's repday-derived caches with the days this round actually solved
        refresh_report_repday_caches!(ALEAF_model_instance, ALEAF_setting, round_network_data, round_ids_y_decision)

        # 6-7) per-round reports into a scratch piece dir (concatenated at finalize)
        round_piece_dir = joinpath(report_parts_dir, string("round_", round_id))
        report_round_pieces!(ALEAF_model_instance, ALEAF_setting, round_solution, round_ids_y_decision, round_piece_dir,
            accumulated_annual_gen_info, accumulated_scarcity_info, report_expansion_enabled, report_scarcity_enabled,
            report_power_flow_enabled, report_dispatch_enabled, report_summary_enabled)

        # drop bulky network data from the checkpoint (already consumed above)
        delete!(GTEP_multi_round_info[string(round_id)], "network_data")

        # Each round builds its own solver instance; the GC has no view of device memory pressure, so
        # without this the freed instances accumulate on the GPU and a later round fails to allocate.
        solver_val = Val(Symbol(ALEAF_setting["ALEAF Master Setup"]["solver_name"]))
        mem_before = gpu_memory_info(solver_val)
        GC.gc()
        gpu_reclaim!(solver_val)
        mem_after = gpu_memory_info(solver_val)
        if mem_before !== nothing && mem_after !== nothing
            gib(x) = round(x / 2^30, digits=2)
            @aleaf_info string("[ALEAF LC_GTEP Expansion Model]: round ", round_id,
                " GPU free ", gib(mem_before[1]), " -> ", gib(mem_after[1]),
                " GiB of ", gib(mem_after[2]), " (reclaimed ", gib(mem_after[1] - mem_before[1]), ")")
        end

        # persist small resume state outside "solution" (survives pruning): decision-year gen slice + objective
        round_scarcity_slice = Dict{String, Any}(
            string(y) => accumulated_scarcity_info[y] for y in round_ids_y_decision if haskey(accumulated_scarcity_info, y))
        GTEP_multi_round_info[string(round_id)]["scarcity_info"] = round_scarcity_slice

        round_annual_gen_slice = Dict{String, Any}()
        for (i, per_year) in accumulated_annual_gen_info
            for y in round_ids_y_decision
                haskey(per_year, y) || continue
                haskey(round_annual_gen_slice, string(i)) || (round_annual_gen_slice[string(i)] = Dict{String, Any}())
                round_annual_gen_slice[string(i)][string(y)] = per_year[y]
            end
        end
        GTEP_multi_round_info[string(round_id)]["annual_gen_info"] = round_annual_gen_slice
        GTEP_multi_round_info[string(round_id)]["objective"] = round_objective

        # persist per-decision-year repday metadata (no bulky hourly "data") so a resume can rebuild
        # multi_round_info["planning_stages"] for the repday-selection report at finalize
        round_repday_meta = Dict{String, Any}()
        for decision_year_id in round_ids_y_decision
            year_meta = Dict{String, Any}()
            for (day_id, rep) in ALEAF_model_instance.ref[:nw][0][:multi_round_info][decision_year_id]["planning_stages"]
                year_meta[day_id] = Dict{String, Any}(k => v for (k, v) in rep if k != "data")
            end
            round_repday_meta[string(decision_year_id)] = year_meta
        end
        GTEP_multi_round_info[string(round_id)]["repday_meta"] = round_repday_meta

        # guard fields let a resumed/downstream run validate the checkpoint against its own config
        GTEP_multi_round_info["schema_version"] = 3
        GTEP_multi_round_info["case_id"] = case_id
        GTEP_multi_round_info["num_stages_per_simulation_round"] = num_stages_per_simulation_round
        GTEP_multi_round_info["num_overlaps_between_simulation_rounds"] = num_overlaps_between_simulation_rounds
        GTEP_multi_round_info["max_planning_year"] = maximum(ids_y)

        if GTEP_multi_round_info[string(round_id)]["round_last_y"] == maximum(ids_y)
            GTEP_multi_round_info["status"] = "completed"
            GTEP_multi_round_info["last_round"] = round_id
            GTEP_multi_round_info["last round"] = round_id
            export_multi_round_info(parameter(ALEAF_model_instance, 0, :output_path), GTEP_multi_round_info)
            break;
        else
            GTEP_multi_round_info["last_round"] = round_id
            GTEP_multi_round_info["last round"] = round_id
            GTEP_multi_round_info["status"] = "in progress"
            export_multi_round_info(parameter(ALEAF_model_instance, 0, :output_path), GTEP_multi_round_info)
            round_id +=1
        end

    end

    # finalize: concatenate per-round pieces into the standard output files
    # checkpoint retention is handled run-wide by cleanup_multi_round_info_files
    finalize_multi_round_reports!(ALEAF_model_instance, ALEAF_setting, output_path, report_parts_dir, round_order,
        accumulated_annual_gen_info, accumulated_scarcity_info, stage_objective_values,
        report_expansion_enabled, report_scarcity_enabled, report_power_flow_enabled,
        report_dispatch_enabled, report_summary_enabled)

    return ALEAF_model_instance

end


# Write per-year report pieces to a scratch dir; accumulate annual-gen-info by decision year only
# (copying just decision-year slices) to avoid zero-year collisions when merging rounds.
function report_round_pieces!(am, ALEAF_setting, round_solution, round_ids_y_decision, round_piece_dir,
        accumulated_annual_gen_info, accumulated_scarcity_info, report_expansion_enabled, report_scarcity_enabled,
        report_power_flow_enabled, report_dispatch_enabled, report_summary_enabled)

    iteration = 1
    reference = am.ref[:nw]
    setting = am.setting   # flattened per-case setting with run_H/run_T, as report_* expect
    result = round_solution[string(iteration)][:nw]
    mkpath(round_piece_dir)

    if report_expansion_enabled
        report_result_expansion_gen_GTEP(result, reference; years=round_ids_y_decision, out_dir=round_piece_dir, setting=setting)
        report_result_expansion_line_GTEP(result, reference, setting; years=round_ids_y_decision, out_dir=round_piece_dir)
    end

    if report_scarcity_enabled
        report_result_scarcity_ens_GTEP(result, reference, setting; years=round_ids_y_decision, out_dir=round_piece_dir)
        report_result_scarcity_reserve_GTEP(result, reference, setting; years=round_ids_y_decision, out_dir=round_piece_dir)
    end

    if report_power_flow_enabled
        report_result_power_flow_expansion_GTEP(result, reference, setting; years=round_ids_y_decision, out_dir=round_piece_dir)
    end

    # dispatch also computes annual_gen_info + annual_scarcity_info (both needed by summaries);
    # write pieces only if dispatch is on, else compute the metrics with all IO suppressed.
    if report_dispatch_enabled
        report_result_dispatch_expansion_GTEP(result, reference, setting; years=round_ids_y_decision, out_dir=round_piece_dir)
    elseif report_summary_enabled
        report_result_dispatch_expansion_GTEP(result, reference, setting; write_csv=false, years=round_ids_y_decision, out_dir=round_piece_dir)
    end

    if report_dispatch_enabled || report_summary_enabled
        round_annual_gen_info = get_annual_gen_info(reference)
        for (i, per_year) in round_annual_gen_info
            haskey(accumulated_annual_gen_info, i) || (accumulated_annual_gen_info[i] = Dict{Int, Dict{String, Float64}}())
            for y in round_ids_y_decision
                accumulated_annual_gen_info[i][y] = per_year[y]
            end
        end

        round_scarcity_info = get_annual_scarcity_info(reference)
        for y in round_ids_y_decision
            haskey(round_scarcity_info, y) && (accumulated_scarcity_info[y] = round_scarcity_info[y])
        end
    end

    # repdays are written once at finalize (its preallocation spans all stages)
    return nothing
end


# Concatenate per-round pieces into standard output files (each gated by its report flag); if summaries
# are on, install accumulated annual-gen-info + stage-keyed objective so they run without re-converting.
function finalize_multi_round_reports!(am, ALEAF_setting, output_path, report_parts_dir, round_order,
        accumulated_annual_gen_info, accumulated_scarcity_info, stage_objective_values,
        report_expansion_enabled, report_scarcity_enabled, report_power_flow_enabled,
        report_dispatch_enabled, report_summary_enabled)

    iteration = 1
    decomp_group = 1
    reference = am.ref[:nw]
    case_name = parameter(reference, 0, :case_name)

    piece_path(round_id, fname) = joinpath(report_parts_dir, string("round_", round_id), fname)

    # concatenate single-file multi-year reports: header from first piece, data rows from all.
    # no-op when no round produced this piece (e.g. its report flag is off) so absent files never crash.
    function concat_pieces(fname)
        any(round_id -> isfile(piece_path(round_id, fname)), round_order) || return
        final_path = joinpath(output_path, fname)
        open(final_path, "w") do out
            wrote_header = false
            for round_id in round_order
                src = piece_path(round_id, fname)
                isfile(src) || continue
                open(src, "r") do inp
                    first_line = true
                    for line in eachline(inp; keep=true)
                        if first_line
                            first_line = false
                            wrote_header || (write(out, line); wrote_header = true)
                        else
                            write(out, line)
                        end
                    end
                end
            end
        end
    end

    # per-year dispatch files: each year is owned by exactly one round, so move each piece
    function move_per_year_pieces(prefix, suffix)
        for round_id in round_order
            piece_round_dir = joinpath(report_parts_dir, string("round_", round_id))
            isdir(piece_round_dir) || continue
            for f in readdir(piece_round_dir)
                (startswith(f, prefix) && endswith(f, suffix)) || continue
                mv(joinpath(piece_round_dir, f), joinpath(output_path, f); force=true)
            end
        end
    end

    if report_expansion_enabled
        concat_pieces(string(case_name, "__gen_expansion_EXP.csv"))
        concat_pieces(string(case_name, "__line_expansion_EXP.csv"))
    end

    if report_scarcity_enabled
        concat_pieces(string(case_name, "__unserved_energy_EXP.csv"))
        concat_pieces(string(case_name, "__reserve_shortfall_EXP.csv"))
    end

    if report_power_flow_enabled
        concat_pieces(string(case_name, "__power_flow_EXP.csv"))
    end

    # market and the bulky pieces exist only when dispatch reporting is on; summaries now take their
    # scarcity totals in memory. concat_pieces/move are defensive, so absent pieces are skipped.
    if report_dispatch_enabled
        concat_pieces(string(case_name, "__market_EXP.csv"))
    end
    if report_dispatch_enabled
        concat_pieces(string(case_name, "__policy_slack_EXP.csv"))
        concat_pieces(string(case_name, "__demand_response_EXP.csv"))
        move_per_year_pieces(string(case_name, "__dispatch_EXP_year_"), ".csv")
    end

    # repdays span all stages (preallocated), so write once over the full horizon
    report_result_repdays_GTEP(reference, am.setting)

    # summaries: read the merged converted expansion + accumulated annual info; do NOT re-apply pu
    if report_summary_enabled
        setting = am.setting   # flattened per-case setting, as report_* expect
        reference[0][:annual_gen_info] = accumulated_annual_gen_info
        reference[0][:annual_scarcity_info] = accumulated_scarcity_info
        merged_result = am.solution[string(iteration)][:nw]
        # stage-keyed objective so system summary finds an objective per stage for any overlap
        merged_result[string(decomp_group)]["Objective Values"] = stage_objective_values
        report_result_tech_summary_expansion_GTEP(merged_result, reference, setting)
        report_result_system_summary_expansion_GTEP(merged_result, reference, setting)
    end

    ispath(report_parts_dir) && rm(report_parts_dir; force=true, recursive=true)

    return nothing
end


function build_solve_LCO_GTEP_expansion_multi_round_direct(ALEAF_setting, case_id, decomp_group, GTEP_multi_round_info, round_id)

    # generate network data (coarse repday selection first; VRE weights re-keyed from precomputed
    # __col_weights recorded with the run reference in record_decision_year_solutions!). Falls back if empty.
    pre_groups, pre_days = (Dict{String,Any}(), Dict{String,Any}())
    coarse = select_repdays_coarse(ALEAF_setting, case_id, "expansion"; GTEP_multi_round_info_data=GTEP_multi_round_info)
    coarse !== nothing && ((pre_groups, pre_days) = coarse)
    network_data = generate_networkdata_LC_GTEP(ALEAF_setting, case_id, "expansion"; print_output_flag=true, GTEP_multi_round_info_data = GTEP_multi_round_info, precomputed_repday_groups=pre_groups, precomputed_repdays=pre_days)

    # build common ALEAF model instance structure. transfer_heavy_data=false: this multi-round builder
    # keeps reading network_data's repday payload after add_ref (LP-relaxation reuse, direct solve, and
    # the repdays returned for the RA / next-round handoff), so preserve the full deepcopy here.
    ALEAF_model_instance = build_model_structure_GTEP(ALEAF_setting, network_data, Abstract_LC_GTEP_Model)
    add_ref_LC_GTEP_model!(ALEAF_model_instance, ALEAF_setting; transfer_heavy_data=false)
    @aleaf_info "[ALEAF LC_GTEP Expansion Model]: Create a GTEP expansion model structure and add reference data."

    # add RA information to the reference
    if ALEAF_setting["Simulation Configuration"][string(case_id)]["update_CAPCRED_in_each_round_of_Expansion_Flag"]== true
        if haskey(GTEP_multi_round_info["RA_Info"], "ELCC_result_applied")
            update_capacity_credits_system_wide!(ALEAF_model_instance, GTEP_multi_round_info)
        end
    end

    # get relaxed investment decisions by solving LP
    LP_expansion_solution = Dict{String,Any}()
    if (ALEAF_model_instance.data["model_type"]["model_type_EXP"] == "MILP") && (ALEAF_setting["Simulation Setting"]["MIP_relaxed_solution_bounds_flag"] == true)
        LP_expansion_solution = get_relaxed_investment_decisions_for_GTEP_expansion_model(ALEAF_setting, decomp_group, network_data, case_id; GTEP_multi_round_info)
    end

    # build LCO GEP model instance (expansion mode)
    build_LCO_GTEP_expansion_mode_instance!(ALEAF_model_instance, decomp_group; GTEP_multi_round_info_data = GTEP_multi_round_info, ref_expansion_result=LP_expansion_solution)
    @aleaf_info "[ALEAF LC_GTEP Expansion Model]: Build GTEP model instance (expansion mode)"

    # Export ALEAF_model_instance model
    if ALEAF_setting["Simulation Setting"]["export_model_lp_expansion_flag"] == true
        output_path = parameter(ALEAF_model_instance, 0, :output_path)
        file_name = ALEAF_setting["Simulation Setting"]["model_lp_file_name_value"]
        file_name = string("Expansion_", round_id, "_", file_name)
        JuMP.write_to_file(ALEAF_model_instance.model[:nw][ALEAF_model_instance.cnw][1], joinpath(output_path, file_name))
        @aleaf_info string("[ALEAF LC_GTEP Expansion Model]: Export ALEAF_model_instance model at ", joinpath(output_path, file_name), ": file saved")
    end

    # Solve
    solve_LC_GTEP_direct_multi_round(ALEAF_model_instance, ALEAF_setting, network_data, decomp_group, GTEP_multi_round_info)

    network_data_to_return = Dict{String, Any}()
    for decision_year_id in GTEP_multi_round_info["round_ids_y_decision"]
        network_data_to_return[string(decision_year_id)] = Dict{String, Any}()
        network_data_to_return[string(decision_year_id)]["planning_stages"] = network_data["planning_stages"][string(GTEP_multi_round_info["round_first_y"])]["repdays"]
    end
        
    return ALEAF_model_instance.solution, network_data_to_return

end


function update_capacity_credits_system_wide!(am, GTEP_multi_round_info)

    @aleaf_info "[ALEAF LC_GTEP Expansion Model]: Update capacity credits based on ELCC values"
    
    ids_i = [(i) for (i) in get_index(am, :gen_index, 0)]
    
    for i in ids_i
        bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
        tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)    
        unit_group = parameter(am, 0, :gen_index, "UNIT_GROUP", i)
        ELCC_Flag = parameter(am, bus_idx, :gen_bus, tech_idx, "ELCC_Flag")

        if ELCC_Flag == true
            try
                new_CAPCRED = GTEP_multi_round_info["RA_Info"]["ELCC_result_applied"][unit_group]
                if new_CAPCRED != "NA"
                    am.ref[:nw][bus_idx][:gen_bus][tech_idx]["CAPCRED"] = new_CAPCRED                
                    @aleaf_info "[ALEAF LC_GTEP Expansion Model]: - Unit Group: $unit_group,\t New CAPCRED: $new_CAPCRED"
                else
                    @aleaf_info "[ALEAF LC_GTEP Expansion Model]: Filed to find a new ELCC value of $unit_group"
                end
            catch
                @aleaf_info "[ALEAF LC_GTEP Expansion Model]: Filed to find a new ELCC value of $unit_group"
            end
        end
    end

end


function initialize_recorded_investment_decisions!(ALEAF_model_instance, recorded_investment_decisions; reset_flag::Bool=false)

    ids_i = [(i) for (i) in get_index(ALEAF_model_instance, :gen_index, 0)]
    ids_i_sto = [(i) for (i) in get_index(ALEAF_model_instance, :gen_index, 0) if parameter(ALEAF_model_instance, 0, :gen_index, "UNIT_CATEGORY", i) == "STORAGE"]
    ids_k = [(k) for (k) in get_index(ALEAF_model_instance, :branch, 0) if parameter(ALEAF_model_instance, 0, :branch, "model_flag", k) == true]
    ids_y = [(y) for (y) in get_index(ALEAF_model_instance, :planning_stages, 0)]
    ids_n = [(k) for (k) in get_index(ALEAF_model_instance, :bus, 0)]

    if reset_flag == false

        recorded_investment_decisions["u_new_G_iy"] = Dict()
        recorded_investment_decisions["u_ret_G_iy"] = Dict()
        recorded_investment_decisions["u_new_ESH_iy"] = Dict()

        recorded_investment_decisions["U_NEW_G_i"] = Dict(string(i) => 0.0 for i in ids_i)
        recorded_investment_decisions["U_RET_G_i"] = Dict(string(i) => 0.0 for i in ids_i)
        recorded_investment_decisions["U_NEW_ESH_i"] = Dict(string(i) => 0.0 for i in ids_i_sto)
        recorded_investment_decisions["U_ESE_i"] = Dict(string(i) => 0.0 for i in ids_i_sto)
        recorded_investment_decisions["U_NEW_T_k"] = Dict(string(k) => 0.0 for k in ids_k)
        recorded_investment_decisions["U_G_i"] = Dict(string(i) => 0.0 for i in ids_i)
        
        recorded_investment_decisions["Wind_TotalMW"] = Dict(string(i) => 0.0 for i in ids_n)   # for scenario selection
        recorded_investment_decisions["PV_TotalMW"] = Dict(string(i) => 0.0 for i in ids_n) # for scenario selection
        recorded_investment_decisions["RTPV_TotalMW"] = Dict(string(i) => 0.0 for i in ids_n) # for scenario selection
        recorded_investment_decisions["4hr_Storage_TotalMW"] = Dict(string(i) => 0.0 for i in ids_n) # for scenario selection
        recorded_investment_decisions["8hr_Storage_TotalMW"] = Dict(string(i) => 0.0 for i in ids_n) # for scenario selection
        recorded_investment_decisions["10hr_Storage_TotalMW"] = Dict(string(i) => 0.0 for i in ids_n) # for scenario selection
        recorded_investment_decisions["20hr_Storage_TotalMW"] = Dict(string(i) => 0.0 for i in ids_n) # for scenario selection

        for i in ids_i
            bus_idx = parameter(ALEAF_model_instance, 0, :gen_index, "bus_idx", i)
            tech_idx = parameter(ALEAF_model_instance, 0, :gen_index, "genco_tech_id", i)    
            EXUNITS = parameter(ALEAF_model_instance, bus_idx, :gen_bus, tech_idx, "EXUNITS")
            profile_type = parameter(ALEAF_model_instance, 0, :gen_index, "Profile_Type", i)
            UNIT_CATEGORY = parameter(ALEAF_model_instance, 0, :gen_index, "UNIT_CATEGORY", i)
            CAP = parameter(ALEAF_model_instance, bus_idx, :gen_bus, tech_idx, "CAP")

            # WIND <- onshore wind only
            if profile_type == "wind_ons"
                recorded_investment_decisions["Wind_TotalMW"][string(bus_idx)] += CAP * (EXUNITS)
            end

            # PV + RTPV, No CSP
            if profile_type == "pv"
                recorded_investment_decisions["PV_TotalMW"][string(bus_idx)] += CAP * (EXUNITS)
            end

            if profile_type == "rtpv"
                recorded_investment_decisions["RTPV_TotalMW"][string(bus_idx)] += CAP * (EXUNITS)
            end

            # 4 hour storage
            if UNIT_CATEGORY == "STORAGE"
                if parameter(ALEAF_model_instance, bus_idx, :gen_bus, tech_idx, "STOHR_MIN") == 4
                    recorded_investment_decisions["4hr_Storage_TotalMW"][string(bus_idx)] += CAP * (EXUNITS)
                elseif parameter(ALEAF_model_instance, bus_idx, :gen_bus, tech_idx, "STOHR_MIN") == 8
                    recorded_investment_decisions["8hr_Storage_TotalMW"][string(bus_idx)] += CAP * (EXUNITS)
                elseif parameter(ALEAF_model_instance, bus_idx, :gen_bus, tech_idx, "STOHR_MIN") == 10
                    recorded_investment_decisions["10hr_Storage_TotalMW"][string(bus_idx)] += CAP * (EXUNITS)
                elseif parameter(ALEAF_model_instance, bus_idx, :gen_bus, tech_idx, "STOHR_MIN") == 20
                    recorded_investment_decisions["20hr_Storage_TotalMW"][string(bus_idx)] += CAP * (EXUNITS)
                end
            end

        end

    else

        recorded_investment_decisions["Wind_TotalMW"] = Dict(string(i) => 0.0 for i in ids_n)   # for scenario selection
        recorded_investment_decisions["PV_TotalMW"] = Dict(string(i) => 0.0 for i in ids_n) # for scenario selection
        recorded_investment_decisions["RTPV_TotalMW"] = Dict(string(i) => 0.0 for i in ids_n) # for scenario selection
        recorded_investment_decisions["4hr_Storage_TotalMW"] = Dict(string(i) => 0.0 for i in ids_n) # for scenario selection
        recorded_investment_decisions["8hr_Storage_TotalMW"] = Dict(string(i) => 0.0 for i in ids_n) # for scenario selection
        recorded_investment_decisions["10hr_Storage_TotalMW"] = Dict(string(i) => 0.0 for i in ids_n) # for scenario selection
        recorded_investment_decisions["20hr_Storage_TotalMW"] = Dict(string(i) => 0.0 for i in ids_n) # for scenario selection
    end


end


function record_decision_year_solutions!(ALEAF_model_instance, result, recorded_investment_decisions, current_round_ids_y_decision, iteration, decomp_group; nw_key=:nw)

    ids_i = [(i) for (i) in get_index(ALEAF_model_instance, :gen_index, 0)]
    ids_i_sto = [(i) for (i) in get_index(ALEAF_model_instance, :gen_index, 0) if parameter(ALEAF_model_instance, 0, :gen_index, "UNIT_CATEGORY", i) == "STORAGE"]
    ids_k = [(k) for (k) in get_index(ALEAF_model_instance, :branch, 0) if parameter(ALEAF_model_instance, 0, :branch, "model_flag", k) == true]
    
    for decision_y in current_round_ids_y_decision
        for i in ids_i
            
            i_y_index = string("(", i, ", ", decision_y, ")")
            u_new_G_iy = result[string(iteration)][nw_key][string(decomp_group)]["solution"]["expansion"][i_y_index]["u_new_G_iy"]
            u_ret_G_iy = result[string(iteration)][nw_key][string(decomp_group)]["solution"]["expansion"][i_y_index]["u_ret_G_iy"]
            u_G_iy = result[string(iteration)][nw_key][string(decomp_group)]["solution"]["expansion"][i_y_index]["u_G_iy"]

            recorded_investment_decisions["u_new_G_iy"][i_y_index] = u_new_G_iy
            if u_new_G_iy > 0.001
                recorded_investment_decisions["U_NEW_G_i"][string(i)] += u_new_G_iy # accumulated investment decisions
            end
            
            recorded_investment_decisions["u_ret_G_iy"][i_y_index] = u_ret_G_iy
            if u_ret_G_iy > 0.001
                recorded_investment_decisions["U_RET_G_i"][string(i)] += u_ret_G_iy # accumulated investment decisions
            end

            if u_G_iy > 0.001
                recorded_investment_decisions["U_G_i"][string(i)] = u_G_iy # NOT accumulated investment decisions
            end
        end

        for i in ids_i_sto

            i_y_index = string("(", i, ", ", decision_y, ")")
            u_new_ESH_iy = result[string(iteration)][nw_key][string(decomp_group)]["solution"]["expansion"][i_y_index]["u_new_ESH_iy"]
            u_ESE_iy = result[string(iteration)][nw_key][string(decomp_group)]["solution"]["expansion"][i_y_index]["u_ESE_iy"]

            recorded_investment_decisions["u_new_ESH_iy"][i_y_index] = u_new_ESH_iy
            if u_new_ESH_iy > 0.001
                recorded_investment_decisions["U_NEW_ESH_i"][string(i)] += u_new_ESH_iy
            end

            if u_ESE_iy > 0.001
                # u_ESE_iy is a stock value for the selected year, not an annual increment.
                recorded_investment_decisions["U_ESE_i"][string(i)] = u_ESE_iy
            end
        end

        for k in ids_k
            if result[string(iteration)][nw_key][string(decomp_group)]["solution"]["expansion_line"][string("(", k, ", ", decision_y, ")")]["u_new_T_ky"] > 0.001
                recorded_investment_decisions["U_NEW_T_k"][string(k)] += result[string(iteration)][nw_key][string(decomp_group)]["solution"]["expansion_line"][string("(", k, ", ", decision_y, ")")]["u_new_T_ky"]
            end
        end
    end

    # reset total wind/pv MW
    initialize_recorded_investment_decisions!(ALEAF_model_instance, recorded_investment_decisions; reset_flag = true)

    # update total wind/pv MW
    for i in ids_i
        bus_idx = parameter(ALEAF_model_instance, 0, :gen_index, "bus_idx", i)
        tech_idx = parameter(ALEAF_model_instance, 0, :gen_index, "genco_tech_id", i)    
        EXUNITS = parameter(ALEAF_model_instance, bus_idx, :gen_bus, tech_idx, "EXUNITS")
        profile_type = parameter(ALEAF_model_instance, 0, :gen_index, "Profile_Type", i)
        UNIT_CATEGORY = parameter(ALEAF_model_instance, 0, :gen_index, "UNIT_CATEGORY", i)
        CAP = parameter(ALEAF_model_instance, bus_idx, :gen_bus, tech_idx, "CAP")

        # WIND <- onshore wind only
        if profile_type == "wind_ons"
            recorded_investment_decisions["Wind_TotalMW"][string(bus_idx)] += CAP * (recorded_investment_decisions["U_G_i"][string(i)])
        end

        # PV + RTPV, No CSP
        if profile_type == "pv"
            recorded_investment_decisions["PV_TotalMW"][string(bus_idx)] += CAP * (recorded_investment_decisions["U_G_i"][string(i)])
        end

        if profile_type == "rtpv"
            recorded_investment_decisions["RTPV_TotalMW"][string(bus_idx)] += CAP * (recorded_investment_decisions["U_G_i"][string(i)])
        end

        # 4 hour storage
        if UNIT_CATEGORY == "STORAGE"
            if parameter(ALEAF_model_instance, bus_idx, :gen_bus, tech_idx, "STOHR_MIN") == 4
                recorded_investment_decisions["4hr_Storage_TotalMW"][string(bus_idx)] += CAP * (recorded_investment_decisions["U_G_i"][string(i)])
            elseif parameter(ALEAF_model_instance, bus_idx, :gen_bus, tech_idx, "STOHR_MIN") == 8
                recorded_investment_decisions["8hr_Storage_TotalMW"][string(bus_idx)] += CAP * (recorded_investment_decisions["U_G_i"][string(i)])
            elseif parameter(ALEAF_model_instance, bus_idx, :gen_bus, tech_idx, "STOHR_MIN") == 10
                recorded_investment_decisions["10hr_Storage_TotalMW"][string(bus_idx)] += CAP * (recorded_investment_decisions["U_G_i"][string(i)])
            elseif parameter(ALEAF_model_instance, bus_idx, :gen_bus, tech_idx, "STOHR_MIN") == 20
                recorded_investment_decisions["20hr_Storage_TotalMW"][string(bus_idx)] += CAP * (recorded_investment_decisions["U_G_i"][string(i)])
            end
        end

    end

    # Precompute per-data-column VRE weights while the run reference (finest-region map) is in scope, so a
    # later coarse repday selection can re-key run-bus-keyed totals exactly without the run network.
    run_region_map = get(ALEAF_model_instance.ref[:nw][0], :profile_data_region_map, nothing)
    recorded_investment_decisions["__col_weights"] = build_repday_selection_col_weights(
        ALEAF_model_instance.ref[:nw][0][:bus], run_region_map,
        recorded_investment_decisions["Wind_TotalMW"],
        recorded_investment_decisions["PV_TotalMW"],
        recorded_investment_decisions["RTPV_TotalMW"])

end


function solve_LC_GTEP_direct_multi_round(ALEAF_model_instance::Abstract_ALEAF_Model, ALEAF_setting, network_data, decomp_group, GTEP_multi_round_info)

    current_round_ids_y = GTEP_multi_round_info["round_ids_y"]
    
    # Define JuMP_model, solution_list, and solver setting
    JuMP_model = ALEAF_model_instance.model[:nw][ALEAF_model_instance.cnw][decomp_group]
    solution_list = ALEAF_model_instance.sol[:nw][ALEAF_model_instance.cnw][decomp_group]
    solver_setting = ALEAF_model_instance.setting["Solver Setting"]
    iteration = 1
    output = Dict{String, Any}()
    
    # Solve ALEAF LC_GEP Expansion Model
    output[string(decomp_group)] = solve_model_GTEP!(JuMP_model, decomp_group, solution_list, solver_setting; iteration, PH_flag=false)
    @aleaf_info "[ALEAF LC_GTEP Expansion Model]: Solved."

    # initialize solution dict
    if !haskey(ALEAF_model_instance.solution, string(iteration))
        ALEAF_model_instance.solution[string(iteration)] = Dict{Symbol,Any}(:nw => Dict{String,Any}())             
        ALEAF_model_instance.solution[string(iteration)][:nw][string(decomp_group)] = Dict{String,Any}()             
    end

    # Save solution
    if output[string(decomp_group)] == "No Solution"
        ALEAF_model_instance.solution[string(iteration)][:nw][string(decomp_group)] = "No Solution"
    else
        ALEAF_model_instance.solution[string(iteration)][:nw][string(decomp_group)] = output[string(decomp_group)]
    end

    # Get Dual Values
    dual_status = JuMP.dual_status(JuMP_model)
    if dual_status == MOI.NO_SOLUTION 
        # fix discrete variables and resolve
        JuMP.fix_discrete_variables(JuMP_model)    # fix discrete variables
        JuMP.optimize!(JuMP_model)
    end
    # collect dual values
    get_dual_LC_GTEP_expansion_model!(ALEAF_model_instance, ALEAF_setting, network_data; current_round_ids_y = current_round_ids_y)
    @aleaf_info "[ALEAF LC_GTEP Expansion Model]: Obtain Dual."

    return ALEAF_model_instance
end


function get_relaxed_investment_decisions_for_GTEP_expansion_model(ALEAF_setting::Dict{String,<:Any}, decomp_group, network_data, case_id; iteration=1, GTEP_multi_round_info::Dict{Any,<:Any} = Dict{Any,Any}())

    @aleaf_info "[ALEAF LC_GTEP Expansion Model]: === MILP bound reconfiguration process"
    
    # build common ALEAF model instance structure. transfer_heavy_data=false: this shares the caller's
    # network_data, which the parent multi-round builder keeps reading after this LP relaxation returns.
    LP_ALEAF_model_instance = build_model_structure_GTEP(ALEAF_setting, network_data, Abstract_LC_GTEP_Model)
    add_ref_LC_GTEP_model!(LP_ALEAF_model_instance, ALEAF_setting; transfer_heavy_data=false)
    @aleaf_info "[ALEAF LC_GTEP Expansion Model]: Create a LP GTEP expansion model structure and add reference data."

    # add RA information to the reference
    if ALEAF_setting["Simulation Configuration"][string(case_id)]["update_CAPCRED_in_each_round_of_Expansion_Flag"]== true
        if haskey(GTEP_multi_round_info["RA_Info"], "ELCC_result_applied")
            update_capacity_credits_system_wide!(LP_ALEAF_model_instance, GTEP_multi_round_info)
        end
    end

    # 2) reset integrality settings to false
    for i in [(i) for (i) in ALEAF.get_index(LP_ALEAF_model_instance, :gen_index, 0)]
        bus_idx = parameter(LP_ALEAF_model_instance, 0, :gen_index, "bus_idx", i)
        tech_idx = parameter(LP_ALEAF_model_instance, 0, :gen_index, "genco_tech_id", i)
        unit_group = parameter(LP_ALEAF_model_instance, 0, :gen_index, "UNIT_GROUP", i)

        if LP_ALEAF_model_instance.ref[:nw][bus_idx][:gen_bus][tech_idx]["Integrality"] == true
            LP_ALEAF_model_instance.ref[:nw][bus_idx][:gen_bus][tech_idx]["Integrality"] = false
        end
    end

    # build LCO GEP model instance (expansion mode)
    build_LCO_GTEP_expansion_mode_instance!(LP_ALEAF_model_instance, decomp_group; GTEP_multi_round_info_data = GTEP_multi_round_info)
    @aleaf_info "[ALEAF LC_GTEP Expansion Model]: Build a LP GTEP expansion model instance"

    # Define JuMP_model, solution_list, and solver setting
    JuMP_model = LP_ALEAF_model_instance.model[:nw][LP_ALEAF_model_instance.cnw][decomp_group]
    solution_list = LP_ALEAF_model_instance.sol[:nw][LP_ALEAF_model_instance.cnw][decomp_group]
    solver_setting = LP_ALEAF_model_instance.setting["Solver Setting"]

    solution = solve_model_GTEP!(JuMP_model, decomp_group, solution_list, solver_setting; iteration, PH_flag=false)
    LP_ALEAF_model_instance = 0.0   # delete
    @aleaf_info "[ALEAF LC_GTEP Expansion Model]: Solved a relaxed LP GTEP model."

    return solution["solution"]["expansion"]

end


function build_and_solve_LC_GTEP_expansion_model(ALEAF_setting::Dict{String,<:Any}, ALEAF_model_instance::Abstract_ALEAF_Model, decomp_group, network_data; iteration=1, GTEP_multi_round_info::Dict{Any,<:Any} = Dict{Any,Any}())

    function get_solution(ALEAF_model_instance, decomp_group; iteration = 1)

        # Define JuMP_model, solution_list, and solver setting
        JuMP_model = ALEAF_model_instance.model[:nw][ALEAF_model_instance.cnw][decomp_group]
        solution_list = ALEAF_model_instance.sol[:nw][ALEAF_model_instance.cnw][decomp_group]
        solver_setting = ALEAF_model_instance.setting["Solver Setting"]

        solution = solve_model_GTEP!(JuMP_model, decomp_group, solution_list, solver_setting; iteration, PH_flag=false)

        return solution

    end

    if ALEAF_model_instance.data["model_type"]["model_type_EXP"] == "LP" 
        
        # build LCO GEP model instance (expansion mode)
        build_LCO_GTEP_expansion_mode_instance!(ALEAF_model_instance, decomp_group)
        @aleaf_info "[ALEAF LC_GTEP Expansion Model]: Build a GTEP model instance (expansion mode)"

        # Export ALEAF_model_instance model
        if ALEAF_setting["Simulation Setting"]["export_model_lp_expansion_flag"] == true
            output_path = parameter(ALEAF_model_instance, 0, :output_path)
            file_name = ALEAF_setting["Simulation Setting"]["model_lp_file_name_value"]
            file_name = string("Expansion_", file_name)
            JuMP.write_to_file(ALEAF_model_instance.model[:nw][ALEAF_model_instance.cnw][1], joinpath(output_path, file_name))
            @aleaf_info "[ALEAF LC_GTEP Expansion Model]: Export ALEAF_model_instance model (lp format)."
        end
        
        # Solve ALEAF LC_GEP Expansion Model
        output = get_solution(ALEAF_model_instance, decomp_group)
        @aleaf_info "[ALEAF LC_GTEP Expansion Model]: Solved."

        return output

    else

        if ALEAF_setting["Simulation Setting"]["MIP_relaxed_solution_bounds_flag"] == false # EXP model is MILP but we will solve it directly

            # build LCO GEP model instance (expansion mode)
            build_LCO_GTEP_expansion_mode_instance!(ALEAF_model_instance, decomp_group)
            @aleaf_info "[ALEAF LC_GTEP Expansion Model]: Build a GTEP model instance (expansion mode)"

            # Export ALEAF_model_instance model
            if ALEAF_setting["Simulation Setting"]["export_model_lp_expansion_flag"] == true
                output_path = parameter(ALEAF_model_instance, 0, :output_path)
                file_name = ALEAF_setting["Simulation Setting"]["model_lp_file_name_value"]
                file_name = string("Expansion_", file_name)
                JuMP.write_to_file(ALEAF_model_instance.model[:nw][ALEAF_model_instance.cnw][1], joinpath(output_path, file_name))
                @aleaf_info "[ALEAF LC_GTEP Expansion Model]: Export ALEAF_model_instance model (lp format)."
            end
                        
            # Solve ALEAF LC_GEP Expansion Model
            output = get_solution(ALEAF_model_instance, decomp_group)
            @aleaf_info "[ALEAF LC_GTEP Expansion Model]: Solved."

            return output

        else     # EXP model is MILP, and we will first solve LP to redefine investment decision bounds

            @aleaf_info "[ALEAF LC_GTEP Expansion Model]: === MILP bound reconfiguration process"
            
            # 1) build LP ALEAF model instance structure
            LP_ALEAF_model_instance = build_model_structure_GTEP(ALEAF_setting, network_data, Abstract_LC_GTEP_Model)
            LP_ALEAF_model_instance.ref = deepcopy(ALEAF_model_instance.ref)
            LP_ALEAF_model_instance.setting = deepcopy(ALEAF_model_instance.setting)
            @aleaf_info "[ALEAF LC_GTEP Expansion Model]: Create a LP GTEP expansion model structure and add reference data."
            
            # 2) reset integrality settings to false
            for i in [(i) for (i) in ALEAF.get_index(LP_ALEAF_model_instance, :gen_index, 0)]
                bus_idx = parameter(LP_ALEAF_model_instance, 0, :gen_index, "bus_idx", i)
                tech_idx = parameter(LP_ALEAF_model_instance, 0, :gen_index, "genco_tech_id", i)
                unit_group = parameter(LP_ALEAF_model_instance, 0, :gen_index, "UNIT_GROUP", i)

                if LP_ALEAF_model_instance.ref[:nw][bus_idx][:gen_bus][tech_idx]["Integrality"] == true
                    LP_ALEAF_model_instance.ref[:nw][bus_idx][:gen_bus][tech_idx]["Integrality"] = false
                end
            end

            # 3) build LP model optimization instance
            build_LCO_GTEP_expansion_mode_instance!(LP_ALEAF_model_instance, decomp_group)
            @aleaf_info "[ALEAF LC_GTEP Expansion Model]: Build a LP GTEP expansion model instance"

            # 4) solve LP expansion model
            LP_output = get_solution(LP_ALEAF_model_instance, decomp_group)
            LP_ALEAF_model_instance = 0.0   # delete
            @aleaf_info "[ALEAF LC_GTEP Expansion Model]: Solved a relaxed LP GTEP model."

            # 3) Collect expansion decisions 
            LP_expansion_solution = LP_output["solution"]["expansion"]

            # 4) Build original MILP model with updated bounds
            build_LCO_GTEP_expansion_mode_instance!(ALEAF_model_instance, decomp_group; ref_expansion_result=LP_expansion_solution)
            @aleaf_info "[ALEAF LC_GTEP Expansion Model]: Build the original GTEP model instance (expansion mode)"

            # Export ALEAF_model_instance model
            if ALEAF_setting["Simulation Setting"]["export_model_lp_expansion_flag"] == true
                output_path = parameter(ALEAF_model_instance, 0, :output_path)
                file_name = ALEAF_setting["Simulation Setting"]["model_lp_file_name_value"]
                file_name = string("Expansion_", file_name)
                JuMP.write_to_file(ALEAF_model_instance.model[:nw][ALEAF_model_instance.cnw][1], joinpath(output_path, file_name))
                @aleaf_info "[ALEAF LC_GTEP Expansion Model]: Export ALEAF_model_instance model (lp format)."
            end
                        
            # Solve ALEAF LC_GEP Expansion Model
            output = get_solution(ALEAF_model_instance, decomp_group)
            @aleaf_info "[ALEAF LC_GTEP Expansion Model]: Solved."

            return output

        end

    end

end


function execute_LC_GTEP_direct(ALEAF_setting::Dict{String,<:Any}, case_id::Int)
    
    # ALEAF_model_instance build-------------------------------------------------
    decomp_group = 1    # Only one group when PH is not being used.
    iteration = 1      # no iteration when PH is not being used.

    # generate network data (pick repdays on a coarse selection network first; fall back if empty)
    pre_groups, pre_days = (Dict{String,Any}(), Dict{String,Any}())
    coarse = select_repdays_coarse(ALEAF_setting, case_id, "expansion")
    coarse !== nothing && ((pre_groups, pre_days) = coarse)
    network_data = generate_networkdata_LC_GTEP(ALEAF_setting, case_id, "expansion"; precomputed_repday_groups=pre_groups, precomputed_repdays=pre_days)
    @aleaf_info "[ALEAF LC_GTEP Expansion Model]: Generate network data."

    # build common ALEAF model instance structure. transfer_heavy_data=false: network_data is reused by
    # build_and_solve below, so preserve exact current deepcopy behavior (OOM fix targets serial nodal OP).
    ALEAF_model_instance = build_model_structure_GTEP(ALEAF_setting, network_data, Abstract_LC_GTEP_Model)
    add_ref_LC_GTEP_model!(ALEAF_model_instance, ALEAF_setting; transfer_heavy_data=false)
    @aleaf_info "[ALEAF LC_GTEP Expansion Model]: Create a model structure and add reference data."

    # ---------------------------------------------------------------------------------
    # build and Solve LC-GTEP model instance (expansion mode)
    output = Dict{String, Any}()
    output[string(decomp_group)] = build_and_solve_LC_GTEP_expansion_model(ALEAF_setting, ALEAF_model_instance, decomp_group, network_data)

    # initialize solution dict
    if !haskey(ALEAF_model_instance.solution, string(iteration))
        ALEAF_model_instance.solution[string(iteration)] = Dict{Symbol,Any}(:nw => Dict{String,Any}())             
        ALEAF_model_instance.solution[string(iteration)][:nw][string(decomp_group)] = Dict{String,Any}()             
    end

    # Save solution
    if output[string(decomp_group)] == "No Solution"
        ALEAF_model_instance.solution[string(iteration)][:nw][string(decomp_group)] = "No Solution"
        ALEAF_model_instance.solution[string(iteration)][:nw][string(decomp_group)]["Objective Values"] = 0.0
    else
        ALEAF_model_instance.solution[string(iteration)][:nw][string(decomp_group)] = output[string(decomp_group)]
        ALEAF_model_instance.solution[string(iteration)][:nw][string(decomp_group)]["Objective Values"] = output[string(decomp_group)]["objective"]
    end


    # Get Dual Values
    JuMP_model = ALEAF_model_instance.model[:nw][ALEAF_model_instance.cnw][decomp_group]
    dual_status = JuMP.dual_status(JuMP_model)
    if dual_status == MOI.NO_SOLUTION 
        # fix discrete variables and resolve
        JuMP.fix_discrete_variables(JuMP_model)    # fix discrete variables
        JuMP.optimize!(JuMP_model)
    end
    # collect dual values
    get_dual_LC_GTEP_expansion_model!(ALEAF_model_instance, ALEAF_setting, network_data)
    @aleaf_info "[ALEAF LC_GTEP Expansion Model]: Obtain Dual."

    return ALEAF_model_instance
end


function generate_rep_days_LCO_GTEP(ALEAF_setting::Dict{String,<:Any}, network_data::Dict{String,<:Any}, rep_days, num_days_in_group, case_id, file_name; LC_GTEP_model_type::String="expansion", print_output_flag::Bool=true, GTEP_multi_round_info_data::Dict{Any,<:Any} = Dict{Any,Any}(), recorded_investment_decisions::Dict{String,<:Any} = Dict{String,Any}(), year_id::Int64=1)

    num_case_ids_list = [rep_days]
    pu_base = ALEAF_setting["Simulation Setting"]["per_unit_base_value"]

    # add fixed_extreme_day_list
    fixed_extreme_day_list = []
    if ALEAF_setting["Scenario Reduction Setting"]["preselected_extreme_days_list"] != "[]"
        fixed_extreme_day_list = parse.(Int, split(strip(ALEAF_setting["Scenario Reduction Setting"]["preselected_extreme_days_list"], ['[', ']']), ","))
    end

    type_of_data_set = []
    num_data_set = 0
    if ALEAF_setting["Scenario Reduction Setting"]["input_type_load_shape_flag"] == true 
        push!(type_of_data_set, "load_shape")
        num_data_set += 1 
    end
    if ALEAF_setting["Scenario Reduction Setting"]["input_type_load_MWh_flag"] == true 
        push!(type_of_data_set, "load_MWh") 
        num_data_set += 1
    end
    if ALEAF_setting["Scenario Reduction Setting"]["input_type_wind_shape_flag"] == true 
        push!(type_of_data_set, "wind_shape") 
        num_data_set += 1
    end
    if ALEAF_setting["Scenario Reduction Setting"]["input_type_wind_MWh_flag"] == true 
        push!(type_of_data_set, "wind_MWh") 
        num_data_set += 1
    end
    if ALEAF_setting["Scenario Reduction Setting"]["input_type_solar_shape_flag"] == true 
        push!(type_of_data_set, "solar_shape") 
        num_data_set += 1
    end
    if ALEAF_setting["Scenario Reduction Setting"]["input_type_solar_MWh_flag"] == true 
        push!(type_of_data_set, "solar_MWh") 
        num_data_set += 1
    end
    if ALEAF_setting["Scenario Reduction Setting"]["input_type_net_load_MWh_flag"] == true 
        push!(type_of_data_set, "net_load_MWh") 
        num_data_set += 1
    end

    # Extreme day set
    extreme_day_set = []
    if ALEAF_setting["Scenario Reduction Setting"]["fix_peak_demand_day_flag"] == true 
        push!(extreme_day_set, "peak_demand")
    end
    if ALEAF_setting["Scenario Reduction Setting"]["fix_peak_net_demand_day_flag"] == true 
        push!(extreme_day_set, "peak_net_demand")
    end
    if ALEAF_setting["Scenario Reduction Setting"]["fix_peak_solar_generation_day_flag"] == true 
        push!(extreme_day_set, "peak_pv_generation")
    end
    if ALEAF_setting["Scenario Reduction Setting"]["fix_peak_wind_generation_day_flag"] == true 
        push!(extreme_day_set, "peak_wind_generation")
    end
    if ALEAF_setting["Scenario Reduction Setting"]["fix_least_solar_generation_day_flag"] == true 
        push!(extreme_day_set, "least_pv_generation")
    end
    if ALEAF_setting["Scenario Reduction Setting"]["fix_least_wind_generation_day_flag"] == true 
        push!(extreme_day_set, "least_wind_generation")
    end

    # data path
    data_output_path = data_input_path = ALEAF_setting["data_location_timeseries"]

    # prepare data for the scenario reduction algorithm
    
    # Aggregate timeseries data for selected regions 

    # 1) list of selected regions
    region_map = get(network_data, "profile_data_region_map", nothing)

    # Precomputed per-data-column VRE capacity weights (Symbol col => MW). Supplied by the coarse
    # selection path (select_repdays_coarse) after re-aggregating run-resolution recorded decisions onto
    # data-columns, so Σ_c is preserved EXACTLY without needing the run-resolution bus keys here.
    col_weights = get(network_data, "__repday_col_weights", nothing)
    wind_col_weights = col_weights === nothing ? nothing : get(col_weights, "wind", nothing)
    pv_col_weights   = col_weights === nothing ? nothing : get(col_weights, "pv", nothing)

    selected_regions = []   # (region_name, original_peak_load_MW)
    for bus_id in keys(network_data["bus"])
        for region_id in keys(network_data["bus"][bus_id]["aggregation_info"]["original_load_(bus_i, MW)"])
            push!(selected_regions, (bus_id, region_id, network_data["bus"][bus_id]["aggregation_info"]["original_load_(bus_i, MW)"][region_id]))
        end
    end
    
    # 2) load 
    df_load = network_data["time_series_data"]["load"]
    if ALEAF_setting["Simulation Configuration"][string(case_id)]["Load_File_ID"] != "Base"
        data_path = joinpath(data_input_path, "0_additional_scenarios", "Load")
        new_file_name = string("timeseries_load_hourly_", ALEAF_setting["Simulation Configuration"][string(case_id)]["Load_File_ID"], ".csv")
        df_load = DataFrame(CSV.File(joinpath(data_path, new_file_name)))
    end
    
    load_value = zeros(8760)    # MWh
    load_shape = zeros(8760)    # Shape
    Peak_demand = 0
    num_regions_with_load = 0

    target_stage = "1"
    if (ALEAF_setting["Planning Design"]["multi_round_solution_process_flag"] == true) && (length(GTEP_multi_round_info_data) > 0)
        target_stage = string(GTEP_multi_round_info_data["round_ids_y"][1])
    elseif (LC_GTEP_model_type == "operation")
        target_stage = string(year_id)
    end

    # Group-by-distinct-column aggregation (distributive equivalence).
    # Tens of thousands of nodal bus-regions broadcast from only a handful of
    # distinct profile columns, so Σ_bus (colvec_bus * w_bus) == Σ_col colvec_col * (Σ_bus w_bus).
    # We accumulate per-distinct-column SCALAR weights first, then do one vector
    # multiply-add per distinct column. This is resolution-agnostic: at coarse
    # resolution distinct_cols≈bus count, so it degenerates to the original with no regression.
    W_load = Dict{Symbol,Float64}()   # value weight per column: Σ region_load*growth/pu
    S_load = Dict{Symbol,Int}()       # shape weight per column: # bus-regions on that col
    for region_id in keys(selected_regions)
        region_name = selected_regions[region_id][2]
        region_load = selected_regions[region_id][3]
        load_growth = get_load_growth_factor(network_data, target_stage, region_name)

        load_col = Symbol(profile_data_region(region_map, "load", region_name))
        W_load[load_col] = get(W_load, load_col, 0.0) + region_load * load_growth / pu_base
        S_load[load_col] = get(S_load, load_col, 0) + 1

        Peak_demand += region_load * load_growth / pu_base

        if region_load > 0
            num_regions_with_load += 1
        end
    end
    for (load_col, w) in W_load
        colvec = df_load[:, load_col]
        load_value += colvec * w
        load_shape += colvec * S_load[load_col]
    end
    load_shape = load_shape / num_regions_with_load

    # 3) WIND <- onshore wind only
    df_wind = network_data["time_series_data"]["wind_ons"]
    if ALEAF_setting["Simulation Configuration"][string(case_id)]["Wind_Ons_File_ID"] != "Base"
        data_path = joinpath(data_input_path, "0_additional_scenarios", "WIND")
        new_file_name = string("timeseries_wind_ons_hourly_", ALEAF_setting["Simulation Configuration"][string(case_id)]["Wind_Ons_File_ID"], ".csv")
        df_wind = DataFrame(CSV.File(joinpath(data_path, new_file_name)))
    end
    
    wind_on_value = zeros(8760) # MWh
    wind_on_shape = zeros(8760) # Shape
    Wind_Capacity = 0   # system-wide capacity (MW)
    num_regions_with_wind_potential = 0
    # Group-by-distinct-column (see load block above for the distributive rationale).
    W_wind = Dict{Symbol,Float64}()   # value weight per column: Σ local_wind_capacity
    S_wind = Dict{Symbol,Int}()       # shape weight per column: # bus-regions on that col
    for region_id in keys(selected_regions)
        region_idx = selected_regions[region_id][1]
        region_name = selected_regions[region_id][2]

        # local wind capacity (skipped when precomputed per-column weights are supplied)
        local_wind_capacity = 0.0
        if wind_col_weights !== nothing
            # W is injected below from the precomputed per-column weights
        elseif (ALEAF_setting["Planning Design"]["multi_round_solution_process_flag"] == true) && (length(GTEP_multi_round_info_data) > 0)
            local_wind_capacity = GTEP_multi_round_info_data["recorded_investment_decisions"]["Wind_TotalMW"][region_idx]
        elseif (LC_GTEP_model_type == "operation") && (length(recorded_investment_decisions) > 0)
            local_wind_capacity = recorded_investment_decisions["Wind_TotalMW"][region_idx]
        else
            for gen_idx in network_data["bus"][region_idx]["aggregation_info"]["local_gen_idx"]
                if network_data["plant"][gen_idx]["Profile_Type"] == "wind_ons"
                    local_wind_capacity += network_data["plant"][gen_idx]["CAP"]    # add existing capacity
                end
            end
        end

        wind_col = Symbol(profile_data_region(region_map, "wind_ons", region_name))
        W_wind[wind_col] = get(W_wind, wind_col, 0.0) + local_wind_capacity
        S_wind[wind_col] = get(S_wind, wind_col, 0) + 1
        Wind_Capacity += local_wind_capacity    # update the system-wide wind capacity
    end
    if wind_col_weights !== nothing
        empty!(W_wind)
        Wind_Capacity = 0.0
        for (col, w) in wind_col_weights
            W_wind[Symbol(col)] = get(W_wind, Symbol(col), 0.0) + w
            Wind_Capacity += w
        end
    end
    for (wind_col, w) in W_wind
        colvec = df_wind[:, wind_col]
        wind_on_value += colvec * w                 # local wind MWh
        wind_on_shape += colvec * S_wind[wind_col]  # system-wide wind shape
        # potential-count originally incremented once per bus-region with sum(col)>0
        if sum(colvec) > 0
            num_regions_with_wind_potential += S_wind[wind_col]
        end
    end
    wind_on_shape = wind_on_shape / num_regions_with_wind_potential   # system-wide average wind shape


    # 4) PV, No RTPV, No CSP
    df_pv = network_data["time_series_data"]["pv"]
    if ALEAF_setting["Simulation Configuration"][string(case_id)]["PV_File_ID"] != "Base"
        data_path = joinpath(data_input_path, "0_additional_scenarios", "PV")
        new_file_name = string("timeseries_pv_hourly_", ALEAF_setting["Simulation Configuration"][string(case_id)]["PV_File_ID"], ".csv")
        df_pv = DataFrame(CSV.File(joinpath(data_path, new_file_name)))
    end

    pv_value = zeros(8760)  # MWh
    pv_shape = zeros(8760)  # Shape
    PV_Capacity = 0 # system-wide capacity (MW)
    num_regions_with_pv_potential = 0
    # Group-by-distinct-column (see load block above for the distributive rationale).
    W_pv = Dict{Symbol,Float64}()   # value weight per column: Σ local_pv_capacity
    S_pv = Dict{Symbol,Int}()       # shape weight per column: # bus-regions on that col
    for region_id in keys(selected_regions)
        region_idx = selected_regions[region_id][1]
        region_name = selected_regions[region_id][2]

        # local pv capacity (skipped when precomputed per-column weights are supplied)
        local_pv_capacity = 0.0
        if pv_col_weights !== nothing
            # W is injected below from the precomputed per-column weights
        elseif (ALEAF_setting["Planning Design"]["multi_round_solution_process_flag"] == true) && (length(GTEP_multi_round_info_data) > 0)
            local_pv_capacity = GTEP_multi_round_info_data["recorded_investment_decisions"]["PV_TotalMW"][region_idx]
        elseif (LC_GTEP_model_type == "operation") && (length(recorded_investment_decisions) > 0)
            local_pv_capacity = recorded_investment_decisions["PV_TotalMW"][region_idx]
        else
            for gen_idx in network_data["bus"][region_idx]["aggregation_info"]["local_gen_idx"]
                if network_data["plant"][gen_idx]["Profile_Type"] == "pv"
                    local_pv_capacity += network_data["plant"][gen_idx]["CAP"]
                end
            end
        end

        pv_col = Symbol(profile_data_region(region_map, "pv", region_name))
        W_pv[pv_col] = get(W_pv, pv_col, 0.0) + local_pv_capacity
        S_pv[pv_col] = get(S_pv, pv_col, 0) + 1
        PV_Capacity += local_pv_capacity    # update the system-wide pv capacity
    end
    if pv_col_weights !== nothing
        empty!(W_pv)
        PV_Capacity = 0.0
        for (col, w) in pv_col_weights
            W_pv[Symbol(col)] = get(W_pv, Symbol(col), 0.0) + w
            PV_Capacity += w
        end
    end
    for (pv_col, w) in W_pv
        colvec = df_pv[:, pv_col]
        pv_value += colvec * w              # system-wide pv output (MWh)
        pv_shape += colvec * S_pv[pv_col]   # system-wide pv shape
        if sum(colvec) > 0
            num_regions_with_pv_potential += S_pv[pv_col]
        end
    end

    pv_shape = pv_shape / num_regions_with_pv_potential

    # Pass the six aggregated time-series to scenario reduction as in-memory DataFrames
    # (no transient CSV round-trip); keys/columns match what generate_input_data_hourly consumes.
    in_memory_timeseries = Dict{String,DataFrame}(
        "load_MW"     => DataFrame(Load = load_value),
        "load_shape"  => DataFrame(Load = load_shape),
        "wind_MW"     => DataFrame(Wind = wind_on_value),
        "wind_shape"  => DataFrame(Wind = wind_on_shape),
        "solar_MW"    => DataFrame(Solar = pv_value),
        "solar_shape" => DataFrame(Solar = pv_shape),
    )


    data_location = joinpath(pwd(), "data", ALEAF_setting["Simulation Setting"]["test_system_name"])

    timeSeries_df = CSV.read(joinpath(data_location, ALEAF_setting["File Path"]["timeseries_data_load_path"]), DataFrame)

    # run scenario reduction; aggregated time-series flow in memory (no transient CSV files)
    # The rep-day result df is returned and threaded back to the caller in memory.
    rep_day_input_df = run_scenario_reduction(ALEAF_setting, num_case_ids_list, type_of_data_set,
                                extreme_day_set, fixed_extreme_day_list,
                                data_location, data_input_path, data_output_path,
                                Wind_Capacity, PV_Capacity, Peak_demand,
                                num_days_in_group, print_output_flag, file_name;
                                in_memory_timeseries = in_memory_timeseries)

    # garbage collection
    GC.gc()

    return rep_day_input_df

end


function run_scenario_reduction(ALEAF_setting, num_case_ids_list, type_of_data_set,
                                extreme_day_set, fixed_extreme_day_list,
                                data_location, data_input_path, data_output_path,
                                Wind_Capacity, PV_Capacity, Peak_demand,
                                num_days_in_group, print_output_flag, file_name;
                                in_memory_timeseries::Union{Nothing,Dict{String,DataFrame}}=nothing)


    # Run the scenario reduction function; the rep-day result df is threaded back in memory.
    rep_day_input_df = run_scenario_reduction(
        in_memory_timeseries = in_memory_timeseries,
        time_resolution = ALEAF_setting["Scenario Reduction Setting"]["time_resolution"],
        num_scenarios_list = num_case_ids_list,
        fixing_extreme_days_flag = ALEAF_setting["Scenario Reduction Setting"]["fixing_extreme_days_flag"],
        generate_input_data_flag = ALEAF_setting["Scenario Reduction Setting"]["generate_input_data_flag"],
        allow_repday_overlap_flag = ALEAF_setting["Scenario Reduction Setting"]["allow_repday_overlap_flag"],
        type_of_data_set = type_of_data_set,
        extreme_day_set = extreme_day_set,
        fixed_extreme_day_list = fixed_extreme_day_list,
        data_location_timeseries = data_location,
        data_input_path = data_input_path,
        data_output_path = data_output_path,
        windCapacity = Wind_Capacity,
        solarCapacity = PV_Capacity,
        peakDemand = Peak_demand,
        num_days_in_group = num_days_in_group,
        print_output = print_output_flag,
        file_name = file_name
    )

    return rep_day_input_df

end


"""
Re-aggregate RUN-resolution VRE capacity totals (Wind_TotalMW / PV_TotalMW / RTPV_TotalMW, keyed by run
bus_idx) onto per-data-column weights, using each run-bus's finest-region set and the finest->data-column
map. This reproduces exactly what generate_rep_days_LCO_GTEP would compute at run resolution: for a run-bus
b with finest regions R(b), each r in R(b) contributes totals[b] to column map(r). Because Σ_c is a pure
sum over data-columns, the coarse selection then yields byte-identical rep days.

`run_bus_dict` maps run bus_idx (String) => bus data with aggregation_info["aggregated_regions_bus_i"].
Returns Dict("wind"=>Dict{String,Float64}, "pv"=>..., "rtpv"=>...) keyed by data-column string.
"""
function build_repday_selection_col_weights(run_bus_dict::AbstractDict, region_map, wind_totals::AbstractDict, pv_totals::AbstractDict, rtpv_totals::AbstractDict)
    out = Dict{String,Any}("wind" => Dict{String,Float64}(), "pv" => Dict{String,Float64}(), "rtpv" => Dict{String,Float64}())
    for (bus_idx, bus_data) in run_bus_dict
        finest_regions = bus_data["aggregation_info"]["aggregated_regions_bus_i"]
        wt = get(wind_totals, string(bus_idx), 0.0)
        pt = get(pv_totals, string(bus_idx), 0.0)
        rt = get(rtpv_totals, string(bus_idx), 0.0)
        (wt == 0.0 && pt == 0.0 && rt == 0.0) && continue
        for r in finest_regions
            if wt != 0.0
                c = string(profile_data_region(region_map, "wind_ons", r))
                out["wind"][c] = get(out["wind"], c, 0.0) + wt
            end
            if pt != 0.0
                c = string(profile_data_region(region_map, "pv", r))
                out["pv"][c] = get(out["pv"], c, 0.0) + pt
            end
            if rt != 0.0
                c = string(profile_data_region(region_map, "rtpv", r))
                out["rtpv"][c] = get(out["rtpv"], c, 0.0) + rt
            end
        end
    end
    return out
end


"""
Ensure the resolution vocabulary (network_resolution_level, Network Setting, sub_area_mapping) is present
in ALEAF_setting without running the full network build, so resolve_repday_selection_resolution is
reliable up front. No-op when already loaded or when no config workbook is configured.
"""
function ensure_resolution_metadata_loaded!(ALEAF_setting::Dict{String,<:Any}, case_id::Int)
    if !isempty(get(ALEAF_setting, "network_resolution_level", Dict{String,Any}())) &&
       !isempty(get(ALEAF_setting, "sub_area_mapping", Dict{String,Any}())) &&
       haskey(get(ALEAF_setting, "Network Setting", Dict{String,Any}()), "network_boundary_type")
        return nothing
    end
    resolve_network_input_files!(ALEAF_setting, case_id)
    config_file_location = get(ALEAF_setting, "network_config_file_location", nothing)
    (config_file_location === nothing || !isfile(config_file_location)) && return nothing
    XLSX.openxlsx(config_file_location) do config_xf
        ALEAF_setting["Network Setting"] = read_setting_xlsx_return_dict_string_any(config_xf, "Network Setting")
        ALEAF_setting["network_resolution_level"] = read_xlsx_return_dict_string_any(config_xf, "network_resolution_level"; first_row_value = 2)
        ALEAF_setting["sub_area_mapping"] = read_xlsx_return_dict_string_any(config_xf, "sub_area_mapping"; first_row_value = 2)
    end
    return nothing
end


"""
Pick representative days on a lightweight coarse SELECTION network (bounded by the repday selection
resolution) instead of the full run network, then return the day assignment for the run build to reuse.

Returns (repday_groups, repdays_meta) via extract_repday_assignment_LC_GTEP, or `nothing` when no coarse
assignment could be built (caller then falls back to the in-line selection). VRE capacity weights from
recorded_investment_decisions/GTEP_multi_round_info are re-aggregated to data-columns via
`run_bus_dict`/`region_map` so results are byte-identical; when those are absent the coarse network's own
(zero, config-path) weights match the run build.
"""
function select_repdays_coarse(ALEAF_setting::Dict{String,<:Any}, case_id::Int, model_type::String; year_id::Int64=1, recorded_investment_decisions::Dict{String,<:Any} = Dict{String,Any}(), GTEP_multi_round_info_data::Dict{Any,<:Any} = Dict{Any,Any}(), run_bus_dict = nothing, region_map = nothing)

    # On a successful coarse selection ("fired"), report the representative-day selection summary
    # (# day groups, # representative days). The run build reuses these precomputed days and would
    # otherwise skip the scenario-reduction count log, so this is the only place it surfaces. On
    # fallback the in-line selection logs its own counts, so nothing is emitted here. The coarse
    # resolution / fire-vs-fallback detail is kept for debugging but only at "detailed" logging level.
    _detailed_log = lowercase(string(get(get(ALEAF_setting, "Simulation Setting", Dict{String,Any}()), "logging_level_value", "simple"))) == "detailed"
    _emit(assignment, res, buses) = begin
        if _detailed_log
            @aleaf_info "[ALEAF LC_GTEP]: coarse repday selection resolution=$(res === nothing ? "none" : res), $(assignment === nothing ? "fell back to in-line selection" : "fired")$(buses === nothing ? "" : ", collapsed buses=$(buses)")"
        end
        if assignment !== nothing
            repday_groups, repdays_meta = assignment
            @aleaf_info "[ALEAF LC_GTEP]: representative day selection: $(length(repday_groups)) day group(s), $(length(repdays_meta)) representative day(s)"
        end
    end

    # Manual mode: the repday file drives selection; skip coarse selection entirely.
    mode_key = model_type == "expansion" ? "Repday_Selection_Mode_EXP" :
               model_type == "operation" ? "Repday_Selection_Mode_OP" : "Repday_Selection_Mode_RA"
    if string(get(ALEAF_setting["Simulation Configuration"][string(case_id)], mode_key, "")) == "Manual"
        _emit(nothing, nothing, nothing)
        return nothing
    end

    # The resolver needs the resolution vocabulary (network_resolution_level / Network Setting /
    # sub_area_mapping), which is normally populated by the full network build. Load it cheaply up front
    # so the resolver is reliable before any build runs.
    ensure_resolution_metadata_loaded!(ALEAF_setting, case_id)

    # nothing => key absent/blank/unrecognized => use today's in-line selection at the run resolution.
    override = resolve_repday_selection_resolution(ALEAF_setting, case_id)
    if override === nothing
        _emit(nothing, nothing, nothing)
        return nothing
    end

    # Re-key VRE capacity totals onto data-columns when decisions are keyed by run bus_idx.
    col_weights = Dict{String,Any}()
    rid = !isempty(GTEP_multi_round_info_data) ? get(GTEP_multi_round_info_data, "recorded_investment_decisions", Dict{String,Any}()) : recorded_investment_decisions
    nonzero_total(d) = any(v -> v != 0.0, values(get(rid, d, Dict{String,Any}())))
    has_vre = nonzero_total("Wind_TotalMW") || nonzero_total("PV_TotalMW") || nonzero_total("RTPV_TotalMW")
    if has_vre
        if haskey(rid, "__col_weights")
            # Precomputed at the recording site (run reference available there) => exact + cheap.
            col_weights = rid["__col_weights"]
        elseif run_bus_dict !== nothing
            col_weights = build_repday_selection_col_weights(run_bus_dict, region_map,
                get(rid, "Wind_TotalMW", Dict{String,Any}()),
                get(rid, "PV_TotalMW", Dict{String,Any}()),
                get(rid, "RTPV_TotalMW", Dict{String,Any}()))
        else
            # Cannot safely re-key run-resolution VRE weights => fall back to in-line selection.
            _emit(nothing, override, nothing)
            return nothing
        end
    end

    selection_network = generate_networkdata_LC_GTEP(ALEAF_setting, case_id, model_type;
        year_id=year_id, GTEP_multi_round_info_data=GTEP_multi_round_info_data,
        selection_only=true, aggregation_resolution_override=override,
        repday_col_weights=col_weights, print_output_flag=false)

    _sel_buses = haskey(selection_network, "bus") ? length(selection_network["bus"]) : nothing
    if !haskey(selection_network, "repday_groups") || isempty(selection_network["repday_groups"])
        _emit(nothing, override, _sel_buses)
        return nothing
    end
    _assignment = extract_repday_assignment_LC_GTEP(selection_network)
    _emit(_assignment, override, _sel_buses)
    return _assignment
end


function generate_networkdata_LC_GTEP(ALEAF_setting::Dict{String,<:Any}, case_id::Int, LC_GTEP_model_type::String; PH_flag::Bool=false, decomp_id::Int=0, print_output_flag::Bool=true, recorded_investment_decisions::Dict{String,<:Any} = Dict{String,Any}(), GTEP_multi_round_info_data::Dict{Any,<:Any} = Dict{Any,Any}(), year_id::Int64=1, target_day_group_id::Int = 0, precomputed_repday_groups::AbstractDict = Dict{String,Any}(), precomputed_repdays::AbstractDict = Dict{String,Any}(), build_repday_hourly_data::Bool = true, selection_only::Bool = false, aggregation_resolution_override = nothing, repday_col_weights::AbstractDict = Dict{String,Any}())

    # Generate network_data for a given case_id
    network_data = Dict{String, Any}()
    network_data["case_name"] = ALEAF_setting["Simulation Configuration"][string(case_id)]["Case_ID"]
    network_data["case_id"] = string(case_id)

    # Selection-only lightweight build: coarsen the regional aggregation and precomputed VRE weights
    # so repday selection runs on a tiny system aggregate (see select_repdays_coarse).
    if aggregation_resolution_override !== nothing
        network_data["__selection_agg_override"] = aggregation_resolution_override
    end
    if !isempty(repday_col_weights)
        network_data["__repday_col_weights"] = repday_col_weights
    end

    resolve_network_input_files!(ALEAF_setting, case_id)
    data_location = ALEAF_setting["data_location"]
    load_file_path_sheet!(ALEAF_setting, ALEAF_setting["network_data_file_location"])
   
    # Get simulation round information 
    network_data["simulation_round_idx"] = 1 # default value is 1
    if !isempty(GTEP_multi_round_info_data)
        network_data["simulation_round_idx"] = GTEP_multi_round_info_data["round_idx"]
    end

    # Get network data: plant, bus, branch
    # Topology is year/round-invariant; cache it so multi-round reruns skip the rebuild.
    get_network_data_cached!(network_data, ALEAF_setting)
   
    network_data["annual_load_growth_by_region"] = Dict{String,Any}()

    # Prepare multi-year optimization
    network_data["planning_stages"] = Dict{String,Any}()
    for id in 1:ALEAF_setting["Planning Design"]["num_stages_value"]
        network_data["planning_stages"][string(id)] = Dict()
        network_data["planning_stages"][string(id)]["stage_id"] = id
        network_data["planning_stages"][string(id)]["year"] = ALEAF_setting["Planning Design"]["first_stage_year_value"] + ((id-1) * ALEAF_setting["Planning Design"]["num_years_per_stage_value"])
        network_data["planning_stages"][string(id)]["stage_length"] = ALEAF_setting["Planning Design"]["num_years_per_stage_value"]
        network_data["planning_stages"][string(id)]["load_growth_by_region"] = Dict{String,Float64}()
        network_data["planning_stages"][string(id)]["repdays"] = Dict{String,Any}()
    end
    build_load_growth_by_region_LCO_GTEP!(ALEAF_setting, network_data, case_id)

    # Selection-only early return: run ONLY repday selection on the coarse aggregate, then stop.
    # Skips zonal data, branch cost, gen/storage/ATB/hybrid tech, per-unit, raw material, and the
    # hourly repday ["data"] materialization. Needs only bus aggregation_info, planning_stages
    # load_growth, profile_data_region_map, and the load/wind_ons/pv time series.
    if selection_only
        build_profile_data_region_maps!(network_data, ALEAF_setting)
        get_timeseries_data!(ALEAF_setting, network_data, case_id)
        get_repday_groups_data!(network_data, ALEAF_setting, case_id, LC_GTEP_model_type, year_id; print_output_flag, GTEP_multi_round_info_data, recorded_investment_decisions, build_hourly_data=false)
        network_data["output_path"] = define_output_path(ALEAF_setting, case_id)
        return network_data
    end

    # Define zonal settings (reserve, policy, supply curve, cost scaling, capacity credits, inertia)
    get_network_zonal_data!(network_data, ALEAF_setting, case_id)
    
    # Get annualized line investment cost
    update_branch_annual_cost_data!(ALEAF_setting, network_data)
    
    # Define Gen Technology
    network_data["gen_technology"] = deepcopy(ALEAF_setting["Gen Technology"])

    # Update Storage Technology with ESGC data
    update_storage_technology_data!(ALEAF_setting, network_data, case_id)

    # Update Gen Technology with ATB data
    update_gen_technology_data!(data_location, ALEAF_setting, network_data, case_id)

    # Update Hybrid Technology data
    update_hybrid_plant_technology_data!(network_data) 

    # Update Hybrid Technology data to LFL 
    add_hybrid_plant_to_LFL!(network_data) 

    # Check model_type (MILP, LP)
    determine_opt_model_type!(network_data, ALEAF_setting, case_id, LC_GTEP_model_type)

    # Get Time Series Data
    get_timeseries_data!(ALEAF_setting, network_data, case_id)

    # Define representative day groups and update the planning stages.
    # When target_day_group_id/precomputed_* are supplied (distributed OP worker), only the target
    # day-group is materialized instead of the full horizon.
    get_repday_groups_data!(network_data, ALEAF_setting, case_id, LC_GTEP_model_type, year_id; print_output_flag, GTEP_multi_round_info_data, recorded_investment_decisions, target_day_group_id, precomputed_repday_groups, precomputed_repdays, build_hourly_data=build_repday_hourly_data)

    # Update system-wide peak demand
    update_PD_value!(ALEAF_setting, network_data, case_id)

    # Apply the per-unit system 
    apply_per_unit_to_network_data!(ALEAF_setting, network_data, case_id)

    # Add raw material data
    if LC_GTEP_model_type == "expansion" 
        add_raw_material_data_to_gen_tech!(ALEAF_setting, network_data)
    end

    # add operation year
    if (LC_GTEP_model_type == "operation")
        network_data["operation_year"] = year_id
    end

    # Define output path
    network_data["output_path"] = define_output_path(ALEAF_setting, case_id)

    #------------------------------------------

    return network_data

end


function get_modeled_load_regions(network_data::Dict{String,<:Any})

    modeled_regions = Set{String}()
    for bus_data in values(network_data["bus"])
        for region_id in keys(bus_data["aggregation_info"]["original_load_(bus_i, MW)"])
            push!(modeled_regions, region_id)
        end
    end

    return sort!(collect(modeled_regions))

end


function resolve_load_growth_file_path(data_location::AbstractString, file_id::AbstractString)

    if isempty(strip(file_id))
        error("Simulation Configuration `load_increase_rate_file_ID` must be provided when `load_increase_rate_mode` is regional.")
    end

    data_path = joinpath(data_location, "timeseries_data_files", "0_additional_scenarios", "Load")
    candidate_paths = [
        joinpath(data_path, string("timeseries_load_growth_", file_id, ".csv")),
        joinpath(data_path, string("timeseries_load_growth_rate_", file_id, ".csv")),
    ]

    for file_path in candidate_paths
        if isfile(file_path)
            return file_path
        end
    end

    error("Load growth file not found for `load_increase_rate_file_ID=$(file_id)`. Expected one of: $(join(candidate_paths, ", ")).")

end


function build_load_growth_by_region_LCO_GTEP!(ALEAF_setting::Dict{String,<:Any}, network_data::Dict{String,<:Any}, case_id)

    sim_config = ALEAF_setting["Simulation Configuration"][string(case_id)]
    raw_mode = string(sim_config["load_increase_rate_mode"])
    mode = lowercase(strip(raw_mode))
    valid_modes = Set(["systemwide", "regional"])
    if !(mode in valid_modes)
        error("Invalid `load_increase_rate_mode=$(raw_mode)`. Allowed values are `systemwide` and `regional`.")
    end

    modeled_regions = get_modeled_load_regions(network_data)
    current_year = ALEAF_setting["Planning Design"]["base_year_value"]
    start_year = ALEAF_setting["Planning Design"]["first_stage_year_value"]
    total_report_years = ALEAF_setting["Planning Design"]["num_stages_value"] * ALEAF_setting["Planning Design"]["num_years_per_stage_value"]
    annual_report_years = collect(start_year:(start_year + total_report_years - 1))

    if mode == "systemwide"
        annual_growth = Dict{String,Any}()
        for year in annual_report_years
            load_growth = (1 + sim_config["load_increase_rate_value"]) ^ (year - current_year)
            annual_growth[string(year)] = Dict(region_id => load_growth for region_id in modeled_regions)
        end
        network_data["annual_load_growth_by_region"] = annual_growth
        for stage_data in values(network_data["planning_stages"])
            year = stage_data["year"]
            stage_data["load_growth_by_region"] = deepcopy(annual_growth[string(year)])
        end
        return
    end

    file_id = string(get(sim_config, "load_increase_rate_file_ID", ""))
    file_path = resolve_load_growth_file_path(ALEAF_setting["data_location"], file_id)
    growth_df = CSV.read(file_path, DataFrame)

    required_columns = ["Year", "Region"]
    growth_column_names = Set(string(col) for col in names(growth_df))
    missing_columns = filter(col -> !(col in growth_column_names), required_columns)
    if !isempty(missing_columns)
        error("Load growth file `$(file_path)` is missing required columns: $(join(missing_columns, ", ")).")
    end

    growth_value_column = if "Load Growth Factor" in growth_column_names
        "Load Growth Factor"
    elseif "Load Growth Rate" in growth_column_names
        "Load Growth Rate"
    else
        error("Load growth file `$(file_path)` is missing the required growth column. Expected `Load Growth Factor` or `Load Growth Rate`.")
    end

    growth_lookup = Dict{Tuple{Int, String}, Float64}()
    for row in eachrow(growth_df)
        year = Int(row["Year"])
        region_id = strip(string(row["Region"]))
        growth_value = Float64(row[growth_value_column])
        lookup_key = (year, region_id)
        if haskey(growth_lookup, lookup_key)
            error("Duplicate load growth entry found in `$(file_path)` for Year=$(year), Region=$(region_id).")
        end
        growth_lookup[lookup_key] = growth_value
    end

    # Resolve each modeled (finest) region to its load-growth data-region before indexing the CSV,
    # so a nodal/county DB reads one growth series per coarser data-region. Absent config => identity.
    build_profile_data_region_maps!(network_data, ALEAF_setting)
    region_map = get(network_data, "profile_data_region_map", nothing)

    annual_growth = Dict{String,Any}()
    for year in annual_report_years
        year_growth = Dict{String,Float64}()
        for region_id in modeled_regions
            data_region_id = strip(string(profile_data_region(region_map, "load_growth", region_id)))
            lookup_key = (year, data_region_id)
            if !haskey(growth_lookup, lookup_key)
                error("Missing load growth entry in `$(file_path)` for Year=$(year), Region=$(data_region_id) (modeled region $(region_id)).")
            end
            base_key = (current_year, data_region_id)
            if !haskey(growth_lookup, base_key)
                error("Missing current-year load growth entry in `$(file_path)` for Year=$(current_year), Region=$(data_region_id) (modeled region $(region_id)).")
            end
            # Normalize CSV factors to current_year so regional matches systemwide semantics (1.0 at run year)
            # and does not double-count growth when the bus load is already current-year data.
            # Store keyed by the finest region id so per-bus callers (get_load_growth_factor) are unchanged.
            # Guard a zero (no-data) base factor: regions with an all-zero growth series (e.g. the zeroed
            # Canadian/Mexican BAs in the shared NA growth file) would otherwise give 0/0 = NaN here and
            # poison every bus demand at that region. Treat a zero base as flat (no) growth = 1.0.
            base_val = growth_lookup[base_key]
            year_growth[region_id] = base_val == 0.0 ? 1.0 : growth_lookup[lookup_key] / base_val
        end
        annual_growth[string(year)] = year_growth
    end
    network_data["annual_load_growth_by_region"] = annual_growth

    for stage_data in values(network_data["planning_stages"])
        year = stage_data["year"]
        stage_data["load_growth_by_region"] = deepcopy(annual_growth[string(year)])
    end

end


function replace_stochastic_time_series_data!(ALEAF_setting, network_data, case_id, data_location)
        
    # Update Stochastic Scenarios
    if ALEAF_setting["Simulation Configuration"][string(case_id)]["Stochastic_Expansion_Flag"] == true

        for sce_id in 1:ALEAF_setting["Simulation Configuration"][string(case_id)]["Num_Sto_Scenarios"]
        
            if ALEAF_setting["Simulation Configuration"][string(case_id)]["Stochastic_Wind_Ons_Flag"] == true
                
                # read new time-series file
                data_path = joinpath(data_location, "timeseries_data_files", "0_additional_scenarios", "WIND")
                file_name = string("timeseries_wind_ons_hourly_", ALEAF_setting["Simulation Configuration"][string(case_id)]["Stochastic_File_ID"], "_", sce_id, ".csv")
                timeSeries_df = CSV.read(joinpath(data_path, file_name), DataFrame)

                # filter and update the network data (use first present repday key, not hardcoded "1")
                if haskey(network_data["repdays"][first(keys(network_data["repdays"]))]["data"]["1"], "wind_ons_BA")
                    select_and_update_stochastic_timeseries_data!(timeSeries_df, ALEAF_setting, network_data, "wind_ons_BA", sce_id)
                else
                    select_and_update_stochastic_timeseries_data!(timeSeries_df, ALEAF_setting, network_data, "wind_ons", sce_id)
                end

            end

            if ALEAF_setting["Simulation Configuration"][string(case_id)]["Stochastic_PV_Flag"] == true
                
                # read new time-series file
                data_path = joinpath(data_location, "timeseries_data_files", "0_additional_scenarios", "PV")
                file_name = string("timeseries_pv_hourly_", ALEAF_setting["Simulation Configuration"][string(case_id)]["Stochastic_File_ID"], "_", sce_id, ".csv")
                timeSeries_df = CSV.read(joinpath(data_path, file_name), DataFrame)

                # filter and update the network data
                select_and_update_stochastic_timeseries_data!(timeSeries_df, ALEAF_setting, network_data, "pv", sce_id)

            end

            if ALEAF_setting["Simulation Configuration"][string(case_id)]["Stochastic_Load_Flag"] == true
                
                # read new time-series file
                data_path = joinpath(data_location, "timeseries_data_files", "0_additional_scenarios", "Load")
                file_name = string("timeseries_load_hourly_", ALEAF_setting["Simulation Configuration"][string(case_id)]["Stochastic_File_ID"], "_", sce_id, ".csv")
                timeSeries_df = CSV.read(joinpath(data_path, file_name), DataFrame)

                # filter and update the network data
                select_and_update_stochastic_timeseries_data!(timeSeries_df, ALEAF_setting, network_data, "load", sce_id)

            end
        end
    end

end


function add_raw_material_data_to_gen_tech!(ALEAF_setting::Dict{String,<:Any}, network_data::Dict{String,<:Any})

    network_data["raw_materials"] = deepcopy(ALEAF_setting["Raw Materials"])

    function find_raw_material_lists(ALEAF_setting, tech_id::String)
        for list_idx in keys(ALEAF_setting["Gen Technology Raw Materials"])
            if ALEAF_setting["Gen Technology Raw Materials"][list_idx]["Tech_ID"] == tech_id
                return ALEAF_setting["Gen Technology Raw Materials"][list_idx]
            end
        end
    end
    
    for tech in keys(network_data["gen_technology"])
        if network_data["gen_technology"][tech]["Material_Flag"] == true
            network_data["gen_technology"][tech]["GenTech_Raw_Materials"] = find_raw_material_lists(ALEAF_setting, network_data["gen_technology"][tech]["Tech_ID"])
        end
    end

end


function apply_per_unit_to_network_data!(ALEAF_setting::Dict{String,<:Any}, network_data::Dict{String,<:Any}, case_id)

    pu_power_base = ALEAF_setting["Simulation Setting"]["per_unit_base_value"]
    pu_econ_base = ALEAF_setting["Simulation Setting"]["per_unit_econ_base_value"] / pu_power_base    # per unit economic base
    
    # bus
    for (bus_id, bus_data) in network_data["bus"]
       
        bus_data["MW load"] /= pu_power_base

        for idx in keys(bus_data["aggregation_info"]["original_load_(bus_i, MW)"])
            bus_data["aggregation_info"]["original_load_(bus_i, MW)"][idx] /= pu_power_base
        end
    end

    # plant
    for (plant_id, plant_data) in network_data["plant"]
        plant_data["CAP"] /= pu_power_base
        plant_data["Charge_CAP"] /= pu_power_base
        plant_data["ES_MWh"] /= pu_power_base

        if get(plant_data, "bypass_parameter_check", false)
            for field in ["NLC", "SUC", "SDC", "DECC"]
                if haskey(plant_data, field)
                    plant_data[field] /= pu_econ_base
                end
            end

            if ALEAF_setting["Planning Design"]["reserve_cost_type_flag"] == "absolute"
                for field in ["reg_cost", "spin_cost", "nspin_cost", "flex_cost"]
                    if haskey(plant_data, field)
                        plant_data[field] /= pu_econ_base
                    end
                end
            end

            if haskey(plant_data, "AET")
                plant_data["AET"] /= pu_power_base
            end

            for field in ["Emission_CO2", "Emission_1", "Emission_2", "Emission_3"]
                if haskey(plant_data, field)
                    plant_data[field] *= pu_power_base
                end
            end

            for field in ["Annual_STO_INV", "Annual_INVC", "Annual_FOM", "Annual_VOM", "Annual_FC", "Annual_MC"]
                if haskey(plant_data, field)
                    for year_id in keys(plant_data[field])
                        plant_data[field][year_id] /= pu_econ_base
                    end
                end
            end
        end
    end
   
    # branch
    for (branch_id, branch_data) in network_data["branch"]
        branch_data["rate_a"] /= pu_power_base
        branch_data["max_rate_a"] /= pu_power_base
        branch_data["transmission_expansion_cost"] /= pu_econ_base
        branch_data["transmission_fom_cost"] /= pu_econ_base
    end

    # repdays
    services = ["flex_up", "flex_down", "reg_up", "reg_down", "spin", "nspin"]

    for (stage, stage_data) in network_data["planning_stages"]
        for (day, day_data) in stage_data["repdays"]
            if haskey(day_data, "data") 
                for (hour, hour_data) in day_data["data"]
                    for (time, time_data) in hour_data
                        for service in services
                            for (key, value) in time_data[service]
                                time_data[service][key] = value / pu_power_base
                            end
                        end
                    end
                end
            end
        end
    end

    # Zone
    for (zone_id, zone_data) in network_data["zone"]["resource_supply_curve"]
        if !(zone_id in ["data_resolution_type", "supply_curve_resource_list", "supply_curve_region_list"])
            for (key, value) in zone_data["data"]
                zone_data["data"][key]["Resource_Limit_MW"] /= pu_power_base
            end
        end
    end

    # Gen_Technology 
    function apply_pu_to_gen_tech(gen_tech_data)
        # Scalar fields
        for field in ["NLC", "SUC", "SDC", "DECC"]
            gen_tech_data[field] /= pu_econ_base
        end

        for field in ["CAP", "Charge_CAP", "AET"]
            gen_tech_data[field] /= pu_power_base
        end

        for field in ["Emission_CO2", "Emission_1", "Emission_2", "Emission_3"]
            gen_tech_data[field] *= pu_power_base
        end

        # Annual cost-related fields
        for year_id in keys(gen_tech_data["Annual_STO_INV"])
            for field in ["Annual_STO_INV", "Annual_INVC", "Annual_FOM", "Annual_VOM", "Annual_FC", "Annual_MC"]
                gen_tech_data[field][year_id] /= pu_econ_base
            end
        end

        # PTC values, excluding specific keys
        for (idx, val) in gen_tech_data["PTC"]
            if idx != "UNITGROUP" && idx != "Tech_ID"
                gen_tech_data["PTC"][idx] = val / pu_econ_base
            end
        end

        # reserve-related fields
        if ALEAF_setting["Planning Design"]["reserve_cost_type_flag"] == "absolute"
            for field in ["reg_cost", "spin_cost", "nspin_cost", "flex_cost"]
                gen_tech_data[field] /= pu_econ_base
            end
        end

    end

    for (gen_tech_id, gen_tech_data) in network_data["gen_technology"]
        apply_pu_to_gen_tech(gen_tech_data)
    end

    # hybrid  
    for (hybrid_id, hybrid_data) in network_data["hybrid"]
        hybrid_data["CAP"] /= pu_power_base
        hybrid_data["INTERCON_LIM"] /= pu_power_base

        if hybrid_data["Hybrid_Gen"] != "NA"
            apply_pu_to_gen_tech(hybrid_data["Component_gen_tech_data"])
        end

        if hybrid_data["Hybrid_ES"] != "NA"
            apply_pu_to_gen_tech(hybrid_data["Component_ES_tech_data"])
        end
    end

    # demand 
    for (demand_id, demand_data) in network_data["demand"]
        demand_data["CAP"] /= pu_power_base
        demand_data["INTERCON_LIM"] /= pu_power_base
        demand_data["Daily_DR_Limit_MWh"] /= pu_power_base

        demand_data["Price_1"] /= pu_econ_base    
        demand_data["Price_2"] /= pu_econ_base
        demand_data["Price_3"] /= pu_econ_base
        demand_data["Price_4"] /= pu_econ_base
        demand_data["Price_5"] /= pu_econ_base    

        if demand_data["Hybrid_Gen"] != "NA"
            apply_pu_to_gen_tech(demand_data["Component_gen_tech_data"])
        end

        if demand_data["Hybrid_ES"] != "NA"
            apply_pu_to_gen_tech(demand_data["Component_ES_tech_data"])
        end
    end
    
end


function update_PD_value!(ALEAF_setting::Dict{String,<:Any}, network_data::Dict{String,<:Any}, case_id)

    PD = 0
    for bus_id in keys(network_data["bus"])
        PD += network_data["bus"][bus_id]["MW load"]
    end

    ALEAF_setting["Simulation Configuration"][string(case_id)]["PD"] = PD

end


function aggregate_vre_timeseries_zdt_GTEP!(am::Abstract_ALEAF_Model)

    num_stages = am.setting["Planning Design"]["num_stages_value"]

    if haskey(am.ref[:nw][0], :operation_year)
        ids_y = [am.ref[:nw][0][:operation_year]]
    else
        ids_y = [i for i in 1:num_stages]
    end

    am.ref[:nw][0][:vre_aggregated_data_zdt] = Dict{Int, Dict{Int, Array{Float64,3}}}()

    ids_zone = [(z) for (z) in get_index(am, :bus, 0)]
    ids_d = [(d) for (d) in get_index(am, :repdays, 0)]
    ids_t = [(t) for (t) in am.setting["run_H"]]

    # Stable position maps: repday-id -> dim1, hour-id -> dim2. sort(ids_d) / ids_t order is the contract.
    d_pos = Dict{Int,Int}(d => i for (i, d) in enumerate(sort(ids_d)))
    h_pos = Dict{Int,Int}(t => i for (i, t) in enumerate(ids_t))
    am.ref[:nw][0][:vre_zdt_index] = (d = d_pos, h = h_pos)
    nd = length(ids_d)
    nh = length(ids_t)

    # Only regions that carry the resource at some hour belong in the zone average; a region that is
    # zero (or absent) across every rep-day hour has no resource and must stay out of the denominator.
    region_masks = Dict{Int64, Dict{Int64, Dict{String, Vector{String}}}}()

    function build_region_masks!(am::Abstract_ALEAF_Model, data_labels::Vector{String}, year::Int64)

        repday_data = am.ref[:nw][0][:planning_stages][year]["repdays"]
        resource_bearing = Dict{String, Set{String}}(label => Set{String}() for label in data_labels)

        for d in ids_d
            hourly_data = repday_data[string(d)]["data"]
            for h_key in keys(hourly_data)
                sub_period_data = hourly_data[h_key]["1"]
                for label in data_labels
                    haskey(sub_period_data, label) || continue
                    bearing_regions = resource_bearing[label]
                    for (region_id, value) in sub_period_data[label]
                        if value != 0
                            push!(bearing_regions, region_id)
                        end
                    end
                end
            end
        end

        region_map = get(am.ref[:nw][0], :profile_data_region_map, nothing)
        year_masks = get!(region_masks, year, Dict{Int64, Dict{String, Vector{String}}}())
        for zone in ids_zone
            region_list = am.ref[:nw][0][:bus][zone]["aggregation_info"]["aggregated_regions_bus_i"]
            zone_masks = get!(year_masks, zone, Dict{String, Vector{String}}())
            for label in data_labels
                bearing_regions = resource_bearing[label]
                # Lever A: keep the finest region id as the mask entry, but test resource-bearing at its data-region.
                profile_type = replace(label, "_BA" => "")
                zone_masks[label] = String[string(region_id) for region_id in region_list if string(profile_data_region(region_map, profile_type, region_id)) in bearing_regions]
            end
        end

    end

    # Initialize data structure for vre_aggregated_data_zdt (dense, zero-filled per [year][bus])
    for year in 1:num_stages
        am.ref[:nw][0][:vre_aggregated_data_zdt][year] = Dict{Int, Array{Float64,3}}()
        for zone in ids_zone
            am.ref[:nw][0][:vre_aggregated_data_zdt][year][zone] = zeros(Float64, nd, nh, 6)
        end
    end


    nw0 = am.ref[:nw][0]

    for year in ids_y

        # Check if BA versions of labels exist once per year. Use the first PRESENT repday key rather
        # than a hardcoded "1": a distributed OP worker materializes only its own day-group's repday,
        # so "1" may not exist. All repdays share the same data columns, so any present one works.
        template_data = nw0[:planning_stages][year]["repdays"][first(keys(nw0[:planning_stages][year]["repdays"]))]["data"]["1"]["1"]
        labels = Dict(
            "csp" => haskey(template_data, "csp_BA") ? "csp_BA" : "csp",
            "wind_ons" => haskey(template_data, "wind_ons_BA") ? "wind_ons_BA" : "wind_ons",
            "wind_ofs" => haskey(template_data, "wind_ofs_BA") ? "wind_ofs_BA" : "wind_ofs"
        )

        build_region_masks!(am, unique(String["hydro", labels["csp"], labels["wind_ons"], labels["wind_ofs"], "rtpv", "pv"]), year)

        # nw_data (the hourly slice) depends only on (d, h, data_label), not on zone; precompute per
        # (data_label, zone) the ORDERED data-region ids to sum so the redundant 7-level fetch drops out.
        region_map = get(nw0, :profile_data_region_map, nothing)

        # tech -> (tech_shape key, resolved data_label)
        tech_specs = [
            ("csp_shape",      labels["csp"]),
            ("pv_shape",       "pv"),
            ("rtpv_shape",     "rtpv"),
            ("wind_ons_shape", labels["wind_ons"]),
            ("wind_ofs_shape", labels["wind_ofs"]),
            ("hydro_shape",    "hydro"),
        ]

        for (shape_key, data_label) in tech_specs
            profile_type = replace(data_label, "_BA" => "")
            tech_idx = VRE_ZDT_TECH_IDX[shape_key]

            # Ordered data-region id vectors per zone, in the SAME order as region_masks (summation order preserved).
            drs_by_zone = Vector{Vector{String}}(undef, length(ids_zone))
            for (zi, zone) in enumerate(ids_zone)
                region_list = region_masks[year][zone][data_label]
                drs_by_zone[zi] = String[string(profile_data_region(region_map, profile_type, region_id)) for region_id in region_list]
            end

            # Group buses by identical ordered data-region list. Every bus in a group draws its value from the
            # same nw_data keys summed in the same order, so acc/n is bit-identical across the group -> compute
            # once per distinct group and fan out. At nodal scale with County-level data resolution this collapses
            # ~25k per-bus dict-sums down to the (much smaller) number of distinct data-region signatures.
            # Empty drs are skipped (dense store stays 0.0), matching the previous n==0 => 0.0 behavior.
            group_drs = Vector{Vector{String}}()
            group_zones = Vector{Vector{Int}}()
            sig_to_group = Dict{Vector{String}, Int}()
            for (zi, zone) in enumerate(ids_zone)
                drs = drs_by_zone[zi]
                isempty(drs) && continue   # dense store already 0.0; matches n==0 => 0.0
                gi = get(sig_to_group, drs, 0)
                if gi == 0
                    push!(group_drs, drs)
                    push!(group_zones, Int[zone])
                    sig_to_group[drs] = length(group_drs)
                else
                    push!(group_zones[gi], zone)
                end
            end

            year_store = nw0[:vre_aggregated_data_zdt][year]
            idx = nw0[:vre_zdt_index]
            for d in ids_d
                d_i = idx.d[d]
                for t in ids_t
                    h_i = idx.h[t]
                    nw_data = nw0[:planning_stages][year]["repdays"][string(d)]["data"][string(t)]["1"][data_label]
                    for gi in eachindex(group_drs)
                        drs = group_drs[gi]
                        n = length(drs)
                        acc = 0.0
                        for dr in drs
                            acc += Float64(get(nw_data, dr, 0.0))
                        end
                        val = acc / n
                        for zone in group_zones[gi]
                            year_store[zone][d_i, h_i, tech_idx] = val
                        end
                    end
                end
            end
        end
    end

end


"""
Re-point the reporting instance at a round's representative days and rebuild the caches derived from them.

The reporting instance is built once, before any round, so its repday-derived caches (VRE shapes,
per-repday fuel prices) come from its own scenario reduction. Each round re-runs scenario reduction and
can select different days, so without this the reports pair a round's dispatch with another day's profile.
"""
function refresh_report_repday_caches!(am::Abstract_ALEAF_Model, ALEAF_setting, round_network_data, round_ids_y_decision)

    round_repdays = nothing
    for decision_year_id in round_ids_y_decision
        year_data = get(round_network_data, string(decision_year_id), nothing)
        year_data === nothing && continue
        repdays = get(year_data, "planning_stages", nothing)
        repdays === nothing && continue

        # a resumed run restores repday metadata with the hourly "data" stripped, so nothing can be rebuilt from it
        if any(!haskey(rep, "data") for rep in values(repdays))
            @aleaf_warn "[ALEAF LC_GTEP]: round repday hourly data unavailable; reported curtailment and per-repday fuel cost may not match the solved representative days"
            return
        end

        am.ref[:nw][0][:planning_stages][decision_year_id]["repdays"] = repdays
        round_repdays === nothing && (round_repdays = repdays)
    end

    round_repdays === nothing && return

    am.ref[:nw][0][:repdays] = Dict{Int64, Any}(parse(Int, day_id) => rep for (day_id, rep) in round_repdays)

    aggregate_vre_timeseries_zdt_GTEP!(am)
    update_regional_fuel_prices_with_pu!(am, ALEAF_setting)

end


function aggregate_vre_timeseries_per_policy_zone_LC_GTEP!(am::Abstract_ALEAF_Model)

    num_stages = am.setting["Planning Design"]["num_stages_value"]
    
    if haskey(am.ref[:nw][0], :operation_year)
        ids_y = [am.ref[:nw][0][:operation_year]]
    else
        ids_y = [i for i in 1:num_stages]
    end

    function get_generation_shape(am::Abstract_ALEAF_Model, data_label::String, policy_zone_id::String, year::Int)

        region_list = am.ref[:nw][0][:zone]["policy"][policy_zone_id]["aggregation_info"]["aggregated_regions_bus_i"]
        num_regions = length(region_list)
        df = am.ref[:nw][0][:time_series_data][data_label]

        # Lever A: read each finest region's shape at its data-region column (identity at finest resolution).
        region_map = get(am.ref[:nw][0], :profile_data_region_map, nothing)
        profile_type = replace(data_label, "_BA" => "")
        total_generation_shape = sum(sum(df[!, Symbol(profile_data_region(region_map, profile_type, col))]) for col in region_list)

        return total_generation_shape / num_regions   # return average generation per zone
    end

    for policy_zone_id in keys(am.ref[:nw][0][:zone]["policy"])
        
        policy_data = am.ref[:nw][0][:zone]["policy"][policy_zone_id]
        policy_data["vre_aggregated_data"] = Dict{Any, Any}()

        for year in ids_y
            year_str = string(year)
            policy_data["vre_aggregated_data"][year_str] = Dict{Any, Any}()

            # Check if BA versions of labels exist once per year
            labels = Dict(
                "csp" => haskey(am.ref[:nw][0][:time_series_data], "csp_BA") ? "csp_BA" : "csp",
                "wind_ons" => haskey(am.ref[:nw][0][:time_series_data], "wind_ons_BA") ? "wind_ons_BA" : "wind_ons",
                "wind_ofs" => haskey(am.ref[:nw][0][:time_series_data], "wind_ofs_BA") ? "wind_ofs_BA" : "wind_ofs"
            )
            
            vre_list = ["hydro", "csp", "wind_ons", "wind_ofs", "rtpv", "pv"]
            vre_data = policy_data["vre_aggregated_data"][year_str]

            for vre_tech in vre_list
                label = get(labels, vre_tech, vre_tech)
                vre_data[vre_tech * "_shape"] = get_generation_shape(am, label, policy_zone_id, year)
            end
            
        end
    end

end


function aggregate_vre_timeseries_per_bus_LC_GTEP!(am::Abstract_ALEAF_Model)

    num_stages = am.setting["Planning Design"]["num_stages_value"]
    
    if haskey(am.ref[:nw][0], :operation_year)
        ids_y = [am.ref[:nw][0][:operation_year]]
    else
        ids_y = [i for i in 1:num_stages]
    end

    function get_generation_shape(am::Abstract_ALEAF_Model, data_label::String, bus_zone_id::Int, year::Int)

        region_list = am.ref[:nw][0][:bus][bus_zone_id]["aggregation_info"]["aggregated_regions_bus_i"]
        num_regions = length(region_list)
        df = am.ref[:nw][0][:time_series_data][data_label]

        # Lever A: read each finest region's shape at its data-region column (identity at finest resolution).
        region_map = get(am.ref[:nw][0], :profile_data_region_map, nothing)
        profile_type = replace(data_label, "_BA" => "")
        total_generation_shape = sum(sum(df[!, Symbol(profile_data_region(region_map, profile_type, col))]) for col in region_list)

        return total_generation_shape / num_regions   # return average generation per zone
    end

    for bus_id in keys(am.ref[:nw][0][:bus])

        bus_data = am.ref[:nw][0][:bus][bus_id]
        bus_data["vre_aggregated_data"] = Dict{Any, Any}()

        for year in ids_y
            year_str = string(year)
            bus_data["vre_aggregated_data"][year_str] = Dict{Any, Any}()

            # Check if BA versions of labels exist once per year
            labels = Dict(
                "csp" => haskey(am.ref[:nw][0][:time_series_data], "csp_BA") ? "csp_BA" : "csp",
                "wind_ons" => haskey(am.ref[:nw][0][:time_series_data], "wind_ons_BA") ? "wind_ons_BA" : "wind_ons",
                "wind_ofs" => haskey(am.ref[:nw][0][:time_series_data], "wind_ofs_BA") ? "wind_ofs_BA" : "wind_ofs"
            )
            
            vre_list = ["hydro", "csp", "wind_ons", "wind_ofs", "rtpv", "pv"]
            vre_data = bus_data["vre_aggregated_data"][year_str]

            for vre_tech in vre_list
                label = get(labels, vre_tech, vre_tech)
                vre_data[vre_tech * "_shape"] = get_generation_shape(am, label, bus_id, year)
            end
            
        end
    end

end


function add_water_management_info_LC_GTEP!(am::Abstract_ALEAF_Model)

    # Sets and data
    ids_i = [(i) for (i) in get_index(am, :gen_index, 0)]
    water_value_table = am.ref[:nw][0][:time_series_data]["hydro_water_value_df"]

    # Add Reference
    am.ref[:nw][0][:water_management] = Dict{String, Any}()

    # Define generators that needs water level variables
    ids_i_with_water_level_variables = []

    for i in ids_i
        bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
        tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)
        if parameter(am, bus_idx, :gen_bus, tech_idx, "Hydro_Flag") in ("ROR", "IMPOUNDMENT")
            fuel_limit = parameter(am, bus_idx, :gen_bus, tech_idx, "FUEL_LIMIT")
            fuel_type = parameter(am, bus_idx, :gen_bus, tech_idx, "FC")

            if (fuel_limit == "Budget (day groups)") && (fuel_type == "Water Value")
                push!(ids_i_with_water_level_variables, i)
            end           
        end
    end
    am.ref[:nw][0][:water_management]["gen_index"] = ids_i_with_water_level_variables

    # Define number of segments
    am.ref[:nw][0][:water_management]["num_segments"] = length(water_value_table[(water_value_table.Day_ID .== 1), "Reservoir_Level"])
    am.ref[:nw][0][:water_management]["segment_index"] = [i for i in 1:length(water_value_table[(water_value_table.Day_ID .== 1), "Reservoir_Level"])]
end


function add_ref_LC_GTEP_model!(am::Abstract_ALEAF_Model, ALEAF_setting::Dict{String, Any}; recorded_investment_decisions::Dict{String, Any}=Dict{String, Any}(), light::Bool=false, transfer_heavy_data::Bool=true)
    
    # model setting
    add_model_setting_GTEP!(am, ALEAF_setting)
    
    # run period
    add_run_period_GTEP!(am)

    # JSON model-reference export needs network_data whole, so it forces the full deepcopy even when the
    # caller opts into transfer-ownership.
    keep_full_data = get(ALEAF_setting["Simulation Setting"], "export_model_reference_json_operation_flag", false) == true ||
                     get(ALEAF_setting["Simulation Setting"], "export_model_reference_json_expansion_flag", false) == true

    # add common reference
    add_common_ref_GTEP!(am; transfer_heavy_data = transfer_heavy_data && !keep_full_data)

    # generation aggregation
    if ALEAF_setting["Simulation Setting"]["generation_representation_option"] == "Aggregated_by_Tech"
        aggregate_generation_GTEP!(am)
    elseif ALEAF_setting["Simulation Setting"]["generation_representation_option"] == "Individual_Plants"
        
        if ALEAF_setting["Simulation Setting"]["generation_parameter_source_option"] == "Individual Plant"
            update_individual_plant_info_GTEP_using_plant_info!(am, ALEAF_setting; recorded_investment_decisions=recorded_investment_decisions)
        elseif ALEAF_setting["Simulation Setting"]["generation_parameter_source_option"] == "Gen Technology"
            update_individual_plant_info_GTEP_using_tech_info!(am; recorded_investment_decisions=recorded_investment_decisions)
        end

    end

    # add gen index
    add_gen_index_ref_GTEP!(am)

    # add local gen index
    add_local_gen_idx_GTEP!(am)

    # define reserve-zone generator groups by technology
    add_reserve_zone_gen_technology_groups_GTEP!(am)

    # adjust gen capacity for UC
    update_CAP_i_noEXP!(am)

    # calculate vre class average generation (per year, bus, tech)
    aggregate_vre_timeseries_per_bus_LC_GTEP!(am)

    # calculate vre class average generation (per year, zone, tech,)
    aggregate_vre_timeseries_per_policy_zone_LC_GTEP!(am)

    # calculate vre class average generation (per year, zone, day, time, tech)
    # light=true (OP-distributed reporting-only reference): skip this. It is the only add_ref step that
    # needs the per-repday hourly ["data"] and it builds the heavy full-horizon vre_aggregated_data_zdt,
    # which the master-side OP summary/repday reporters never read (per-day dispatch/power-flow reporting
    # runs on the workers). The cheap per-bus/per-policy annual aggregations above are kept.
    if !light
        aggregate_vre_timeseries_zdt_GTEP!(am)
    end

    # water level management data
    add_water_management_info_LC_GTEP!(am)

    # add per unit to the market parameters
    apply_per_unit_to_market_parameters!(am, ALEAF_setting)

    # update water budget values. This reads vre_aggregated_data_zdt (built by aggregate_vre_timeseries_zdt_GTEP!,
    # skipped above when light=true) and produces :hydro_budget, a solve-only input. In the light reporting-only
    # reference (never solved; no master OP reporter reads :hydro_budget) skip it and leave :hydro_budget empty.
    if !light
        update_water_budget_values_with_pu!(am)
    else
        am.ref[:nw][0][:hydro_budget] = Dict{Int64, Any}()
    end

    # update regional fuel prices and marginal costs
    update_regional_fuel_prices_with_pu!(am, ALEAF_setting)

    # Free the redundant original network_data. `add_common_ref_GTEP!` already deepcopied it into the
    # independent `ref[:nw]` (the single keeper used by every solve + report path), so `am.data` now
    # only duplicates it. Rebind am.data to a small scalar-only dict rather than deleting in place:
    # am.data is the SAME object as the caller's network_data, which the serial OP driver reuses across
    # its per-day-group loop. Rebinding drops this instance's reference (memory freed once the caller
    # releases network_data) without corrupting the shared object. Keep full data when JSON export asks.
    # (keep_full_data computed above, before add_common_ref_GTEP!, and reused here.)
    if !keep_full_data
        am.data = Dict{String, Any}(k => am.data[k] for k in ("model_type", "output_path", "case_id", "case_name") if haskey(am.data, k))
    end

end


function apply_per_unit_to_market_parameters!(am::Abstract_ALEAF_Model, ALEAF_setting::Dict{String, Any})

    pu_power_base = am.setting["Simulation Setting"]["per_unit_base_value"]
    pu_econ_base = am.setting["Simulation Setting"]["per_unit_econ_base_value"] / pu_power_base    # per unit economic base

    # Apply per unit to market parameters
    for (key, value) in am.setting["Simulation Configuration"]
        if key in ["VOLL", "RegRSP", "SRSP", "NSRSP", "FLEXRSP", "Clean_Energy_Generation_Penalty", "RPS_Penalty"]
            am.setting["Simulation Configuration"][key] /= pu_econ_base
        end

        if key in ["CTAX"]
            am.setting["Simulation Configuration"][key] /= am.setting["Simulation Setting"]["per_unit_econ_base_value"]
        end
    end

end


function update_regional_fuel_prices_with_pu!(am::Abstract_ALEAF_Model, ALEAF_setting)

                                                                                                                                                    
    pu_power_base = am.setting["Simulation Setting"]["per_unit_base_value"]                                                                       
    pu_econ_base = am.setting["Simulation Setting"]["per_unit_econ_base_value"] / pu_power_base    # per unit economic base                       
    region_data_identifier = get(ALEAF_setting["Network Setting"], "regional_fuel_zone_resolution_type", "system")
    fuel_price_df = am.ref[:nw][0][:time_series_data]["fuel_price"]                                                                               
                                                                                                                                                
    has_monthly_prices = hasproperty(fuel_price_df, :Month)   # annual-price files have no Month column
    fuel_price_map = Dict{Tuple{String, String, Int, Int}, Any}()                                                                                 
    for row in eachrow(fuel_price_df)                                                                                                             
        key = (row.Region, row.Fuel, row.Year, has_monthly_prices ? row.Month : 0)   # month 0 = annual price, used for every month                                                                                         
        if !haskey(fuel_price_map, key)                                                                                                           
            fuel_price_map[key] = row.Price                                                                                                       
        end                                                                                                                                       
    end                                                                                                                                           
                                                                                                                                                
    function find_monthly_fuel_cost(region::String, fuel::String, year::Int, month::Int)                                                          
        key = (region, fuel, year, has_monthly_prices ? month : 0)
        if haskey(fuel_price_map, key)                                                                                                            
            return fuel_price_map[key]                                                                                                            
        end                                                                                                                                       
                                                                                                                                                
        if lowercase(string(get(am.setting["Simulation Configuration"], "logging_level_value", "simple"))) == "detailed"
            @aleaf_warn "Failed to load fuel price for type $fuel in year $year in region $region. System-wide price will be used"
        end                                        
        fallback_key = ("System-wide", fuel, year, has_monthly_prices ? month : 0)
        if haskey(fuel_price_map, fallback_key)                                                                                                   
            return fuel_price_map[fallback_key]                                                                                                   
        end                                                                                                                                       
                                                                                                                                                
        terminate_with_error(; msg = "Failed to load fuel price for type $fuel in year $year in region $region. Terminate the program")                         
    end                                                                                                                                           
                                                                                                                                                
    planning_year_map = Dict{Int, Int}()                                                                                                          
    for year_id in keys(am.ref[:nw][0][:planning_stages])                                                                                         
        planning_year_map[year_id] = am.ref[:nw][0][:planning_stages][year_id]["year"]                                                            
    end                                                                                                                                           
                                                                                                                                                
    repday_month_map = Dict{Int, Int}()                                                                                                           
    for day_id in keys(am.ref[:nw][0][:repdays])                                                                                                  
        repday_month_map[day_id] = am.ref[:nw][0][:repdays][day_id]["Month"]                                                                      
    end

    for bus_id in keys(am.ref[:nw])
        if bus_id != 0

            for local_gen_idx in keys(am.ref[:nw][bus_id][:gen_bus])

                gen_ref = am.ref[:nw][bus_id][:gen_bus][local_gen_idx] 

                # check region based on data_identifier 
                fuel_region = "System-wide"
                if gen_ref["FC"] == "Fuel-Regional"                    
                    fuel_region = am.ref[:nw][0][:bus][bus_id]["region_mapping_info"][region_data_identifier]
                end

                fuel_type = gen_ref["FUEL"]
                
                for year_id in keys(am.ref[:nw][0][:planning_stages])
                    
                    year_string = string(year_id)
                    gen_ref["Annual_FC"][year_string] = Dict{Int64, Any}() # initilize
                    gen_ref["Annual_MC"][year_string] = Dict{Int64, Any}() # initilize
                    
                    vom = gen_ref["Annual_VOM"][year_string] * pu_econ_base # $/MWh

                    for day_id in keys(am.ref[:nw][0][:repdays])

                        month = am.ref[:nw][0][:repdays][day_id]["Month"]

                        fuel_cost = gen_ref["FC"]
                        if fuel_cost in ["Fuel", "Fuel-Regional"]
                            fuel_cost = find_monthly_fuel_cost(fuel_region, fuel_type, am.ref[:nw][0][:planning_stages][year_id]["year"], month) 
                        end

                        # add fuel cost
                        gen_ref["Annual_FC"][year_string][day_id] = fuel_cost 
                        gen_ref["Annual_MC"][year_string][day_id] = gen_ref["HR"] * fuel_cost + vom
                                                                
                        # apply per unit
                        gen_ref["Annual_FC"][year_string][day_id] /= pu_econ_base
                        gen_ref["Annual_MC"][year_string][day_id] /= pu_econ_base

                    end
                    
                end
            end
        end
    end    
                                                                                                                                             
end


function update_water_budget_values_with_pu!(am::Abstract_ALEAF_Model)

    pu_power_base = am.setting["Simulation Setting"]["per_unit_base_value"]
    pu_econ_base = am.setting["Simulation Setting"]["per_unit_econ_base_value"] / pu_power_base    # per unit economic base

    # set reference year
    ref_year_idx = 1
    if haskey(am.ref[:nw][0], :operation_year)
        ref_year_idx = am.ref[:nw][0][:operation_year]
    end
    
    am.ref[:nw][0][:hydro_budget] = Dict{Int64, Any}()
    # calculate hydro budget for each bus
    for bus_idx in keys(am.ref[:nw][0][:bus])

        am.ref[:nw][0][:hydro_budget][bus_idx] = Dict{Int64, Any}() # hydro budget for each day group
        
        for day_group_id in keys(am.ref[:nw][0][:repday_groups])

            start_day = am.ref[:nw][0][:repday_groups][day_group_id]["Start_Day_Id"] 
            end_day = am.ref[:nw][0][:repday_groups][day_group_id]["End_Day_Id"] 

            # total budget
            # coarse hydro data resolution can fold many finest regions onto one budget column; summing per
            # finest region N-counts it, so count each distinct hydro column once (unique() no-ops on nodal/NA identity maps)
            total_budget_of_bus_in_day_group = 0.0
            region_map = get(am.ref[:nw][0], :profile_data_region_map, nothing)
            hydro_cols = unique(profile_data_region(region_map, "hydro", ba_ids)
                for ba_ids in am.ref[:nw][0][:bus][bus_idx]["aggregation_info"]["aggregated_regions_bus_i"])
            for col in hydro_cols
                total_budget_of_bus_in_day_group += sum(am.ref[:nw][0][:time_series_data]["hydro_budget_df"][start_day:end_day, col]) / pu_power_base
            end

            regional_total_impoundment_capacity = 0.0
            for local_gen_idx in am.ref[:nw][0][:bus][bus_idx]["aggregation_info"]["new_local_gen_idx"]
                if parameter(am, bus_idx, :gen_bus, parameter(am, 0, :gen_index, "genco_tech_id", local_gen_idx), "Hydro_Flag") == "IMPOUNDMENT"
                    CAP = parameter(am, bus_idx, :gen_bus, parameter(am, 0, :gen_index, "genco_tech_id", local_gen_idx), "CAP")
                    EXUNITS = parameter(am, bus_idx, :gen_bus, parameter(am, 0, :gen_index, "genco_tech_id", local_gen_idx), "EXUNITS")
                    regional_total_impoundment_capacity += CAP * EXUNITS
                end
            end

            if regional_total_impoundment_capacity > 0
                # check budget feasibility
                total_fixed_output = 0.0

                for d in am.ref[:nw][0][:repday_groups][day_group_id]["Day_Idx_List"]
                    for (h, t) in [(h,t) for h in am.setting["run_H"] for t in am.setting["run_T"]]
                        total_fixed_output += get_vre_zdt_shape(am.ref[:nw][0], ref_year_idx, bus_idx, d, h, "hydro_shape") * regional_total_impoundment_capacity
                    end
                end

                if total_fixed_output > total_budget_of_bus_in_day_group
                    # raise the day-group budget to the fixed hydro output when it would be infeasible
                    total_budget_of_bus_in_day_group = total_fixed_output
                end

                # convert MWh to MWh/MW
                am.ref[:nw][0][:hydro_budget][bus_idx][day_group_id] = total_budget_of_bus_in_day_group / regional_total_impoundment_capacity
            else
                am.ref[:nw][0][:hydro_budget][bus_idx][day_group_id] = 0.0
            end
        end
    end

end


function add_gen_index_ref_GTEP!(am::Abstract_ALEAF_Model)

    # add gen_index
    am.ref[:nw][0][:gen_index] = Dict{Int64, Any}()
    gen_idx = 1

    for sub_area_idx in keys(am.ref[:nw][0][:bus])
        for genco_tech_id in sort!(collect(keys(am.ref[:nw][sub_area_idx][:gen_bus])))
            UNIT_REPORT_LABEL_1 = am.ref[:nw][sub_area_idx][:gen_bus][genco_tech_id]["UNIT_REPORT_LABEL_1"]
            UNIT_REPORT_LABEL_2 = am.ref[:nw][sub_area_idx][:gen_bus][genco_tech_id]["UNIT_REPORT_LABEL_2"]
            UNIT_GROUP = am.ref[:nw][sub_area_idx][:gen_bus][genco_tech_id]["UNITGROUP"]
            UNIT_CATEGORY = am.ref[:nw][sub_area_idx][:gen_bus][genco_tech_id]["UNIT_CATEGORY"]
            Profile_Type = am.ref[:nw][sub_area_idx][:gen_bus][genco_tech_id]["Profile_Type"]
            hybrid_type = am.ref[:nw][sub_area_idx][:gen_bus][genco_tech_id]["hybrid_type"]
            am.ref[:nw][0][:gen_index][gen_idx] = Dict{String, Any}("bus_idx" => sub_area_idx,
                                                                    "genco_tech_id" => genco_tech_id,
                                                                    "genco_tech_UNIT_REPORT_LABEL_1" => UNIT_REPORT_LABEL_1,
                                                                    "genco_tech_UNIT_REPORT_LABEL_2" => UNIT_REPORT_LABEL_2,
                                                                    "UNIT_GROUP" => UNIT_GROUP,
                                                                    "UNIT_CATEGORY" => UNIT_CATEGORY,
                                                                    "Profile_Type" => Profile_Type,
                                                                    "hybrid_type" => hybrid_type
                                                                    )
            am.ref[:nw][sub_area_idx][:gen_bus][genco_tech_id]["gen_idx"] = gen_idx
            gen_idx += 1
        end
    end

    # update hybrid gen index 
    for gen_idx in keys(am.ref[:nw][0][:gen_index])
        hybrid_type = am.ref[:nw][0][:gen_index][gen_idx]["hybrid_type"]
        if hybrid_type == "GEN" # this is the main hybrid unit
            bus_idx = am.ref[:nw][0][:gen_index][gen_idx]["bus_idx"]
            genco_tech_id = am.ref[:nw][0][:gen_index][gen_idx]["genco_tech_id"]
            Hybrid_ID = am.ref[:nw][bus_idx][:gen_bus][genco_tech_id]["PLANT_NAME"]

            # find other hybrid members
            for other_genco_tech_id in keys(am.ref[:nw][bus_idx][:gen_bus])
                other_gen_idx = am.ref[:nw][bus_idx][:gen_bus][other_genco_tech_id]["gen_idx"]
                hybrid_type = am.ref[:nw][0][:gen_index][other_gen_idx]["hybrid_type"]
                if (hybrid_type == "ES") && (other_genco_tech_id != genco_tech_id)
                    other_Hybrid_ID = am.ref[:nw][bus_idx][:gen_bus][other_genco_tech_id]["PLANT_NAME"]                    
                    if (other_Hybrid_ID == Hybrid_ID) 
                        am.ref[:nw][0][:gen_index][other_gen_idx]["hybrid_main_gen_idx"] = gen_idx
                        am.ref[:nw][0][:gen_index][gen_idx]["hybrid_ES_gen_idx"] = other_gen_idx
                    end
                end
            end
        end
    end

    am.ref[:nw][0][:gen_index] = sort(am.ref[:nw][0][:gen_index])

    # free per-nw dicts no longer needed; repdays/bus/branch/annualData/output_path kept for downstream use
    for idx in keys(am.ref[:nw])
        if idx != 0
            delete!(am.ref[:nw][idx], :sub_area_mapping)
            delete!(am.ref[:nw][idx], :gen)
            delete!(am.ref[:nw][idx], :plant)
            delete!(am.ref[:nw][idx], :gen_technology)
        end
    end

end


function add_common_ref_GTEP!(am::Abstract_ALEAF_Model; transfer_heavy_data::Bool=true)

    am.ref[:nw][0] = Dict{Symbol,Any}()

    # Transfer-ownership: at nodal scale a blanket deepcopy(am.data) 2x-peaks the heavy hourly payload
    # (per-repday ["data"] tensors + raw time_series_data) and OOM-kills the serial master. Detach the
    # heavy data (by reference, no copy) before the light deepcopy, then reattach it into nw[0] so nw[0]
    # is the sole keeper and network_data is left light. transfer_heavy_data=false preserves the full
    # deepcopy for the one caller that add_refs twice on a shared network_data (distributed OP worker).
    detached = transfer_heavy_data ? detach_heavy_network_data!(am.data) : nothing

    nws_data = deepcopy(am.data)
    for (key, item) in nws_data
        if isa(item, Dict{String,Any})
            try
                item_lookup = Dict{Int,Any}([(parse(Int, k), v) for (k,v) in item])
                am.ref[:nw][0][Symbol(key)] = item_lookup
            catch
                am.ref[:nw][0][Symbol(key)] = item
            end
        elseif isa(item, Dict{Int64,Any})
            item_lookup = Dict{Int,Any}([(k, v) for (k,v) in item])
            am.ref[:nw][0][Symbol(key)] = item_lookup
        else
            am.ref[:nw][0][Symbol(key)] = item
        end
    end

    detached === nothing || reattach_heavy_network_data!(am.ref[:nw][0], detached)

end


# Pop (by reference, no copy) the heavy hourly payload out of a network_data dict so a subsequent
# deepcopy stays light. Returns the detached objects for reattachment into nw[0]; leaves data light.
function detach_heavy_network_data!(data::Dict)
    tsd = pop!(data, "time_series_data", nothing)

    # per-repday ["data"] tensors: top-level repdays and each planning_stages[s]["repdays"] are separate
    # dicts (generate_network deepcopies stage-1 repdays), so detach both, keyed by identity path.
    rep_top = Dict{Any,Any}()
    if haskey(data, "repdays")
        for (d, rd) in data["repdays"]
            (isa(rd, Dict) && haskey(rd, "data")) || continue
            rep_top[d] = pop!(rd, "data")
        end
    end
    rep_stage = Dict{Any,Any}()
    if haskey(data, "planning_stages")
        for (s, stage) in data["planning_stages"]
            (isa(stage, Dict) && haskey(stage, "repdays")) || continue
            inner = Dict{Any,Any}()
            for (d, rd) in stage["repdays"]
                (isa(rd, Dict) && haskey(rd, "data")) || continue
                inner[d] = pop!(rd, "data")
            end
            isempty(inner) || (rep_stage[s] = inner)
        end
    end

    return (time_series_data=tsd, rep_top=rep_top, rep_stage=rep_stage)
end


# Attach the detached heavy payload (same objects) into nw[0]'s light containers, making nw[0] whole.
# nw[0] top-level keys were rekeyed String->Int in add_common_ref; nested "data"/"repdays" keys are not.
function reattach_heavy_network_data!(nw0::Dict{Symbol,Any}, detached)
    detached.time_series_data === nothing || (nw0[:time_series_data] = detached.time_series_data)

    if haskey(nw0, :repdays)
        for (d, dat) in detached.rep_top
            key = _rekeyed_lookup(nw0[:repdays], d)
            key === nothing && continue
            isa(nw0[:repdays][key], Dict) && (nw0[:repdays][key]["data"] = dat)
        end
    end
    if haskey(nw0, :planning_stages)
        for (s, inner) in detached.rep_stage
            skey = _rekeyed_lookup(nw0[:planning_stages], s)
            skey === nothing && continue
            stage = nw0[:planning_stages][skey]
            (isa(stage, Dict) && haskey(stage, "repdays")) || continue
            for (d, dat) in inner
                haskey(stage["repdays"], d) || continue
                isa(stage["repdays"][d], Dict) && (stage["repdays"][d]["data"] = dat)
            end
        end
    end
end


# add_common_ref rekeys top-level String keys to Int when parseable; map an original key to its nw[0] key.
_rekeyed_lookup(d::AbstractDict, k) = haskey(d, k) ? k :
    (isa(k, AbstractString) ? (tryparse(Int, k) !== nothing && haskey(d, parse(Int, k)) ? parse(Int, k) : nothing) : nothing)


function aggregate_generation_GTEP!(am::Abstract_ALEAF_Model; recorded_investment_decisions::Dict{String, Any}=Dict{String, Any}())

    function find_gen_tech_GTEP(am::Abstract_ALEAF_Model, gen_dict::Dict{String, Any}, sub_area_id::Int64)

        gen_tech_id = 0
        for idx in keys(am.ref[:nw][sub_area_id][:gen_technology])
            if am.ref[:nw][sub_area_id][:gen_technology][idx]["UNITGROUP"] == gen_dict["UNITGROUP"]
                gen_tech_id = idx
            end
        end
    
        if gen_tech_id == 0
            gen_id = gen_dict["GEN UID"]
            @aleaf_info "ERROR: Failed to find gen_technology type; gen uid: $gen_id; temporal id is given"
            gen_tech_id = 1
        end
    
        return gen_tech_id        
    end

    function merge_missing!(plant_dict::Dict, gen_tech_info)
        
        for (k, v) in gen_tech_info
            get!(plant_dict, k, v)
        end
    end

    function find_gen_tech_info(am::Abstract_ALEAF_Model, gen_dict::Dict{String, Any}, sub_area_id::Int64)

        gen_tech_id = 0
        for idx in keys(am.ref[:nw][sub_area_id][:gen_technology])
            if (am.ref[:nw][sub_area_id][:gen_technology][idx]["UNITGROUP"] == gen_dict["UNITGROUP"]) 
                gen_tech_id = idx
            end
        end
    
        if gen_tech_id == 0
            gen_id = gen_dict["GEN UID"]
            @aleaf_info "ERROR: Failed to find gen_technology type; gen uid: $gen_id; temporal id is given"
            gen_tech_id = 1
        end
    
        return am.ref[:nw][sub_area_id][:gen_technology][gen_tech_id]        
    end

    add_investment_option_flag = false 
    if (am.setting["Simulation Configuration"]["Run_expansion_flag"] == true) 
        add_investment_option_flag = true
    elseif (am.setting["Simulation Configuration"]["Run_RA_flag"] == true) && (am.setting["RA Setting"]["generate_gen_index_with_investment_options_flag"] == true)
        add_investment_option_flag = true
    elseif length(recorded_investment_decisions) > 0    # we need to add investment options when we have recorded investment data
        add_investment_option_flag = true
    end

    # index plants by their bus region id once, so per-sub_area gathering is O(regions+plants)
    # instead of re-scanning every plant per sub_area (was O(N_sub_areas x N_plants) at nodal)
    plants_by_region_bus_i = Dict{String, Vector{Any}}()
    for (plant_id, plant_data) in am.ref[:nw][0][:plant]
        push!(get!(plants_by_region_bus_i, string(plant_data["bus_ID"]), Vector{Any}()), plant_id)
    end

    for sub_area_id in keys(am.ref[:nw][0][:bus])    # the "bus" dict already completed the network reduction

        # add sub_area dict first
        am.ref[:nw][sub_area_id][:gen_bus] = Dict{Int64, Any}()

        # Prepare aggregation
        am.ref[:nw][sub_area_id][:gen_bus] = deepcopy(am.ref[:nw][sub_area_id][:gen_technology])
        for gen_tech_id in keys(am.ref[:nw][sub_area_id][:gen_bus])
            
            # aggregated plant list
            am.ref[:nw][sub_area_id][:gen_bus][gen_tech_id]["Agg_Gen_UID"] = [] 

            # add new Gen uid (bus_id + tech)
            gen_uid =  am.ref[:nw][sub_area_id][:gen_bus][gen_tech_id]["UNITGROUP"]                    
            new_gen_uid = string(sub_area_id, "_", gen_uid)                    
            am.ref[:nw][sub_area_id][:gen_bus][gen_tech_id]["NEW_Gen_UID"] = new_gen_uid
            
            PLANT_NAME = string("agg_", gen_uid, "_", sub_area_id)        
            am.ref[:nw][sub_area_id][:gen_bus][gen_tech_id]["PLANT_NAME"] = PLANT_NAME
            
            am.ref[:nw][sub_area_id][:gen_bus][gen_tech_id]["EXCAPS"] = 0.0
            am.ref[:nw][sub_area_id][:gen_bus][gen_tech_id]["EXUNITS"] = 0.0
            am.ref[:nw][sub_area_id][:gen_bus][gen_tech_id]["ES_MWh"] = 0.0

            # add timeseries data tag
            am.ref[:nw][sub_area_id][:gen_bus][gen_tech_id]["Timeseries_Tag"] = "LOCAL"

            # add planned retirement tag
            am.ref[:nw][sub_area_id][:gen_bus][gen_tech_id]["Planned_Retirement"] = Dict()
            for stage_id in keys(am.ref[:nw][0][:planning_stages])
                year = am.ref[:nw][0][:planning_stages][stage_id]["year"]
                tag = string("Ret_", year)

                am.ref[:nw][sub_area_id][:gen_bus][gen_tech_id]["Planned_Retirement"][tag] = 0.0
            end

            am.ref[:nw][sub_area_id][:gen_bus][gen_tech_id]["bypass_parameter_check"] = false # false => perform check_and_update_plant_data!
            am.ref[:nw][sub_area_id][:gen_bus][gen_tech_id]["hybrid_type"] = "NA" # will be updated in the update_hybrid_plant_technology_data! function
            am.ref[:nw][sub_area_id][:gen_bus][gen_tech_id]["hybrid_ID"] = 0  # will be updated in the update_hybrid_plant_technology_data! function
        end

        # update CAPAX based on the scaling factors
        if am.setting["Simulation Configuration"]["Regional_CAPAX_scaling_flag"] == true
            for gen_tech_id in keys(am.ref[:nw][sub_area_id][:gen_bus])
                
                if am.ref[:nw][sub_area_id][:gen_bus][gen_tech_id]["Locational_Scaling_Flag"] == true

                    gen_tech_idx = am.ref[:nw][sub_area_id][:gen_bus][plant_id]["Tech_ID"]
                    data_identifier = am.ref[:nw][0][:zone]["cost_scaling"]["data_identifier"]
                    data_region = am.ref[:nw][0][:bus][sub_area_id]["region_mapping_info"][data_identifier]

                    regional_scaling = am.ref[:nw][0][:zone]["cost_scaling"][data_region][gen_tech_idx]["CAPAX_scale"]

                    for year_idx in keys(am.ref[:nw][sub_area_id][:gen_bus][gen_tech_id]["Annual_INVC"])
                        am.ref[:nw][sub_area_id][:gen_bus][gen_tech_id]["Annual_INVC"][year_idx] *= regional_scaling
                    end

                    if haskey(am.ref[:nw][sub_area_id][:gen_bus][gen_tech_id], "STO_CAPEX")
                        for year_idx in keys(am.ref[:nw][sub_area_id][:gen_bus][gen_tech_id]["Annual_STO_INV"])
                            am.ref[:nw][sub_area_id][:gen_bus][gen_tech_id]["Annual_STO_INV"][year_idx] *= regional_scaling
                        end
                    end
                end
            end
        end

        # update gen dict based on existing unit_type; gather from the region->plant index
        # instead of scanning every plant per sub_area (set-identical to the prior filter)
        regs = am.ref[:nw][0][:bus][sub_area_id]["aggregation_info"]["aggregated_regions_bus_i"]
        plants_at_this_bus = Dict(pid => am.ref[:nw][0][:plant][pid] for r in regs for pid in get(plants_by_region_bus_i, string(r), Any[]))

        for gen_id in keys(plants_at_this_bus)

             # find existing hybrid plants
            hybrid_type = plants_at_this_bus[gen_id]["hybrid_type"]

            # find gen technology id
            gen_tech_id = find_gen_tech_GTEP(am, am.ref[:nw][0][:plant][gen_id], sub_area_id)

            if hybrid_type == "NA"  # this is not a hybrid plant so we can aggregate this plant

                # update gen dict existing capacity
                am.ref[:nw][sub_area_id][:gen_bus][gen_tech_id]["EXCAPS"] += am.ref[:nw][0][:plant][gen_id]["CAP"] 

                # update ES MWh capacity
                am.ref[:nw][sub_area_id][:gen_bus][gen_tech_id]["ES_MWh"] += am.ref[:nw][0][:plant][gen_id]["ES_MWh"] 

                # update gen dict existing number of units
                tech_cap = am.ref[:nw][sub_area_id][:gen_bus][gen_tech_id]["CAP"] # this value is pu already
                new_existing_units = am.ref[:nw][sub_area_id][:gen_bus][gen_tech_id]["EXCAPS"] / tech_cap
                am.ref[:nw][sub_area_id][:gen_bus][gen_tech_id]["EXUNITS"] = new_existing_units

                # ES storage hour
                if am.ref[:nw][sub_area_id][:gen_bus][gen_tech_id]["ES_MWh"] != 0
                    ES_MWh_hour = am.ref[:nw][sub_area_id][:gen_bus][gen_tech_id]["ES_MWh"] / am.ref[:nw][sub_area_id][:gen_bus][gen_tech_id]["EXCAPS"]
                else
                    ES_MWh_hour = 0.0
                end

                # update retirement year if not available
                if am.ref[:nw][0][:plant][gen_id]["RetireYear"] == "NA"
                    online_year = am.ref[:nw][0][:plant][gen_id]["Online_Year"]
                    if online_year != "NA"
                        am.ref[:nw][0][:plant][gen_id]["RetireYear"] = online_year + am.ref[:nw][sub_area_id][:gen_bus][gen_tech_id]["Life"] - 1
                    else
                        am.ref[:nw][0][:plant][gen_id]["RetireYear"] = 9999
                    end
                end

                # add planned retirement tag
                for stage_id in keys(am.ref[:nw][0][:planning_stages])
                    year_at_y = am.ref[:nw][0][:planning_stages][stage_id]["year"]
                    stage_length = am.ref[:nw][0][:planning_stages][stage_id]["stage_length"]

                    # add previous planned retirement to the first planning stage
                    retire_year = am.ref[:nw][0][:plant][gen_id]["RetireYear"]
                    
                    if stage_id == 1
                        if retire_year <= year_at_y
                            am.ref[:nw][sub_area_id][:gen_bus][gen_tech_id]["Planned_Retirement"][string("Ret_", year_at_y)] += am.ref[:nw][0][:plant][gen_id]["CAP"] 
                        elseif retire_year > year_at_y && retire_year <= (year_at_y + stage_length - 1)
                            am.ref[:nw][sub_area_id][:gen_bus][gen_tech_id]["Planned_Retirement"][string("Ret_", year_at_y)] += am.ref[:nw][0][:plant][gen_id]["CAP"] 
                        end
                    else
                        if retire_year > year_at_y && retire_year <= (year_at_y + stage_length - 1)
                            am.ref[:nw][sub_area_id][:gen_bus][gen_tech_id]["Planned_Retirement"][string("Ret_", year_at_y)] += am.ref[:nw][0][:plant][gen_id]["CAP"] 
                        end
                        
                    end
                end
            
            elseif hybrid_type == "GEN" || hybrid_type == "ES" # this is hybrid plant so we cannot aggregate this plant

                plant_id = parse(Int, "$(sub_area_id)010$(gen_id)010$(gen_tech_id)")  # generate new plant id

                plant_data = deepcopy(plants_at_this_bus[gen_id])

                gen_tech_info = find_gen_tech_info(am, plant_data, sub_area_id)
                merge_missing!(plant_data, gen_tech_info)

                # copy plant data to each bus
                original_cap = deepcopy(plant_data["CAP"])
                am.ref[:nw][sub_area_id][:gen_bus][plant_id] = plant_data

                # add new Gen uid (bus_id + tech)
                gen_uid =  gen_tech_info["UNITGROUP"]                    
                new_gen_uid = string(sub_area_id, "_", gen_uid)                    
                am.ref[:nw][sub_area_id][:gen_bus][plant_id]["NEW_Gen_UID"] = new_gen_uid
    
                # update plant information
                am.ref[:nw][sub_area_id][:gen_bus][plant_id]["CAP"] = original_cap
                am.ref[:nw][sub_area_id][:gen_bus][plant_id]["EXCAPS"] = original_cap
                am.ref[:nw][sub_area_id][:gen_bus][plant_id]["EXUNITS"] = 1.0

                # add timeseries data tag
                am.ref[:nw][sub_area_id][:gen_bus][plant_id]["Timeseries_Tag"] = "LOCAL"
                
                # update retirement year if not available
                if am.ref[:nw][sub_area_id][:gen_bus][plant_id]["RetireYear"] == "NA"
                    online_year = am.ref[:nw][sub_area_id][:gen_bus][plant_id]["Online_Year"]
                    if online_year != "NA"
                        am.ref[:nw][sub_area_id][:gen_bus][plant_id]["RetireYear"] = online_year + gen_tech_info["Life"] - 1
                    else
                        am.ref[:nw][sub_area_id][:gen_bus][plant_id]["RetireYear"] = 9999
                    end
                end

                # add planned retirement tag
                am.ref[:nw][sub_area_id][:gen_bus][plant_id]["Planned_Retirement"] = Dict()
                for stage_id in keys(am.ref[:nw][0][:planning_stages])
                    year_at_y = am.ref[:nw][0][:planning_stages][stage_id]["year"]
                    stage_length = am.ref[:nw][0][:planning_stages][stage_id]["stage_length"]

                    am.ref[:nw][sub_area_id][:gen_bus][plant_id]["Planned_Retirement"][string("Ret_", year_at_y)] = 0.0

                    # add previous planned retirement to the first planning stage
                    retire_year = am.ref[:nw][sub_area_id][:gen_bus][plant_id]["RetireYear"]
                    if stage_id == 1
                        if retire_year <= year_at_y
                            am.ref[:nw][sub_area_id][:gen_bus][plant_id]["Planned_Retirement"][string("Ret_", year_at_y)] += original_cap
                        elseif retire_year > year_at_y && retire_year <= (year_at_y + stage_length - 1)
                            am.ref[:nw][sub_area_id][:gen_bus][plant_id]["Planned_Retirement"][string("Ret_", year_at_y)] += original_cap
                        end
                    else
                        if retire_year > year_at_y && retire_year <= (year_at_y + stage_length - 1)
                            am.ref[:nw][sub_area_id][:gen_bus][plant_id]["Planned_Retirement"][string("Ret_", year_at_y)] += original_cap
                        end
                    end
                end

                # update plant data
                if plant_data["bypass_parameter_check"] == false
                    check_and_update_plant_data!(am, ALEAF_setting, plant_data, gen_tech_info, parse(Int, am.ref[:nw][0][:case_id]))
                end
            end
        end

        # delete sub-dicts that won't be used (no existing unit & no investment option)
        for gen_tech_id in keys(am.ref[:nw][sub_area_id][:gen_bus])
            if am.ref[:nw][sub_area_id][:gen_bus][gen_tech_id]["EXUNITS"] == 0.0
                if am.ref[:nw][sub_area_id][:gen_bus][gen_tech_id]["INVEST_FLAG"] == false
                    delete!(am.ref[:nw][sub_area_id][:gen_bus], gen_tech_id)
                end
            end
        end
    end
end


function update_individual_plant_info_GTEP_using_tech_info!(am::Abstract_ALEAF_Model; recorded_investment_decisions=[])

    function find_gen_tech_info(am::Abstract_ALEAF_Model, gen_dict::Dict{String, Any}, sub_area_id::Int64)

        gen_tech_id = 0
        for idx in keys(am.ref[:nw][sub_area_id][:gen_technology])
            if (am.ref[:nw][sub_area_id][:gen_technology][idx]["UNITGROUP"] == gen_dict["UNITGROUP"]) 
                gen_tech_id = idx
            end
        end
    
        if gen_tech_id == 0
            gen_id = gen_dict["GEN UID"]
            @aleaf_info "ERROR: Failed to find gen_technology type; gen uid: $gen_id; temporal id is given"
            gen_tech_id = 1
        end
    
        return am.ref[:nw][sub_area_id][:gen_technology][gen_tech_id]        
    end

    add_investment_option_flag = false 
    if (am.setting["Simulation Configuration"]["Run_expansion_flag"] == true) 
        add_investment_option_flag = true
    elseif (am.setting["Simulation Configuration"]["Run_RA_flag"] == true) && (am.setting["RA Setting"]["generate_gen_index_with_investment_options_flag"] == true)
        add_investment_option_flag = true
    elseif length(recorded_investment_decisions) > 0    # we need to add investment options when we have recorded investment data
        add_investment_option_flag = true
    end

    # index plants by their bus region id once, so per-sub_area gathering is O(regions+plants)
    # instead of re-scanning every plant per sub_area (was O(N_sub_areas x N_plants) at nodal)
    plants_by_region_bus_i = Dict{String, Vector{Any}}()
    for (plant_id, plant_data) in am.ref[:nw][0][:plant]
        push!(get!(plants_by_region_bus_i, string(plant_data["bus_ID"]), Vector{Any}()), plant_id)
    end

    for sub_area_id in keys(am.ref[:nw][0][:bus])    # the "bus" dict already completed the network reduction

        # add sub_area dict first
        am.ref[:nw][sub_area_id][:gen_bus] = Dict{Int64, Any}()

        # allocate plant data to each bus
        if length(am.ref[:nw][0][:plant]) > 0

            # copy plant data to each bus; gather from the region->plant index (set-identical to prior filter)
            regs = am.ref[:nw][0][:bus][sub_area_id]["aggregation_info"]["aggregated_regions_bus_i"]
            src = am.ref[:nw][0][:plant]
            am.ref[:nw][sub_area_id][:gen_bus] = Dict{keytype(src), valtype(src)}(pid => src[pid] for r in regs for pid in get(plants_by_region_bus_i, string(r), Any[]))

            for plant_id in keys(am.ref[:nw][sub_area_id][:gen_bus])

                original_cap = deepcopy(am.ref[:nw][sub_area_id][:gen_bus][plant_id]["CAP"])
                original_ES_MWh = deepcopy(am.ref[:nw][sub_area_id][:gen_bus][plant_id]["ES_MWh"])
                original_Charge_CAP = deepcopy(am.ref[:nw][sub_area_id][:gen_bus][plant_id]["Charge_CAP"])
                
                # find and merge tech info
                gen_tech_info = copy(find_gen_tech_info(am, am.ref[:nw][sub_area_id][:gen_bus][plant_id], sub_area_id))
                am.ref[:nw][sub_area_id][:gen_bus][plant_id] = merge(am.ref[:nw][sub_area_id][:gen_bus][plant_id], gen_tech_info)
                
                # add new Gen uid (bus_id + tech)
                gen_uid =  gen_tech_info["UNITGROUP"]                    
                new_gen_uid = string(sub_area_id, "_", gen_uid)                    
                am.ref[:nw][sub_area_id][:gen_bus][plant_id]["NEW_Gen_UID"] = new_gen_uid
    
                # update plant information
                am.ref[:nw][sub_area_id][:gen_bus][plant_id]["CAP"] = original_cap
                am.ref[:nw][sub_area_id][:gen_bus][plant_id]["EXCAPS"] = original_cap
                am.ref[:nw][sub_area_id][:gen_bus][plant_id]["EXUNITS"] = 1.0

                am.ref[:nw][sub_area_id][:gen_bus][plant_id]["ES_MWh"] = original_ES_MWh 
                am.ref[:nw][sub_area_id][:gen_bus][plant_id]["Charge_CAP"] = original_Charge_CAP
                
                # add timeseries data tag
                am.ref[:nw][sub_area_id][:gen_bus][plant_id]["Timeseries_Tag"] = "LOCAL"
                
                # update retirement year if not available
                if am.ref[:nw][sub_area_id][:gen_bus][plant_id]["RetireYear"] == "NA"
                    online_year = am.ref[:nw][sub_area_id][:gen_bus][plant_id]["Online_Year"]
                    if online_year != "NA"
                        am.ref[:nw][sub_area_id][:gen_bus][plant_id]["RetireYear"] = online_year + gen_tech_info["Life"] - 1
                    else
                        am.ref[:nw][sub_area_id][:gen_bus][plant_id]["RetireYear"] = 9999
                    end
                end

                # add planned retirement tag
                am.ref[:nw][sub_area_id][:gen_bus][plant_id]["Planned_Retirement"] = Dict()
                for stage_id in keys(am.ref[:nw][0][:planning_stages])
                    year_at_y = am.ref[:nw][0][:planning_stages][stage_id]["year"]
                    stage_length = am.ref[:nw][0][:planning_stages][stage_id]["stage_length"]

                    am.ref[:nw][sub_area_id][:gen_bus][plant_id]["Planned_Retirement"][string("Ret_", year_at_y)] = 0.0

                    # add previous planned retirement to the first planning stage
                    retire_year = am.ref[:nw][sub_area_id][:gen_bus][plant_id]["RetireYear"]
                    if stage_id == 1
                        if retire_year <= year_at_y
                            am.ref[:nw][sub_area_id][:gen_bus][plant_id]["Planned_Retirement"][string("Ret_", year_at_y)] += original_cap
                        elseif retire_year > year_at_y && retire_year <= (year_at_y + stage_length - 1)
                            am.ref[:nw][sub_area_id][:gen_bus][plant_id]["Planned_Retirement"][string("Ret_", year_at_y)] += original_cap
                        end
                    else
                        if retire_year > year_at_y && retire_year <= (year_at_y + stage_length - 1)
                            am.ref[:nw][sub_area_id][:gen_bus][plant_id]["Planned_Retirement"][string("Ret_", year_at_y)] += original_cap
                        end
                    end
                end
            end
        else
            am.ref[:nw][sub_area_id][:gen_bus] = Dict()
        end

        # Add new investment options in the expansion run 
        if add_investment_option_flag == true
            
            for idx in keys(am.ref[:nw][sub_area_id][:gen_technology])
                if am.ref[:nw][sub_area_id][:gen_technology][idx]["Tech_Type"] == "New"

                    if am.ref[:nw][sub_area_id][:gen_technology][idx]["INVEST_FLAG"] == true
                    
                        # new_idx (bus_id + techid)
                        new_genco_tech_id = parse(Int, string(sub_area_id) * "000" *string(idx))
                        
                        am.ref[:nw][sub_area_id][:gen_bus][new_genco_tech_id] = copy(am.ref[:nw][sub_area_id][:gen_technology][idx])

                        # add new Gen uid (bus_id + tech)
                        gen_uid =  am.ref[:nw][sub_area_id][:gen_bus][new_genco_tech_id]["UNITGROUP"]                    
                        new_gen_uid = string(sub_area_id, "_", gen_uid)                    
                        am.ref[:nw][sub_area_id][:gen_bus][new_genco_tech_id]["NEW_Gen_UID"] = new_gen_uid

                        # update plant information
                        am.ref[:nw][sub_area_id][:gen_bus][new_genco_tech_id]["EXCAPS"] = 0.0
                        am.ref[:nw][sub_area_id][:gen_bus][new_genco_tech_id]["EXUNITS"] = 0.0
                        am.ref[:nw][sub_area_id][:gen_bus][new_genco_tech_id]["ES_MWh"] = 0.0
                        
                        # add timeseries data tag
                        am.ref[:nw][sub_area_id][:gen_bus][new_genco_tech_id]["Timeseries_Tag"] = "LOCAL"

                        # add planned retirement tag
                        am.ref[:nw][sub_area_id][:gen_bus][new_genco_tech_id]["Planned_Retirement"] = Dict()
                        for stage_id in keys(am.ref[:nw][0][:planning_stages])
                            year = am.ref[:nw][0][:planning_stages][stage_id]["year"]
                            tag = string("Ret_", year)
                            am.ref[:nw][sub_area_id][:gen_bus][new_genco_tech_id]["Planned_Retirement"][tag] = 0.0
                        end

                        # add hybrid parameters
                        am.ref[:nw][sub_area_id][:gen_bus][new_genco_tech_id]["bypass_parameter_check"] = false # false => perform check_and_update_plant_data!
                        am.ref[:nw][sub_area_id][:gen_bus][new_genco_tech_id]["hybrid_type"] = "NA"
                        am.ref[:nw][sub_area_id][:gen_bus][new_genco_tech_id]["hybrid_ID"] = 0
                    end
                end
            end
        end

        # update CAPAX based on the scaling factors
        if am.setting["Simulation Configuration"]["Regional_CAPAX_scaling_flag"] == true
            for plant_id in keys(am.ref[:nw][sub_area_id][:gen_bus])
                
                if am.ref[:nw][sub_area_id][:gen_bus][plant_id]["Locational_Scaling_Flag"] == true

                    gen_tech_id = am.ref[:nw][sub_area_id][:gen_bus][plant_id]["Tech_ID"]
                    data_identifier = am.ref[:nw][0][:zone]["cost_scaling"]["data_identifier"]
                    data_region = am.ref[:nw][0][:bus][sub_area_id]["region_mapping_info"][data_identifier]
                    

                    regional_scaling = am.ref[:nw][0][:zone]["cost_scaling"][data_region][gen_tech_id]["CAPAX_scale"]

                    for year_idx in keys(am.ref[:nw][sub_area_id][:gen_bus][plant_id]["Annual_INVC"])
                        am.ref[:nw][sub_area_id][:gen_bus][plant_id]["Annual_INVC"][year_idx] *= regional_scaling
                    end

                    if haskey(am.ref[:nw][sub_area_id][:gen_bus][plant_id], "STO_CAPEX")
                        for year_idx in keys(am.ref[:nw][sub_area_id][:gen_bus][plant_id]["Annual_STO_INV"])
                            am.ref[:nw][sub_area_id][:gen_bus][plant_id]["Annual_STO_INV"][year_idx] *= regional_scaling
                        end
                    end
                end
            end
        end

        # delete sub-dicts that won't be used (no existing unit & no investment option)
        for gen_tech_id in keys(am.ref[:nw][sub_area_id][:gen_bus])
            if am.ref[:nw][sub_area_id][:gen_bus][gen_tech_id]["EXUNITS"] == 0.0
                if am.ref[:nw][sub_area_id][:gen_bus][gen_tech_id]["INVEST_FLAG"] == false
                    delete!(am.ref[:nw][sub_area_id][:gen_bus], gen_tech_id)
                end
            end
        end
        
    end
end


function update_individual_plant_info_GTEP_using_plant_info!(am::Abstract_ALEAF_Model, ALEAF_setting; recorded_investment_decisions=[])

    function find_gen_tech_info(am::Abstract_ALEAF_Model, gen_dict::Dict{String, Any}, sub_area_id::Int64)

        gen_tech_id = 0
        for idx in keys(am.ref[:nw][sub_area_id][:gen_technology])
            if (am.ref[:nw][sub_area_id][:gen_technology][idx]["UNITGROUP"] == gen_dict["UNITGROUP"]) 
                gen_tech_id = idx
            end
        end
    
        if gen_tech_id == 0
            gen_id = gen_dict["GEN UID"]
            @aleaf_info "ERROR: Failed to find gen_technology type; gen uid: $gen_id; temporal id is given"
            gen_tech_id = 1
        end
    
        return am.ref[:nw][sub_area_id][:gen_technology][gen_tech_id]        
    end

    function merge_missing!(plant_dict::Dict, gen_tech_info)
        
        for (k, v) in gen_tech_info
            get!(plant_dict, k, v)
        end
    end

    function find_ESGC_data(unit_group_id::String, data_label::String)
        for idx in keys(ALEAF_setting["Storage Cost and Performance"])
            if ALEAF_setting["Storage Cost and Performance"][idx]["ESGC_Setting_ID"] == ESGC_ID
                if ALEAF_setting["Storage Cost and Performance"][idx]["UNITGROUP"] == unit_group_id
                    return ALEAF_setting["Storage Cost and Performance"][idx][data_label]
                end
            end
        end
    end

    add_investment_option_flag = false 
    if (am.setting["Simulation Configuration"]["Run_expansion_flag"] == true) 
        add_investment_option_flag = true
    elseif (am.setting["Simulation Configuration"]["Run_RA_flag"] == true) && (am.setting["RA Setting"]["generate_gen_index_with_investment_options_flag"] == true)
        add_investment_option_flag = true
    elseif length(recorded_investment_decisions) > 0    # we need to add investment options when we have recorded investment data
        add_investment_option_flag = true
    end

    # index plants by their bus region id once, so per-sub_area gathering is O(regions+plants)
    # instead of re-scanning every plant per sub_area (was O(N_sub_areas x N_plants) at nodal)
    plants_by_region_bus_i = Dict{String, Vector{Any}}()
    for (plant_id, plant_data) in am.ref[:nw][0][:plant]
        push!(get!(plants_by_region_bus_i, string(plant_data["bus_ID"]), Vector{Any}()), plant_id)
    end

    for sub_area_id in keys(am.ref[:nw][0][:bus])    # the "bus" dict already completed the network reduction

        # add sub_area dict first
        am.ref[:nw][sub_area_id][:gen_bus] = Dict{Int64, Any}()

        # allocate plant data to each bus
        if length(am.ref[:nw][0][:plant]) > 0

            # copy plant data to each bus; gather from the region->plant index (set-identical to prior filter)
            regs = am.ref[:nw][0][:bus][sub_area_id]["aggregation_info"]["aggregated_regions_bus_i"]
            src = am.ref[:nw][0][:plant]
            am.ref[:nw][sub_area_id][:gen_bus] = Dict{keytype(src), valtype(src)}(pid => src[pid] for r in regs for pid in get(plants_by_region_bus_i, string(r), Any[]))

            for plant_id in keys(am.ref[:nw][sub_area_id][:gen_bus])

                plant_data = am.ref[:nw][sub_area_id][:gen_bus][plant_id]
           
                gen_tech_info = copy(find_gen_tech_info(am, plant_data, sub_area_id))
                merge_missing!(plant_data, gen_tech_info)

                original_cap = deepcopy(am.ref[:nw][sub_area_id][:gen_bus][plant_id]["CAP"])

                # add new Gen uid (bus_id + tech)
                gen_uid =  gen_tech_info["UNITGROUP"]                    
                new_gen_uid = string(sub_area_id, "_", gen_uid)                    
                am.ref[:nw][sub_area_id][:gen_bus][plant_id]["NEW_Gen_UID"] = new_gen_uid
    
                # update plant information
                am.ref[:nw][sub_area_id][:gen_bus][plant_id]["CAP"] = original_cap
                am.ref[:nw][sub_area_id][:gen_bus][plant_id]["EXCAPS"] = original_cap
                am.ref[:nw][sub_area_id][:gen_bus][plant_id]["EXUNITS"] = 1.0

                # add timeseries data tag
                am.ref[:nw][sub_area_id][:gen_bus][plant_id]["Timeseries_Tag"] = "LOCAL"
                
                # update retirement year if not available
                if am.ref[:nw][sub_area_id][:gen_bus][plant_id]["RetireYear"] == "NA"
                    online_year = am.ref[:nw][sub_area_id][:gen_bus][plant_id]["Online_Year"]
                    if online_year != "NA"
                        am.ref[:nw][sub_area_id][:gen_bus][plant_id]["RetireYear"] = online_year + gen_tech_info["Life"] - 1
                    else
                        am.ref[:nw][sub_area_id][:gen_bus][plant_id]["RetireYear"] = 9999
                    end
                end

                # add planned retirement tag
                am.ref[:nw][sub_area_id][:gen_bus][plant_id]["Planned_Retirement"] = Dict()
                for stage_id in keys(am.ref[:nw][0][:planning_stages])
                    year_at_y = am.ref[:nw][0][:planning_stages][stage_id]["year"]
                    stage_length = am.ref[:nw][0][:planning_stages][stage_id]["stage_length"]

                    am.ref[:nw][sub_area_id][:gen_bus][plant_id]["Planned_Retirement"][string("Ret_", year_at_y)] = 0.0

                    # add previous planned retirement to the first planning stage
                    retire_year = am.ref[:nw][sub_area_id][:gen_bus][plant_id]["RetireYear"]
                    if stage_id == 1
                        if retire_year <= year_at_y
                            am.ref[:nw][sub_area_id][:gen_bus][plant_id]["Planned_Retirement"][string("Ret_", year_at_y)] += original_cap
                        elseif retire_year > year_at_y && retire_year <= (year_at_y + stage_length - 1)
                            am.ref[:nw][sub_area_id][:gen_bus][plant_id]["Planned_Retirement"][string("Ret_", year_at_y)] += original_cap
                        end
                    else
                        if retire_year > year_at_y && retire_year <= (year_at_y + stage_length - 1)
                            am.ref[:nw][sub_area_id][:gen_bus][plant_id]["Planned_Retirement"][string("Ret_", year_at_y)] += original_cap
                        end
                    end
                end

                # update plant data
                if plant_data["bypass_parameter_check"] == false
                    check_and_update_plant_data!(am, ALEAF_setting, plant_data, gen_tech_info, parse(Int, am.ref[:nw][0][:case_id]))
                end

            end
        else
            am.ref[:nw][sub_area_id][:gen_bus] = Dict()
        end

        # Add new investment options in the expansion run 
        if add_investment_option_flag == true
            
            for idx in keys(am.ref[:nw][sub_area_id][:gen_technology])
                if am.ref[:nw][sub_area_id][:gen_technology][idx]["Tech_Type"] == "New"

                    if am.ref[:nw][sub_area_id][:gen_technology][idx]["INVEST_FLAG"] == true
                    
                        # new_idx (bus_id + techid)
                        new_genco_tech_id = parse(Int, string(sub_area_id) * "000" *string(idx))
                        
                        am.ref[:nw][sub_area_id][:gen_bus][new_genco_tech_id] = copy(am.ref[:nw][sub_area_id][:gen_technology][idx])

                        # add new Gen uid (bus_id + tech)
                        gen_uid =  am.ref[:nw][sub_area_id][:gen_bus][new_genco_tech_id]["UNITGROUP"]                    
                        new_gen_uid = string(sub_area_id, "_", gen_uid)                    
                        am.ref[:nw][sub_area_id][:gen_bus][new_genco_tech_id]["NEW_Gen_UID"] = new_gen_uid

                        # update plant information
                        am.ref[:nw][sub_area_id][:gen_bus][new_genco_tech_id]["EXCAPS"] = 0.0
                        am.ref[:nw][sub_area_id][:gen_bus][new_genco_tech_id]["EXUNITS"] = 0.0
                        am.ref[:nw][sub_area_id][:gen_bus][new_genco_tech_id]["ES_MWh"] = 0.0
                        
                        # add timeseries data tag
                        am.ref[:nw][sub_area_id][:gen_bus][new_genco_tech_id]["Timeseries_Tag"] = "LOCAL"

                        # add planned retirement tag
                        am.ref[:nw][sub_area_id][:gen_bus][new_genco_tech_id]["Planned_Retirement"] = Dict()
                        for stage_id in keys(am.ref[:nw][0][:planning_stages])
                            year = am.ref[:nw][0][:planning_stages][stage_id]["year"]
                            tag = string("Ret_", year)
                            am.ref[:nw][sub_area_id][:gen_bus][new_genco_tech_id]["Planned_Retirement"][tag] = 0.0
                        end

                        # add hybrid parameters
                        am.ref[:nw][sub_area_id][:gen_bus][new_genco_tech_id]["bypass_parameter_check"] = false # false => perform check_and_update_plant_data!
                        am.ref[:nw][sub_area_id][:gen_bus][new_genco_tech_id]["hybrid_type"] = "NA"
                        am.ref[:nw][sub_area_id][:gen_bus][new_genco_tech_id]["hybrid_ID"] = 0
                    end
                end
            end
        end

        # update CAPAX based on the scaling factors
        if am.setting["Simulation Configuration"]["Regional_CAPAX_scaling_flag"] == true
            for plant_id in keys(am.ref[:nw][sub_area_id][:gen_bus])
                
                if am.ref[:nw][sub_area_id][:gen_bus][plant_id]["Locational_Scaling_Flag"] == true

                    gen_tech_id = am.ref[:nw][sub_area_id][:gen_bus][plant_id]["Tech_ID"]
                    data_identifier = am.ref[:nw][0][:zone]["cost_scaling"]["data_identifier"]
                    data_region = am.ref[:nw][0][:bus][sub_area_id]["region_mapping_info"][data_identifier]
                    

                    regional_scaling = am.ref[:nw][0][:zone]["cost_scaling"][data_region][gen_tech_id]["CAPAX_scale"]

                    for year_idx in keys(am.ref[:nw][sub_area_id][:gen_bus][plant_id]["Annual_INVC"])
                        am.ref[:nw][sub_area_id][:gen_bus][plant_id]["Annual_INVC"][year_idx] *= regional_scaling
                    end

                    if haskey(am.ref[:nw][sub_area_id][:gen_bus][plant_id], "STO_CAPEX")
                        for year_idx in keys(am.ref[:nw][sub_area_id][:gen_bus][plant_id]["Annual_STO_INV"])
                            am.ref[:nw][sub_area_id][:gen_bus][plant_id]["Annual_STO_INV"][year_idx] *= regional_scaling
                        end
                    end
                end
            end
        end

        # delete sub-dicts that won't be used (no existing unit & no investment option)
        for gen_tech_id in keys(am.ref[:nw][sub_area_id][:gen_bus])
            if am.ref[:nw][sub_area_id][:gen_bus][gen_tech_id]["EXUNITS"] == 0.0
                if am.ref[:nw][sub_area_id][:gen_bus][gen_tech_id]["INVEST_FLAG"] == false
                    delete!(am.ref[:nw][sub_area_id][:gen_bus], gen_tech_id)
                end
            end
        end
        
    end
end


function add_local_gen_idx_GTEP!(am::Abstract_ALEAF_Model)
    
    ids_i = [(i) for (i) in get_index(am, :gen_index, 0)]

    # The no-config path still needs reserve/policy data-region lists for timeseries lookups.
    for n in keys(am.ref[:nw][0][:zone]["reserve"])
        aggregation_info = am.ref[:nw][0][:zone]["reserve"][n]["aggregation_info"]
        if !haskey(aggregation_info, "aggregated_data_regions_id")
            aggregation_info["aggregated_data_regions_id"] = deepcopy(get(am.ref[:nw][0][:zone]["reserve"][n], "data_area_list", String[]))
        end
    end
    for n in keys(am.ref[:nw][0][:zone]["policy"])
        aggregation_info = am.ref[:nw][0][:zone]["policy"][n]["aggregation_info"]
        if !haskey(aggregation_info, "aggregated_data_regions_id")
            aggregation_info["aggregated_data_regions_id"] = deepcopy(get(am.ref[:nw][0][:zone]["policy"][n], "data_area_list", String[]))
        end
    end

    for i in ids_i
        push!(am.ref[:nw][0][:bus][parameter(am, 0, :gen_index, "bus_idx", i)]["aggregation_info"]["new_local_gen_idx"], i)
    end

    
    # planning reserve zone
    for n in keys(am.ref[:nw][0][:zone]["planning_reserve"])
        
        for local_n in am.ref[:nw][0][:zone]["planning_reserve"][n]["aggregation_info"]["zone_bus_idx"]
            
            append!(am.ref[:nw][0][:zone]["planning_reserve"][n]["aggregation_info"]["new_local_gen_idx"], 
                    am.ref[:nw][0][:bus][local_n]["aggregation_info"]["new_local_gen_idx"])     
            
        end
    end
    
    
    # reserve zone
    for n in keys(am.ref[:nw][0][:zone]["reserve"])
        
        for local_n in am.ref[:nw][0][:zone]["reserve"][n]["aggregation_info"]["zone_bus_idx"]
            
            append!(am.ref[:nw][0][:zone]["reserve"][n]["aggregation_info"]["new_local_gen_idx"], 
                    am.ref[:nw][0][:bus][local_n]["aggregation_info"]["new_local_gen_idx"])     
            
        end
    end

    # policy zone
    for n in keys(am.ref[:nw][0][:zone]["policy"])
        
        for local_n in am.ref[:nw][0][:zone]["policy"][n]["aggregation_info"]["zone_bus_idx"]
            
            append!(am.ref[:nw][0][:zone]["policy"][n]["aggregation_info"]["new_local_gen_idx"], 
                    am.ref[:nw][0][:bus][local_n]["aggregation_info"]["new_local_gen_idx"])     
            
        end
    end

    # resource supply zone
    for n in am.ref[:nw][0][:zone]["resource_supply_curve"]["supply_curve_region_list"]

        for local_n in am.ref[:nw][0][:zone]["resource_supply_curve"][string(n)]["aggregation_info"]["zone_bus_idx"]
            
            append!(am.ref[:nw][0][:zone]["resource_supply_curve"][string(n)]["aggregation_info"]["new_local_gen_idx"], 
                    am.ref[:nw][0][:bus][local_n]["aggregation_info"]["new_local_gen_idx"])     
            
        end
    end

end


function add_reserve_zone_gen_technology_groups_GTEP!(am::Abstract_ALEAF_Model)

    reserve_group_lookup = Dict{Int64, Any}()
    reserve_group_reverse_lookup = Dict{Tuple{String, Int64}, Int64}()
    reserve_group_idx = 1

    for zone_idx in keys(am.ref[:nw][0][:zone]["reserve"])
        aggregation_info = am.ref[:nw][0][:zone]["reserve"][zone_idx]["aggregation_info"]

        technology_list = String[]
        gen_technology_groups = Dict{String, Any}()

        for gen_idx in aggregation_info["new_local_gen_idx"]
            unit_group = parameter(am, 0, :gen_index, "UNIT_GROUP", gen_idx)
            bus_idx = parameter(am, 0, :gen_index, "bus_idx", gen_idx)
            tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", gen_idx)
            PMAX = parameter(am, bus_idx, :gen_bus, tech_idx, "PMAX")
            MAXR = parameter(am, bus_idx, :gen_bus, tech_idx, "MAXR")
            RUL = parameter(am, bus_idx, :gen_bus, tech_idx, "RUL")
            RDL = parameter(am, bus_idx, :gen_bus, tech_idx, "RDL")
            MAXC = parameter(am, bus_idx, :gen_bus, tech_idx, "MAXC")
            CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")
            reserve_up_fraction = min(RUL, MAXC)
            reserve_down_fraction = min(RDL, MAXC)

            if !haskey(gen_technology_groups, unit_group)
                push!(technology_list, unit_group)
                gen_technology_groups[unit_group] = Dict{String, Any}(
                    "unit_group" => unit_group,
                    "UNIT_CATEGORY" => parameter(am, 0, :gen_index, "UNIT_CATEGORY", gen_idx),
                    "gen_idx" => Int[],
                    "bus_idx" => Int[],
                    "genco_tech_id" => Int[],
                    "total_capacity" => 0.0,
                    "regulation_fraction_numerator" => 0.0,
                    "regulation_fraction_denominator" => 0.0,
                    "regulation_fraction" => 0.0,
                    "reserve_up_fraction_numerator" => 0.0,
                    "reserve_up_fraction_denominator" => 0.0,
                    "reserve_up_fraction" => 0.0,
                    "reserve_down_fraction_numerator" => 0.0,
                    "reserve_down_fraction_denominator" => 0.0,
                    "reserve_down_fraction" => 0.0,
                )
            end

            group_data = gen_technology_groups[unit_group]
            push!(group_data["gen_idx"], gen_idx)
            push!(group_data["bus_idx"], bus_idx)
            push!(group_data["genco_tech_id"], tech_idx)
            group_data["total_capacity"] += CAP * PMAX
            group_data["regulation_fraction_numerator"] += CAP * PMAX * MAXR
            group_data["regulation_fraction_denominator"] += CAP * PMAX
            group_data["reserve_up_fraction_numerator"] += CAP * PMAX * reserve_up_fraction
            group_data["reserve_up_fraction_denominator"] += CAP * PMAX
            group_data["reserve_down_fraction_numerator"] += CAP * PMAX * reserve_down_fraction
            group_data["reserve_down_fraction_denominator"] += CAP * PMAX
        end

        for unit_group in keys(gen_technology_groups)
            group_data = gen_technology_groups[unit_group]
            denominator = group_data["regulation_fraction_denominator"]
            if denominator > 0.0
                group_data["regulation_fraction"] = group_data["regulation_fraction_numerator"] / denominator
            end

            denominator = group_data["reserve_up_fraction_denominator"]
            if denominator > 0.0
                group_data["reserve_up_fraction"] = group_data["reserve_up_fraction_numerator"] / denominator
            end

            denominator = group_data["reserve_down_fraction_denominator"]
            if denominator > 0.0
                group_data["reserve_down_fraction"] = group_data["reserve_down_fraction_numerator"] / denominator
            end
        end

        aggregation_info["technology_list"] = technology_list
        aggregation_info["gen_technology_groups"] = gen_technology_groups

        zone_idx_int = parse(Int64, zone_idx)
        for unit_group in technology_list
            reserve_group_lookup[reserve_group_idx] = Dict{String, Any}(
                "unit_group" => unit_group,
                "zone_idx" => zone_idx_int,
                "UNIT_CATEGORY" => aggregation_info["gen_technology_groups"][unit_group]["UNIT_CATEGORY"],
            )
            reserve_group_reverse_lookup[(unit_group, zone_idx_int)] = reserve_group_idx
            aggregation_info["gen_technology_groups"][unit_group]["reserve_group_idx"] = reserve_group_idx
            reserve_group_idx += 1
        end
    end

    am.ref[:nw][0][:reserve_group_lookup] = reserve_group_lookup
    am.ref[:nw][0][:reserve_group_reverse_lookup] = reserve_group_reverse_lookup

end


function add_model_setting_GTEP!(am::Abstract_ALEAF_Model, ALEAF_setting::Dict{String, Any}; nw::Int=am.cnw)

    case_id = am.data["case_id"]

    # Add LC_GEP setting to the model instance
    am.setting["Planning Design"] = deepcopy(ALEAF_setting["Planning Design"])
    am.setting["Simulation Configuration"] = deepcopy(ALEAF_setting["Simulation Configuration"][case_id])
    am.setting["Simulation Setting"] = deepcopy(ALEAF_setting["Simulation Setting"])
    am.setting["RA Setting"] = deepcopy(ALEAF_setting["RA Setting"])
    am.setting["File Path"] = deepcopy(ALEAF_setting["File Path"])
    am.setting["ALEAF Master Setup"] = deepcopy(ALEAF_setting["ALEAF Master Setup"])
    am.setting["Storage Cost and Performance"] = deepcopy(ALEAF_setting["Storage Cost and Performance"])
    am.setting["Network Setting"] = deepcopy(ALEAF_setting["Network Setting"])

    # External Constraints
    am.setting["External Constraints"] = Dict{String, Any}()
    EC_ID = ALEAF_setting["Simulation Configuration"][case_id]["External_Constraints_ID"]
    for idx in keys(ALEAF_setting["External Constraints"])
        if ALEAF_setting["External Constraints"][idx]["External_Constraints_ID"] == EC_ID
            am.setting["External Constraints"][idx] = ALEAF_setting["External Constraints"][idx]
        end
    end

    # Solver setting
    if ALEAF_setting["ALEAF Master Setup"]["solver_name"] == "CPLEX"
        am.setting["Solver Setting"] = ALEAF_setting["CPLEX Setting"]
        am.setting["Solver Setting"]["solver_name"] = "CPLEX"
        am.setting["Solver Setting"]["optimizer"] = cplex_optimizer_type(Val(:CPLEX))
    elseif ALEAF_setting["ALEAF Master Setup"]["solver_name"] == "HiGHS"
        am.setting["Solver Setting"] = ALEAF_setting["HiGHS Setting"]
        am.setting["Solver Setting"]["solver_name"] = "HiGHS"
        am.setting["Solver Setting"]["optimizer"] = HiGHS.Optimizer
    elseif ALEAF_setting["ALEAF Master Setup"]["solver_name"] == "MadNLP"
        am.setting["Solver Setting"] = ALEAF_setting["MadNLP Setting"]
        am.setting["Solver Setting"]["solver_name"] = "MadNLP"
        am.setting["Solver Setting"]["optimizer"] = gpu_setting_optimizer(Val(:MadNLP))
    elseif ALEAF_setting["ALEAF Master Setup"]["solver_name"] == "cuOpt"
        am.setting["Solver Setting"] = ALEAF_setting["cuOpt Setting"]
        am.setting["Solver Setting"]["solver_name"] = "cuOpt"
        am.setting["Solver Setting"]["optimizer"] = gpu_setting_optimizer(Val(:cuOpt))
    else
        error("Unsupported solver_name \"$(ALEAF_setting["ALEAF Master Setup"]["solver_name"])\"; supported solvers are HiGHS, CPLEX (optional), MadNLP, cuOpt")
    end
end


function add_run_period_GTEP!(am::Abstract_ALEAF_Model)

    "Hour"
    num_hours_per_day = am.setting["Planning Design"]["num_hours_per_day_value"]
    index_list = []
    for h = 1:num_hours_per_day
        if h%am.setting["Simulation Configuration"]["HFREQ"] == 0
            append!(index_list, h)
        end
    end
    am.setting["run_H"] = index_list

    "Sub-Hour"
    num_sub_period = am.setting["Planning Design"]["num_sub_period_value"]
    index_list = []
    if (am.setting["Simulation Configuration"]["FIVEMIN"] == 1) || (am.setting["Simulation Configuration"]["FIVEMIN_OP"] == 1)   
        indicator = 1
        for h = 1:num_sub_period
            if h%indicator == 0
                append!(index_list, h)
            end
        end
    else
        append!(index_list, 1)
    end
    am.setting["run_T"] = index_list
end


# Demand-side (Load Resource) contingency reserve provision variable; per-index upper bound is
# set in constraint_R_Cont_zdhty (0 when demand_reserve_provision_fraction=0 -> feature off).
function variable_demand_reserve_zdhty_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, category::Symbol, variable_name::String, decomp_group, ids_z, ids_d, ids_h, ids_t, ids_y; nw::Int=am.cnw, report::Bool=true)
    ids = [(z, d, h, t, y) for z in ids_z for d in ids_d for h in ids_h for t in ids_t for y in ids_y]
    var = am.var[:nw][nw][decomp_group][Symbol(variable_name)] = JuMP.@variable(JuMP_model,
        [idx in ids],
        base_name=variable_name,
        lower_bound=0.0,
        integer=false,
        binary=false
    )
    report && add_sol_component(am, nw, category, Symbol(variable_name), ids, var, decomp_group)
end


function build_LCO_GTEP_operational_mode_instance!(am::Abstract_ALEAF_Model, result_LC_GTEP_expansion::Dict{String,<:Any}, day_group_id::Int, y::Int; nw::Int=am.cnw, expansion_fix::Bool=true, operation_fix::Bool=false, result_operation::Dict{String,<:Any} = Dict{String,Any}())

    decomp_group = day_group_id  

    # Reset a JuMP model
    if (am.setting["Solver Setting"]["solver_name"] == "CPLEX") && (am.setting["Solver Setting"]["1"]["Value"] == true)
        
        JuMP_model = cplex_direct_model(Val(:CPLEX); pass_names = true)
    else
        JuMP_model = JuMP.Model()
    end

    const_name_flag = am.setting["Simulation Setting"]["const_name_flag"]
    if const_name_flag == false
        JuMP.set_string_names_on_creation(JuMP_model, false)
    end

    _t_build = time()   # detailed build timing: instance build start
    _detailed_log = lowercase(string(get(am.setting["Simulation Setting"], "logging_level_value", "simple"))) == "detailed"

    am.var = Dict{Symbol,Any}(:nw => Dict{Int,Any}())
    am.con = Dict{Symbol,Any}(:nw => Dict{Int,Any}())
    am.sol = Dict{Symbol,Any}(:nw => Dict{Int,Any}())

    for (nw_id, nw) in am.ref[:nw]
        am.var[:nw][nw_id] = Dict{Int,Any}()
        am.con[:nw][nw_id] = Dict{Int,Any}()
        am.sol[:nw][nw_id] = Dict{Int,Any}()
        for decomp_id in keys(am.ref[:nw][0][:repday_groups])
            am.var[:nw][nw_id][decomp_id] = Dict{Symbol,Any}()
            am.con[:nw][nw_id][decomp_id] = Dict{Symbol,Any}()
            am.sol[:nw][nw_id][decomp_id] = Dict{Symbol,Any}()
        end
    end

    #------ run for each day group
    # save the original day group information
    original_repday_group_info = am.ref[:nw][0][:repday_groups][day_group_id]

    original_Day_Idx_List = am.ref[:nw][0][:repday_groups][day_group_id]["Day_Idx_List"]

    # Adjust the day group with the given look-ahead days
    num_seq_day = 0

    start_day_id = original_Day_Idx_List[1]
    end_day_id = original_Day_Idx_List[end]
    
    if day_group_id != 1
        start_day_id = max(1, original_Day_Idx_List[1] - num_seq_day)
    else
        start_day_id = original_Day_Idx_List[1]
    end

    # Use the GLOBAL day-group count (config), not length(repday_groups): a distributed worker's
    # instance holds only its own day-group, so length()==1 would misidentify the last group once
    # look-ahead (num_seq_day) is ever nonzero. day_group_id is the global id, so this stays correct.
    total_num_day_groups = am.setting["Simulation Configuration"]["NDAY_Groups_OP"]
    if day_group_id != total_num_day_groups
        end_day_id = min(365, original_Day_Idx_List[end] + num_seq_day)
    else
        end_day_id = original_Day_Idx_List[end]
    end
   
    ids_d = [(d) for (d) in start_day_id:end_day_id]

    # update day list, start_day_id, and end_day_id 
    new_day_list = []
    for d in ids_d 
        push!(new_day_list, am.ref[:nw][0][:repdays][d]["Day"])
    end

    am.ref[:nw][0][:repday_groups][day_group_id]["Start_Day_Id"] = new_day_list[1]
    am.ref[:nw][0][:repday_groups][day_group_id]["End_Day_Id"] = new_day_list[end]
    am.ref[:nw][0][:repday_groups][day_group_id]["Day_Idx_List"] = ids_d
    am.ref[:nw][0][:repday_groups][day_group_id]["Day_List"] = new_day_list

    ids_y = [(y)]

    # Enhanced-hybrid OP physics: fold the realized expansion u* into susceptance so the
    # B-theta base flow sees the expanded admittance b0·(1+u*), consistent with the
    # rate_a·(1+u*) thermal caps. Exact (u* is fixed) — no relaxation downstream. Idempotent
    # via a stored base, since am.ref is shared across day groups within a year.
    if hybrid_tx_enabled(am) && expansion_fix && haskey(result_LC_GTEP_expansion, "solution")
        for k in get_index(am, :branch, 0)
            parameter(am, 0, :branch, "model_flag", k) == true || continue
            get(am.ref[:nw][0][:branch][k], "dc_line", false) != true || continue
            br = am.ref[:nw][0][:branch][k]
            haskey(br, "br_x_pu_hybrid_base") || (br["br_x_pu_hybrid_base"] = br["br_x_pu"])
            key = string("(", k, ", ", y, ")")
            u = haskey(result_LC_GTEP_expansion["solution"]["expansion"], key) ? result_LC_GTEP_expansion["solution"]["expansion"][key]["u_T_ky"] : 0.0
            br["br_x_pu"] = br["br_x_pu_hybrid_base"] / (1 + u)
        end
    end

    ###################################
    #------ Pre-processing
    ###################################

    ids_i = [(i) for (i) in get_index(am, :gen_index, 0)]
    ids_i_sto = [(i) for (i) in get_index(am, :gen_index, 0) if parameter(am, 0, :gen_index, "UNIT_CATEGORY", i) == "STORAGE"]
    
    ids_i_water = am.ref[:nw][0][:water_management]["gen_index"]

    ids_i_commit = Int64[]   # Resources with commitment decisions
    ids_i_hybrid_gen = Int64[]  # Hybrid: main gen (non ES)
    for i in ids_i
        bus_idx = am.ref[:nw][0][:gen_index][i]["bus_idx"]
        tech_idx = am.ref[:nw][0][:gen_index][i]["genco_tech_id"]
        hybrid_type = am.ref[:nw][0][:gen_index][i]["hybrid_type"]
        if parameter(am, 0, :gen_index, "UNIT_CATEGORY", i) != "STORAGE" # Exclude Storage
            if (parameter(am, bus_idx, :gen_bus, tech_idx, "Commitment") == true) && (parameter(am, bus_idx, :gen_bus, tech_idx, "Must_Run_Flag") == false)
                push!(ids_i_commit, i)
            end
        end
        if hybrid_type == "GEN"
            push!(ids_i_hybrid_gen, i)
        end
    end

    # Up-reserve eligible units (Dispatchable only): single source of truth for creating
    # reg_up/flex_up/spin variables, so non-dispatchable units never reference an unbuilt one.
    ids_i_up_res = [i for i in ids_i if up_reserve_eligible(am, i)]

    ids_h = [(h) for (h) in am.setting["run_H"]]
    ids_t = [(h) for (h) in am.setting["run_T"]]
    ids_k = [(k) for (k) in get_index(am, :branch, 0) if parameter(am, 0, :branch, "model_flag", k) == true]
    ids_n = [(k) for (k) in get_index(am, :bus, 0)]
    ids_z = [(z) for (z) in get_index(am, :zone, 0, "reserve")]
    ids_p = [(z) for (z) in get_index(am, :zone, 0, "policy")]
    ids_lfl = [(lfl) for (lfl) in get_index(am, :demand, 0)]    # Large Flexible Load

    # Aggregated reserve mode requires GTEP-embedded OP; standalone OP only supports "individual"
    if am.setting["Simulation Configuration"]["operating_reserve_modeling_option"] != "individual"
        error("aggregated reserve mode is not supported in standalone Operation; set operating_reserve_modeling_option = \"individual\"")
    end
    
    ###################################
    #------ DEFINE DECISION VARIABLES
    ###################################
    # Fixed Expansion variables (line expansion omitted: OP does not re-optimize transmission)
    variable_u_G_iy_real(JuMP_model, am, :expansion, "u_G_iy", decomp_group, ids_i, ids_y; fix_flag = expansion_fix, result = result_LC_GTEP_expansion)                             # total number of units in the system
    variable_u_ESE_iy_real(JuMP_model, am, :expansion, "u_ESE_iy", decomp_group, ids_i_sto, ids_y; fix_flag = expansion_fix, result = result_LC_GTEP_expansion)                     # total number of units in the system

    # Power flow variables
    if am.setting["Simulation Setting"]["power_flow_mode_flag"] in ["PTDF", "B-theta"]
        variable_f_kdhty_real(JuMP_model, am, :powerflow, "f_kdhty", decomp_group, ids_k, ids_d, ids_h, ids_t, ids_y, bounded=false)    # Power flow (unbounded for DC power flow)
    else # power_flow_mode_flag == "Network_Flow"
        variable_f_kdhty_real(JuMP_model, am, :powerflow, "f_kdhty", decomp_group, ids_k, ids_d, ids_h, ids_t, ids_y, bounded=false)                   # Power flow (bounded >=0 for network flow)
    end

    # Dispatch Variables
    variable_g_idhty_real(JuMP_model, am, :dispatch, "g_idhty", decomp_group, ids_i, ids_d, ids_h, ids_t, ids_y)                  # Unit Generation (MW)
    reserve_enabled(am, :reg_up)  && variable_reg_up_idhty_real(JuMP_model, am, :reserve, "reg_up_idhty", decomp_group, ids_i_up_res, ids_d, ids_h, ids_t, ids_y)        # Regulation up Reserves (MW) - dispatchable units only
    reserve_enabled(am, :reg_dn)  && variable_reg_dn_idhty_real(JuMP_model, am, :reserve, "reg_dn_idhty", decomp_group, ids_i, ids_d, ids_h, ids_t, ids_y)        # Regulation down Reserves (MW)
    reserve_enabled(am, :flex_up) && variable_flex_up_idhty_real(JuMP_model, am, :reserve, "flex_up_idhty", decomp_group, ids_i_up_res, ids_d, ids_h, ids_t, ids_y)      # Flex up Reserves (MW) - dispatchable units only
    reserve_enabled(am, :flex_dn) && variable_flex_dn_idhty_real(JuMP_model, am, :reserve, "flex_dn_idhty", decomp_group, ids_i, ids_d, ids_h, ids_t, ids_y)      # Flex down Reserves (MW)
    reserve_enabled(am, :spin)    && variable_spin_idhty_real(JuMP_model, am, :reserve, "spin_idhty", decomp_group, ids_i_up_res, ids_d, ids_h, ids_t, ids_y)            # Spinning Reserves (MW) - dispatchable units only
    reserve_enabled(am, :nonspin) && variable_nonspin_idhty_real(JuMP_model, am, :reserve, "nonspin_idhty", decomp_group, ids_i, ids_d, ids_h, ids_t, ids_y)            # Non-spinning Reserves (MW)

    # Scarcity variables
    ens_flag = am.setting["Simulation Configuration"]["Load_Shed_in_OP_Flag"]  # true = ens is allowed
    rns_flag = am.setting["Simulation Configuration"]["OR_Shortage_in_OP_Flag"]    # tru = rns in allowed

    variable_ens_ndhty_real(JuMP_model, am, :scarcity, "ens_ndhty", decomp_group, ids_n, ids_d, ids_h, ids_t, ids_y; bounded=ens_flag)                              # Energy Not Served (MW).
    reserve_enabled(am, :reg_up)  && variable_rns_reg_up_zdhty_real(JuMP_model, am, :scarcity_reserve, "rns_reg_up_zdhty", decomp_group, ids_z, ids_d, ids_h, ids_t, ids_y; bounded=rns_flag)            # Reg-up Reserves Not Served (MW)
    reserve_enabled(am, :reg_dn)  && variable_rns_reg_dn_zdhty_real(JuMP_model, am, :scarcity_reserve, "rns_reg_dn_zdhty", decomp_group, ids_z, ids_d, ids_h, ids_t, ids_y; bounded=rns_flag)            # Reg-down Reserves Not Served (MW)
    reserve_enabled(am, :spin)    && variable_rns_cont_zdhty_real(JuMP_model, am, :scarcity_reserve, "rns_cont_zdhty", decomp_group, ids_z, ids_d, ids_h, ids_t, ids_y; bounded=rns_flag)            # Spinning/Contingency Reserves Not Served (MW)
    reserve_enabled(am, :spin)    && variable_demand_reserve_zdhty_real(JuMP_model, am, :scarcity_reserve, "demand_reserve_zdhty", decomp_group, ids_z, ids_d, ids_h, ids_t, ids_y)            # Load-Resource contingency reserve provision (MW)
    reserve_enabled(am, :nonspin) && variable_rns_nonspin_zdhty_real(JuMP_model, am, :scarcity_reserve, "rns_nonspin_zdhty", decomp_group, ids_z, ids_d, ids_h, ids_t, ids_y; bounded=rns_flag)            # Non-spinning Reserves Not Served (MW)
    reserve_enabled(am, :flex_up) && variable_rns_flex_up_zdhty_real(JuMP_model, am, :scarcity_reserve, "rns_flex_up_zdhty", decomp_group, ids_z, ids_d, ids_h, ids_t, ids_y; bounded=rns_flag)      # Flex-up Reserves Not Served (MW)
    reserve_enabled(am, :flex_dn) && variable_rns_flex_dn_zdhty_real(JuMP_model, am, :scarcity_reserve, "rns_flex_dn_zdhty", decomp_group, ids_z, ids_d, ids_h, ids_t, ids_y; bounded=rns_flag)      # Flex-down Reserves Not Served (MW)

    # Storage
    variable_soc_idhty_real(JuMP_model, am, :storage, "soc_idhty", decomp_group, ids_i_sto, ids_d, ids_h, ids_t, ids_y)            # Storage charge level at end of the period
    variable_chg_idhty_real(JuMP_model, am, :storage, "chg_idhty", decomp_group, ids_i_sto, ids_d, ids_h, ids_t, ids_y)            # Storage charging MW

    # Commitment Variables
    if am.setting["Simulation Configuration"]["Dispatch_Mode_in_OP"] == "Unit Commitment"
        variable_c_idhty_integer(JuMP_model, am, :commitment, "c_idhty", decomp_group, ids_i_commit, ids_d, ids_h, ids_t, ids_y; fix_flag=operation_fix, result=result_operation)                     # Number of Units Committed 
        variable_su_idhty_real(JuMP_model, am, :commitment, "su_idhty", decomp_group, ids_i_commit, ids_d, ids_h, ids_t, ids_y)               # Number of Units Started Up 
    end

    variable_sto_c_idhty_integer(JuMP_model, am, :storage_commitment, "sto_c_idhty", decomp_group, ids_i, ids_d, ids_h, ids_t, ids_y; fix_flag=operation_fix, result=result_operation)               
    variable_water_use_jidhty_real(JuMP_model, am, :water_management, "water_use_ijdhty", decomp_group, ids_d, ids_h, ids_t, ids_y; fix_flag=operation_fix, result=result_operation)
    
    # Large Flexible Load
    define_variable_idhty_real(JuMP_model, am, :demand, "lfl_ldhty", decomp_group, ids_lfl, ids_d, ids_h, ids_t, ids_y; bounded_lower=true, lower_bound=0)  
    define_variable_idhty_real(JuMP_model, am, :demand, "lfl_DR_ldhty", decomp_group, ids_lfl, ids_d, ids_h, ids_t, ids_y; bounded_lower=true, lower_bound=0)  
    variable_lfl_seg_lsdhty_real(JuMP_model, am, :demand, "lfl_seg_lsdhty", decomp_group, ids_lfl, ids_d, ids_h, ids_t, ids_y; bounded_lower=true, lower_bound=0)  
    variable_lfl_ind_lsdhty_binary(JuMP_model, am, :demand, "lfl_ind_lsdhty", decomp_group, ids_lfl, ids_d, ids_h, ids_t, ids_y)  

    # Define LFL onsite hybrid gen/storage variables only for LFL nodes that have them; the
    # per-lfl constraints access them under the same Hybrid_Gen/Hybrid_ES != "NA" guard.
    ids_lfl_hybrid_gen = [lfl for lfl in ids_lfl if am.ref[:nw][0][:demand][lfl]["Hybrid_Gen"] != "NA"]
    ids_lfl_hybrid_es  = [lfl for lfl in ids_lfl if am.ref[:nw][0][:demand][lfl]["Hybrid_ES"]  != "NA"]

    if !isempty(ids_lfl_hybrid_gen)
        define_variable_idhty_real(JuMP_model, am, :demand, "lfl_g_G_LFL_ldhty", decomp_group, ids_lfl_hybrid_gen, ids_d, ids_h, ids_t, ids_y; bounded_lower=true, lower_bound=0)
        define_variable_idhty_real(JuMP_model, am, :demand, "lfl_g_G_Grid_ldhty", decomp_group, ids_lfl_hybrid_gen, ids_d, ids_h, ids_t, ids_y; bounded_lower=true, lower_bound=0)
        define_variable_idhty_real(JuMP_model, am, :demand, "lfl_g_G_ES_ldhty", decomp_group, ids_lfl_hybrid_gen, ids_d, ids_h, ids_t, ids_y; bounded_lower=true, lower_bound=0)
    end

    if !isempty(ids_lfl_hybrid_es)
        define_variable_idhty_real(JuMP_model, am, :demand, "lfl_g_ES_LFL_ldhty", decomp_group, ids_lfl_hybrid_es, ids_d, ids_h, ids_t, ids_y; bounded_lower=true, lower_bound=0)
        define_variable_idhty_real(JuMP_model, am, :demand, "lfl_g_ES_Grid_ldhty", decomp_group, ids_lfl_hybrid_es, ids_d, ids_h, ids_t, ids_y; bounded_lower=true, lower_bound=0)
        define_variable_idhty_real(JuMP_model, am, :demand, "lfl_chg_Grid_ES_ldhty", decomp_group, ids_lfl_hybrid_es, ids_d, ids_h, ids_t, ids_y; bounded_lower=true, lower_bound=0)
        define_variable_idhty_real(JuMP_model, am, :demand, "lfl_soc_ldhty", decomp_group, ids_lfl_hybrid_es, ids_d, ids_h, ids_t, ids_y; bounded_lower=true, lower_bound=0)
    end

    # Hybrid Generation to ES variable 
    define_variable_idhty_real(JuMP_model, am, :hybrid, "g_G_ES_idhty", decomp_group, ids_i_hybrid_gen, ids_d, ids_h, ids_t, ids_y; bounded_lower=true, lower_bound=0)  
    
    ###################################
    #------ DEFINE Constraints
    ###################################
    _detailed_log && @aleaf_info "[OP build] variables defined: nvar=$(JuMP.num_variables(JuMP_model)) elapsed=$(round(time()-_t_build, digits=1))s"   # detailed build timing

    # Unit constraints
    if expansion_fix == false   # we don't have expansion decisions in this case (i.e., operational run only)
    
        for i in ids_i
            constraint_u_G_balance_iy_OP_real(JuMP_model, am, "constraint_u_G_balance_iy_OP_real", decomp_group, i, y; const_name_flag)  # Balance number of units
        end

        for i in ids_i_sto
            constraint_u_ESE_balance_iy_OP_real(JuMP_model, am, "constraint_u_ESE_balance_iy_OP_real", decomp_group, i, y; const_name_flag)
        end
    end

    # System Balance & power flow
    if am.setting["Simulation Setting"]["power_flow_mode_flag"] == "PTDF"

        variable_p_inj_ndhty_real(JuMP_model, am, :powerflow, "p_inj_ndhty", decomp_group, ids_n, ids_d, ids_h, ids_t, ids_y)                              # power injection at node n
        variable_demand_ndhty_real(JuMP_model, am, :powerflow, "demand_ndhty", decomp_group, ids_n, ids_d, ids_h, ids_t, ids_y)                              # delivered demand

        for y in ids_y
            for d in ids_d, h in ids_h, t in ids_t
                constraint_sum_p_injdhty(JuMP_model, am, "constraint_sum_p_injdhty", decomp_group, d, h, t, y)
                
                for n in ids_n
                    constraint_PTDF_Power_injection_ndhty_real(JuMP_model, am, "constraint_PTDF_Power_injection_ndhty_real", decomp_group, n, d, h, t, y)
                    constraint_PTDF_LoadBalance_ndhty_real(JuMP_model, am, "constraint_PTDF_LoadBalance_ndhty_real", decomp_group, n, d, h, t, y; ens_flag)
                end
                for k in ids_k
                    constraint_dc_power_flow_kdhty(JuMP_model, am, "constraint_dc_power_flow_kdhty", decomp_group, k, d, h, t, y; const_name_flag)
                    constraint_dc_power_flow_max_kdhty_OP(JuMP_model, am, "constraint_dc_power_flow_max_kdhty_OP", decomp_group, k, d, h, t, y; const_name_flag, expansion_record_flag = expansion_fix, result = result_LC_GTEP_expansion)
                    constraint_dc_power_flow_min_kdhty_OP(JuMP_model, am, "constraint_dc_power_flow_min_kdhty_OP", decomp_group, k, d, h, t, y; const_name_flag, expansion_record_flag = expansion_fix, result = result_LC_GTEP_expansion)
                end
            end
        end
    elseif am.setting["Simulation Setting"]["power_flow_mode_flag"] == "Network_Flow"
        for d in ids_d, h in ids_h, t in ids_t
            for n in ids_n
                constraint_LoadBalance_ndhty_real(JuMP_model, am, "constraint_LoadBalance_ndhty_real", decomp_group, n, d, h, t, y; ens_flag)
            end
            for k in ids_k
                constraint_power_flow_max_kdhty_OP(JuMP_model, am, "constraint_power_flow_max_kdhty_OP", decomp_group, k, d, h, t, y; const_name_flag, expansion_record_flag = expansion_fix, result = result_LC_GTEP_expansion)
                constraint_power_flow_min_kdhty_OP(JuMP_model, am, "constraint_power_flow_min_kdhty_OP", decomp_group, k, d, h, t, y; const_name_flag, expansion_record_flag = expansion_fix, result = result_LC_GTEP_expansion)
            end
        end

    elseif am.setting["Simulation Setting"]["power_flow_mode_flag"] == "B-theta"

        # per-island slack: fix angle=0 at one reference bus per AC-connected component
        ac_ref_bus_ids = get_ac_reference_buses(am; nw=0)
        variable_bus_angle_ndhty_real(JuMP_model, am, :powerflow, "bus_angle_ndhty", decomp_group, ids_n, ids_d, ids_h, ids_t, ids_y; ref_bus_ids = ac_ref_bus_ids)                              # power injection at node n

        for d in ids_d, h in ids_h, t in ids_t
            for n in ids_n
                constraint_LoadBalance_ndhty_real(JuMP_model, am, "constraint_LoadBalance_ndhty_real", decomp_group, n, d, h, t, y; ens_flag)
            end

            for k in ids_k
                # DC ties do not couple bus angles; skip only the angle equation, keep flow limits
                if get(am.ref[:nw][0][:branch][k], "dc_line", false) != true
                    constraint_b_theta_power_flow_kdhty(JuMP_model, am, "constraint_b_theta_power_flow_kdhty", decomp_group, k, d, h, t, y; const_name_flag)
                end
                constraint_dc_power_flow_max_kdhty_OP(JuMP_model, am, "constraint_dc_power_flow_max_kdhty_OP", decomp_group, k, d, h, t, y; const_name_flag, expansion_record_flag = expansion_fix, result = result_LC_GTEP_expansion)
                constraint_dc_power_flow_min_kdhty_OP(JuMP_model, am, "constraint_dc_power_flow_min_kdhty_OP", decomp_group, k, d, h, t, y; const_name_flag, expansion_record_flag = expansion_fix, result = result_LC_GTEP_expansion)
            end
        end

    end

    _detailed_log && @aleaf_info "[OP build] power-flow + load-balance constraints done: ncon=$(JuMP.num_constraints(JuMP_model; count_variable_in_set_constraints=false)) elapsed=$(round(time()-_t_build, digits=1))s"   # detailed build timing

    # Operating reserve balance
    for d in ids_d, h in ids_h, t in ids_t
        for z in ids_z
            reserve_enabled(am, :reg_up) && constraint_reg_Up_zdhty(JuMP_model, am, "constraint_reg_Up_zdhty", decomp_group, z, d, h, t, y)
            reserve_enabled(am, :reg_dn) && constraint_reg_Dn_zdhty(JuMP_model, am, "constraint_reg_Dn_zdhty", decomp_group, z, d, h, t, y)

            reserve_enabled(am, :flex_up) && constraint_flex_Up_zdhty(JuMP_model, am, "constraint_flex_Up_zdhty", decomp_group, z, d, h, t, y; rns_flag)
            reserve_enabled(am, :flex_dn) && constraint_flex_Dn_zdhty(JuMP_model, am, "constraint_flex_Dn_zdhty", decomp_group, z, d, h, t, y; rns_flag)

            if am.setting["Simulation Configuration"]["Dispatch_Mode_in_OP"] == "Unit Commitment"
                reserve_enabled(am, :spin)    && constraint_R_Spin_zdhty(JuMP_model, am, "constraint_R_Spin_zdhty", decomp_group, z, d, h, t, y; rns_flag)
                reserve_enabled(am, :nonspin) && constraint_R_NonSpin_zdhty(JuMP_model, am, "constraint_R_NonSpin_zdhty", decomp_group, z, d, h, t, y; rns_flag)
            else
                reserve_enabled(am, :spin)    && constraint_R_Cont_zdhty(JuMP_model, am, "constraint_R_Cont_zdhty", decomp_group, z, d, h, t, y; rns_flag)
            end
        end
    end

    # Unit Dispatch and commitment
    for d in ids_d, h in ids_h, t in ids_t
        for i in ids_i

            unit_category = parameter(am, 0, :gen_index, "UNIT_CATEGORY", i)

            if (i in ids_i_commit) && (am.setting["Simulation Configuration"]["Dispatch_Mode_in_OP"] == "Unit Commitment") # if unit i requires a commitment decision
                
                # Unit commitment constraints
                constraint_Commit_Limit_idhty(JuMP_model, am, "constraint_Commit_Limit_idhty", decomp_group, i, d, h, t, y; const_name_flag)
                constraint_Start_Up_Limit_idhty(JuMP_model, am, "constraint_Start_Up_Limit_idhty", decomp_group, i, d, h, t, y; const_name_flag)
                constraint_Start_Up_Status_Dn_idhty(JuMP_model, am, "constraint_Start_Up_Status_Dn_idhty", decomp_group, i, d, h, t, y; const_name_flag)
                constraint_Start_Up_Status_Up_idhty(JuMP_model, am, "constraint_Start_Up_Status_Up_idhty", decomp_group, i, d, h, t, y; const_name_flag)

                # Unit dispatch constraints
                constraint_TherMin_Dispatch_UC_idhty(JuMP_model, am, "constraint_TherMin_Dispatch_UC_idhty", decomp_group, i, d, h, t, y; const_name_flag)
                constraint_TherMax_Dispatch_UC_idhty(JuMP_model, am, "constraint_TherMax_Dispatch_UC_idhty", decomp_group, i, d, h, t, y; const_name_flag)
                
                # Ramping constraints
                constraint_RampUpMax_UC_idhty(JuMP_model, am, "constraint_RampUpMax_UC_idhty", decomp_group, i, d, h, t, y; const_name_flag)
                constraint_RampDnMax_UC_idhty(JuMP_model, am, "constraint_RampDnMax_UC_idhty", decomp_group, i, d, h, t, y; const_name_flag)
                reserve_enabled(am, :reg_up) && up_reserve_eligible(am, i) && constraint_RegUpMax_UC_idhty(JuMP_model, am, "constraint_RegUpMax_UC_idhty", decomp_group, i, d, h, t, y; const_name_flag)
                reserve_enabled(am, :reg_dn)  && constraint_RegDnMax_UC_idhty(JuMP_model, am, "constraint_RegDnMax_UC_idhty", decomp_group, i, d, h, t, y; const_name_flag)
                reserve_enabled(am, :nonspin) && constraint_NonSpinMax_UC_idhty(JuMP_model, am, "constraint_NonSpinMax_UC_idhty", decomp_group, i, d, h, t, y; const_name_flag)

                # Inter-temporal Ramping Constraints (with commitment)
                if unit_category in ["THERMAL", "NUCLEAR", "OTHER"]
                    constraint_RampUp_InterTemporal_Hour_UC_idhy(JuMP_model, am, "constraint_RampUp_InterTemporal_Hour_UC_idhy", decomp_group, i, d, h, t, y; const_name_flag)
                    constraint_RampDn_InterTemporal_Hour_UC_idhy(JuMP_model, am, "constraint_RampDn_InterTemporal_Hour_UC_idhy", decomp_group, i, d, h, t, y; const_name_flag)
                end

            else
                
                # Unit dispatch constraint for thermal resources
                constraint_TherMinDispatch_idhty(JuMP_model, am, "constraint_TherMinDispatch_idhty", decomp_group, i, d, h, t, y; const_name_flag)
                constraint_TherMaxDispatch_idhty(JuMP_model, am, "constraint_TherMaxDispatch_idhty", decomp_group, i, d, h, t, y; const_name_flag)
                    
                # Inter-temporal Ramping Constraints for thermal resources (w/o commitment)
                if unit_category in ["THERMAL", "NUCLEAR", "OTHER"]
                    constraint_RampUp_InterTemporal_Hour_idhy_OP(JuMP_model, am, "constraint_RampUp_InterTemporal_Hour_idhy", decomp_group, i, d, h, t, y, day_group_id; const_name_flag)
                    constraint_RampDn_InterTemporal_Hour_idhy_OP(JuMP_model, am, "constraint_RampDn_InterTemporal_Hour_idhy", decomp_group, i, d, h, t, y, day_group_id; const_name_flag)
                end

                # Ramping constraints for all resources (including VRE and ES)
                constraint_RampUpMax_idhty(JuMP_model, am, "constraint_RampUpMax_idhty", decomp_group, i, d, h, t, y; const_name_flag)
                constraint_RampDnMax_idhty(JuMP_model, am, "constraint_RampDnMax_idhty", decomp_group, i, d, h, t, y; const_name_flag)
                reserve_enabled(am, :reg_up) && up_reserve_eligible(am, i) && constraint_RegUp_Max_idhty(JuMP_model, am, "constraint_RegUp_Max_idhty", decomp_group, i, d, h, t, y; const_name_flag)
                reserve_enabled(am, :reg_dn)  && constraint_RegDn_Max_idhty(JuMP_model, am, "constraint_RegDn_Max_idhty", decomp_group, i, d, h, t, y; const_name_flag)
                reserve_enabled(am, :nonspin) && constraint_NonSpinMax_idhty(JuMP_model, am, "constraint_NonSpinMax_idhty", decomp_group, i, d, h, t, y; const_name_flag)

            end
        end
    end


    _detailed_log && @aleaf_info "[OP build] reserve + unit dispatch/ramping constraints done: ncon=$(JuMP.num_constraints(JuMP_model; count_variable_in_set_constraints=false)) elapsed=$(round(time()-_t_build, digits=1))s"   # detailed build timing

    # VRE Balance
    for i in ids_i

        if !(i in ids_i_hybrid_gen) # Exclude hybrid main gen units
            bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
            tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)

            if parameter(am, bus_idx, :gen_bus, tech_idx, "VRE_Flag") == true
                unit_group = parameter(am, 0, :gen_index, "UNIT_GROUP", i)
                fuel_limit = parameter(am, bus_idx, :gen_bus, tech_idx, "FUEL_LIMIT")

                if (fuel_limit == "Fixed Profile") 
                    for d in ids_d, h in ids_h, t in ids_t
                        constraint_VREBalance_Fixed_Profile_idhty(JuMP_model, am, "constraint_VREBalance_Fixed_Profile_idhty", decomp_group, i, d, h, t, y, unit_group; const_name_flag)
                    end
                
                elseif (fuel_limit == "Budget (day groups)") || (fuel_limit == "Budget (annual)")

                    if am.setting["Simulation Configuration"]["Hydro_Flexibility_Flag"] == true
                        for d in ids_d, h in ids_h, t in ids_t
                            constraint_VREBalance_Fixed_Profile_wifh_Flexibilty_idhty(JuMP_model, am, "constraint_VREBalance_Fixed_Profile_wifh_Flexibilty_idhty", decomp_group, i, d, h, t, y, unit_group; const_name_flag)
                        end
                    end

                    if am.setting["Simulation Configuration"]["Hydro_Budget_Flag"] == true
                        constraint_VREBalance_Budget_iy(JuMP_model, am, "constraint_VREBalance_Budget_iy", decomp_group, i, day_group_id, y, unit_group; const_name_flag)
                    end

                    # When neither the flexibility nor the budget energy-limit mechanism is active, treat impoundment as ROR.
                    if (am.setting["Simulation Configuration"]["Hydro_Flexibility_Flag"] == false) && (am.setting["Simulation Configuration"]["Hydro_Budget_Flag"] == false)
                        for d in ids_d, h in ids_h, t in ids_t
                            constraint_VREBalance_Fixed_Profile_idhty(JuMP_model, am, "constraint_VREBalance_Fixed_Profile_idhty", decomp_group, i, d, h, t, y, unit_group; const_name_flag)
                        end
                    end

                end

            end
        end
    end

    # Storage Balance and Power Constraints
    for i in ids_i_sto
        hybrid_type = am.ref[:nw][0][:gen_index][i]["hybrid_type"] 

        if hybrid_type == "ES"      # limit hybrid ES grid charge
            constraint_hybrid_ES_Charge_from_Grid_idhty(JuMP_model, am, "constraint_hybrid_ES_Charge_from_Grid_idhty", decomp_group, i, ids_d, ids_h, ids_t, ids_y; const_name_flag)
        end
    
        for d in ids_d, h in ids_h, t in ids_t
            
            # charge limit
            if hybrid_type == "ES"
                constraint_hybrid_ES_Charge_Max_idhty(JuMP_model, am, "constraint_hybrid_ES_Charge_Max_idhty", decomp_group, i, d, h, t, y; const_name_flag)
            else
                constraint_ES_Charge_Max_idhty(JuMP_model, am, "constraint_ES_Charge_Max_idhty", decomp_group, i, d, h, t, y; const_name_flag)
            end
            
            # SOC bounds
            constraint_ES_SOC_Max_idhty(JuMP_model, am, "constraint_ES_SOC_Max_idhty", decomp_group, i, d, h, t, y; const_name_flag)
            constraint_ES_SOC_Min_idhty(JuMP_model, am, "constraint_ES_SOC_Min_idhty", decomp_group, i, d, h, t, y; const_name_flag)

            # Storage commitment constraints
            constraint_ES_Sto_UC_Limit_idhty(JuMP_model, am, "constraint_ES_Sto_UC_Limit_idhty", decomp_group, i, d, h, t, y; const_name_flag)
            constraint_ES_DisCharge_Max_Sto_UC_idhty(JuMP_model, am, "constraint_ES_DisCharge_Max_Sto_UC_idhty", decomp_group, i, d, h, t, y; const_name_flag)

            if hybrid_type == "ES"
                constraint_hybrid_ES_Charge_Max_Sto_UC_idhty(JuMP_model, am, "constraint_hybrid_ES_Charge_Max_Sto_UC_idhty", decomp_group, i, d, h, t, y; const_name_flag)
            else
                constraint_ES_Charge_Max_Sto_UC_idhty(JuMP_model, am, "constraint_ES_Charge_Max_Sto_UC_idhty", decomp_group, i, d, h, t, y; const_name_flag)
            end
            

        end
    end
    
    # Storage SOC Balance
    for d in am.ref[:nw][0][:repday_groups][day_group_id]["Day_Idx_List"]
        for i in ids_i_sto
            hybrid_type = am.ref[:nw][0][:gen_index][i]["hybrid_type"] 

            if hybrid_type == "ES"
                for h in ids_h, t in ids_t
                    constraint_hybrid_ES_SOC_Balance_Inter_Hour_idhty(JuMP_model, am, "constraint_hybrid_ES_SOC_Balance_Inter_Hour_idhty", decomp_group, i, d, h, t, y, day_group_id; const_name_flag)
                end
            else
                for h in ids_h, t in ids_t
                    constraint_ES_SOC_Balance_Inter_Hour_idhty(JuMP_model, am, "constraint_ES_SOC_Balance_Inter_Hour_idhty", decomp_group, i, d, h, t, y, day_group_id; const_name_flag)
                end 
            end
        end
    end
       
    # Storage SOC Neutral
    for d in am.ref[:nw][0][:repday_groups][day_group_id]["Day_Idx_List"]
        for i in ids_i_sto
            hybrid_type = am.ref[:nw][0][:gen_index][i]["hybrid_type"] 
            end_day = am.ref[:nw][0][:repday_groups][day_group_id]["End_Day_Id"]
            start_day_idx = am.ref[:nw][0][:repday_groups][day_group_id]["Day_Idx_List"][1]
            if am.ref[:nw][0][:repdays][d]["Day"] == end_day
                if hybrid_type == "ES"
                    constraint_hybrid_ES_SOC_Neutral_idy(JuMP_model, am, "constraint_hybrid_ES_SOC_Neutral_idy", decomp_group, i, d, y, start_day_idx; const_name_flag)
                else
                    constraint_ES_SOC_Neutral_idy(JuMP_model, am, "constraint_ES_SOC_Neutral_idy", decomp_group, i, d, y, start_day_idx; const_name_flag)
                end
            end
        end
    end

    # Water Management Constraints
    for i in ids_i_water
        for d in ids_d, h in ids_h, t in ids_t
            constraint_water_use_to_generation_idhty(JuMP_model, am, "constraint_water_use_to_generation_idhty", decomp_group, i, d, h, t, y; const_name_flag)
        end
        
        constraint_water_use_limit_per_daygroup_idy(JuMP_model, am, "constraint_water_use_limit_per_daygroup_idy", decomp_group, i, day_group_id, y; const_name_flag)
        constraint_water_use_limit_per_segment_idy(JuMP_model, am, "constraint_water_use_limit_per_segment_idy", decomp_group, i, day_group_id, y; const_name_flag)
        constraint_water_use_segment_order_idy(JuMP_model, am, "constraint_water_use_segment_order_idy", decomp_group, i, day_group_id, y; const_name_flag)
    end

    # Must run
    for i in ids_i, d in ids_d, h in ids_h, t in ids_t
        constraint_MustRun_idhty(JuMP_model, am, "constraint_MustRun_idhty", decomp_group, i, d, h, t, y; const_name_flag)
    end

    # Large Flexible Load
    for d in ids_d, h in ids_h, t in ids_t
        for lfl in ids_lfl

            constraint_LFL_power_balance_ldhty(JuMP_model, am, "constraint_LFL_power_balance_ldhty", decomp_group, lfl, d, h, t, y; const_name_flag)
            constraint_LFL_inter_connection_limit_injection_ldhty(JuMP_model, am, "constraint_LFL_inter_connection_limit_injection_ldhty", decomp_group, lfl, d, h, t, y; const_name_flag)
            constraint_LFL_inter_connection_limit_withdraw_ldhty(JuMP_model, am, "constraint_LFL_inter_connection_limit_withdraw_ldhty", decomp_group, lfl, d, h, t, y; const_name_flag)
            
            constraint_LFL_segment_bound_lsdhty(JuMP_model, am, "constraint_LFL_segment_bound_lsdhty", decomp_group, lfl, d, h, t, y; const_name_flag)
            constraint_LFL_segment_relation_lsdhty(JuMP_model, am, "constraint_LFL_segment_relation_lsdhty", decomp_group, lfl, d, h, t, y; const_name_flag)

            constraint_lfl_DR_balance_ldhty(JuMP_model, am, "constraint_lfl_DR_balance_ldhty", decomp_group, lfl, d, h, t, y; const_name_flag)

            if am.ref[:nw][0][:demand][lfl]["Hybrid_Gen"] != "NA"
                constraint_lfl_onsite_gen_thermal_cap_ldhty(JuMP_model, am, "constraint_lfl_onsite_gen_thermal_cap_ldhty", decomp_group, lfl, d, h, t, y; const_name_flag)
            end

            if am.ref[:nw][0][:demand][lfl]["Hybrid_ES"] != "NA"
                constraint_lfl_onsite_ES_charge_cap_ldhty(JuMP_model, am, "constraint_lfl_onsite_ES_charge_cap_ldhty", decomp_group, lfl, d, h, t, y; const_name_flag)
                constraint_lfl_onsite_ES_discharge_cap_ldhty(JuMP_model, am, "constraint_lfl_onsite_ES_discharge_cap_ldhty", decomp_group, lfl, d, h, t, y; const_name_flag)
                constraint_lfl_onsite_ES_SOC_cap_ldhty(JuMP_model, am, "constraint_lfl_onsite_ES_SOC_cap_ldhty", decomp_group, lfl, d, h, t, y; const_name_flag)
            end
        end
    end

    # Large Flexible Load DR Daily Limit
    for d in ids_d
        for lfl in ids_lfl
            constraint_lfl_DR_daily_limit_ldhty(JuMP_model, am, "constraint_lfl_DR_daily_limit_ldhty", decomp_group, lfl, d, y; const_name_flag)

        end
    end

    # Large Flexible Load Onsite Storage SOC Balance
    for lfl in ids_lfl
        if am.ref[:nw][0][:demand][lfl]["Hybrid_ES"] != "NA"
            for d in am.ref[:nw][0][:repday_groups][day_group_id]["Day_Idx_List"]
            
                for h in ids_h, t in ids_t
                    constraint_lfl_onsite_ES_SOC_balance_ldhty(JuMP_model, am, "constraint_lfl_onsite_ES_SOC_balance_ldhty", decomp_group, lfl, d, h, t, y, day_group_id; const_name_flag)
                end

                if am.ref[:nw][0][:repdays][d]["Day"] == am.ref[:nw][0][:repday_groups][day_group_id]["End_Day_Id"]
                    start_day_idx = am.ref[:nw][0][:repday_groups][day_group_id]["Day_Idx_List"][1]
                    constraint_lfl_onsite_ES_SOC_neutral_ldhty(JuMP_model, am, "constraint_lfl_onsite_ES_SOC_neutral_ldhty", decomp_group, lfl, d, y, start_day_idx; const_name_flag)
                end 

            end
        end
    end

    # Hybrid Plant Generation Constraint
    for i in ids_i_hybrid_gen
        for d in ids_d, h in ids_h, t in ids_t
            constraint_hybrid_onsite_gen_thermal_cap_ldhty(JuMP_model, am, "constraint_hybrid_onsite_gen_thermal_cap_ldhty", decomp_group, i, d, h, t, y; const_name_flag)
            constraint_hybrid_inter_connection_limit_injection_ldhty(JuMP_model, am, "constraint_hybrid_inter_connection_limit_injection_ldhty", decomp_group, i, d, h, t, y; const_name_flag)
            constraint_hybrid_inter_connection_limit_withdraw_ldhty(JuMP_model, am, "constraint_hybrid_inter_connection_limit_withdraw_ldhty", decomp_group, i, d, h, t, y; const_name_flag)
        end
    end

    # Policy Constraints

    # Storage annual energy throughput limits
    if am.setting["Simulation Configuration"]["Energy_Storage_AET_Limit_Flag"] == true
        for i in ids_i_sto
            constraint_ES_AET_OP_y(JuMP_model, am, "constraint_ES_AET_y", decomp_group, i, y, ids_d; const_name_flag)
        end
    end

    # Clean Energy Generation Constraint
    if am.setting["Simulation Configuration"]["Clean_Energy_Generation_Target_OP_Flag"] == true

        variable_slack_CEG_ndy_real(JuMP_model, am, :slack, "slack_CEG_ndy", decomp_group, ids_p, ids_d, ids_y)                              # Energy Not Served (MW).

        if am.ref[:nw][0][:planning_stages][y]["year"] >= am.setting["Simulation Configuration"]["Clean_Energy_Generation_Target_Start_Year"]
            for n in ids_p
                for d in ids_d
                    constraint_clean_energy_generation_ndy(JuMP_model, am, "constraint_clean_energy_generation_ndy", decomp_group, n, d, y)
                end
            end
        end
    end

    # System inertia
    if (am.setting["Planning Design"]["enforce_rotational_inertia_constraints_flag"] == true) && (am.setting["Simulation Configuration"]["Dispatch_Mode_in_OP"] == "Unit Commitment")
        for n in ids_p
            for d in ids_d, h in ids_h, t in ids_t
                constraint_inertia_ndhty(JuMP_model, am, "constraint_inertia_ndhty", decomp_group, n, d, h, t, y; const_name_flag)
            end
        end
    end


    # Ramping: Inter-hour Ramping Constraints (place holder for now)
    if am.setting["Simulation Configuration"]["FIVEMIN"] == 1      # Intra Hour (only applies with sub_hourly time resolution)
        for d in am.ref[:nw][0][:repday_groups][day_group_id]["Day_Idx_List"]
            for i in ids_i, h in ids_h
                unit_category = parameter(am, 0, :gen_index, "UNIT_CATEGORY", i)
                if unit_category in ["THERMAL", "NUCLEAR", "OTHER"]
                    for t in ids_t
                        constraint_RampUp_Inter_SubHour_idhy(JuMP_model, am, "constraint_RampUp_Inter_SubHour_idhy", decomp_group, i, d, h, t, y, day_group_id; const_name_flag) 
                        constraint_RampDn_Inter_SubHour_idhy(JuMP_model, am, "constraint_RampDn_Inter_SubHour_idhy", decomp_group, i, d, h, t, y, day_group_id; const_name_flag) 
                        constraint_ES_SOC_Balance_Inter_SubHour_idhty(JuMP_model, am, "constraint_ES_SOC_Balance_Inter_SubHour_idhty", decomp_group, i, d, h, t, y, day_group_id; const_name_flag)
                    end
                end
            end
        end
    end

    ###################################
    #------ DEFINE Objective
    ###################################
    _detailed_log && @aleaf_info "[OP build] all constraints done: ncon=$(JuMP.num_constraints(JuMP_model; count_variable_in_set_constraints=false)) elapsed=$(round(time()-_t_build, digits=1))s"   # detailed build timing
    objective_function_operation(JuMP_model, am, ids_i, ids_d, ids_h, ids_t, ids_y, ids_z, ids_p, ids_i_commit, ids_k, ids_i_water, ids_lfl, ids_n, day_group_id)

    am.model[:nw][nw][day_group_id] = JuMP_model
    _detailed_log && @aleaf_info "[OP build] objective built, model complete: nvar=$(JuMP.num_variables(JuMP_model)) ncon=$(JuMP.num_constraints(JuMP_model; count_variable_in_set_constraints=false)) total_build=$(round(time()-_t_build, digits=1))s"   # detailed build timing

    # return the original day group values
    am.ref[:nw][0][:repday_groups][day_group_id]= original_repday_group_info
    
end


function expansion_constraint_family_enabled(setting::Dict{String,<:Any}, family::String)

    diagnostics = get(setting, "Diagnostics", nothing)
    if !(diagnostics isa AbstractDict)
        return true
    end

    disabled_families = get(diagnostics, "disabled_constraint_families", String[])
    return !(family in string.(collect(disabled_families)))

end


function build_LCO_GTEP_expansion_mode_instance!(am::Abstract_ALEAF_Model, decomp_group::Int; nw::Int=am.cnw, fix_decision::Bool=false, result_data::Dict{String,<:Any} = Dict{String,Any}(), GTEP_multi_round_info_data::Dict{Any,<:Any} = Dict{Any,Any}(), ref_expansion_result::Dict{String,<:Any} = Dict{String,Any}())

    # Reset a JuMP model
    if (am.setting["Solver Setting"]["solver_name"] == "CPLEX") && (am.setting["Solver Setting"]["1"]["Value"] == true)
        
        JuMP_model = cplex_direct_model(Val(:CPLEX); pass_names = true)
    else
        JuMP_model = JuMP.Model()
    end

    # Define ids_d
    ids_d = [(d) for (d) in get_index(am, :repdays, 0)]
    ids_day_groups = [(c) for (c) in get_index(am, :repday_groups, 0)]

    # Define ids_y
    if am.setting["Planning Design"]["multi_round_solution_process_flag"] == true
        ids_y = GTEP_multi_round_info_data["round_ids_y"]
    else
        ids_y = [(y) for (y) in get_index(am, :planning_stages, 0)]
    end

    # Define set_string_names_on_creation
    const_name_flag = am.setting["Simulation Setting"]["const_name_flag"]
    if const_name_flag == false
        JuMP.set_string_names_on_creation(JuMP_model, false)
    end
       
    
    ###################################
    #------ Pre-processing
    ###################################
    ids_new_tech = [(i) for (i) in get_index(am, :gen_technology, 0) if parameter(am, 0, :gen_technology, "INVEST_FLAG", i) == true]

    ids_i = [(i) for (i) in get_index(am, :gen_index, 0)]
    ids_i_sto = [(i) for (i) in get_index(am, :gen_index, 0) if parameter(am, 0, :gen_index, "UNIT_CATEGORY", i) == "STORAGE"]
    ids_i_commit = Int64[]   # Resources with commitment decisions
    ids_i_hybrid_gen = Int64[]  # Hybrid: main gen (non ES)
    for i in ids_i
        bus_idx = am.ref[:nw][0][:gen_index][i]["bus_idx"]
        tech_idx = am.ref[:nw][0][:gen_index][i]["genco_tech_id"]
        hybrid_type = am.ref[:nw][0][:gen_index][i]["hybrid_type"]
        if parameter(am, 0, :gen_index, "UNIT_CATEGORY", i) != "STORAGE" # Exclude Storage
            if (parameter(am, bus_idx, :gen_bus, tech_idx, "Commitment") == true) && (parameter(am, bus_idx, :gen_bus, tech_idx, "Must_Run_Flag") == false)
                push!(ids_i_commit, i)
            end
        end
        if hybrid_type == "GEN"
            push!(ids_i_hybrid_gen, i)
        end
    end

    # Up-reserve eligible units (Dispatchable only): single source of truth for creating
    # reg_up/flex_up/spin variables, so non-dispatchable units never reference an unbuilt one.
    ids_i_up_res = [i for i in ids_i if up_reserve_eligible(am, i)]

    ids_h = [(h) for (h) in am.setting["run_H"]]
    ids_t = [(h) for (h) in am.setting["run_T"]]
    ids_k = [(k) for (k) in get_index(am, :branch, 0) if parameter(am, 0, :branch, "model_flag", k) == true]
    ids_n = [(k) for (k) in get_index(am, :bus, 0)]
    ids_z = [(z) for (z) in get_index(am, :zone, 0, "reserve")]
    ids_p = [(z) for (z) in get_index(am, :zone, 0, "policy")]
    ids_pr = [(z) for (z) in get_index(am, :zone, 0, "planning_reserve")]
    ids_m = [(m) for (m) in get_index(am, :raw_materials, 0)]
    ids_lfl = [(lfl) for (lfl) in get_index(am, :demand, 0)]    # Large Flexible Load

    PH_flag = false
    record_solution = true  # default is true
    
    # Reserve modeling options 
    reserve_modeling_option = am.setting["Simulation Configuration"]["operating_reserve_modeling_option"] # [aggregated, individual]
    ids_r = Int64[]
    if reserve_modeling_option == "aggregated" 
        ids_r = sort(collect(keys(am.ref[:nw][0][:reserve_group_lookup])))
    end
    
    ###################################
    #------ DEFINE DECISION VARIABLES
    ###################################
    # Expansion variables
    variable_u_newG_iy_integer_real(JuMP_model, am, :expansion, "u_new_G_iy", decomp_group, ids_i, ids_y; fix_flag=fix_decision, result=result_data, PH_val=PH_flag, ref_expansion_result)          # expansion decision to build or not build for a unit i
    variable_u_retG_iy_integer_real(JuMP_model, am, :expansion, "u_ret_G_iy", decomp_group, ids_i, ids_y; fix_flag=fix_decision, result=result_data, PH_val=false)               # retirement decision to retire or not for a unit i    
    variable_u_newESH_iy_integer_real(JuMP_model, am, :expansion, "u_new_ESH_iy", decomp_group, ids_i_sto, ids_y; fix_flag=false, result=result_data, PH_val=false)       # expansion decision for storage duration
    variable_u_newT_ky_real(JuMP_model, am, :expansion_line, "u_new_T_ky", decomp_group, ids_k, ids_y; fix_flag=fix_decision, result=result_data, PH_val=PH_flag, cycling_flag=true)                  # expansion decision for line k
    
    # Unit status variables
    variable_u_T_ky_real(JuMP_model, am, :expansion, "u_T_ky", decomp_group, ids_k, ids_y; report=record_solution)                             # total amount of new line expansion 
    variable_u_G_iy_real(JuMP_model, am, :expansion, "u_G_iy", decomp_group, ids_i, ids_y; report=record_solution)                             # total number of units in the system    
    variable_u_ESE_iy_real(JuMP_model, am, :expansion, "u_ESE_iy", decomp_group, ids_i_sto, ids_y; report=record_solution)                     # total number of units in the system
    
    # Power flow variables
    if am.setting["Simulation Setting"]["power_flow_mode_flag"] in ["PTDF", "B-theta"]
        variable_f_kdhty_real(JuMP_model, am, :powerflow, "f_kdhty", decomp_group, ids_k, ids_d, ids_h, ids_t, ids_y; report=record_solution, bounded=false)    # Power flow (unbounded for DC power flow)
    else # power_flow_mode_flag == "Network_Flow"
        variable_f_kdhty_real(JuMP_model, am, :powerflow, "f_kdhty", decomp_group, ids_k, ids_d, ids_h, ids_t, ids_y; report=record_solution, bounded=false)                   # Power flow (bounded >=0 for network flow)
    end

    # Enhanced-hybrid expansion-increment flow (B-theta only), registered for eligible AC corridors
    if hybrid_tx_enabled(am)
        ids_k_hybrid = [k for k in ids_k if hybrid_exp_branch(am, k)]
        variable_f_exp_kdhty_real(JuMP_model, am, :powerflow, "f_exp_kdhty", decomp_group, ids_k_hybrid, ids_d, ids_h, ids_t, ids_y; report=record_solution)
    end

    # Dispatch Variables
    variable_g_idhty_real(JuMP_model, am, :dispatch, "g_idhty", decomp_group, ids_i, ids_d, ids_h, ids_t, ids_y; report=record_solution)                  # Unit Generation (MW)

    # Reserves
    if reserve_modeling_option == "individual"
        reserve_enabled(am, :reg_up)  && variable_reg_up_idhty_real(JuMP_model, am, :reserve, "reg_up_idhty", decomp_group, ids_i_up_res, ids_d, ids_h, ids_t, ids_y; report=record_solution)        # Regulation up Reserves (MW) - dispatchable units only
        reserve_enabled(am, :reg_dn)  && variable_reg_dn_idhty_real(JuMP_model, am, :reserve, "reg_dn_idhty", decomp_group, ids_i, ids_d, ids_h, ids_t, ids_y; report=record_solution)        # Regulation down Reserves (MW)
        reserve_enabled(am, :flex_up) && variable_flex_up_idhty_real(JuMP_model, am, :reserve, "flex_up_idhty", decomp_group, ids_i_up_res, ids_d, ids_h, ids_t, ids_y; report=record_solution)      # Flex up Reserves (MW) - dispatchable units only
        reserve_enabled(am, :flex_dn) && variable_flex_dn_idhty_real(JuMP_model, am, :reserve, "flex_dn_idhty", decomp_group, ids_i, ids_d, ids_h, ids_t, ids_y; report=record_solution)      # Flex down Reserves (MW)
        reserve_enabled(am, :spin)    && variable_spin_idhty_real(JuMP_model, am, :reserve, "spin_idhty", decomp_group, ids_i_up_res, ids_d, ids_h, ids_t, ids_y; report=record_solution)            # Spinning Reserves (MW) - dispatchable units only
        if am.setting["Simulation Configuration"]["Dispatch_Mode_in_EXP"] == "Unit Commitment"
            reserve_enabled(am, :nonspin) && variable_nonspin_idhty_real_EXP(JuMP_model, am, :reserve, "nonspin_idhty", decomp_group, ids_i, ids_d, ids_h, ids_t, ids_y; report=record_solution)            # Non Spinning Reserves (MW)
        end
    elseif reserve_modeling_option == "aggregated"
        reserve_enabled(am, :reg_up)  && variable_reg_up_idhty_real(JuMP_model, am, :reserve, "reg_up_idhty", decomp_group, ids_r, ids_d, ids_h, ids_t, ids_y; report=record_solution)        # Regulation up Reserves (MW)
        reserve_enabled(am, :reg_dn)  && variable_reg_dn_idhty_real(JuMP_model, am, :reserve, "reg_dn_idhty", decomp_group, ids_r, ids_d, ids_h, ids_t, ids_y; report=record_solution)        # Regulation down Reserves (MW)
        reserve_enabled(am, :flex_up) && variable_flex_up_idhty_real(JuMP_model, am, :reserve, "flex_up_idhty", decomp_group, ids_r, ids_d, ids_h, ids_t, ids_y; report=record_solution)      # Flex up Reserves (MW)
        reserve_enabled(am, :flex_dn) && variable_flex_dn_idhty_real(JuMP_model, am, :reserve, "flex_dn_idhty", decomp_group, ids_r, ids_d, ids_h, ids_t, ids_y; report=record_solution)      # Flex down Reserves (MW)
        reserve_enabled(am, :spin)    && variable_spin_idhty_real(JuMP_model, am, :reserve, "spin_idhty", decomp_group, ids_r, ids_d, ids_h, ids_t, ids_y; report=record_solution)            # Spinning Reserves (MW)
        if am.setting["Simulation Configuration"]["Dispatch_Mode_in_EXP"] == "Unit Commitment"
            reserve_enabled(am, :nonspin) && variable_nonspin_idhty_real_EXP(JuMP_model, am, :reserve, "nonspin_idhty", decomp_group, ids_r, ids_d, ids_h, ids_t, ids_y; report=record_solution)            # Non Spinning Reserves (MW)
        end
    end

    # Scarcity variables
    ens_flag = am.setting["Simulation Configuration"]["Load_Shed_in_EXP_Flag"]  # true = ens is allowed
    rns_flag = am.setting["Simulation Configuration"]["OR_Shortage_in_EXP_Flag"]    # true = rns in allowed
    variable_ens_ndhty_real(JuMP_model, am, :scarcity, "ens_ndhty", decomp_group, ids_n, ids_d, ids_h, ids_t, ids_y; report=record_solution, bounded=ens_flag)           # Energy Not Served (MW).
    reserve_enabled(am, :reg_up)  && variable_rns_reg_up_zdhty_real(JuMP_model, am, :scarcity_reserve, "rns_reg_up_zdhty", decomp_group, ids_z, ids_d, ids_h, ids_t, ids_y; report=record_solution, bounded=rns_flag)            # Reg-up Reserves Not Served (MW)
    reserve_enabled(am, :reg_dn)  && variable_rns_reg_dn_zdhty_real(JuMP_model, am, :scarcity_reserve, "rns_reg_dn_zdhty", decomp_group, ids_z, ids_d, ids_h, ids_t, ids_y; report=record_solution, bounded=rns_flag)            # Reg-down Reserves Not Served (MW)
    reserve_enabled(am, :spin)    && variable_rns_cont_zdhty_real(JuMP_model, am, :scarcity_reserve, "rns_cont_zdhty", decomp_group, ids_z, ids_d, ids_h, ids_t, ids_y; report=record_solution, bounded=rns_flag)            # Spinning/Contingency Reserves Not Served (MW)
    reserve_enabled(am, :spin)    && variable_demand_reserve_zdhty_real(JuMP_model, am, :scarcity_reserve, "demand_reserve_zdhty", decomp_group, ids_z, ids_d, ids_h, ids_t, ids_y; report=record_solution)            # Load-Resource contingency reserve provision (MW)
    reserve_enabled(am, :flex_up) && variable_rns_flex_up_zdhty_real(JuMP_model, am, :scarcity_reserve, "rns_flex_up_zdhty", decomp_group, ids_z, ids_d, ids_h, ids_t, ids_y; report=record_solution, bounded=rns_flag)      # Flex-up Reserves Not Served (MW)
    reserve_enabled(am, :flex_dn) && variable_rns_flex_dn_zdhty_real(JuMP_model, am, :scarcity_reserve, "rns_flex_dn_zdhty", decomp_group, ids_z, ids_d, ids_h, ids_t, ids_y; report=record_solution, bounded=rns_flag)      # Flex-down Reserves Not Served (MW)
    # Non-spinning reserve shortfall exists only where non-spinning reserve is modeled in the expansion model (Unit Commitment dispatch)
    (am.setting["Simulation Configuration"]["Dispatch_Mode_in_EXP"] == "Unit Commitment") && reserve_enabled(am, :nonspin) && variable_rns_nonspin_zdhty_real(JuMP_model, am, :scarcity_reserve, "rns_nonspin_zdhty", decomp_group, ids_z, ids_d, ids_h, ids_t, ids_y; report=record_solution, bounded=rns_flag)

    # Storage
    variable_soc_idhty_real(JuMP_model, am, :storage, "soc_idhty", decomp_group, ids_i_sto, ids_d, ids_h, ids_t, ids_y; report=record_solution)            # Storage charge level at end of the period
    variable_chg_idhty_real(JuMP_model, am, :storage, "chg_idhty", decomp_group, ids_i_sto, ids_d, ids_h, ids_t, ids_y; report=record_solution)            # Storage charging MW
    variable_sto_c_idhty_integer(JuMP_model, am, :storage_commitment, "sto_c_idhty", decomp_group, ids_i, ids_d, ids_h, ids_t, ids_y; report=record_solution, fix_flag=fix_decision, result=result_data)      
    if am.setting["Simulation Configuration"]["Dispatch_Mode_in_EXP"] == "Unit Commitment"
        variable_c_idhty_integer(JuMP_model, am, :commitment, "c_idhty", decomp_group, ids_i_commit, ids_d, ids_h, ids_t, ids_y; report=record_solution, fix_flag=fix_decision, result=result_data)                     # Number of Units Committed 
        variable_su_idhty_real(JuMP_model, am, :commitment, "su_idhty", decomp_group, ids_i_commit, ids_d, ids_h, ids_t, ids_y; report=record_solution)               # Number of Units Started Up 
    end

    # Large Flexible Load (LFL)
    define_variable_idhty_real(JuMP_model, am, :demand, "lfl_ldhty", decomp_group, ids_lfl, ids_d, ids_h, ids_t, ids_y; bounded_lower=true, lower_bound=0)  
    define_variable_idhty_real(JuMP_model, am, :demand, "lfl_DR_ldhty", decomp_group, ids_lfl, ids_d, ids_h, ids_t, ids_y; bounded_lower=true, lower_bound=0)  
    variable_lfl_seg_lsdhty_real(JuMP_model, am, :demand, "lfl_seg_lsdhty", decomp_group, ids_lfl, ids_d, ids_h, ids_t, ids_y; bounded_lower=true, lower_bound=0)  
    variable_lfl_ind_lsdhty_binary(JuMP_model, am, :demand, "lfl_ind_lsdhty", decomp_group, ids_lfl, ids_d, ids_h, ids_t, ids_y)  

    # Define hybrid gen/storage variables only for LFL nodes that have them; the per-lfl
    # constraints access them under the same Hybrid_Gen/Hybrid_ES != "NA" guard.
    ids_lfl_hybrid_gen = [lfl for lfl in ids_lfl if am.ref[:nw][0][:demand][lfl]["Hybrid_Gen"] != "NA"]
    ids_lfl_hybrid_es  = [lfl for lfl in ids_lfl if am.ref[:nw][0][:demand][lfl]["Hybrid_ES"]  != "NA"]

    if !isempty(ids_lfl_hybrid_gen)
        define_variable_idhty_real(JuMP_model, am, :demand, "lfl_g_G_LFL_ldhty", decomp_group, ids_lfl_hybrid_gen, ids_d, ids_h, ids_t, ids_y; bounded_lower=true, lower_bound=0)
        define_variable_idhty_real(JuMP_model, am, :demand, "lfl_g_G_Grid_ldhty", decomp_group, ids_lfl_hybrid_gen, ids_d, ids_h, ids_t, ids_y; bounded_lower=true, lower_bound=0)
        define_variable_idhty_real(JuMP_model, am, :demand, "lfl_g_G_ES_ldhty", decomp_group, ids_lfl_hybrid_gen, ids_d, ids_h, ids_t, ids_y; bounded_lower=true, lower_bound=0)
    end

    if !isempty(ids_lfl_hybrid_es)
        define_variable_idhty_real(JuMP_model, am, :demand, "lfl_g_ES_LFL_ldhty", decomp_group, ids_lfl_hybrid_es, ids_d, ids_h, ids_t, ids_y; bounded_lower=true, lower_bound=0)
        define_variable_idhty_real(JuMP_model, am, :demand, "lfl_g_ES_Grid_ldhty", decomp_group, ids_lfl_hybrid_es, ids_d, ids_h, ids_t, ids_y; bounded_lower=true, lower_bound=0)
        define_variable_idhty_real(JuMP_model, am, :demand, "lfl_chg_Grid_ES_ldhty", decomp_group, ids_lfl_hybrid_es, ids_d, ids_h, ids_t, ids_y; bounded_lower=true, lower_bound=0)
        define_variable_idhty_real(JuMP_model, am, :demand, "lfl_soc_ldhty", decomp_group, ids_lfl_hybrid_es, ids_d, ids_h, ids_t, ids_y; bounded_lower=true, lower_bound=0)
    end

    # Hybrid Generation to ES variable 
    define_variable_idhty_real(JuMP_model, am, :hybrid, "g_G_ES_idhty", decomp_group, ids_i_hybrid_gen, ids_d, ids_h, ids_t, ids_y; bounded_lower=true, lower_bound=0)

    ###################################
    #------ DEFINE Constraints
    ###################################

    # External Constraints
    for external_const_id in keys(am.setting["External Constraints"])
        if (am.setting["External Constraints"][external_const_id]["Model"] == "Expansion") && (am.setting["External Constraints"][external_const_id]["Apply_Flag"] == true)
            const_name = string("Expansion_External_Constraint_", external_const_id)
            external_constraint_generator(JuMP_model, am, const_name, decomp_group, external_const_id, ids_y; const_name_flag)
        end
    end
    
    # System-wide investment constraints per technology
    for tech in ids_new_tech
        constraint_system_total_investment_tech_real(JuMP_model, am, "constraint_system_total_investment_tech_real", decomp_group, tech, ids_i, ids_y; const_name_flag, GTEP_multi_round_info = GTEP_multi_round_info_data)  
    end

    # Expansion Constraints
    for i in ids_i
        for y in ids_y
            constraint_planned_ret_G_balance_iy_real(JuMP_model, am, "constraint_planned_ret_G_balance_i_real", decomp_group, i, y; const_name_flag, GTEP_multi_round_info = GTEP_multi_round_info_data)
        end
    end

    # Storage duration constraints
    for y in ids_y
        for i in ids_i_sto
            constraint_ESH_investment_min_iy_real(JuMP_model, am, "constraint_ESH_investment_min_iy_real", decomp_group, i, y; const_name_flag)
        end
    end

    # Transmission expansion constraints
    for y in ids_y
        for k in ids_k
            constraint_u_newT_ky_real(JuMP_model, am, "constraint_u_newT_ky_real", decomp_group, k, y; const_name_flag, fix_flag=fix_decision, GTEP_multi_round_info = GTEP_multi_round_info_data)
        end
    end

    # Resource balance
    for y in ids_y
        for i in ids_i
            constraint_u_G_balance_iy_real(JuMP_model, am, "constraint_u_G_balance_iy_real", decomp_group, i, y; const_name_flag, GTEP_multi_round_info = GTEP_multi_round_info_data)  # Balance number of units
        end
    
        for i in ids_i_sto
            constraint_u_ESE_balance_iy_real(JuMP_model, am, "constraint_u_ESE_balance_iy_real", decomp_group, i, y; const_name_flag, GTEP_multi_round_info = GTEP_multi_round_info_data)
        end
    end

    # Resource Supply Curve (Limits)
    if am.setting["Simulation Configuration"]["Regional_resource_limits_flag"] == true
        
        resource_zone = [(z) for (z) in am.ref[:nw][0][:zone]["resource_supply_curve"]["supply_curve_region_list"]]
        resource_list = am.ref[:nw][0][:zone]["resource_supply_curve"]["supply_curve_resource_list"]

        for resource_zone_id in resource_zone
            for resource_idx in eachindex(resource_list) 
                for y in ids_y
                    constraint_resource_limit_per_tech(JuMP_model, am, "constraint_resource_limit_per_tech", decomp_group, y, resource_zone_id, resource_idx, resource_list[resource_idx]; const_name_flag, GTEP_multi_round_info = GTEP_multi_round_info_data)
                end
            end
        end

    end

    # System Balance & power flow
    if am.setting["Simulation Setting"]["power_flow_mode_flag"] == "PTDF"
        
        variable_p_inj_ndhty_real(JuMP_model, am, :powerflow, "p_inj_ndhty", decomp_group, ids_n, ids_d, ids_h, ids_t, ids_y)                              # power injection at node n
        variable_demand_ndhty_real(JuMP_model, am, :powerflow, "demand_ndhty", decomp_group, ids_n, ids_d, ids_h, ids_t, ids_y)                              # delivered demand

        for y in ids_y
            for d in ids_d, h in ids_h, t in ids_t

                constraint_sum_p_injdhty(JuMP_model, am, "constraint_sum_p_injdhty", decomp_group, d, h, t, y)

                for n in ids_n
                    constraint_PTDF_Power_injection_ndhty_real(JuMP_model, am, "constraint_PTDF_Power_injection_ndhty_real", decomp_group, n, d, h, t, y)
                    constraint_PTDF_LoadBalance_ndhty_real(JuMP_model, am, "constraint_PTDF_LoadBalance_ndhty_real", decomp_group, n, d, h, t, y; ens_flag)
                end

                for k in ids_k
                    constraint_dc_power_flow_kdhty(JuMP_model, am, "constraint_dc_power_flow_kdhty", decomp_group, k, d, h, t, y)
                    constraint_dc_power_flow_max_kdhty(JuMP_model, am, "constraint_dc_power_flow_max_kdhty", decomp_group, k, d, h, t, y; const_name_flag, GTEP_multi_round_info = GTEP_multi_round_info_data)
                    constraint_dc_power_flow_min_kdhty(JuMP_model, am, "constraint_dc_power_flow_min_kdhty", decomp_group, k, d, h, t, y; const_name_flag, GTEP_multi_round_info = GTEP_multi_round_info_data)
                end
            end
        end
    elseif am.setting["Simulation Setting"]["power_flow_mode_flag"] == "Network_Flow"
        for y in ids_y
            for d in ids_d, h in ids_h, t in ids_t
                for n in ids_n
                    constraint_LoadBalance_ndhty_real(JuMP_model, am, "constraint_LoadBalance_ndhty_real", decomp_group, n, d, h, t, y; ens_flag)
                end
            
                for k in ids_k
                    constraint_power_flow_max_kdhty(JuMP_model, am, "constraint_power_flow_max_kdhty", decomp_group, k, d, h, t, y; const_name_flag, GTEP_multi_round_info = GTEP_multi_round_info_data)
                    constraint_power_flow_min_kdhty(JuMP_model, am, "constraint_power_flow_min_kdhty", decomp_group, k, d, h, t, y; const_name_flag, GTEP_multi_round_info = GTEP_multi_round_info_data)
                end
            end
        end
    elseif am.setting["Simulation Setting"]["power_flow_mode_flag"] == "B-theta"

        # per-island slack: fix angle=0 at one reference bus per AC-connected component
        ac_ref_bus_ids = get_ac_reference_buses(am; nw=0)
        variable_bus_angle_ndhty_real(JuMP_model, am, :powerflow, "bus_angle_ndhty", decomp_group, ids_n, ids_d, ids_h, ids_t, ids_y; ref_bus_ids = ac_ref_bus_ids)                              # power injection at node n

        for y in ids_y
            for d in ids_d, h in ids_h, t in ids_t
                for n in ids_n
                    constraint_LoadBalance_ndhty_real(JuMP_model, am, "constraint_LoadBalance_ndhty_real", decomp_group, n, d, h, t, y; ens_flag)
                end

                for k in ids_k
                    # DC ties do not couple bus angles; skip only the angle equation, keep flow limits
                    if get(am.ref[:nw][0][:branch][k], "dc_line", false) != true
                        constraint_b_theta_power_flow_kdhty(JuMP_model, am, "constraint_b_theta_power_flow_kdhty", decomp_group, k, d, h, t, y; const_name_flag)
                    end
                    if hybrid_exp_branch(am, k)
                        # Split flow: KVL base cap at rate_a + increment fenced by the McCormick hull
                        constraint_hybrid_base_cap_kdhty(JuMP_model, am, "constraint_hybrid_base_cap_kdhty", decomp_group, k, d, h, t, y; const_name_flag)
                        constraint_hybrid_exp_capacity_kdhty(JuMP_model, am, "constraint_hybrid_exp_capacity_kdhty", decomp_group, k, d, h, t, y; const_name_flag)
                        constraint_hybrid_exp_coupling_kdhty(JuMP_model, am, "constraint_hybrid_exp_coupling_kdhty", decomp_group, k, d, h, t, y; const_name_flag)
                    else
                        constraint_dc_power_flow_max_kdhty(JuMP_model, am, "constraint_dc_power_flow_max_kdhty", decomp_group, k, d, h, t, y; const_name_flag, GTEP_multi_round_info = GTEP_multi_round_info_data)
                        constraint_dc_power_flow_min_kdhty(JuMP_model, am, "constraint_dc_power_flow_min_kdhty", decomp_group, k, d, h, t, y; const_name_flag, GTEP_multi_round_info = GTEP_multi_round_info_data)
                    end
                end
            end
        end

    end
    
    # Operating reserve balance
    for y in ids_y
        for d in ids_d, h in ids_h, t in ids_t
            for z in ids_z
                reserve_enabled(am, :reg_up)  && constraint_reg_Up_zdhty(JuMP_model, am, "constraint_reg_Up_zdhty", decomp_group, z, d, h, t, y; rns_flag)
                reserve_enabled(am, :reg_dn)  && constraint_reg_Dn_zdhty(JuMP_model, am, "constraint_reg_Dn_zdhty", decomp_group, z, d, h, t, y; rns_flag)
                reserve_enabled(am, :flex_up) && constraint_flex_Up_zdhty(JuMP_model, am, "constraint_flex_Up_zdhty", decomp_group, z, d, h, t, y; rns_flag)
                reserve_enabled(am, :flex_dn) && constraint_flex_Dn_zdhty(JuMP_model, am, "constraint_flex_Dn_zdhty", decomp_group, z, d, h, t, y; rns_flag)

                if am.setting["Simulation Configuration"]["Dispatch_Mode_in_EXP"] == "Unit Commitment"
                    reserve_enabled(am, :spin)    && constraint_R_Spin_zdhty(JuMP_model, am, "constraint_R_Spin_zdhty", decomp_group, z, d, h, t, y; rns_flag)
                    reserve_enabled(am, :nonspin) && constraint_R_NonSpin_zdhty(JuMP_model, am, "constraint_R_NonSpin_zdhty", decomp_group, z, d, h, t, y; rns_flag)
                else
                    reserve_enabled(am, :spin)    && constraint_R_Cont_zdhty(JuMP_model, am, "constraint_R_Cont_zdhty", decomp_group, z, d, h, t, y; rns_flag)
                end
            end

            if reserve_modeling_option == "aggregated"
                for r in ids_r  # aggregated reserve headroom
                    constraint_TherMax_Agg_Res_idhty(JuMP_model, am, "constraint_TherMax_Agg_Res_idhty", decomp_group, r, d, h, t, y)
                    constraint_Agg_RampUpMax_idhty(JuMP_model, am, "constraint_Agg_RampUpMax_idhty", decomp_group, r, d, h, t, y)
                    constraint_Agg_RampDnMax_idhty(JuMP_model, am, "constraint_Agg_RampDnMax_idhty", decomp_group, r, d, h, t, y)
                    reserve_enabled(am, :reg_up) && constraint_Agg_RegUp_Max_idhty(JuMP_model, am, "constraint_Agg_RegUp_Max_idhty", decomp_group, r, d, h, t, y)
                    reserve_enabled(am, :reg_dn) && constraint_Agg_RegDn_Max_idhty(JuMP_model, am, "constraint_Agg_RegDn_Max_idhty", decomp_group, r, d, h, t, y)

                    if am.ref[:nw][0][:reserve_group_lookup][r]["UNIT_CATEGORY"] == "STORAGE" 
                        constraint_agg_ES_Charge_Max_idhty(JuMP_model, am, "constraint_agg_ES_Charge_Max_idhty", decomp_group, r, d, h, t, y; const_name_flag)
                        constraint_agg_ES_SOC_Max_idhty(JuMP_model, am, "constraint_agg_ES_SOC_Max_idhty", decomp_group, r, d, h, t, y; const_name_flag)
                        constraint_agg_ES_SOC_Min_idhty(JuMP_model, am, "constraint_agg_ES_SOC_Min_idhty", decomp_group, r, d, h, t, y; const_name_flag)
                    end
                end
            end
        end
    end

    # Unit Dispatch and commitment
    for y in ids_y
        for d in ids_d, h in ids_h, t in ids_t
            for i in ids_i

                if (i in ids_i_commit) && (am.setting["Simulation Configuration"]["Dispatch_Mode_in_EXP"] == "Unit Commitment") # if unit i requires a commitment decision
                    
                    # Unit commitment constraints
                    constraint_Commit_Limit_idhty(JuMP_model, am, "constraint_Commit_Limit_idhty", decomp_group, i, d, h, t, y; const_name_flag)
                    constraint_Start_Up_Limit_idhty(JuMP_model, am, "constraint_Start_Up_Limit_idhty", decomp_group, i, d, h, t, y; const_name_flag)
                    constraint_Start_Up_Status_Dn_idhty(JuMP_model, am, "constraint_Start_Up_Status_Dn_idhty", decomp_group, i, d, h, t, y; const_name_flag)
                    constraint_Start_Up_Status_Up_idhty(JuMP_model, am, "constraint_Start_Up_Status_Up_idhty", decomp_group, i, d, h, t, y; const_name_flag)

                    # Unit dispatch constraints
                    constraint_TherMin_Dispatch_UC_idhty(JuMP_model, am, "constraint_TherMin_Dispatch_UC_idhty", decomp_group, i, d, h, t, y; const_name_flag)
                    constraint_TherMax_Dispatch_UC_idhty(JuMP_model, am, "constraint_TherMax_Dispatch_UC_idhty", decomp_group, i, d, h, t, y; const_name_flag)
                    
                    # Ramping constraints
                    constraint_RampUpMax_UC_idhty(JuMP_model, am, "constraint_RampUpMax_UC_idhty", decomp_group, i, d, h, t, y; const_name_flag)
                    constraint_RampDnMax_UC_idhty(JuMP_model, am, "constraint_RampDnMax_UC_idhty", decomp_group, i, d, h, t, y; const_name_flag)
                    reserve_enabled(am, :reg_up) && up_reserve_eligible(am, i) && constraint_RegUpMax_UC_idhty(JuMP_model, am, "constraint_RegUpMax_UC_idhty", decomp_group, i, d, h, t, y; const_name_flag)
                    reserve_enabled(am, :reg_dn)  && constraint_RegDnMax_UC_idhty(JuMP_model, am, "constraint_RegDnMax_UC_idhty", decomp_group, i, d, h, t, y; const_name_flag)
                    reserve_enabled(am, :nonspin) && constraint_NonSpinMax_UC_idhty(JuMP_model, am, "constraint_NonSpinMax_UC_idhty", decomp_group, i, d, h, t, y; const_name_flag)
                else
                    
                    # Unit dispatch constraint for thermal resources
                    constraint_TherMaxDispatch_idhty(JuMP_model, am, "constraint_TherMaxDispatch_idhty", decomp_group, i, d, h, t, y; const_name_flag)
                    if reserve_modeling_option == "individual"
                        constraint_TherMinDispatch_idhty(JuMP_model, am, "constraint_TherMinDispatch_idhty", decomp_group, i, d, h, t, y; const_name_flag)
                    end
                        
                    # Ramping constraints for all resources (including VRE and ES)
                    if reserve_modeling_option == "individual"
                        constraint_RampUpMax_idhty(JuMP_model, am, "constraint_RampUpMax_idhty", decomp_group, i, d, h, t, y; const_name_flag)
                        constraint_RampDnMax_idhty(JuMP_model, am, "constraint_RampDnMax_idhty", decomp_group, i, d, h, t, y; const_name_flag)
                        reserve_enabled(am, :reg_up) && up_reserve_eligible(am, i) && constraint_RegUp_Max_idhty(JuMP_model, am, "constraint_RegUp_Max_idhty", decomp_group, i, d, h, t, y; const_name_flag)
                        reserve_enabled(am, :reg_dn) && constraint_RegDn_Max_idhty(JuMP_model, am, "constraint_RegDn_Max_idhty", decomp_group, i, d, h, t, y; const_name_flag)
                    end
                end 

                # Must run
                constraint_MustRun_idhty(JuMP_model, am, "constraint_MustRun_idhty", decomp_group, i, d, h, t, y; const_name_flag)
            end
        end  
    end


    # Inter-temporal Ramping Constraints
    if am.setting["Simulation Configuration"]["include_dispatch_ramping_flag"] == true
        for y in ids_y
            for d in ids_d, h in ids_h, t in ids_t
                
                if am.setting["Simulation Configuration"]["Dispatch_Mode_in_EXP"] == "Unit Commitment"  # no ramping aggregation

                    for i in ids_i
                        unit_category = parameter(am, 0, :gen_index, "UNIT_CATEGORY", i)
                        if unit_category in ["THERMAL", "NUCLEAR", "OTHER"]
                            if i in ids_i_commit
                                constraint_RampUp_InterTemporal_Hour_UC_idhy(JuMP_model, am, "constraint_RampUp_InterTemporal_Hour_UC_idhy", decomp_group, i, d, h, t, y; const_name_flag)
                                constraint_RampDn_InterTemporal_Hour_UC_idhy(JuMP_model, am, "constraint_RampDn_InterTemporal_Hour_UC_idhy", decomp_group, i, d, h, t, y; const_name_flag)
                            end
                        end 
                    end
                
                elseif am.setting["Simulation Configuration"]["dispatch_ramping_modeling_option"] == "individual"
                    for i in ids_i
                        unit_category = parameter(am, 0, :gen_index, "UNIT_CATEGORY", i)
                        if unit_category in ["THERMAL", "NUCLEAR", "OTHER"]
                            constraint_RampUp_InterTemporal_Hour_idhy(JuMP_model, am, "constraint_RampUp_InterTemporal_Hour_idhy", decomp_group, i, d, h, t, y; const_name_flag)
                            constraint_RampDn_InterTemporal_Hour_idhy(JuMP_model, am, "constraint_RampDn_InterTemporal_Hour_idhy", decomp_group, i, d, h, t, y; const_name_flag)
                        end 
                    end
                else
                    for r in ids_r
                        constraint_Agg_Ramp_InterTemporal_Hour_idhy(JuMP_model, am, "constraint_Agg_Ramp_InterTemporal_Hour_idhy", decomp_group, r, d, h, t, y; const_name_flag)    # both up and down
                    end
                end
            end  
        end
    end

    # VRE Balance
    for y in ids_y
        for i in ids_i
            
            if !(i in ids_i_hybrid_gen) # Exclude hybrid main gen units
                bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
                tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)

                if parameter(am, bus_idx, :gen_bus, tech_idx, "VRE_Flag") == true
                    unit_group = parameter(am, 0, :gen_index, "UNIT_GROUP", i)
                    fuel_limit = parameter(am, bus_idx, :gen_bus, tech_idx, "FUEL_LIMIT")

                    if (fuel_limit == "Fixed Profile") 
                        for d in ids_d, h in ids_h, t in ids_t
                            constraint_VREBalance_Fixed_Profile_idhty(JuMP_model, am, "constraint_VREBalance_Fixed_Profile_idhty", decomp_group, i, d, h, t, y, unit_group; const_name_flag)
                        end


                    elseif fuel_limit != "NA"

                        if am.setting["Simulation Configuration"]["Hydro_Flexibility_Flag"] == true
                            for d in ids_d, h in ids_h, t in ids_t
                                constraint_VREBalance_Fixed_Profile_wifh_Flexibilty_idhty(JuMP_model, am, "constraint_VREBalance_Fixed_Profile_wifh_Flexibilty_idhty", decomp_group, i, d, h, t, y, unit_group; const_name_flag)
                            end
                        end

                        if am.setting["Simulation Configuration"]["Hydro_Budget_Flag"] == true

                            if am.setting["Simulation Configuration"]["Reservoir_Hydro_Operation_Option"] == "Budget (day groups)"
                                for day_group_id in ids_day_groups
                                    constraint_VREBalance_Budget_iy(JuMP_model, am, "constraint_VREBalance_Budget_iy", decomp_group, i, day_group_id, y, unit_group; const_name_flag)
                                end
                            
                            elseif am.setting["Simulation Configuration"]["Reservoir_Hydro_Operation_Option"] == "Budget (annual)"
                                constraint_VREBalance_Budget_Annual_iy(JuMP_model, am, "constraint_VREBalance_Budget_Annual_iy", decomp_group, i, y, ids_d, unit_group; const_name_flag)
                            end

                        end

                        # When neither the flexibility nor the budget energy-limit mechanism is active, treat impoundment as ROR.
                        if (am.setting["Simulation Configuration"]["Hydro_Flexibility_Flag"] == false) && (am.setting["Simulation Configuration"]["Hydro_Budget_Flag"] == false)
                            for d in ids_d, h in ids_h, t in ids_t
                                constraint_VREBalance_Fixed_Profile_idhty(JuMP_model, am, "constraint_VREBalance_Fixed_Profile_idhty", decomp_group, i, d, h, t, y, unit_group; const_name_flag)
                            end
                        end
                    end
                end
            end
        end
    end

    # Storage Balance and Power Constraints
    for y in ids_y
        for i in ids_i_sto
            hybrid_type = am.ref[:nw][0][:gen_index][i]["hybrid_type"] 

            if hybrid_type == "ES"      # limit hybrid ES grid charge
                constraint_hybrid_ES_Charge_from_Grid_idhty(JuMP_model, am, "constraint_hybrid_ES_Charge_from_Grid_idhty", decomp_group, i, ids_d, ids_h, ids_t, ids_y; const_name_flag)
            end

            for d in ids_d, h in ids_h, t in ids_t
                
                # charge limit
                if hybrid_type == "ES"
                    constraint_hybrid_ES_Charge_Max_idhty(JuMP_model, am, "constraint_hybrid_ES_Charge_Max_idhty", decomp_group, i, d, h, t, y; const_name_flag)
                else
                    constraint_ES_Charge_Max_idhty(JuMP_model, am, "constraint_ES_Charge_Max_idhty", decomp_group, i, d, h, t, y; const_name_flag)
                end
                
                # SOC bounds
                constraint_ES_SOC_Max_idhty(JuMP_model, am, "constraint_ES_SOC_Max_idhty", decomp_group, i, d, h, t, y; const_name_flag)
                constraint_ES_SOC_Min_idhty(JuMP_model, am, "constraint_ES_SOC_Min_idhty", decomp_group, i, d, h, t, y; const_name_flag)
    
                # Storage commitment constraints
                constraint_ES_DisCharge_Max_Sto_UC_idhty(JuMP_model, am, "constraint_ES_DisCharge_Max_Sto_UC_idhty", decomp_group, i, d, h, t, y; const_name_flag)
                if hybrid_type == "ES"
                    constraint_hybrid_ES_Charge_Max_Sto_UC_idhty(JuMP_model, am, "constraint_hybrid_ES_Charge_Max_Sto_UC_idhty", decomp_group, i, d, h, t, y; const_name_flag)
                else
                    constraint_ES_Charge_Max_Sto_UC_idhty(JuMP_model, am, "constraint_ES_Charge_Max_Sto_UC_idhty", decomp_group, i, d, h, t, y; const_name_flag)
                end
    
            end
        end
    end

    # Storage SOC Balance
    for y in ids_y
        for day_group in ids_day_groups
            for d in am.ref[:nw][0][:repday_groups][day_group]["Day_Idx_List"]
                
                for i in ids_i_sto
                    hybrid_type = am.ref[:nw][0][:gen_index][i]["hybrid_type"] 
                    if hybrid_type == "ES"
                        for h in ids_h, t in ids_t
                            constraint_hybrid_ES_SOC_Balance_Inter_Hour_idhty(JuMP_model, am, "constraint_hybrid_ES_SOC_Balance_Inter_Hour_idhty", decomp_group, i, d, h, t, y, day_group; const_name_flag)
                        end
                    else
                        for h in ids_h, t in ids_t
                            constraint_ES_SOC_Balance_Inter_Hour_idhty(JuMP_model, am, "constraint_ES_SOC_Balance_Inter_Hour_idhty", decomp_group, i, d, h, t, y, day_group; const_name_flag)
                        end
                    end
                end
            end
        end
    end

    # Storage SOC Neutral
    for y in ids_y
        for day_group in ids_day_groups
            for d in am.ref[:nw][0][:repday_groups][day_group]["Day_Idx_List"]
                for i in ids_i_sto
                    hybrid_type = am.ref[:nw][0][:gen_index][i]["hybrid_type"] 
                    end_day = am.ref[:nw][0][:repday_groups][day_group]["End_Day_Id"]
                    start_day_idx = am.ref[:nw][0][:repday_groups][day_group]["Day_Idx_List"][1]
                    if am.ref[:nw][0][:repdays][d]["Day"] == end_day
                        if hybrid_type == "ES"
                            constraint_hybrid_ES_SOC_Neutral_idy(JuMP_model, am, "constraint_hybrid_ES_SOC_Neutral_idy", decomp_group, i, d, y, start_day_idx; const_name_flag)
                        else
                            constraint_ES_SOC_Neutral_idy(JuMP_model, am, "constraint_ES_SOC_Neutral_idy", decomp_group, i, d, y, start_day_idx; const_name_flag)
                        end
                    end
                end
            end
        end
    end


    # Ramping: Inter-hour Ramping Constraints (place holder for now)
    if am.setting["Simulation Configuration"]["FIVEMIN"] == 1      # Intra Hour (only applies with sub_hourly time resolution)
        for y in ids_y
            for day_group in ids_day_groups
                for d in am.ref[:nw][0][:repday_groups][day_group]["Day_Idx_List"]
                    for i in ids_i, h in ids_h
                        unit_category = parameter(am, 0, :gen_index, "UNIT_CATEGORY", i)
                        if unit_category in ["THERMAL", "NUCLEAR", "OTHER"]
                            for t in ids_t
                                constraint_RampUp_Inter_SubHour_idhy(JuMP_model, am, "constraint_RampUp_Inter_SubHour_idhy", decomp_group, i, d, h, t, y, day_group; const_name_flag) 
                                constraint_RampDn_Inter_SubHour_idhy(JuMP_model, am, "constraint_RampDn_Inter_SubHour_idhy", decomp_group, i, d, h, t, y, day_group; const_name_flag) 
                            end
                        end
                    end
                end
            end
        end
    end

    # Large Flexible Load
    for y in ids_y
        for d in ids_d, h in ids_h, t in ids_t
            for lfl in ids_lfl

                constraint_LFL_power_balance_ldhty(JuMP_model, am, "constraint_LFL_power_balance_ldhty", decomp_group, lfl, d, h, t, y; const_name_flag)
                constraint_LFL_inter_connection_limit_injection_ldhty(JuMP_model, am, "constraint_LFL_inter_connection_limit_injection_ldhty", decomp_group, lfl, d, h, t, y; const_name_flag)
                constraint_LFL_inter_connection_limit_withdraw_ldhty(JuMP_model, am, "constraint_LFL_inter_connection_limit_withdraw_ldhty", decomp_group, lfl, d, h, t, y; const_name_flag)
                
                constraint_LFL_segment_bound_lsdhty(JuMP_model, am, "constraint_LFL_segment_bound_lsdhty", decomp_group, lfl, d, h, t, y; const_name_flag)
                constraint_LFL_segment_relation_lsdhty(JuMP_model, am, "constraint_LFL_segment_relation_lsdhty", decomp_group, lfl, d, h, t, y; const_name_flag)

                constraint_lfl_DR_balance_ldhty(JuMP_model, am, "constraint_lfl_DR_balance_ldhty", decomp_group, lfl, d, h, t, y; const_name_flag)

                if am.ref[:nw][0][:demand][lfl]["Hybrid_Gen"] != "NA"
                    constraint_lfl_onsite_gen_thermal_cap_ldhty(JuMP_model, am, "constraint_lfl_onsite_gen_thermal_cap_ldhty", decomp_group, lfl, d, h, t, y; const_name_flag)
                end

                if am.ref[:nw][0][:demand][lfl]["Hybrid_ES"] != "NA"
                    constraint_lfl_onsite_ES_charge_cap_ldhty(JuMP_model, am, "constraint_lfl_onsite_ES_charge_cap_ldhty", decomp_group, lfl, d, h, t, y; const_name_flag)
                    constraint_lfl_onsite_ES_discharge_cap_ldhty(JuMP_model, am, "constraint_lfl_onsite_ES_discharge_cap_ldhty", decomp_group, lfl, d, h, t, y; const_name_flag)
                    constraint_lfl_onsite_ES_SOC_cap_ldhty(JuMP_model, am, "constraint_lfl_onsite_ES_SOC_cap_ldhty", decomp_group, lfl, d, h, t, y; const_name_flag)
                end
            end
        end

        for d in ids_d
            for lfl in ids_lfl
                constraint_lfl_DR_daily_limit_ldhty(JuMP_model, am, "constraint_lfl_DR_daily_limit_ldhty", decomp_group, lfl, d, y; const_name_flag)

            end
        end

        # Large Flexible Load Onsite Storage SOC Balance
        for lfl in ids_lfl
            if am.ref[:nw][0][:demand][lfl]["Hybrid_ES"] != "NA"
                for day_group in ids_day_groups
                    for d in am.ref[:nw][0][:repday_groups][day_group]["Day_Idx_List"]
                
                        for h in ids_h, t in ids_t
                            constraint_lfl_onsite_ES_SOC_balance_ldhty(JuMP_model, am, "constraint_lfl_onsite_ES_SOC_balance_ldhty", decomp_group, lfl, d, h, t, y, day_group; const_name_flag)
                        end

                        if am.ref[:nw][0][:repdays][d]["Day"] == am.ref[:nw][0][:repday_groups][day_group]["End_Day_Id"]
                            start_day_idx = am.ref[:nw][0][:repday_groups][day_group]["Day_Idx_List"][1]
                            constraint_lfl_onsite_ES_SOC_neutral_ldhty(JuMP_model, am, "constraint_lfl_onsite_ES_SOC_neutral_ldhty", decomp_group, lfl, d, y, start_day_idx; const_name_flag)
                        end 
                    end
                end
            end
        end
    end

    # Hybrid Plant Generation Constraint
    for y in ids_y
        for i in ids_i_hybrid_gen
            for d in ids_d, h in ids_h, t in ids_t
                constraint_hybrid_onsite_gen_thermal_cap_ldhty(JuMP_model, am, "constraint_hybrid_onsite_gen_thermal_cap_ldhty", decomp_group, i, d, h, t, y; const_name_flag)
                constraint_hybrid_inter_connection_limit_injection_ldhty(JuMP_model, am, "constraint_hybrid_inter_connection_limit_injection_ldhty", decomp_group, i, d, h, t, y; const_name_flag)
                constraint_hybrid_inter_connection_limit_withdraw_ldhty(JuMP_model, am, "constraint_hybrid_inter_connection_limit_withdraw_ldhty", decomp_group, i, d, h, t, y; const_name_flag)
            end
        end
    end

    # Reliability Constraints

    # 1) Planning reserve margin
    if am.setting["Simulation Configuration"]["enforce_min_reserve_margin_flag"] == true
        for y in ids_y
            for n in ids_pr
                constraint_PRM_y_p_real(JuMP_model, am, "constraint_PRM_y_p_real", decomp_group, y, n)                # Planning Reserve Margin
            end
        end
    end

    # 2) Total ENS MWh cap
    if am.setting["Simulation Configuration"]["Total_ENS_MWh_Cap_Flag"] == true
        for y in ids_y 
            constraint_total_ENS_MWh_cap_y(JuMP_model, am, "constraint_total_ENS_MWh_cap_y", decomp_group, y, ids_n, ids_d)
        end
    end

    # 3) Maximum ENS cap
    if am.setting["Simulation Configuration"]["Max_ENS_MWh_Cap_Flag"] == true
        
        variable_max_ENS_MWh_y_real(JuMP_model, am, :slack, "max_ENS_MWh_y", decomp_group, ids_y)       
                
        for y in ids_y 
            constraint_max_ENS_MWh_y(JuMP_model, am, "constraint_max_ENS_MWh_y", decomp_group, y, ids_n, ids_d) # define max_ENS_MWh_y
            constraint_max_ENS_MWh_cap_y(JuMP_model, am, "constraint_max_ENS_MWh_cap_y", decomp_group, y) # define max_ENS_MWh_y
        end
    end

    # 4) Maximum number of ENS hours cap
    if am.setting["Simulation Configuration"]["ENS_Hours_Cap_Flag"] == true

        for y in ids_y
            constraint_ENS_hours_approx_cap_y_real(JuMP_model, am, "constraint_ENS_hours_approx_cap_y_real", decomp_group, y, ids_n, ids_d)
        end
    end


    # Policy Constraints

    # Storage annual energy throughput limits
    if am.setting["Simulation Configuration"]["Energy_Storage_AET_Limit_Flag"] == true
        for y in ids_y 
            for i in ids_i_sto
                constraint_ES_AET_y(JuMP_model, am, "constraint_ES_AET_y", decomp_group, i, y, ids_d; const_name_flag)
            end
        end
    end

    # RPS Constraint
    if am.setting["Simulation Configuration"]["RPS_Flag"] == true

        variable_slack_RPS_ny_real(JuMP_model, am, :slack, "slack_RPS_ny", decomp_group, ids_p, ids_y)   

        for y in ids_y 
            for n in ids_p
                constraint_RPS_regional_ny(JuMP_model, am, "constraint_RPS_regional_ny", decomp_group, n, y)
            end
        end
    end

    # Clean Energy Generation Constraint
    if am.setting["Simulation Configuration"]["Clean_Energy_Generation_Target_EXP_Flag"] == true

        if am.setting["Simulation Configuration"]["Clean_Energy_Generation_Target_EXP_Type"] == "Daygroup"

            variable_slack_CEG_ndy_real(JuMP_model, am, :slack, "slack_CEG_ndy", decomp_group, ids_p, ids_d, ids_y)       

            for y in ids_y 
                if am.ref[:nw][0][:planning_stages][y]["year"] >= am.setting["Simulation Configuration"]["Clean_Energy_Generation_Target_Start_Year"]
                    for n in ids_p, d in ids_d
                        constraint_clean_energy_generation_ndy(JuMP_model, am, "constraint_clean_energy_generation_ndy", decomp_group, n, d, y)
                    end
                end
            end
        
        elseif am.setting["Simulation Configuration"]["Clean_Energy_Generation_Target_EXP_Type"] == "Annual"

            variable_slack_CEG_ny_real(JuMP_model, am, :slack, "slack_CEG_ny", decomp_group, ids_p, ids_y)       

            for y in ids_y 
                if am.ref[:nw][0][:planning_stages][y]["year"] >= am.setting["Simulation Configuration"]["Clean_Energy_Generation_Target_Start_Year"]
                    for n in ids_p
                        constraint_clean_energy_generation_ny(JuMP_model, am, "constraint_clean_energy_generation_ny", decomp_group, n, y, ids_d)
                    end
                end
            end
        end

    end

    # Carbon emission reduction target constraint
    if am.setting["Simulation Configuration"]["Carbon_Emission_Reduction_Target_Flag"] == true
        for y in ids_y 
            if am.ref[:nw][0][:planning_stages][y]["year"] >= am.setting["Simulation Configuration"]["Carbon_Emission_Reduction_Target_Start_Year"]
                for n in ids_p
                    constraint_carbon_emission_reduction_target_ny(JuMP_model, am, "constraint_carbon_emission_reduction_target_ny", decomp_group, n, y, ids_d)
                end
            end
        end
    end

    # Material Constraints
    if am.setting["Planning Design"]["enforce_material_constraints_flag"] == true
        for y in ids_y
            for m in ids_m
                constraint_annual_raw_mmaterial_limits_my_real(JuMP_model, am, "constraint_annual_raw_mmaterial_limits_my_real", decomp_group, m, y; const_name_flag)
            end
        end
    end

    # System inertia
    if (am.setting["Planning Design"]["enforce_rotational_inertia_constraints_flag"] == true) && (am.setting["Simulation Configuration"]["Dispatch_Mode_in_EXP"] == "Unit Commitment")
        for n in ids_p
            for y in ids_y
                for d in ids_d, h in ids_h, t in ids_t
                    constraint_inertia_ndhty(JuMP_model, am, "constraint_inertia_ndhty", decomp_group, n, d, h, t, y; const_name_flag)
                end
            end
        end
    end
       
    ###################################
    #------ DEFINE Objective
    ###################################
    objective_function_expansion(JuMP_model, am, ids_i, ids_d, ids_h, ids_t, ids_y, ids_k, ids_z, ids_p, ids_n, ids_i_commit, ids_lfl, decomp_group)

    am.model[:nw][nw][decomp_group] = JuMP_model

end


function update_CAP_i_noEXP!(am::Abstract_ALEAF_Model; nw::Int=am.cnw)

    for i in [(i) for (i) in get_index(am, :gen_index, 0)]
        bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
        tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)    
        EXUNITS = parameter(am, bus_idx, :gen_bus, tech_idx, "EXUNITS")
        CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP") 
        EXCAP = EXUNITS * CAP

        if EXCAP > 0
            if (EXCAP <= CAP) 
                am.ref[:nw][bus_idx][:gen_bus][tech_idx]["CAP"] = EXCAP 
                am.ref[:nw][bus_idx][:gen_bus][tech_idx]["EXUNITS"] = 1.0
                am.ref[:nw][bus_idx][:gen_bus][tech_idx]["U_G_iy"] = 1.0

                if am.ref[:nw][bus_idx][:gen_bus][tech_idx]["Charge_CAP"] > 0
                    am.ref[:nw][bus_idx][:gen_bus][tech_idx]["Charge_CAP"] = EXCAP 
                end

            else
                if (isinteger(EXUNITS) == false) #&& (EXUNITS > 1)
                    reminder = mod(EXCAP, CAP)
                    num_units = round(EXUNITS, RoundDown)
                    capacity_to_add = reminder/num_units
                    final_capacity = (CAP + capacity_to_add) 
                    am.ref[:nw][bus_idx][:gen_bus][tech_idx]["CAP"] = final_capacity
                    am.ref[:nw][bus_idx][:gen_bus][tech_idx]["EXUNITS"] = num_units
                    am.ref[:nw][bus_idx][:gen_bus][tech_idx]["U_G_iy"] = num_units

                    if am.ref[:nw][bus_idx][:gen_bus][tech_idx]["Charge_CAP"] > 0
                        am.ref[:nw][bus_idx][:gen_bus][tech_idx]["Charge_CAP"] = final_capacity
                    end

                else 

                    am.ref[:nw][bus_idx][:gen_bus][tech_idx]["CAP"] = CAP
                    am.ref[:nw][bus_idx][:gen_bus][tech_idx]["EXUNITS"] = EXUNITS
                    am.ref[:nw][bus_idx][:gen_bus][tech_idx]["U_G_iy"] = EXUNITS

                    if am.ref[:nw][bus_idx][:gen_bus][tech_idx]["Charge_CAP"] > 0
                        am.ref[:nw][bus_idx][:gen_bus][tech_idx]["Charge_CAP"] = CAP
                    end
                end
            end
        else
            am.ref[:nw][bus_idx][:gen_bus][tech_idx]["EXUNITS"] = 0.0
            am.ref[:nw][bus_idx][:gen_bus][tech_idx]["U_G_iy"] = 0.0
        end
    end

end


function external_constraint_generator(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, external_const_id::String, modeled_year_list; nw::Int=am.cnw, const_name_flag::Bool=false)

    pu_base = am.setting["Simulation Setting"]["per_unit_base_value"]
    
    # constraint information
    const_unit_group = am.setting["External Constraints"][external_const_id]["Resource Unit Group"]
    const_parameter = am.setting["External Constraints"][external_const_id]["Parameter"]
    const_region_id = am.setting["External Constraints"][external_const_id]["Region ID"]
    const_operator = am.setting["External Constraints"][external_const_id]["Operator"]
    const_value = am.setting["External Constraints"][external_const_id]["Value"]
    const_value_unit = am.setting["External Constraints"][external_const_id]["Value_Unit"]
    const_compared_value = am.setting["External Constraints"][external_const_id]["Compared To"]
    const_years = am.setting["External Constraints"][external_const_id]["Planning Year"]
    const_description = am.setting["External Constraints"][external_const_id]["Description"]

    # index
    ids_i = [(i) for (i) in get_index(am, :gen_index, 0)]
    
    if const_years in modeled_year_list
    
        ids_y = [const_years]
        
        # Retirement
        if const_parameter == "Econ Retirement" # Retirement

            for y in ids_y 
                for i in ids_i
                    
                    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
                    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)    
                    CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP") 

                    if parameter(am, 0, :bus, "region_id", bus_idx) == const_region_id
                        if parameter(am, bus_idx, :gen_bus, tech_idx, "UNITGROUP") == const_unit_group
                            constraint_ret_G_external_iy_real(JuMP_model, am, "constraint_ret_G_external_iy_real", decomp_group, i, y, const_value/pu_base)
                        end
                    end

                end
            end
        end

        # Investment
        if const_parameter == "Investment" # Investment

            for y in ids_y 
                for i in ids_i
                    
                    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
                    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)    
                    CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP") 

                    if parameter(am, 0, :bus, "region_id", bus_idx) == const_region_id
                        if parameter(am, bus_idx, :gen_bus, tech_idx, "UNITGROUP") == const_unit_group
                            constraint_investment_external_iy_real(JuMP_model, am, "constraint_investment_external_iy_real", decomp_group, i, y, const_value/pu_base)
                        end
                    end

                end
            end

        end


        for y in ids_y  # in the selected year
                    
            sum_of_variables = JuMP.AffExpr(0.0)   
            if const_parameter == "Total ICAP" # Existing + New Capacity
                for i in ids_i
                    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
                    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)    
                    CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP") 

                    if parameter(am, 0, :bus, "region_id", bus_idx) == const_region_id
                        if parameter(am, bus_idx, :gen_bus, tech_idx, "UNITGROUP") == const_unit_group
                            JuMP.add_to_expression!(sum_of_variables, CAP, variable(am, nw, decomp_group, :u_new_G_iy, (i, y)))      
                            JuMP.add_to_expression!(sum_of_variables, -CAP, variable(am, nw, decomp_group, :u_ret_G_iy, (i, y))) 
                        end
                    end
                end       
            end

            # sum_of_variables now includes the investment and retirement of const_unit_group resources in const_region_id in year y
            RHS = 0.0
            if const_value_unit == "%"
                if const_compared_value == "EXCAP"
                    for i in ids_i
                        bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
                        tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)
                        CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")
                        EXUNITS = parameter(am, bus_idx, :gen_bus, tech_idx, "EXUNITS")

                        if parameter(am, 0, :bus, "region_id", bus_idx) == const_region_id
                            if parameter(am, bus_idx, :gen_bus, tech_idx, "UNITGROUP") == const_unit_group
                                RHS += parameter(am, bus_idx, :gen_bus, tech_idx, "EXCAPS") * const_value * 0.01
                                RHS += - CAP * EXUNITS
                            end
                        end
                    end
                end
            elseif const_value_unit == "MW"
                # absolute capacity bound in MW; divide by pu_base to match per-unit sum_of_variables
                RHS = const_value / pu_base
            end

            # constraint expression
            expr = JuMP.@expression(JuMP_model,  
                sum_of_variables - RHS
                )

            if const_operator == "less than"
                constraint = JuMP.@constraint(JuMP_model, 
                    expr <= 0
                    )  
            elseif const_operator == "greater than"
                constraint = JuMP.@constraint(JuMP_model, 
                    expr >= 0
                    )  
            elseif const_operator == "equal to"
                constraint = JuMP.@constraint(JuMP_model, 
                    expr == 0
                    )  
            end

            # Set constraint name and add to the aleaf model instance

            if const_name_flag JuMP.set_name(constraint, string(const_name, "_($external_const_id)")) end    
            
        end
    end
    
end


###################################
#------ Objectives
###################################

function objective_function_expansion(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, ids_i, ids_d, ids_h, ids_t, ids_y, ids_k, ids_z, ids_p, ids_n, ids_i_commit, ids_lfl, decomp_group; nw::Int=am.cnw, report::Bool=true)

    card_T = length(am.setting["run_T"])
    inv_cart_T = 1.0 / card_T

    # Anti-degeneracy tie-breaker penalty (per-unit). Negligible vs. the smallest real
    # per-unit marginal cost (~3e-3), so dispatch is unchanged; lifts the objective floor.
    tie_breaker_penalty = 1e-6

    # define the reference year for discount
    current_year = am.setting["Planning Design"]["dollar_year_value"] # Discount factor relative to the current dollor year
    discount_rate = am.setting["Planning Design"]["discount_rate_value"]

    # Minimum regulation-reserve cost floor ($/MWh) so regulation is never free, removing the
    # degeneracy where reg-down is procured arbitrarily when its cost (0.3*MC) is ~0; default 0.0.
    pu_econ_base = am.setting["Simulation Setting"]["per_unit_econ_base_value"] / am.setting["Simulation Setting"]["per_unit_base_value"]
    # Divide the raw $/MWh value by pu_econ_base to match Annual_MC's per-unit basis,
    # so the floor compares like-for-like against reg_cost in the objective.
    min_reg_cost = get(am.setting["Planning Design"], "min_regulation_cost_value", 0.0) / pu_econ_base

    # Define the objective function
    objective = JuMP.AffExpr(0.0)

    # Scale cost functions
    total_num_days = 0.0
    for d in ids_d
        total_num_days += parameter(am, 0, :repdays, "NumDays", d)
    end
    num_day_scale = total_num_days / 365    # this will be 1 if PH is not used.
    num_day_scale_OP =  1.0    # this will be 1 if PH is not used.

    # Investment costs for generators
    for y in ids_y

        stage_length = am.ref[:nw][0][:planning_stages][y]["stage_length"]
        if am.setting["Planning Design"]["multi_round_solution_process_flag"] == true
            num_stage = (am.setting["Planning Design"]["num_decision_stages_per_round_value"] + am.setting["Planning Design"]["num_lookahead_stages_per_round_value"])
            remaining_years = (num_stage + ids_y[1] - y) * stage_length
        else
            num_stage = am.setting["Planning Design"]["num_stages_value"]
            remaining_years = (num_stage + 1 - y) * stage_length
        end

        for i in ids_i
            bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
            tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)
            capacity = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")
            CRP_i = parameter(am, bus_idx, :gen_bus, tech_idx, "crpyears")

            # Determine the cost duration (investment payments)
            payment_duration = min(CRP_i, remaining_years)

            # Get the investment costs
            investment_cost = parameter(am, bus_idx, :gen_bus, tech_idx, "Annual_INVC")[string(y)]
            storage_investment_cost = parameter(am, bus_idx, :gen_bus, tech_idx, "Annual_STO_INV")[string(y)]
            
            for future_y in 1:payment_duration

                future_year = am.ref[:nw][0][:planning_stages][y]["year"] + future_y - 1    # Actual future year 
                discount_factor = (1 + discount_rate)^-(future_year - current_year) # Discount factor relative to the current dollor year

                # Expansion Cost (new units only) 
                JuMP.add_to_expression!(objective, num_day_scale * 1000 * discount_factor * investment_cost * capacity, variable(am, nw, decomp_group, :u_new_G_iy, (i, y)))  

                # Storage Duration Cost ($/kWh)
                if parameter(am, 0, :gen_index, "UNIT_CATEGORY", i) == "STORAGE"
                    JuMP.add_to_expression!(objective, num_day_scale * 1000 * discount_factor * storage_investment_cost * capacity, variable(am, nw, decomp_group, :u_new_ESH_iy, (i, y)))  
                end
                
            end
        end
    end

    # Investment costs for transmission lines
    for y in ids_y

        stage_length = am.ref[:nw][0][:planning_stages][y]["stage_length"]
        if am.setting["Planning Design"]["multi_round_solution_process_flag"] == true
            num_stage = (am.setting["Planning Design"]["num_decision_stages_per_round_value"] + am.setting["Planning Design"]["num_lookahead_stages_per_round_value"])
            remaining_years = (num_stage + ids_y[1] - y) * stage_length
        else
            num_stage = am.setting["Planning Design"]["num_stages_value"]
            remaining_years = (num_stage + 1 - y) * stage_length
        end

        trans_crpyears = am.setting["Planning Design"]["transmission_investment_CRP_value"]
        payment_duration = min(trans_crpyears, remaining_years)

        for k in ids_k
            
            if parameter(am, 0, :branch, "expansion_flag", k) == true
                branch_length = get(am.ref[:nw][0][:branch][k], "length", get(am.ref[:nw][0][:branch][k], "Length", 0.0))
                # DC-tie expansion cost is per-MW (no mile factor); AC lines are per-MW-mile
                len_factor = get(am.ref[:nw][0][:branch][k], "dc_line", false) == true ? 1.0 : branch_length

                for future_y in 1:payment_duration

                    future_year = am.ref[:nw][0][:planning_stages][y]["year"] + future_y - 1
                    discount_factor = (1 + discount_rate)^-(future_year - current_year)

                    # Transmission expansion cost ($/mile)
                    JuMP.add_to_expression!(objective, num_day_scale * discount_factor * parameter(am, 0, :branch, "transmission_expansion_cost", k) * parameter(am, 0, :branch, "rate_a", k) * len_factor, variable(am, nw, decomp_group, :u_new_T_ky, (k, y)))

                end
            end
        end
        
    end

    # Retirement costs for generators
    for y in ids_y

        future_year = am.ref[:nw][0][:planning_stages][y]["year"]
        discount_factor = (1 + discount_rate)^-(future_year - current_year)

        for i in ids_i
            bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
            tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)

            # Decommissioning Cost: k$/MW
            DECC = parameter(am, bus_idx, :gen_bus, tech_idx, "DECC") 
            if DECC > 0
                JuMP.add_to_expression!(objective, num_day_scale * 1000 * discount_factor * DECC * parameter(am, bus_idx, :gen_bus, tech_idx, "CAP"), variable(am, nw, decomp_group, :u_ret_G_iy, (i, y)))  
            end
        end

    end


    # Investment Credit: already considered in the annual investment costs

    # FOM
    for y in ids_y
        
        stage_length = am.ref[:nw][0][:planning_stages][y]["stage_length"]

        for i in ids_i
            bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
            tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)

            capacity = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")
            FOM = parameter(am, bus_idx, :gen_bus, tech_idx, "Annual_FOM")[string(y)]

            if (FOM > 0) && (capacity > 0)
        
                for future_y in 1:stage_length
                
                    future_year = am.ref[:nw][0][:planning_stages][y]["year"] + future_y - 1 
                    discount_factor = (1 + discount_rate)^-(future_year - current_year)

                    # Fixed OM Cost ($/kW-Year)
                    JuMP.add_to_expression!(objective, num_day_scale * 1000 * discount_factor * parameter(am, bus_idx, :gen_bus, tech_idx, "Annual_FOM")[string(y)] * parameter(am, bus_idx, :gen_bus, tech_idx, "CAP"), variable(am, nw, decomp_group, :u_G_iy, (i, y)))  
                    
                end
        
            end
        end
    end

    # Fixed O&M for transmission (expansion only). Existing-grid FOM is a fixed constant -- transmission
    # cannot retire, so no decision depends on it -- and is added in reporting, not the objective.
    # Charged on cumulative in-service expansion (u_T_ky), every operating year of the stage, like gen FOM.
    for y in ids_y
        stage_length = am.ref[:nw][0][:planning_stages][y]["stage_length"]
        for k in ids_k
            if parameter(am, 0, :branch, "expansion_flag", k) == true
                branch_length = get(am.ref[:nw][0][:branch][k], "length", get(am.ref[:nw][0][:branch][k], "Length", 0.0))
                len_factor = get(am.ref[:nw][0][:branch][k], "dc_line", false) == true ? 1.0 : branch_length
                fom_cost = parameter(am, 0, :branch, "transmission_fom_cost", k)
                rate_a = parameter(am, 0, :branch, "rate_a", k)

                for future_y in 1:stage_length
                    future_year = am.ref[:nw][0][:planning_stages][y]["year"] + future_y - 1
                    discount_factor = (1 + discount_rate)^-(future_year - current_year)

                    JuMP.add_to_expression!(objective, num_day_scale * discount_factor * fom_cost * rate_a * len_factor, variable(am, nw, decomp_group, :u_T_ky, (k, y)))
                end
            end
        end
    end


    # Hourly Operating Costs
    ctax = am.setting["Simulation Configuration"]["CTAX"]

    for y in ids_y

        stage_length = am.ref[:nw][0][:planning_stages][y]["stage_length"]

        # LFL demand response cost
        for lfl in ids_lfl

            local_lfl = am.ref[:nw][0][:demand][lfl]
            num_seg = local_lfl["Num_DR_Segments"]

            for future_y in 1:stage_length

                future_year = am.ref[:nw][0][:planning_stages][y]["year"] + future_y - 1 
                discount_factor = (1 + discount_rate)^-(future_year - current_year)

                # skip s==1 by construction
                for s in 2:num_seg
                    dr_cost = local_lfl["Price_$(s)"]  # computed once per s

                    for d in ids_d 
                        num_days = parameter(am, 0, :repdays, "NumDays", d)
                        coef = -(num_day_scale_OP * num_days * discount_factor * dr_cost * inv_cart_T)
                        for h in ids_h, t in ids_t
                            JuMP.add_to_expression!(objective, coef, variable(am, nw, decomp_group, :lfl_seg_lsdhty, (lfl,s,d,h,t,y)))  
                        end
                    end
                end

                # Ancillary penalties for hybrid ES grid charging
                if am.ref[:nw][0][:demand][lfl]["Hybrid_ES"] != "NA"
                    for d in ids_d, h in ids_h, t in ids_t
                        JuMP.add_to_expression!(objective, tie_breaker_penalty, variable(am, nw, decomp_group, :lfl_chg_Grid_ES_ldhty, (lfl,d,h,t,y)))  
                    end
                end
            end
        end
        
        for i in ids_i

            bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
            tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)
            unit_group = parameter(am, 0, :gen_index, "UNIT_GROUP", i)
            
            emission_rate = parameter(am, bus_idx, :gen_bus, tech_idx, "Emission_CO2")

            reg_cost = parameter(am, bus_idx, :gen_bus, tech_idx, "reg_cost")
            spin_cost = parameter(am, bus_idx, :gen_bus, tech_idx, "spin_cost")
            nspin_cost = parameter(am, bus_idx, :gen_bus, tech_idx, "nspin_cost")
            flex_cost = parameter(am, bus_idx, :gen_bus, tech_idx, "flex_cost")
            start_up_cost = parameter(am, bus_idx, :gen_bus, tech_idx, "SUC")
            no_load_cost = parameter(am, bus_idx, :gen_bus, tech_idx, "NLC")

            for future_y in 1:stage_length

                future_year = am.ref[:nw][0][:planning_stages][y]["year"] + future_y - 1 
                discount_factor = (1 + discount_rate)^-(future_year - current_year)

                for d in ids_d 

                    num_days = parameter(am, 0, :repdays, "NumDays", d)
                    Annual_MC = parameter(am, bus_idx, :gen_bus, tech_idx, "Annual_MC")[string(y)][d]
                
                    if am.setting["Planning Design"]["reserve_cost_type_flag"] == "percentage"
                        reg_cost = parameter(am, bus_idx, :gen_bus, tech_idx, "reg_cost") * Annual_MC
                        spin_cost = parameter(am, bus_idx, :gen_bus, tech_idx, "spin_cost") * Annual_MC
                        nspin_cost = parameter(am, bus_idx, :gen_bus, tech_idx, "nspin_cost") * Annual_MC
                        flex_cost = parameter(am, bus_idx, :gen_bus, tech_idx, "flex_cost") * Annual_MC
                    end

                    # Apply the minimum regulation cost floor ($/MWh) so regulation is never free.
                    reg_cost = max(reg_cost, min_reg_cost)

                    coef = num_day_scale_OP * num_days * discount_factor * inv_cart_T

                    for h in ids_h, t in ids_t

                        # Generation Cost
                        JuMP.add_to_expression!(objective, coef * Annual_MC, variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y)))

                        if am.setting["Simulation Configuration"]["operating_reserve_modeling_option"] == "individual"
                            up_res_eligible = up_reserve_eligible(am, i)  # up-reserves exist only for dispatchable units
                            # Regulation Cost = x% of generation cost
                            reserve_enabled(am, :reg_up) && up_res_eligible && JuMP.add_to_expression!(objective, coef * reg_cost, variable(am, nw, decomp_group, :reg_up_idhty, (i,d,h,t,y)))
                            reserve_enabled(am, :reg_dn) && JuMP.add_to_expression!(objective, coef * reg_cost, variable(am, nw, decomp_group, :reg_dn_idhty, (i,d,h,t,y)))

                            # spin reserve Cost = 80% of regulation cost
                            reserve_enabled(am, :spin) && up_res_eligible && JuMP.add_to_expression!(objective, coef * spin_cost, variable(am, nw, decomp_group, :spin_idhty, (i,d,h,t,y)))

                            if am.setting["Simulation Configuration"]["Dispatch_Mode_in_EXP"] == "Unit Commitment"
                                reserve_enabled(am, :nonspin) && JuMP.add_to_expression!(objective, coef * nspin_cost, variable(am, nw, decomp_group, :nonspin_idhty, (i,d,h,t,y)))
                            end

                            # Flex reserve Cost = 60% of regulation cost
                            reserve_enabled(am, :flex_up) && up_res_eligible && JuMP.add_to_expression!(objective, coef * flex_cost, variable(am, nw, decomp_group, :flex_up_idhty, (i,d,h,t,y)))
                            reserve_enabled(am, :flex_dn) && JuMP.add_to_expression!(objective, coef * flex_cost, variable(am, nw, decomp_group, :flex_dn_idhty, (i,d,h,t,y)))
                        end

                        # Carbon Tax 
                        if (emission_rate > 0) && (ctax > 0)
                            # Emissions Cost ($/MW) = $/tonne * kg/MMbtu * MMbtu/MWh * MWh * 1 tonne/1000kg = $/tonne * metric tonnes CO2 per timestep
                            JuMP.add_to_expression!(objective, coef * ctax * emission_rate, variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y)))  
                        end
                    end

                    # Start up and no load cost
                    if (i in ids_i_commit) && (am.setting["Simulation Configuration"]["Dispatch_Mode_in_EXP"] == "Unit Commitment")
                        for h in ids_h

                            # Startup Cost ($/MW)
                            JuMP.add_to_expression!(objective, coef * start_up_cost, variable(am, nw, decomp_group, :su_idhty, (i,d,h,1,y)))  

                            # No Load Cost ($/MW)
                            JuMP.add_to_expression!(objective, coef * no_load_cost, variable(am, nw, decomp_group, :c_idhty, (i,d,h,1,y)))  
                        end
                    end
                end
            end
        end

        if am.setting["Simulation Configuration"]["operating_reserve_modeling_option"] == "aggregated"
            ids_r = sort(collect(keys(am.ref[:nw][0][:reserve_group_lookup])))
            for r in ids_r

                group_info = am.ref[:nw][0][:reserve_group_lookup][r]                                         
                unit_group = group_info["unit_group"]                                                         
                z = group_info["zone_idx"]                                                                    
                gen_idx = am.ref[:nw][0][:zone]["reserve"][string(z)]["aggregation_info"]["gen_technology_groups"][unit_group]["gen_idx"]

                i = gen_idx[1]  

                bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
                tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)

                reg_cost = parameter(am, bus_idx, :gen_bus, tech_idx, "reg_cost")
                spin_cost = parameter(am, bus_idx, :gen_bus, tech_idx, "spin_cost")
                nspin_cost = parameter(am, bus_idx, :gen_bus, tech_idx, "nspin_cost")
                flex_cost = parameter(am, bus_idx, :gen_bus, tech_idx, "flex_cost")

                for future_y in 1:stage_length

                    future_year = am.ref[:nw][0][:planning_stages][y]["year"] + future_y - 1 
                    discount_factor = (1 + discount_rate)^-(future_year - current_year)

                    for d in ids_d 

                        num_days = parameter(am, 0, :repdays, "NumDays", d)
                        Annual_MC = parameter(am, bus_idx, :gen_bus, tech_idx, "Annual_MC")[string(y)][d]
                    
                        if am.setting["Planning Design"]["reserve_cost_type_flag"] == "percentage"
                            reg_cost = parameter(am, bus_idx, :gen_bus, tech_idx, "reg_cost") * Annual_MC
                            spin_cost = parameter(am, bus_idx, :gen_bus, tech_idx, "spin_cost") * Annual_MC
                            nspin_cost = parameter(am, bus_idx, :gen_bus, tech_idx, "nspin_cost") * Annual_MC
                            flex_cost = parameter(am, bus_idx, :gen_bus, tech_idx, "flex_cost") * Annual_MC
                        end

                        # Apply the minimum regulation cost floor ($/MWh) so regulation is never free.
                        reg_cost = max(reg_cost, min_reg_cost)

                        coef = num_day_scale_OP * num_days * discount_factor * inv_cart_T

                        for h in ids_h, t in ids_t

                            # Regulation Cost = x% of generation cost
                            reserve_enabled(am, :reg_up) && JuMP.add_to_expression!(objective, coef * reg_cost, variable(am, nw, decomp_group, :reg_up_idhty, (r,d,h,t,y)))
                            reserve_enabled(am, :reg_dn) && JuMP.add_to_expression!(objective, coef * reg_cost, variable(am, nw, decomp_group, :reg_dn_idhty, (r,d,h,t,y)))

                            # spin reserve Cost = 80% of regulation cost
                            reserve_enabled(am, :spin) && JuMP.add_to_expression!(objective, coef * spin_cost, variable(am, nw, decomp_group, :spin_idhty, (r,d,h,t,y)))

                            if am.setting["Simulation Configuration"]["Dispatch_Mode_in_EXP"] == "Unit Commitment"
                                reserve_enabled(am, :nonspin) && JuMP.add_to_expression!(objective, coef * nspin_cost, variable(am, nw, decomp_group, :nonspin_idhty, (r,d,h,t,y)))
                            end

                            # Flex reserve Cost = 60% of regulation cost
                            reserve_enabled(am, :flex_up) && JuMP.add_to_expression!(objective, coef * flex_cost, variable(am, nw, decomp_group, :flex_up_idhty, (r,d,h,t,y)))
                            reserve_enabled(am, :flex_dn) && JuMP.add_to_expression!(objective, coef * flex_cost, variable(am, nw, decomp_group, :flex_dn_idhty, (r,d,h,t,y)))
                        end
                    end
                end
            end
        end
    end

    # Production Tax Credit
    if am.setting["Simulation Configuration"]["PTC_Flag"] == true

        PTC_final_year_for_existing_assets = 5
        PTC_final_year_for_new_assets = 10

        for i in ids_i

            bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
            tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)
            Tech_Type = parameter(am, bus_idx, :gen_bus, tech_idx, "Tech_Type")
            unit_group = parameter(am, 0, :gen_index, "UNIT_GROUP", i)
            profile_type = parameter(am, 0, :gen_index, "Profile_Type", i)
            CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")
            PMAX = parameter(am, bus_idx, :gen_bus, tech_idx, "PMAX")

            if parameter(am, bus_idx, :gen_bus, tech_idx, "PTC Flag") == true
                
                if (Tech_Type == "Existing") && (unit_group != "nuclear")
                    count = 1
                    
                    for y in ids_y
                        
                        stage_length = am.ref[:nw][0][:planning_stages][y]["stage_length"]
        
                        for future_y in 1:stage_length

                            future_year = am.ref[:nw][0][:planning_stages][y]["year"] + future_y - 1 
                            discount_factor = (1 + discount_rate)^-(future_year - current_year)
        
                            if count > PTC_final_year_for_existing_assets
                                break
                            end
        
                            ptc_year = min(2050, future_year)
                            ptc = parameter(am, bus_idx, :gen_bus, tech_idx, "PTC")[string(ptc_year)] * 10 

                            for d in ids_d

                                coef = -1 * num_day_scale_OP * parameter(am, 0, :repdays, "NumDays", d) * discount_factor * inv_cart_T

                                for h in ids_h, t in ids_t
                                    # ptc = cents/kWh -> dollar/MWh
                                    JuMP.add_to_expression!(objective, coef * ptc, variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y)))  
                                end

                            end

                            count += 1
                        end
                    end

                elseif (Tech_Type == "Existing") && (unit_group == "nuclear")

                    for y in ids_y
                        
                        stage_length = am.ref[:nw][0][:planning_stages][y]["stage_length"]
        
                        for future_y in 1:stage_length

                            future_year = am.ref[:nw][0][:planning_stages][y]["year"] + future_y - 1 
                            discount_factor = (1 + discount_rate)^-(future_year - current_year)
        
                            ptc_year = min(2050, future_year)
                            ptc = parameter(am, bus_idx, :gen_bus, tech_idx, "PTC")[string(ptc_year)] * 10 

                            for d in ids_d

                                coef = -1 * num_day_scale_OP * parameter(am, 0, :repdays, "NumDays", d) * discount_factor * inv_cart_T

                                for h in ids_h, t in ids_t
                                    # ptc = cents/kWh -> dollar/MWh
                                    JuMP.add_to_expression!(objective, coef * ptc, variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y)))  

                                end
                            end
                        end
                    end

                else # Tech_Type == "New" 

                    for y in ids_y

                        ptc_year = min(2050, am.ref[:nw][0][:planning_stages][y]["year"])   # ptc is fixed to the investment year ptc
                        ptc = parameter(am, bus_idx, :gen_bus, tech_idx, "PTC")[string(ptc_year)] * 10 

                        stage_length = am.ref[:nw][0][:planning_stages][y]["stage_length"]
                        if am.setting["Planning Design"]["multi_round_solution_process_flag"] == true
                            num_stage = (am.setting["Planning Design"]["num_decision_stages_per_round_value"] + am.setting["Planning Design"]["num_lookahead_stages_per_round_value"])
                            remaining_years = (num_stage + ids_y[1] - y) * stage_length
                        else
                            num_stage = am.setting["Planning Design"]["num_stages_value"]
                            remaining_years = (num_stage + 1 - y) * stage_length
                        end

                        ptc_duration = min(PTC_final_year_for_new_assets, remaining_years)
        
                        for future_y in 1:ptc_duration # for the (next 10 years or remaining years)

                            future_year = am.ref[:nw][0][:planning_stages][y]["year"] + future_y - 1 
                            discount_factor = (1 + discount_rate)^-(future_year - current_year)

                            # Estimated annual generation in year y 
                            vre_data = am.ref[:nw][0][:bus][bus_idx]["vre_aggregated_data"][string(y)]

                            # Estimated generation per MW from the shape class (Profile_Type)
                            annual_generation_per_MW = 0.0
                            if profile_type != "NA"
                                annual_generation_per_MW = vre_data[profile_type * "_shape"]
                            else    # nuclear
                                annual_generation_per_MW = 8760 * PMAX # derated by FOR
                            end

                            # estimated annual generation = annual_generation_per_MW * capacity * (1-FOR)
                            JuMP.add_to_expression!(objective, -1 * discount_factor * ptc * CAP * PMAX * annual_generation_per_MW, variable(am, nw, decomp_group, :u_new_G_iy, (i, y)))  
                        end
                    end
                
                
                end
            end
        end
    end


    # Scarcity
    for y in ids_y

        stage_length = am.ref[:nw][0][:planning_stages][y]["stage_length"]
        VOLL = am.setting["Simulation Configuration"]["VOLL"]
        RegRSP = am.setting["Simulation Configuration"]["RegRSP"]
        SRSP = am.setting["Simulation Configuration"]["SRSP"]
        NSRSP = am.setting["Simulation Configuration"]["NSRSP"]
        FLEXRSP = am.setting["Simulation Configuration"]["FLEXRSP"]
        # Load-Resource reserve priced above generation reserve, far below the RNS penalty.
        DRESCOST = get(am.setting["Simulation Configuration"], "demand_reserve_cost", 5.0)

        for future_y in 1:stage_length

            future_year = am.ref[:nw][0][:planning_stages][y]["year"] + future_y - 1 
            discount_factor = (1 + discount_rate)^-(future_year - current_year)

            for d in ids_d

                coef = num_day_scale_OP * parameter(am, 0, :repdays, "NumDays", d) * discount_factor * inv_cart_T

                for h in ids_h, t in ids_t 

                    for n in ids_n
                        JuMP.add_to_expression!(objective, coef * VOLL, variable(am, nw, decomp_group, :ens_ndhty, (n,d,h,t,y)))  
                    end

                    for z in ids_z
                        reserve_enabled(am, :reg_up)  && JuMP.add_to_expression!(objective, coef * RegRSP, variable(am, nw, decomp_group, :rns_reg_up_zdhty, (z,d,h,t,y)))
                        reserve_enabled(am, :reg_dn)  && JuMP.add_to_expression!(objective, coef * RegRSP, variable(am, nw, decomp_group, :rns_reg_dn_zdhty, (z,d,h,t,y)))
                        reserve_enabled(am, :spin)    && JuMP.add_to_expression!(objective, coef * SRSP, variable(am, nw, decomp_group, :rns_cont_zdhty, (z,d,h,t,y)))
                        reserve_enabled(am, :spin)    && JuMP.add_to_expression!(objective, coef * DRESCOST, variable(am, nw, decomp_group, :demand_reserve_zdhty, (z,d,h,t,y)))
                        (am.setting["Simulation Configuration"]["Dispatch_Mode_in_EXP"] == "Unit Commitment") && reserve_enabled(am, :nonspin) && JuMP.add_to_expression!(objective, coef * NSRSP, variable(am, nw, decomp_group, :rns_nonspin_zdhty, (z,d,h,t,y)))
                        reserve_enabled(am, :flex_up) && JuMP.add_to_expression!(objective, coef * NSRSP, variable(am, nw, decomp_group, :rns_flex_up_zdhty, (z,d,h,t,y)))
                        reserve_enabled(am, :flex_dn) && JuMP.add_to_expression!(objective, coef * FLEXRSP, variable(am, nw, decomp_group, :rns_flex_dn_zdhty, (z,d,h,t,y)))

                    end
                end
            end
        end
    end

    # Ancillary penalties for storage commitment and hybrid ES grid charging
    for y in ids_y
        for i in ids_i
            bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
            tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)
            hybrid_type = am.ref[:nw][0][:gen_index][i]["hybrid_type"] 

            if parameter(am, 0, :gen_index, "UNIT_CATEGORY", i) == "STORAGE"
                if parameter(am, bus_idx, :gen_bus, tech_idx, "Storage Commitment") == true
                    for d in ids_d, h in ids_h, t in ids_t
                        JuMP.add_to_expression!(objective, tie_breaker_penalty, variable(am, nw, decomp_group, :sto_c_idhty, (i,d,h,t,y)))    
                    end 
                end

                if hybrid_type == "ES"
                    for d in ids_d, h in ids_h, t in ids_t
                        JuMP.add_to_expression!(objective, tie_breaker_penalty, variable(am, nw, decomp_group, :chg_idhty, (i,d,h,t,y)))    
                    end 
                end
            end
        end
    end

    # Clean energy generation slack
    clean_energy_generation_penalty = am.setting["Simulation Configuration"]["Clean_Energy_Generation_Penalty"]
    if am.setting["Simulation Configuration"]["Clean_Energy_Generation_Target_EXP_Flag"] == true
        for y in ids_y

            stage_length = am.ref[:nw][0][:planning_stages][y]["stage_length"]
            for future_y in 1:stage_length
    
                future_year = am.ref[:nw][0][:planning_stages][y]["year"] + future_y - 1 
                discount_factor = (1 + discount_rate)^-(future_year - current_year)

                coef = num_day_scale_OP * discount_factor * inv_cart_T
                
                if am.setting["Simulation Configuration"]["Clean_Energy_Generation_Target_EXP_Type"] == "Annual"
                    for n in ids_p
                        JuMP.add_to_expression!(objective, coef * clean_energy_generation_penalty, variable(am, nw, decomp_group, :slack_CEG_ny, (n,y)))  
                    end

                elseif am.setting["Simulation Configuration"]["Clean_Energy_Generation_Target_EXP_Type"] == "Daygroup"
                    for d in ids_d
                        num_day = parameter(am, 0, :repdays, "NumDays", d)
                        for n in ids_p
                            JuMP.add_to_expression!(objective, coef * num_day * clean_energy_generation_penalty, variable(am, nw, decomp_group, :slack_CEG_ndy, (n,d,y)))  
                        end
                    end
                end

            end
        end
    end


    # RPS slack
    if (am.setting["Simulation Configuration"]["RPS_Flag"] == true) && (am.setting["Simulation Configuration"]["Allow_Alternative_RPS_Compliance_Flag"] == true)
        
        penalty = am.setting["Simulation Configuration"]["RPS_Penalty"]

        for y in ids_y

            stage_length = am.ref[:nw][0][:planning_stages][y]["stage_length"]
            for future_y in 1:stage_length
    
                future_year = am.ref[:nw][0][:planning_stages][y]["year"] + future_y - 1 
                discount_factor = (1 + discount_rate)^-(future_year - current_year)

                for n in ids_p
                    JuMP.add_to_expression!(objective, discount_factor * penalty, variable(am, nw, decomp_group, :slack_RPS_ny, (n,y)))  
                end
            end
        end
    end

    # DEFINE OBJECTIVE FUNCTION
    return JuMP.@objective(JuMP_model, Min, objective     
    )

end


function objective_function_operation(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, ids_i, ids_d, ids_h, ids_t, ids_y, ids_z, ids_p, ids_i_commit, ids_k, ids_i_water, ids_lfl, ids_n, day_group_id; nw::Int=am.cnw, report::Bool=true)

    decomp_group = day_group_id  # Only one group when PH is not being used.
    objective = JuMP.AffExpr(0.0)

    # Anti-degeneracy tie-breaker penalty (per-unit). Negligible vs. the smallest real
    # per-unit marginal cost (~3e-3), so dispatch is unchanged; lifts the objective floor.
    tie_breaker_penalty = 1e-6

    y = ids_y[1]
    future_year = am.ref[:nw][0][:planning_stages][y]["year"]
    card_T = length(am.setting["run_T"])
    inv_cart_T = 1.0 / card_T

    # Minimum regulation-reserve cost floor ($/MWh). Ensures regulation (reg-up + reg-down)
    # is never free. Safe default 0.0 keeps existing workbooks unchanged.
    pu_econ_base = am.setting["Simulation Setting"]["per_unit_econ_base_value"] / am.setting["Simulation Setting"]["per_unit_base_value"]
    # Divide the raw $/MWh value by pu_econ_base to match Annual_MC's per-unit basis,
    # so the floor compares like-for-like against reg_cost in the objective.
    min_reg_cost = get(am.setting["Planning Design"], "min_regulation_cost_value", 0.0) / pu_econ_base

    # LFL demand response cost
    for lfl in ids_lfl
        local_lfl = am.ref[:nw][0][:demand][lfl]
        num_seg = local_lfl["Num_DR_Segments"]

        # skip s==1 by construction
        for s in 2:num_seg
            dr_cost = local_lfl["Price_$(s)"]  # computed once per s

            for d in ids_d 
                num_days = parameter(am, 0, :repdays, "NumDays", d)
                coef = -(num_days * dr_cost * inv_cart_T)
                for h in ids_h, t in ids_t
                    JuMP.add_to_expression!(objective, coef, variable(am, nw, decomp_group, :lfl_seg_lsdhty, (lfl,s,d,h,t,y)))  
                end
            end
        end

        # Ancillary penalties for hybrid ES grid charging
        if am.ref[:nw][0][:demand][lfl]["Hybrid_ES"] != "NA"
            for d in ids_d, h in ids_h, t in ids_t
                JuMP.add_to_expression!(objective, tie_breaker_penalty, variable(am, nw, decomp_group, :lfl_chg_Grid_ES_ldhty, (lfl,d,h,t,y)))  
            end
        end
    end

    # Hourly Costs, divide by number of timesteps in an hour (card(T))
    
    for i in ids_i

        bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
        tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)
        CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")

        CTAX = am.setting["Simulation Configuration"]["CTAX"]
        Emission_rate = parameter(am, bus_idx, :gen_bus, tech_idx, "Emission_CO2")
        # Emissions Cost ($/MW) = $/tonne * kg/MMbtu * MMbtu/MWh * MWh * 1 tonne/1000kg = $/tonne * metric tonnes CO2 per timestep

        reg_cost = parameter(am, bus_idx, :gen_bus, tech_idx, "reg_cost")
        spin_cost = parameter(am, bus_idx, :gen_bus, tech_idx, "spin_cost")
        nspin_cost = parameter(am, bus_idx, :gen_bus, tech_idx, "nspin_cost")
        flex_cost = parameter(am, bus_idx, :gen_bus, tech_idx, "flex_cost")
        start_up_cost = parameter(am, bus_idx, :gen_bus, tech_idx, "SUC")
        no_load_cost = parameter(am, bus_idx, :gen_bus, tech_idx, "NLC")

        for d in ids_d 
            
            num_days = parameter(am, 0, :repdays, "NumDays", d)
            Annual_MC = parameter(am, bus_idx, :gen_bus, tech_idx, "Annual_MC")[string(y)][d]

            if am.setting["Planning Design"]["reserve_cost_type_flag"] == "percentage"
                reg_cost = parameter(am, bus_idx, :gen_bus, tech_idx, "reg_cost") * Annual_MC
                spin_cost = parameter(am, bus_idx, :gen_bus, tech_idx, "spin_cost") * Annual_MC
                nspin_cost = parameter(am, bus_idx, :gen_bus, tech_idx, "nspin_cost") * Annual_MC
                flex_cost = parameter(am, bus_idx, :gen_bus, tech_idx, "flex_cost") * Annual_MC
            end

            # Apply the minimum regulation cost floor ($/MWh) so regulation is never free.
            reg_cost = max(reg_cost, min_reg_cost)

            coef = num_days * inv_cart_T

            for h in ids_h, t in ids_t

                # Generation Cost
                JuMP.add_to_expression!(objective, coef * Annual_MC, variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y)))  

                # Production Tax Credit
                if am.setting["Simulation Configuration"]["PTC_Flag"] == true
                    if parameter(am, bus_idx, :gen_bus, tech_idx, "PTC Flag") == true
                        # ptc = cents/kWh -> dollar/MWh
                        ptc_year = min(2050, future_year)
                        ptc = parameter(am, bus_idx, :gen_bus, tech_idx, "PTC")[string(ptc_year)] * 10
                        JuMP.add_to_expression!(objective, -1 * coef * ptc, variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y)))  
                    end
                end

                up_res_eligible = up_reserve_eligible(am, i)  # up-reserves exist only for dispatchable units
                # Regulation Cost = x% of generation cost
                reserve_enabled(am, :reg_up) && up_res_eligible && JuMP.add_to_expression!(objective, coef * reg_cost, variable(am, nw, decomp_group, :reg_up_idhty, (i,d,h,t,y)))
                reserve_enabled(am, :reg_dn) && JuMP.add_to_expression!(objective, coef * reg_cost, variable(am, nw, decomp_group, :reg_dn_idhty, (i,d,h,t,y)))

                reserve_enabled(am, :spin) && up_res_eligible && JuMP.add_to_expression!(objective, coef * spin_cost, variable(am, nw, decomp_group, :spin_idhty, (i,d,h,t,y)))

                if am.setting["Simulation Configuration"]["Dispatch_Mode_in_OP"] == "Unit Commitment"
                    reserve_enabled(am, :nonspin) && JuMP.add_to_expression!(objective, num_days * nspin_cost / card_T, variable(am, nw, decomp_group, :nonspin_idhty, (i,d,h,t,y)))
                end

                # Flex reserve Cost = 60% of regulation cost
                reserve_enabled(am, :flex_up) && up_res_eligible && JuMP.add_to_expression!(objective, coef * flex_cost, variable(am, nw, decomp_group, :flex_up_idhty, (i,d,h,t,y)))
                reserve_enabled(am, :flex_dn) && JuMP.add_to_expression!(objective, coef * flex_cost, variable(am, nw, decomp_group, :flex_dn_idhty, (i,d,h,t,y)))

                # Emission
                JuMP.add_to_expression!(objective, coef * CTAX * Emission_rate, variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y)))  

            end

            # Start up and no load cost
            if (i in ids_i_commit) && (am.setting["Simulation Configuration"]["Dispatch_Mode_in_OP"] == "Unit Commitment")

                # Startup Cost ($/MW)
                JuMP.add_to_expression!(objective, coef * start_up_cost, variable(am, nw, decomp_group, :su_idhty, (i,d,h,1,y)))  

                # No Load Cost ($/MW)
                JuMP.add_to_expression!(objective, coef * no_load_cost, variable(am, nw, decomp_group, :c_idhty, (i,d,h,1,y)))  


            end

        end
    end

    # Scarcity
    for d in ids_d
        VOLL = am.setting["Simulation Configuration"]["VOLL"]
        RegRSP = am.setting["Simulation Configuration"]["RegRSP"]
        SRSP = am.setting["Simulation Configuration"]["SRSP"]
        NSRSP = am.setting["Simulation Configuration"]["NSRSP"]
        FLEXRSP = am.setting["Simulation Configuration"]["FLEXRSP"]
        # Load-Resource reserve priced above generation reserve, far below the RNS penalty.
        DRESCOST = get(am.setting["Simulation Configuration"], "demand_reserve_cost", 5.0)

        coef = parameter(am, 0, :repdays, "NumDays", d) * inv_cart_T
        
        for h in ids_h, t in ids_t 

            # Load Shedding
            for n in ids_n
                JuMP.add_to_expression!(objective, coef * VOLL, variable(am, nw, decomp_group, :ens_ndhty, (n,d,h,t,y)))  
            end

            # Reserve Shortage
            for z in ids_z
                reserve_enabled(am, :reg_up)  && JuMP.add_to_expression!(objective, coef * RegRSP, variable(am, nw, decomp_group, :rns_reg_up_zdhty, (z,d,h,t,y)))
                reserve_enabled(am, :reg_dn)  && JuMP.add_to_expression!(objective, coef * RegRSP, variable(am, nw, decomp_group, :rns_reg_dn_zdhty, (z,d,h,t,y)))
                reserve_enabled(am, :spin)    && JuMP.add_to_expression!(objective, coef * SRSP, variable(am, nw, decomp_group, :rns_cont_zdhty, (z,d,h,t,y)))
                reserve_enabled(am, :spin)    && JuMP.add_to_expression!(objective, coef * DRESCOST, variable(am, nw, decomp_group, :demand_reserve_zdhty, (z,d,h,t,y)))
                reserve_enabled(am, :nonspin) && JuMP.add_to_expression!(objective, coef* NSRSP, variable(am, nw, decomp_group, :rns_nonspin_zdhty, (z,d,h,t,y)))
                reserve_enabled(am, :flex_up) && JuMP.add_to_expression!(objective, coef * FLEXRSP, variable(am, nw, decomp_group, :rns_flex_up_zdhty, (z,d,h,t,y)))
                reserve_enabled(am, :flex_dn) && JuMP.add_to_expression!(objective, coef * FLEXRSP, variable(am, nw, decomp_group, :rns_flex_dn_zdhty, (z,d,h,t,y)))
            end
        end
    end

    # Clean energy generation slack
    if am.setting["Simulation Configuration"]["Clean_Energy_Generation_Target_OP_Flag"] == true
        Clean_Energy_Generation_Penalty = am.setting["Simulation Configuration"]["Clean_Energy_Generation_Penalty"]
        
        for d in ids_d
            coef = parameter(am, 0, :repdays, "NumDays", d) * Clean_Energy_Generation_Penalty * inv_cart_T
            for n in ids_p
                JuMP.add_to_expression!(objective, coef, variable(am, nw, decomp_group, :slack_CEG_ndy, (n,d,y)))  
            end
        end
    end

    # Ancillary penalties for storage commitment and hybrid ES grid charging
    for i in ids_i
        bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
        tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)
        hybrid_type = am.ref[:nw][0][:gen_index][i]["hybrid_type"] 

        if parameter(am, 0, :gen_index, "UNIT_CATEGORY", i) == "STORAGE"
            if parameter(am, bus_idx, :gen_bus, tech_idx, "Storage Commitment") == true
                for d in ids_d, h in ids_h, t in ids_t
                    JuMP.add_to_expression!(objective, tie_breaker_penalty, variable(am, nw, decomp_group, :sto_c_idhty, (i,d,h,t,y)))    
                end 
            end

            if hybrid_type == "ES"
                for d in ids_d, h in ids_h, t in ids_t
                    JuMP.add_to_expression!(objective, tie_breaker_penalty, variable(am, nw, decomp_group, :chg_idhty, (i,d,h,t,y)))    
                end 
            end
        end
    end

    # DEFINE OBJECTIVE FUNCTION
    return JuMP.@objective(JuMP_model, Min, objective     
    )

end
