# Copyright © 2026, UChicago Argonne, LLC. All Rights Reserved.
# SPDX-License-Identifier: BSD-3-Clause  (see LICENSE)

# Top-level ALEAF run orchestration and model dispatch.

using Logging
using LoggingExtras
using Dates

const _ALEAF_LOG_MODEL = Ref("ALEAF")

function _aleaf_logger(; format::String="simple", level::Logging.LogLevel=Logging.Info)
    detailed = lowercase(format) == "detailed"
    logger = FormatLogger(stdout) do io, args
        t = Dates.now()
        ts = Dates.format(t, "yyyy-mm-dd HH:MM:SS")
        lvl = uppercase(string(args.level))
        model = _ALEAF_LOG_MODEL[]
        msg = string(args.message)
        color = if lvl == "ERROR"
            :red
        elseif lvl == "WARN"
            :yellow
        elseif lvl == "DEBUG"
            :light_black
        else
            :green
        end

        if detailed
            ms = lpad(string(Dates.value(Dates.Millisecond(t)) % 1000), 3, '0')
            mod = hasproperty(args, :_module) ? getfield(args, :_module) : getfield(args, :module)
            mod_str = mod === nothing ? "" : string(mod)
            func = ""
            if hasproperty(args, :kwargs)
                func_val = get(args.kwargs, :func, nothing)
                if func_val !== nothing
                    func = string(func_val)
                end
            end
            loc = ""
            if hasproperty(args, :file) && hasproperty(args, :line) && args.file !== nothing && args.line !== nothing
                loc = string(Base.basename(String(args.file)), ":", args.line)
            end
            detail = ""
            if mod_str != "" && func != "" && loc != ""
                detail = string(mod_str, ":", func, ":", loc)
            elseif mod_str != "" && loc != ""
                detail = string(mod_str, ":", loc)
            elseif mod_str != ""
                detail = mod_str
            end

            if detail != ""
                print(io, "$ts,$ms | ")
                printstyled(io, lvl, color=color)
                println(io, " | $model | $detail | $msg")
            else
                print(io, "$ts,$ms | ")
                printstyled(io, lvl, color=color)
                println(io, " | $model | $msg")
            end
        else
            print(io, "$ts | ")
            printstyled(io, lvl, color=color)
            println(io, " | $msg")
        end
    end

    return MinLevelLogger(logger, level)
end

function normalize_master_setting_file_name(master_setting_file_name::String)

    if isempty(strip(master_setting_file_name))
        error("master_setting_file_name is required. Specify a setting workbook in the setting folder.")
    end

    if !endswith(master_setting_file_name, ".xlsx")
        master_setting_file_name *= ".xlsx"
    end

    return master_setting_file_name
end


function read_ALEAF_setting(ALEAF_setting::Dict{String,<:Any})

    setting_file_path = string("setting/", ALEAF_setting["ALEAF Master Setup"]["Master_setting_file_name"])

    ALEAF_setting["Simulation Configuration"] = read_xlsx_return_dict_string_any_for_sim_config(setting_file_path, "Simulation Configuration"; first_row_value = 2)
    setting_category_list = ["Simulation Setting", "Planning Design", "Scenario Reduction Setting", "RA Setting"]
            
    for setting_category in setting_category_list
        ALEAF_setting[setting_category] = Dict{String, Any}()
        data = DataFrame(XLSX.readtable(setting_file_path, setting_category))
        for row in eachrow(data)
            ALEAF_setting[setting_category][values(row)[1]] = values(row)[2]
        end
    end

    tab_names_list = ["External Constraints", "ATB Setting", "Storage Cost and Performance", "Raw Materials", "Gen Technology Raw Materials", "Gen Technology", "RA Scenarios"]
    for data_category in tab_names_list
        ALEAF_setting[data_category] = read_xlsx_return_dict_string_any(setting_file_path, data_category; first_row_value = 2)
    end

    for data_category in ["ITC", "PTC"]
        ALEAF_setting[data_category] = read_xlsx_return_dict_string_any(setting_file_path, data_category)
    end

    solver_setting_type = string(ALEAF_setting["Simulation Setting"]["solver_name"], " Setting")
    ALEAF_setting[solver_setting_type] = read_xlsx_return_dict_string_any(setting_file_path, solver_setting_type)
    
end


function initiate_ALEAF(; worker_info::Dict{Any, Any}=Dict(), master_setting_file_name="")

    ALEAF_setting = Dict{String, Any}()
    ALEAF_setting["ALEAF Master Setup"] = Dict{String, Any}()
    ALEAF_setting["ALEAF Master Setup"]["Master_setting_file_name"] = normalize_master_setting_file_name(master_setting_file_name)

    @aleaf_info string("Simulation setting file: ", ALEAF_setting["ALEAF Master Setup"]["Master_setting_file_name"])

    "read ALEAF module setting and data"
    read_ALEAF_setting(ALEAF_setting)

    "check model type"
    model_type = "LC_GTEP"
    ALEAF_setting["ALEAF Master Setup"]["model_type"] = model_type
    ALEAF_setting["ALEAF Master Setup"]["solver_name"] = ALEAF_setting["Simulation Setting"]["solver_name"]
    if model_type == "LC_GTEP"
        ALEAF_setting["ALEAF Master Setup"]["ALEAF_model_type"] = Abstract_LC_GTEP_Model
    end

    return ALEAF_setting
end


function build_model_structure_GTEP(ALEAF_setting::Dict{String,<:Any}, data::Dict{String,<:Any}, model_type::Type; build_jump_models::Bool=true, kwargs...)

    direct_mode_flag = 0
    if (ALEAF_setting["Simulation Setting"]["solver_name"] == "CPLEX")
        if ALEAF_setting["CPLEX Setting"]["1"]["Value"] == true
            direct_mode_flag = 1 # CPLEX misbehaves in auto mode here; force direct mode.
        end
    end

    if direct_mode_flag == 1
        imo = initialize_model_instance_direct_mode_GTEP(model_type, data; build_jump_models=build_jump_models)
    else
        imo = initialize_model_instance_auto_mode_GTEP(model_type, data; build_jump_models=build_jump_models)
    end
    @aleaf_debug "Initialize model time:"

    return imo
end


function initialize_model_instance_RA(model_type::Type, outage_scenario_list)

    @assert model_type <: Abstract_ALEAF_Model

    setting = Dict{String,Any}()

    ref = Dict{Symbol,Any}()

    var = Dict{Symbol,Any}(:nw => Dict{Int,Any}())
    con = Dict{Symbol,Any}(:nw => Dict{Int,Any}())
    sol = Dict{Symbol,Any}(:nw => Dict{Int,Any}())
    jump_model = Dict{Symbol,Any}(:nw => Dict{Int,Any}())

    var[:nw][1] = Dict{Int,Any}()
    con[:nw][1] = Dict{Int,Any}()
    sol[:nw][1] = Dict{Int,Any}()
    jump_model[:nw][1] = Dict{Int,Any}()
    for scenario_id in outage_scenario_list
        var[:nw][1][scenario_id] = Dict{Symbol,Any}()
        con[:nw][1][scenario_id] = Dict{Symbol,Any}()
        sol[:nw][1][scenario_id]= Dict{Symbol,Any}()
        jump_model[:nw][1][scenario_id] = JuMP.Model()
    end

    solution = Dict{String,Any}()
    cnw = 1
    
    imo = ALEAF_Model_Structure_RA(
        jump_model,
        string(model_type),
        setting,
        solution, 
        ref,
        var,
        con,
        sol,
        cnw
    )

    return imo
end


function initialize_model_instance_direct_mode_GTEP(model_type::Type, data::Dict{String,<:Any}; optimizer=nothing, build_jump_models::Bool=true)

    @assert model_type <: Abstract_ALEAF_Model

    setting = Dict{String,Any}()

    ref = network_ref_initialize(model_type, data)

    var = Dict{Symbol,Any}(:nw => Dict{Int,Any}())
    con = Dict{Symbol,Any}(:nw => Dict{Int,Any}())
    sol = Dict{Symbol,Any}(:nw => Dict{Int,Any}())

    jump_model = Dict{Symbol,Any}(:nw => Dict{Int,Any}())

    # build_jump_models=false skips the empty JuMP.Model() allocations for an instance that is never
    # solved (e.g. the reporting-only reference built on the master). With NDAY_Groups_OP=365 that
    # avoids creating 365 unused JuMP models per subarea.
    # Even when build_jump_models=true, only the cnw network (== length(ref[:nw]), the single network
    # every build/solve/export uses via am.cnw) needs a real JuMP.Model(). The other sub-area nws keep
    # empty scaffolding and are never optimized; at nodal scale allocating a model for each was
    # O(buses × repday_groups) pure waste.
    cnw_id = length(ref[:nw])
    for (nw_id, nw) in ref[:nw]
        var[:nw][nw_id] = Dict{Int,Any}()
        con[:nw][nw_id] = Dict{Int,Any}()
        sol[:nw][nw_id] = Dict{Int,Any}()
        jump_model[:nw][nw_id] = Dict{Int,Any}()
        for repdays in keys(data["repday_groups"])
            var[:nw][nw_id][parse(Int, repdays)] = Dict{Symbol,Any}()
            con[:nw][nw_id][parse(Int, repdays)] = Dict{Symbol,Any}()
            sol[:nw][nw_id][parse(Int, repdays)] = Dict{Symbol,Any}()
            if build_jump_models && nw_id == cnw_id
                jump_model[:nw][nw_id][parse(Int, repdays)] = JuMP.Model()
            end
        end
    end

    solution = Dict{String,Any}()
    cnw = length(var[:nw])
    PH = Dict{Symbol,Any}(
        :PH_var => Dict{String, Any}(),
        :PH_var_res => Dict{String, Any}(),
        :PH_var_hat_res => Dict{String, Any}(),
        :PH_status => Dict{String, Any}(),
        :PH_result => Dict{String, Any}(),
        :PH_hash => Dict{String, Any}()
    )

    imo = ALEAF_Model_Structure_PH(
        jump_model,
        string(model_type),
        data,
        setting,
        solution, 
        ref,
        var,
        con,
        sol,
        cnw,
        PH
    )

    return imo
end


function initialize_model_instance_auto_mode_GTEP(model_type::Type, data::Dict{String,<:Any}; build_jump_models::Bool=true)
    @assert model_type <: Abstract_ALEAF_Model

    setting = Dict{String,Any}()

    ref = network_ref_initialize(model_type, data)

    var = Dict{Symbol,Any}(:nw => Dict{Int,Any}())
    con = Dict{Symbol,Any}(:nw => Dict{Int,Any}())
    sol = Dict{Symbol,Any}(:nw => Dict{Int,Any}())

    jump_model = Dict{Symbol,Any}(:nw => Dict{Int,Any}())

    # Only the cnw network (== length(ref[:nw]), the single network every build/solve/export uses via
    # am.cnw) needs a real JuMP.Model(). Other sub-area nws keep empty scaffolding and are never
    # optimized; at nodal scale allocating a model for each nw × repday_group was O(buses) pure waste.
    # build_jump_models=false skips all of them (reporting-only instance that is never solved).
    cnw_id = length(ref[:nw])
    for nw_id in 1:length(ref[:nw])
        var[:nw][nw_id] = Dict{Int,Any}()
        con[:nw][nw_id] = Dict{Int,Any}()
        sol[:nw][nw_id] = Dict{Int,Any}()
        jump_model[:nw][nw_id] = Dict{Int,Any}()

        for repdays in keys(data["repday_groups"])
            var[:nw][nw_id][parse(Int, repdays)] = Dict{Symbol,Any}()
            con[:nw][nw_id][parse(Int, repdays)] = Dict{Symbol,Any}()
            sol[:nw][nw_id][parse(Int, repdays)] = Dict{Symbol,Any}()
            if build_jump_models && nw_id == cnw_id
                jump_model[:nw][nw_id][parse(Int, repdays)] = JuMP.Model()
            end
        end

    end

    solution = Dict{String,Any}()
    cnw = length(var[:nw])
    PH = Dict{Symbol,Any}(
        :PH_var => Dict{String, Any}(),
        :PH_var_res => Dict{String, Any}(),
        :PH_var_hat_res => Dict{String, Any}(),
        :PH_status => Dict{String, Any}(),
        :PH_result => Dict{String, Any}(),
        :PH_hash => Dict{String, Any}()
    )

    imo = ALEAF_Model_Structure_PH(
        jump_model,
        string(model_type),
        data,
        setting,
        solution, 
        ref,
        var,
        con,
        sol,
        cnw,
        PH
    )

    return imo
end


function network_ref_initialize(model_type::Type, data::Dict{String,<:Any})

    refs = Dict{Symbol,Any}()

    if model_type == Abstract_LC_GTEP_Model
        # create multiple nws data for each sub area.
        # gen_technology is identical for every sub area and downstream it is only READ (then deepcopied
        # into a per-subarea :gen_bus at aggregate_generation, and its key deleted). So share ONE copy
        # across all nws instead of deepcopying per bus — at nodal scale (25k+ buses) the per-bus
        # deepcopy was O(buses) redundant. One deepcopy keeps it isolated from the source data dict.
        nws_data = Dict{String,Any}()
        shared_gen_technology = deepcopy(data["gen_technology"])
        for id in keys(data["bus"])
            nws_data[id] = Dict{String,Any}()
            nws_data[id]["gen_technology"] = shared_gen_technology
        end
    else
        nws_data = Dict("0" => data)
    end
    
    nws = refs[:nw] = Dict{Int,Any}()

    for (n, nw_data) in nws_data
        nw_id = parse(Int, n)
        ref = nws[nw_id] = Dict{Symbol,Any}()

        for (key, item) in nw_data
            if isa(item, Dict{String,Any})
                item_lookup = Dict{Int,Any}([(parse(Int, k), v) for (k,v) in item])
                ref[Symbol(key)] = item_lookup
            elseif isa(item, Dict{Int64,Any})
                item_lookup = Dict{Int,Any}([(k, v) for (k,v) in item])
                ref[Symbol(key)] = item_lookup
            else
                ref[Symbol(key)] = item
            end
        end
    end

    return refs
end


function run_ALEAF(; worker_info::Dict{Any, Any}=Dict(), master_setting_file_name="")

    global_logger(_aleaf_logger(; format="simple", level=Logging.Info))
    @aleaf_info "Initiate ALEAF"
    ALEAF_setting = initiate_ALEAF(; worker_info, master_setting_file_name)

    "Check ALEAF model type"
    model_type = ALEAF_setting["ALEAF Master Setup"]["ALEAF_model_type"]
    @assert model_type <: Abstract_ALEAF_Model
    _ALEAF_LOG_MODEL[] = string(model_type)
    sim_setting = get(ALEAF_setting, "Simulation Setting", Dict{String, Any}())
    log_setting = get(sim_setting, "logging_level_value", "simple")
    log_setting_str = lowercase(string(log_setting))
    effective_format = log_setting_str == "detailed" ? "detailed" : "simple"
    log_label = effective_format == "detailed" ? "Detailed" : "Simple"
    global_logger(_aleaf_logger(; format=effective_format, level=Logging.Info))
    @aleaf_info string("Logging format: ", log_label)
    @aleaf_info string("Model type: ", model_type)

    "Execute ALEAF"
    @aleaf_info "Execute ALEAF"
    execution_status = "success"
    if model_type == Abstract_LC_GTEP_Model
        execution_status = run_LC_GTEP_Model(ALEAF_setting)
    end
    
    if execution_status == "success" 
        @aleaf_info "Execute ALEAF: DONE!"
    else
        @aleaf_error "ALEAF execution ended with ERRORS"
    end
   

end


function run_LC_GTEP_Model(ALEAF_setting::Dict{String,<:Any}, kwargs...)

    execution_status =build_run_LC_GTEP_model(ALEAF_setting)

    return execution_status

end




