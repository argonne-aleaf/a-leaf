# Copyright © 2026, UChicago Argonne, LLC. All Rights Reserved.
# SPDX-License-Identifier: BSD-3-Clause  (see LICENSE)

# Purpose: Reliability assessment (RA) workflow and simulation logic.
# Scope: Build outage scenarios, run ED-based RA, and compute metrics.
# Notes: Uses shared components plus RA-specific functions in this file.


using XLSX
using DataFrames
using JSON
using Statistics
using CSV
using Random
using Distributions
using Serialization
using JLD2
using Base.Threads: @spawn
using Distributed: workers, nprocs, myid, remotecall_fetch, @spawnat, RemoteChannel
using Dates

    

function execute_ALEAF_RA_model_using_external_data(ALEAF_setting::Dict{String,<:Any}, case_id; RA_input::Dict{String,<:Any} = Dict{String,Any}(), external_GTEP_multi_round_info::Dict{String,<:Any} = Dict{String,Any}())

    output_path = define_output_path(ALEAF_setting, case_id)

    round_id = ALEAF_setting["RA Setting"]["round_id_to_run_RA"]

    if length(external_GTEP_multi_round_info) != 0

        @aleaf_info "[ALEAF RA Model]\tRunning RA using given generation mix"

        info = is_single_round_expansion_record(external_GTEP_multi_round_info) ?
            convert_single_round_record_to_multi_round_info(external_GTEP_multi_round_info, round_id) :
            external_GTEP_multi_round_info

        updated_investment_decisions = info[string(round_id)]["updated_investment_decisions"]
        current_year = info[string(round_id)]["round_ids_y_decision"][1]

        RA_Info = execute_ALEAF_RA_model(ALEAF_setting, case_id; RA_input, current_year, recorded_investment_decisions=updated_investment_decisions)

    else
        if ispath(joinpath(output_path, "GTEP_multi_round_info.json"))

            prior_GTEP_multi_round_info = JSON.parse(open(joinpath(output_path, "GTEP_multi_round_info.json")))
            updated_investment_decisions = prior_GTEP_multi_round_info[string(round_id)]["updated_investment_decisions"]
            current_year = prior_GTEP_multi_round_info[string(round_id)]["round_ids_y_decision"][1]
            @aleaf_info "[ALEAF RA Model]\tRunning RA using saved generation mix"

            RA_Info = execute_ALEAF_RA_model(ALEAF_setting, case_id; RA_input, current_year, recorded_investment_decisions=updated_investment_decisions)

        else

            @aleaf_info "[ALEAF RA Model]\tERROR: Recorded GTEP_multi_round_info file is not available"

        end
    end


end


function execute_ALEAF_RA_model(ALEAF_setting, case_id; RA_input::Dict{String,<:Any} = Dict{String,Any}(), current_year=1, recorded_investment_decisions::Dict{String,<:Any} = Dict{String,Any}(), predefined_network_data::Dict{String,<:Any} = Dict{String,Any}())

    # RA_input optional overrides: setting_type, risk_scenario_dict, filtered_joint_scenario_map,
    # ELCC_existing_asset_target_list (each: target_PLANT_NAME, target_UNIT_GROUP; opt target_name).

    @aleaf_info "[ALEAF RA Model]: Perform Probabilistic Reliability Assessment Simulation"

    run_expansion_flag = false
    gen_disaggregation_flag = false

    RA_ALEAF_setting = deepcopy(ALEAF_setting)

    update_solver_setting!(RA_ALEAF_setting)

    RA_Info = update_RA_setting(RA_ALEAF_setting, case_id, RA_input, current_year)

    if length(predefined_network_data) > 0
        network_data = deepcopy(predefined_network_data)
    else
        network_data = generate_networkdata_LC_GTEP(RA_ALEAF_setting, case_id, "RA"; print_output_flag=true, year_id=current_year)
    end
    @aleaf_info "[ALEAF RA Model]: Network Generation"

    network_data["output_path"] = define_output_path(ALEAF_setting, case_id)

    ALEAF_model_instance = ALEAF.build_ALEAF_model_instance_for_RA(RA_ALEAF_setting, case_id, RA_Info["setting"]["current_year"], network_data, recorded_investment_decisions);
    @aleaf_info "[ALEAF RA Model]: Build Reference System"

    if length(recorded_investment_decisions) > 0
        update_expansion_results!(ALEAF_model_instance, recorded_investment_decisions)
    end

    # Disaggregation happens here.
    prepare_RA_gen!(ALEAF_model_instance)
    @aleaf_info "[ALEAF RA Model]: Update Network Data for RA"

    risk_scenario_dict = Dict{Int, Any}()
    risk_scenario_dict, RA_Info["RA_reference_risk_data"]["filtered_joint_scenario_map"] = generate_and_filter_risk_scenarios(RA_ALEAF_setting, RA_Info, ALEAF_model_instance; RA_input)
    @aleaf_info "[ALEAF RA Model]: Risk Sampling Completed"

    output_path = ALEAF_model_instance.ref[:nw][0][:output_path]
    outage_cache_path = joinpath(output_path, "RA_outage_cache")
    _write_risk_scenario_cache!(outage_cache_path, risk_scenario_dict; prefix="day")
    RA_Info["setting"]["outage_cache_path"] = outage_cache_path
    RA_Info["setting"]["outage_cache_format"] = "jld2"
    @aleaf_info "[ALEAF RA Model]: Outage scenario cache written to $outage_cache_path"
    
    total_sce = 0
    total_num_days = 0
    for day_idx in keys(RA_Info["RA_reference_risk_data"]["filtered_joint_scenario_map"])
        num_of_risk = length(RA_Info["RA_reference_risk_data"]["filtered_joint_scenario_map"][day_idx])
        if num_of_risk > 0
            total_sce += num_of_risk
            total_num_days += 1
        end
    end
    total_num_risk = RA_ALEAF_setting["Simulation Configuration"][string(case_id)]["NDAY_Groups_RA"] * RA_Info["setting"]["num_risk_scenario"] * RA_Info["setting"]["num_renewable_scenarios"]
    identified_risk_ratio = total_num_risk > 0 ? (total_sce / total_num_risk) * 100 : 0.0
    @aleaf_info "[ALEAF RA Model]: Number of joint scenarios: $total_sce ($(round(identified_risk_ratio, digits=3)) %, in $total_num_days day groups)"

    # Run RA simulations
    RA_Info = ALEAF.execute_RA_model_with_ED!(RA_Info, ALEAF_model_instance, risk_scenario_dict)    

    # Assess capacity credits
    if RA_Info["setting"]["calculate_capacity_credit_flag"] == true
        @aleaf_info "[ALEAF RA Model]: Starting Capacity Credit Analysis"

        if RA_Info["setting"]["capacity_credit_type"] == "DLOL"
            @aleaf_info "[ALEAF RA Model]: Capacity Credit Type: DLOL (Direct Loss of Load)"
            RA_Info["capacity_credit_result"] = Dict{String, Any}()
            RA_Info["capacity_credit_result"]["capacity_credit_type"] = "DLOL"

            if haskey(RA_Info, "RA_solutions") && haskey(RA_Info["RA_solutions"], "DLOL")
                RA_Info["capacity_credit_result"]["DLOL"] = deepcopy(RA_Info["RA_solutions"]["DLOL"])
                @aleaf_info "[ALEAF RA Model]: DLOL result wired to capacity_credit_result."
            else
                @aleaf_warn "[ALEAF RA Model]: DLOL requested but RA_solutions[\"DLOL\"] is missing. Check RA method and DLOL flags."
            end

        elseif RA_Info["setting"]["capacity_credit_type"] == "ELCC"
            @aleaf_info "[ALEAF RA Model]: Capacity Credit Type: ELCC"

            # Drop ENS-free risk scenarios up front; they cannot bind ELCC and slow it down.
            source_filtered_joint_scenario_map = get(
                RA_Info["RA_reference_risk_data"],
                "filtered_joint_scenario_map",
                Dict{Int64, Dict{Int, Dict{String, Any}}}()
            )
            filtered_joint_scenario_map = Dict{Int64, Dict{Int, Dict{String, Any}}}()
            for day_idx in keys(source_filtered_joint_scenario_map)
                retained_joint_map = Dict{Int, Dict{String, Any}}()
                for (joint_id, joint_scenario) in source_filtered_joint_scenario_map[day_idx]
                    if get(RA_Info["RA_solutions"][string(day_idx)][joint_id], "total ENS", 0.0) > 0.0
                        retained_joint_map[joint_id] = deepcopy(joint_scenario)
                    end
                end
                filtered_joint_scenario_map[day_idx] = retained_joint_map
            end
            RA_Info["RA_reference_risk_data"]["filtered_joint_scenario_map"] = filtered_joint_scenario_map

            RA_Info["RA_solutions"] = Dict{String, Any}()
            GC.gc(true)

            total_sce = 0
            total_num_days = 0
            for day_idx in keys(filtered_joint_scenario_map)
                num_of_risk = length(filtered_joint_scenario_map[day_idx])
                if num_of_risk > 0
                    total_sce += length(filtered_joint_scenario_map[day_idx])
                    total_num_days += 1
                end
            end
            @aleaf_info "[ALEAF RA Model]: Remaining joint risk scenarios for ELCC: $total_sce in $total_num_days days"
            
            elcc_assessment_mode = RA_Info["setting"]["capacity_credit_assessment_mode"]

            if (elcc_assessment_mode != "add_new") && (elcc_assessment_mode != "deactivate_existing")
                @aleaf_warn "[ALEAF RA Model]: Unknown ELCC_assessment_mode=$(elcc_assessment_mode); defaulting to add_new."
                elcc_assessment_mode = "add_new"
            end
            
            if elcc_assessment_mode == "deactivate_existing"
                RA_Info["capacity_credit_result"] = Dict{String, Any}()
                RA_Info["capacity_credit_result"]["capacity_credit_type"] = "ELCC"
                RA_Info["capacity_credit_result"]["capacity_credit_assessment_mode"] = elcc_assessment_mode
                RA_Info["capacity_credit_result"]["ELCC_method"] = RA_Info["setting"]["capacity_credit_RA_simulation_method"]
                RA_Info["capacity_credit_result"]["reference_RA_metric"] = RA_Info["setting"]["capacity_credit_reference_RA_metric"]
                RA_Info["capacity_credit_result"]["reference_RA_metric_spatial_resolution"] = RA_Info["setting"]["capacity_credit_reference_RA_metric_spatial_resolution"]
                RA_Info["capacity_credit_result"]["targets"] = Dict{String, Any}()

                target_list = get(RA_input, "ELCC_existing_asset_target_list", Any[])
                if isempty(target_list)
                    @aleaf_warn "[ALEAF RA Model]: ELCC_existing_asset_target_list is empty for deactivate_existing mode; skipping ELCC deactivation runs."
                end

                for (idx, target_cfg_any) in enumerate(target_list)
                    if !(target_cfg_any isa AbstractDict)
                        @aleaf_warn "[ALEAF RA Model]: ELCC_existing_asset_target_list[$idx] is not a Dict; skipping."
                        continue
                    end
                    target_cfg = target_cfg_any
                    target_plant_name = string(get(target_cfg, "target_PLANT_NAME", ""))
                    target_unit_group = string(get(target_cfg, "target_UNIT_GROUP", ""))

                    if isempty(target_plant_name) || isempty(target_unit_group)
                        @aleaf_warn "[ALEAF RA Model]: ELCC_existing_asset_target_list[$idx] missing target_PLANT_NAME or target_UNIT_GROUP; skipping."
                        continue
                    end

                    target_key_default = string(target_plant_name, "|", lowercase(target_unit_group))
                    target_key = string(get(target_cfg, "target_name", target_key_default))

                    RA_Info["capacity_credit_result"]["targets"][target_key] = Dict{String, Any}()
                    RA_Info["capacity_credit_result"]["targets"][target_key]["target_PLANT_NAME"] = target_plant_name
                    RA_Info["capacity_credit_result"]["targets"][target_key]["target_UNIT_GROUP"] = target_unit_group

                    simulation_setting = Dict{String, Any}(
                        "ELCC_method" => RA_Info["setting"]["capacity_credit_RA_simulation_method"],
                        "reference_RA_metric" => RA_Info["setting"]["capacity_credit_reference_RA_metric"],
                        "reference_RA_metric_spatial_resolution" => RA_Info["setting"]["capacity_credit_reference_RA_metric_spatial_resolution"],
                        "deactivated_asset_plant_name" => target_plant_name,
                        "deactivated_asset_unit_group" => target_unit_group
                    )

                    RA_Info["capacity_credit_result"]["targets"][target_key]["ELCC"] =
                        ALEAF.perform_ELCC_analysis_existing_asset_ED(
                            ALEAF_model_instance,
                            RA_Info,
                            simulation_setting,
                            risk_scenario_dict
                        )

                end
             
            else    # add_new resource mode

                if isempty(RA_Info["setting"]["ELCC_resource_list"])
                    @aleaf_warn "[ALEAF RA Model]: ELCC_resource_list is empty for add_new mode; skipping ELCC runs."
                end

                RA_Info["capacity_credit_result"] = Dict{String, Any}()
                RA_Info["capacity_credit_result"]["capacity_credit_type"] = "ELCC"
                RA_Info["capacity_credit_result"]["capacity_credit_assessment_mode"] = elcc_assessment_mode
                RA_Info["capacity_credit_result"]["ELCC_method"] = RA_Info["setting"]["capacity_credit_RA_simulation_method"]
                RA_Info["capacity_credit_result"]["reference_RA_metric"] = RA_Info["setting"]["capacity_credit_reference_RA_metric"]
                RA_Info["capacity_credit_result"]["reference_RA_metric_spatial_resolution"] = RA_Info["setting"]["capacity_credit_reference_RA_metric_spatial_resolution"]
                RA_Info["capacity_credit_result"]["targets"] = Dict{String, Any}()

                for resource_id in RA_Info["setting"]["ELCC_resource_list"]

                    for bus_idx in keys(ALEAF_model_instance.ref[:nw][0][:bus])

                        if ALEAF_model_instance.ref[:nw][0][:bus][bus_idx]["RA_ELCC_Calculation_Flag"] == true

                            target_key = string(resource_id, "|", string(bus_idx))
                            
                            RA_Info["capacity_credit_result"]["targets"][target_key] = Dict{String, Any}()
                            RA_Info["capacity_credit_result"]["targets"][target_key]["target_UNIT_ID"] = resource_id
                            RA_Info["capacity_credit_result"]["targets"][target_key]["target_UNIT_Location"] = bus_idx

                            simulation_setting = Dict{String, Any}(
                                "new_resource_unit_group" => resource_id,
                                "reference_RA_metric" => RA_Info["setting"]["capacity_credit_reference_RA_metric"],
                                "reference_RA_metric_spatial_resolution" => RA_Info["setting"]["capacity_credit_reference_RA_metric_spatial_resolution"],
                                "new_resource_location" => bus_idx                                
                            )
                            RA_Info["capacity_credit_result"]["targets"][target_key]["ELCC"] = ALEAF.perform_ELCC_analysis_ED(ALEAF_model_instance, RA_Info, simulation_setting, risk_scenario_dict)
                            
                        end
                    end
        
                end
            end
        end
    end

    delete!(RA_Info, "RA_risk_data")
    output_verbose_level = lowercase(string(get(RA_Info["setting"], "output_verbose_level", "compact")))
    if output_verbose_level == "detailed"
        if haskey(RA_Info, "RA_reference_risk_data")
            filtered = get(RA_Info["RA_reference_risk_data"], "filtered_joint_scenario_map", Dict{String, Any}())
            RA_Info["RA_reference_risk_data"] = Dict{String, Any}(
                "filtered_joint_scenario_map" => filtered
            )
        end
    else
        delete!(RA_Info, "RA_solutions")
        delete!(RA_Info, "RA_reference_risk_data")
        delete!(RA_Info, "renewable_scenario_data")
        delete!(RA_Info, "renewable_scenario_daygroup_data")        
    end

    # One result per planning stage (multi-round runs call RA once per round; a single file would be overwritten)
    ra_stage = current_year
    pd_ra = ALEAF_setting["Planning Design"]
    ra_year = pd_ra["first_stage_year_value"] + (ra_stage - 1) * pd_ra["num_years_per_stage_value"]
    RA_Info["stage"] = ra_stage
    RA_Info["year"] = ra_year

    stringdata = JSON.json(RA_Info);
    file_name = string(RA_Info["scenario"], "__RA_result_stage_", ra_stage, ".json")

    open(string(ALEAF_model_instance.ref[:nw][0][:output_path], file_name), "w") do f
        write(f, stringdata)
    end

    write_RA_metrics_csv(string(ALEAF_model_instance.ref[:nw][0][:output_path], RA_Info["scenario"], "__RA_metrics_stage_", ra_stage, ".csv"),
                         RA_Info, ra_stage, ra_year)

    # Optional post-run cleanup: remove temporary RA cache directories once
    # all analyses and result export are complete.
    cleanup_cache_flag = get(RA_Info["setting"], "cleanup_RA_cache_after_run_flag", true)
    if cleanup_cache_flag
        cache_paths = String[]
        if haskey(RA_Info["setting"], "outage_cache_path")
            push!(cache_paths, string(RA_Info["setting"]["outage_cache_path"]))
        end
        if haskey(RA_Info["setting"], "reference_cache_path")
            push!(cache_paths, string(RA_Info["setting"]["reference_cache_path"]))
        end
        for cache_path in unique(cache_paths)
            if ispath(cache_path)
                try
                    rm(cache_path; recursive=true, force=true)                    
                catch err
                    @aleaf_warn "[ALEAF RA Model]: Failed to remove RA cache path: $cache_path | $(sprint(showerror, err))"
                end
            end
        end
        @aleaf_info "[ALEAF RA Model]: Removed RA cache"
    end

    @aleaf_info "[ALEAF RA Model]: Completed"

    return RA_Info

end


# Long-format RA metrics table: one row per (metric, scope, day group, region).
const RA_METRIC_UNITS = Dict("EUE" => "MWh", "NEUE" => "ppm", "LOLH" => "hours/year", "LOLE" => "days/year",
    "Max_Consecutive_Outage_Hours" => "hours", "Max_MW_Loss" => "MW", "Max_MWh_Loss" => "MWh")

function write_RA_metrics_csv(file_path::AbstractString, RA_Info::Dict, stage::Integer, year::Integer)
    rows = NamedTuple[]
    scenario = string(RA_Info["scenario"])
    add!(metric, unit, scope, day_group, region, value) = push!(rows,
        (Scenario=scenario, Stage=stage, Year=year, Metric=metric, Unit=unit, Scope=scope, Day_Group=day_group, Region=region, Value=value))

    for (metric, unit) in RA_METRIC_UNITS
        haskey(RA_Info["RA_metrics"], metric) || continue
        m = RA_Info["RA_metrics"][metric]
        if haskey(m, "Day Group")
            for (dg, v) in sort(collect(m["Day Group"]); by = x -> parse(Int, string(x[1])))
                add!(metric, unit, "Day Group", parse(Int, string(dg)), "Systemwide", v["Systemwide"])
                for (reg, val) in sort(collect(v["Regional"]); by = x -> string(x[1]))
                    add!(metric, unit, "Day Group", parse(Int, string(dg)), string(reg), val)
                end
            end
        end
        if haskey(m, "Annual")
            a = m["Annual"]
            add!(metric, unit, "Annual", missing, "Systemwide", a["Systemwide"])
            for (reg, val) in sort(collect(a["Regional"]); by = x -> string(x[1]))
                add!(metric, unit, "Annual", missing, string(reg), val)
            end
        end
    end

    # Sample sizes and stress level, so the metrics can be judged next to how they were produced
    demand = get(get(RA_Info, "system_info", Dict{String,Any}()), "demand", Dict{String,Any}())
    sampled = get(demand, "num_risk_scenario_setting", missing)
    for (dg, n_solved) in sort(collect(get(demand, "solved_scenarios_by_daygroup", Dict{String,Any}())); by = x -> parse(Int, string(x[1])))
        add!("Scenarios_Sampled", "count", "Day Group", parse(Int, string(dg)), "Systemwide", sampled)
        add!("Scenarios_Solved", "count", "Day Group", parse(Int, string(dg)), "Systemwide", n_solved)
    end
    add!("System_Peak_Scale", "multiplier", "Setting", missing, "Systemwide", get(RA_Info["setting"], "system_peak_scale", missing))

    CSV.write(file_path, DataFrame(rows))
end


function update_expansion_results!(ALEAF_model_instance, recorded_investment_decisions)
    
    for i in [(i) for (i) in get_index(ALEAF_model_instance, :gen_index, 0)]
        bus_idx = parameter(ALEAF_model_instance, 0, :gen_index, "bus_idx", i)
        tech_idx = parameter(ALEAF_model_instance, 0, :gen_index, "genco_tech_id", i)    
        EXUNITS = parameter(ALEAF_model_instance, bus_idx, :gen_bus, tech_idx, "EXUNITS")
        CAP = parameter(ALEAF_model_instance, bus_idx, :gen_bus, tech_idx, "CAP")

        if length(recorded_investment_decisions) > 0    
            ALEAF_model_instance.ref[:nw][bus_idx][:gen_bus][tech_idx]["U_G_iy"] = recorded_investment_decisions["U_G_i"][string(i)] 
            ALEAF_model_instance.ref[:nw][bus_idx][:gen_bus][tech_idx]["Original_EXUNITS"] = ALEAF_model_instance.ref[:nw][bus_idx][:gen_bus][tech_idx]["EXUNITS"]
            ALEAF_model_instance.ref[:nw][bus_idx][:gen_bus][tech_idx]["EXUNITS"] = recorded_investment_decisions["U_G_i"][string(i)] 
            if parameter(ALEAF_model_instance, 0, :gen_index, "UNIT_CATEGORY", i) == "STORAGE"
                # Storage units without a u_ESE_iy decision (e.g. hybrid storage with fixed energy)
                # are absent from U_ESE_i — fall back to the network's ES_MWh in that case.
                network_ES_MWh = parameter(ALEAF_model_instance, bus_idx, :gen_bus, tech_idx, "ES_MWh")
                ES_MWh_value = get(recorded_investment_decisions["U_ESE_i"], string(i), network_ES_MWh)
                ALEAF_model_instance.ref[:nw][bus_idx][:gen_bus][tech_idx]["ES_MWh"] = ES_MWh_value
                ALEAF_model_instance.ref[:nw][bus_idx][:gen_bus][tech_idx]["Existing_ES_MWh"] = ES_MWh_value
            end
        else
            ALEAF_model_instance.ref[:nw][bus_idx][:gen_bus][tech_idx]["U_G_iy"] = EXUNITS
            if parameter(ALEAF_model_instance, 0, :gen_index, "UNIT_CATEGORY", i) == "STORAGE"
                ALEAF_model_instance.ref[:nw][bus_idx][:gen_bus][tech_idx]["ES_MWh"] = parameter(ALEAF_model_instance, bus_idx, :gen_bus, tech_idx, "ES_MWh")
                ALEAF_model_instance.ref[:nw][bus_idx][:gen_bus][tech_idx]["Existing_ES_MWh"] = parameter(ALEAF_model_instance, bus_idx, :gen_bus, tech_idx, "ES_MWh")
            end
        end
    end

    if length(recorded_investment_decisions) > 0
        hybrid = hybrid_tx_enabled(ALEAF_model_instance)
        for k in [(k) for (k) in get_index(ALEAF_model_instance, :branch, 0) if parameter(ALEAF_model_instance, 0, :branch, "model_flag", k) == true]
            rate_a = parameter(ALEAF_model_instance, 0, :branch, "rate_a", k)
            u_new_T = recorded_investment_decisions["U_NEW_T_k"][string(k)]
            ALEAF_model_instance.ref[:nw][0][:branch][k]["rate_a"] = rate_a * (1 + u_new_T)
            # Enhanced-hybrid: fold susceptance for AC branches so the B-theta base flow sees
            # the expanded admittance b0·(1+u*), consistent with the scaled rate_a. Idempotent
            # via a stored base to guard against a repeated ref update.
            if hybrid && get(ALEAF_model_instance.ref[:nw][0][:branch][k], "dc_line", false) != true
                br = ALEAF_model_instance.ref[:nw][0][:branch][k]
                haskey(br, "br_x_pu_hybrid_base") || (br["br_x_pu_hybrid_base"] = br["br_x_pu"])
                br["br_x_pu"] = br["br_x_pu_hybrid_base"] / (1 + u_new_T)
            end
        end
    end

    # Diagnostic: report post-update state for all STORAGE units to confirm expansion record application.
    # CAP / Charge_CAP come from the network file (not overwritten); EXUNITS and ES_MWh are restored from the expansion record.
    for i in [(i) for (i) in get_index(ALEAF_model_instance, :gen_index, 0)]
        if parameter(ALEAF_model_instance, 0, :gen_index, "UNIT_CATEGORY", i) == "STORAGE"
            bus_idx  = parameter(ALEAF_model_instance, 0, :gen_index, "bus_idx", i)
            tech_idx = parameter(ALEAF_model_instance, 0, :gen_index, "genco_tech_id", i)
            gb = ALEAF_model_instance.ref[:nw][bus_idx][:gen_bus][tech_idx]
            plant_name = get(gb, "PLANT_NAME", string(tech_idx))
            @aleaf_info "[ALEAF RA update_expansion_results!]: STORAGE i=$i PLANT_NAME=$plant_name " *
                "EXUNITS=$(gb["EXUNITS"]) CAP=$(gb["CAP"]) Charge_CAP=$(gb["Charge_CAP"]) " *
                "ES_MWh=$(gb["ES_MWh"]) Existing_ES_MWh=$(gb["Existing_ES_MWh"])"
        end
    end
end


"""
    couple_hybrid_outages!(risk_scenario_dict, ref)

Treat each hybrid asset as one unit for outages: for every co-located
(renewable `r`, storage `b`) pair, mark both down at any hour where either is down
(per-hour elementwise AND on the availability BitMatrices). This keeps a downed hybrid
battery from being force-charged by the reference onsite-generation term in the ED SOC
balance, and forces the paired renewable off with it. Idempotent; mutates in place, so it
must run once (single-threaded) before scenario filtering and the cache write.
"""
function couple_hybrid_outages!(risk_scenario_dict, ref)
    gen_index = ref[:nw][0][:gen_index]
    hybrid_pairs = Tuple{Int,Int}[]
    for (i, meta) in gen_index
        (i isa Integer && meta isa AbstractDict) || continue
        if get(meta, "hybrid_type", "NA") == "ES"
            push!(hybrid_pairs, (Int(meta["hybrid_main_gen_idx"]), Int(i)))
        end
    end
    isempty(hybrid_pairs) && return risk_scenario_dict

    for day_dict in values(risk_scenario_dict)
        day_dict isa AbstractDict || continue
        for mat in values(day_dict)
            mat isa AbstractMatrix || continue
            ncol = size(mat, 2)
            for (r, b) in hybrid_pairs
                (1 <= r <= ncol && 1 <= b <= ncol) || continue
                both = mat[:, r] .& mat[:, b]
                mat[:, r] .= both
                mat[:, b] .= both
            end
        end
    end
    return risk_scenario_dict
end


function generate_and_filter_risk_scenarios(ALEAF_setting, RA_Info, ALEAF_model_instance; RA_input::Dict{String,<:Any} = Dict{String,Any}())

    # read temperature data
    data_location = ALEAF_setting["data_location"]
    file_path = joinpath(data_location, ALEAF_setting["File Path"]["timeseries_data_temperature_path"])
    if ALEAF_model_instance.setting["Simulation Configuration"]["Temperature_File_ID"] != "Base" 
        data_path = joinpath(data_location, "timeseries_data_files", "0_additional_scenarios", "Weather")
        file_name = string("timeseries_temp_hourly_", ALEAF_model_instance.setting["Simulation Configuration"]["Temperature_File_ID"], ".csv")
        file_path = joinpath(data_path, file_name)
    end
    ALEAF_model_instance.ref[:nw][0][:temperatureData] = DataFrame(CSV.File(file_path))

    # read RA data
    ALEAF_model_instance.ref[:nw][0][:regressionAD] = DataFrame(XLSX.readtable(joinpath(data_location, "RA_data.xlsx"), "RegressionParametersAD"))
    ALEAF_model_instance.ref[:nw][0][:regressionDD] = DataFrame(XLSX.readtable(joinpath(data_location, "RA_data.xlsx"), "RegressionParametersDD"))

    # reference temp
    reference_temp = RA_Info["setting"]["reference_temp"]

    # day groups
    simulation_list = [i for i in 1:ALEAF_model_instance.setting["Simulation Configuration"]["NDAY_Groups_RA"]]
    num_simulations = ALEAF_model_instance.setting["Simulation Configuration"]["NDAY_Groups_RA"]
    build_ra_renewable_daygroup_cache!(RA_Info, ALEAF_model_instance)   # map renewable scenario timeseries to daygroups
    build_ra_joint_scenario_map!(RA_Info, simulation_list)

    # get risk scenarios
    if haskey(RA_input, "risk_scenario_dict")
        risk_scenario_dict = RA_input["risk_scenario_dict"]        
        @aleaf_info "[ALEAF RA Model]: Using pre-determined risk sample data"
    else
        risk_scenario_dict = generate_outage_samples(ALEAF_model_instance.ref, ALEAF_model_instance.setting, simulation_list, reference_temp, RA_Info["setting"]["num_risk_scenario"], RA_Info["setting"]; output_option="Status")
    end

    # Hybrid assets are one unit for outages: either component down takes the whole asset down.
    # Must run before filtering and caching so every downstream consumer sees the coupled matrices.
    couple_hybrid_outages!(risk_scenario_dict, ALEAF_model_instance.ref)

    # filter risk scenarios
    if haskey(RA_input, "risk_scenario_dict")
        if haskey(RA_input, "filtered_joint_scenario_map")
            filtered_joint_scenario_map = RA_input["filtered_joint_scenario_map"]
        else
            error("RA_input includes risk_scenario_dict but missing filtered_joint_scenario_map.")
        end
        RA_Info["RA_reference_risk_data"]["filtered_joint_scenario_map"] = filtered_joint_scenario_map
        @aleaf_info "[ALEAF RA Model]: Using pre-determined risk filtering data"
    else
        filtered_joint_scenario_map = filter_generated_outage_scenarios!(RA_Info, ALEAF_model_instance.ref, ALEAF_model_instance.setting, simulation_list, RA_Info["setting"], risk_scenario_dict)
    end

    return risk_scenario_dict, filtered_joint_scenario_map

end

function filter_generated_outage_scenarios!(RA_Info, am_reference_data, am_setting_data, simulation_list, RA_setting, risk_scenario_dict)
    filtered_joint_scenario_list = Dict{Int64, Any}()
    filtered_joint_scenario_map = Dict{Int64, Dict{Int, Dict{String, Any}}}()
    for day_group_id in keys(risk_scenario_dict)
        filtered_joint_scenario_list[day_group_id] = []
        filtered_joint_scenario_map[day_group_id] = Dict{Int, Dict{String, Any}}()
    end

    for day_group_idx in eachindex(simulation_list)
        day_group_id = simulation_list[day_group_idx]

        if RA_setting["risk_filtering_flag"] == false
            if length(RA_setting["preselected_days_list"]) > 0
                if day_group_id in RA_setting["preselected_days_list"]
                    filtered_joint_scenario_list[day_group_id] = collect(keys(RA_Info["RA_reference_risk_data"]["joint_scenario_map"][day_group_id]))
                end
            else
                filtered_joint_scenario_list[day_group_id] = collect(keys(RA_Info["RA_reference_risk_data"]["joint_scenario_map"][day_group_id]))
            end
        else
            if length(RA_setting["preselected_days_list"]) > 0
                if day_group_id in RA_setting["preselected_days_list"]
                    filtered_joint_scenario_list[day_group_id] = filter_joint_risk_scenarios(RA_Info, am_reference_data, am_setting_data, day_group_id, RA_setting, risk_scenario_dict[day_group_id], RA_Info["RA_reference_risk_data"]["joint_scenario_map"][day_group_id])
                end
            else
                filtered_joint_scenario_list[day_group_id] = filter_joint_risk_scenarios(RA_Info, am_reference_data, am_setting_data, day_group_id, RA_setting, risk_scenario_dict[day_group_id], RA_Info["RA_reference_risk_data"]["joint_scenario_map"][day_group_id])
            end
        end
        for joint_id in filtered_joint_scenario_list[day_group_id]
            filtered_joint_scenario_map[day_group_id][joint_id] = deepcopy(RA_Info["RA_reference_risk_data"]["joint_scenario_map"][day_group_id][joint_id])
        end
    end

    RA_Info["RA_reference_risk_data"]["filtered_joint_scenario_list"] = filtered_joint_scenario_list
    RA_Info["RA_reference_risk_data"]["filtered_joint_scenario_map"] = filtered_joint_scenario_map

    return filtered_joint_scenario_map
end


function filter_risk_scenarios(am_reference_data, am_setting_data, day_group_id, RA_setting, risk_scenario_dict)

    # Total gen capacity
    num_gen = length(am_reference_data[:nw][0][:RA_gen_index])
    num_workers= length(workers())

    ids_d = [(d) for (d) in am_reference_data[:nw][0][:repday_groups][day_group_id]["Day_Idx_List"]]
    ids_h = [(h) for (h) in am_setting_data["run_H"]]
    ids_dh = [(d,h) for d in ids_d for h in ids_h]
    ids_n = [(k) for (k) in keys(am_reference_data[:nw][0][:bus])]
    num_hours = length(ids_h)*length(ids_d)
    y = RA_setting["current_year"]

    # Demand
    Demand = zeros(num_hours)
    idx = 1
    for d in ids_d
        for h in ids_h
            for n in ids_n
                Demand[idx] += get_bus_demand_with_growth(am_reference_data[:nw], n, d, h, 1, y; nw=0, system_peak_scale=RA_setting["system_peak_scale"])
            end
            idx += 1
        end
    end

    # Rank hours by demand (highest = rank 1).
    sorted_demand = sortperm(Demand; rev=true)

    ranks = zeros(Int, length(sorted_demand))
    for (rank, idx) in enumerate(sorted_demand)
        ranks[idx] = rank
    end

    battery_discharge_vector_dict = Dict()
    for tech_id in keys(am_reference_data[:nw][0][:gen_technology]) 
        storage_hour = am_reference_data[:nw][0][:gen_technology][tech_id]["STOHR_MAX"]
        if storage_hour > 0
            if !haskey(battery_discharge_vector_dict, storage_hour)
                battery_discharge_vector_dict[storage_hour] = zeros(length(sorted_demand))
                battery_discharge_vector_dict[storage_hour][ranks .<= storage_hour] .= 1.0  # discharge in the storage_hour peak-demand hours
            end
        end
    end
    
    # Gen Firm Capacity
    gen_firm_capacity = zeros(num_hours, num_gen)

    for i in 1:num_gen
        
        bus_idx = am_reference_data[:nw][0][:RA_gen_index][i]["bus_idx"]
        tech_idx = am_reference_data[:nw][0][:RA_gen_index][i]["genco_tech_id"]
        unit_category = am_reference_data[:nw][0][:RA_gen_index][i]["UNIT_CATEGORY"]
        vre_flag = am_reference_data[:nw][bus_idx][:RA_gen_bus][tech_idx]["VRE_Flag"]
        fuel_limit = am_reference_data[:nw][bus_idx][:RA_gen_bus][tech_idx]["FUEL_LIMIT"]
        ICAP = am_reference_data[:nw][0][:RA_gen_index][i]["ICAP"]

        if (vre_flag == true) && (fuel_limit == "Fixed Profile")
            # VRE shape (available shapes: wind_ons, wind_ofs, pv, rtpv, hydro, csp)
            profile_type = am_reference_data[:nw][0][:RA_gen_index][i]["Profile_Type"]
            Timeseries_Tag =  am_reference_data[:nw][bus_idx][:RA_gen_bus][tech_idx]["Timeseries_Tag"]

            vre_firm_output = zeros(num_hours)  # Local to each thread
            for (idx, (d, h)) in enumerate(ids_dh)
                shape = 0.0
                if Timeseries_Tag == "LOCAL"
                    type_key = profile_type == "NA" ? "NA" : string(profile_type, "_shape")

                    shape = get_vre_zdt_shape(am_reference_data[:nw][0], y, bus_idx, d, h, type_key)

                else
                    if profile_type in ["wind_ons", "wind_ofs", "csp"]
                        shape = am_reference_data[:nw][0][:planning_stages][y]["repdays"][string(d)]["data"][string(h)]["1"][profile_type][Timeseries_Tag]
                    end
                end

                vre_firm_output[idx] = shape * ICAP
            end

            gen_firm_capacity[:, i] = vre_firm_output

        elseif unit_category in ["STORAGE"]

            storage_duration = am_reference_data[:nw][bus_idx][:RA_gen_bus][tech_idx]["STOHR_MAX"]
            if storage_duration > 0
                gen_firm_capacity[:, i] = ICAP .* battery_discharge_vector_dict[storage_duration]
            else
                plant_name = am_reference_data[:nw][bus_idx][:RA_gen_bus][tech_idx]["PLANT_NAME"]
                @aleaf_error "[ALEAF RA Model]: Storage unit with non-positive storage duration found. Please check data for bus_idx=$bus_idx, tech_idx=$tech_idx, plant_name=$plant_name."
                terminate_with_error()
            end
        
        else
            gen_firm_capacity[:, i] = [ICAP for t in 1:num_hours]
        end
    end

    tol = RA_setting["risk_tol"]
    tol_total_gen_loss = 000.0
    tol_max_hourly_gen_loss = 500.0

    num_risk_scenarios = length(risk_scenario_dict)
    risk_indicator_list = falses(num_risk_scenarios)
    risk_list = zeros(num_risk_scenarios)
    selected_risk_data = []
    risk_key_list = collect(keys(risk_scenario_dict))
    gen_margin_vector = (sum(gen_firm_capacity, dims=2) .- Demand)
    gen_margin_computed = any(x -> x <= 0, gen_margin_vector) ? max.(gen_margin_vector, 100.0) : gen_margin_vector

    ################### FILTER scenarios
    Threads.@threads for scenario_id in 1:num_risk_scenarios
        risk_id = risk_key_list[scenario_id]

        # Screening metric: max hourly (gen loss MW / gen margin). Alternatives tried:
        # loss/firm-cap, loss/demand, margin/demand.
        risk_list[scenario_id] = maximum(sum(.!risk_scenario_dict[risk_id] .* gen_firm_capacity, dims=2) ./ gen_margin_computed)

        risk_list[scenario_id] > tol ? risk_indicator_list[scenario_id] = true : nothing
    end

    # Enforce a per-day-group floor on retained scenarios by rank when the tol filter keeps too few.
    if RA_setting["min_num_risk_in_each_day_value"] > 0
        min_num_of_risk_scenarios = max(10, RA_setting["min_num_risk_in_each_day_value"] * length(risk_scenario_dict))
        if count(risk_indicator_list) < min_num_of_risk_scenarios

            sorted_indices = sortperm(risk_list)
            risk_ranks = zeros(Int, length(sorted_indices))

            for (rank, idx) in enumerate(sorted_indices)
                risk_ranks[idx] = rank
            end

            for scenario_id in 1:num_risk_scenarios
                risk_ranks[scenario_id] <= min_num_of_risk_scenarios ? risk_indicator_list[scenario_id] = true : nothing
            end

        end
    end

    for scenario_id in 1:num_risk_scenarios
        if risk_indicator_list[scenario_id] == true
            push!(selected_risk_data, risk_key_list[scenario_id])
        end
    end

    return selected_risk_data

end

function filter_joint_risk_scenarios(RA_Info, am_reference_data, am_setting_data, day_group_id, RA_setting, risk_scenario_dict, joint_scenario_map)

    num_gen = length(am_reference_data[:nw][0][:RA_gen_index])

    ids_d = [(d) for (d) in am_reference_data[:nw][0][:repday_groups][day_group_id]["Day_Idx_List"]]
    ids_h = [(h) for (h) in am_setting_data["run_H"]]
    ids_dh = [(d,h) for d in ids_d for h in ids_h]
    ids_n = [(k) for (k) in keys(am_reference_data[:nw][0][:bus])]
    num_hours = length(ids_h) * length(ids_d)
    y = RA_setting["current_year"]

    Demand = zeros(num_hours)
    idx = 1
    for d in ids_d
        for h in ids_h
            for n in ids_n
                Demand[idx] += get_bus_demand_with_growth(am_reference_data[:nw], n, d, h, 1, y; nw=0, system_peak_scale=RA_setting["system_peak_scale"])
            end
            idx += 1
        end
    end

    sorted_demand = sortperm(Demand; rev=true)
    ranks = zeros(Int, length(sorted_demand))
    for (rank, idx) in enumerate(sorted_demand)
        ranks[idx] = rank
    end

    battery_discharge_vector_dict = Dict()
    for tech_id in keys(am_reference_data[:nw][0][:gen_technology])
        storage_hour = am_reference_data[:nw][0][:gen_technology][tech_id]["STOHR_MAX"]
        if storage_hour > 0
            if !haskey(battery_discharge_vector_dict, storage_hour)
                battery_discharge_vector_dict[storage_hour] = zeros(length(sorted_demand))
                battery_discharge_vector_dict[storage_hour][ranks .<= storage_hour] .= 1.0
            end
        end
    end

    constant_gen_capacity = zeros(num_hours, num_gen)
    vre_unit_meta = Dict{Int, Dict{String, Any}}()

    for i in 1:num_gen
        bus_idx = am_reference_data[:nw][0][:RA_gen_index][i]["bus_idx"]
        tech_idx = am_reference_data[:nw][0][:RA_gen_index][i]["genco_tech_id"]
        unit_category = am_reference_data[:nw][0][:RA_gen_index][i]["UNIT_CATEGORY"]
        vre_flag = am_reference_data[:nw][bus_idx][:RA_gen_bus][tech_idx]["VRE_Flag"]
        fuel_limit = am_reference_data[:nw][bus_idx][:RA_gen_bus][tech_idx]["FUEL_LIMIT"]
        ICAP = am_reference_data[:nw][0][:RA_gen_index][i]["ICAP"]

        if (vre_flag == true) && (fuel_limit == "Fixed Profile")
            vre_unit_meta[i] = Dict(
                "bus_idx" => bus_idx,
                "profile_type" => am_reference_data[:nw][0][:RA_gen_index][i]["Profile_Type"],
                "timeseries_tag" => am_reference_data[:nw][bus_idx][:RA_gen_bus][tech_idx]["Timeseries_Tag"],
                "icap" => ICAP,
            )
        elseif unit_category in ["STORAGE"]
            storage_duration = am_reference_data[:nw][bus_idx][:RA_gen_bus][tech_idx]["STOHR_MAX"]
            if storage_duration > 0
                constant_gen_capacity[:, i] = ICAP .* battery_discharge_vector_dict[storage_duration]
            else
                plant_name = am_reference_data[:nw][bus_idx][:RA_gen_bus][tech_idx]["PLANT_NAME"]
                @aleaf_error "[ALEAF RA Model]: Storage unit with non-positive storage duration found. Please check data for bus_idx=$bus_idx, tech_idx=$tech_idx, plant_name=$plant_name."
                terminate_with_error()
            end
        else
            constant_gen_capacity[:, i] = [ICAP for t in 1:num_hours]
        end
    end

    function get_vre_output_for_joint_scenario(unit_idx, renewable_scenario_id)
        meta = vre_unit_meta[unit_idx]
        bus_idx = meta["bus_idx"]
        profile_type = meta["profile_type"]
        timeseries_tag = meta["timeseries_tag"]
        ICAP = meta["icap"]

        if timeseries_tag == "LOCAL"
            type_key = profile_type == "NA" ? "NA" : string(profile_type, "_shape")
            vre_firm_output = zeros(num_hours)
            for (idx_local, (d, h)) in enumerate(ids_dh)
                vre_firm_output[idx_local] = get_vre_zdt_shape(am_reference_data[:nw][0], y, bus_idx, d, h, type_key) * ICAP
            end
            return vre_firm_output
        end

        # Renewable-scenario perturbations only carry wind_ons/pv classes.
        if profile_type in ["wind_ons", "pv"]
            if haskey(RA_Info["renewable_scenario_daygroup_data"][day_group_id][renewable_scenario_id][profile_type], timeseries_tag)
                return RA_Info["renewable_scenario_daygroup_data"][day_group_id][renewable_scenario_id][profile_type][timeseries_tag] .* ICAP
            end
        end

        vre_firm_output = zeros(num_hours)
        if profile_type in ["wind_ons", "wind_ofs", "csp"]
            for (idx_local, (d, h)) in enumerate(ids_dh)
                vre_firm_output[idx_local] = am_reference_data[:nw][0][:planning_stages][y]["repdays"][string(d)]["data"][string(h)]["1"][profile_type][timeseries_tag] * ICAP
            end
        end
        return vre_firm_output
    end

    tol = RA_setting["risk_tol"]
    num_joint_scenarios = length(joint_scenario_map)
    risk_indicator_list = falses(num_joint_scenarios)
    risk_list = zeros(num_joint_scenarios)
    selected_joint_data = Int[]
    joint_key_list = collect(keys(joint_scenario_map))

    Threads.@threads for scenario_idx in 1:num_joint_scenarios
        joint_id = joint_key_list[scenario_idx]
        outage_risk_id = joint_scenario_map[joint_id]["outage_risk_id"]
        renewable_scenario_id = joint_scenario_map[joint_id]["renewable_scenario_id"]

        gen_firm_capacity = copy(constant_gen_capacity)
        for unit_idx in keys(vre_unit_meta)
            gen_firm_capacity[:, unit_idx] = get_vre_output_for_joint_scenario(unit_idx, renewable_scenario_id)
        end

        gen_margin_vector = vec(sum(gen_firm_capacity, dims=2) .- Demand)
        gen_margin_computed = any(x -> x <= 0, gen_margin_vector) ? max.(gen_margin_vector, 100.0) : gen_margin_vector
        risk_list[scenario_idx] = maximum(sum(.!risk_scenario_dict[outage_risk_id] .* gen_firm_capacity, dims=2) ./ gen_margin_computed)
        risk_list[scenario_idx] > tol ? risk_indicator_list[scenario_idx] = true : nothing
    end

    if RA_setting["min_num_risk_in_each_day_value"] > 0
        min_num_of_risk_scenarios = max(10, RA_setting["min_num_risk_in_each_day_value"] * length(joint_scenario_map))
        if count(risk_indicator_list) < min_num_of_risk_scenarios
            sorted_indices = sortperm(risk_list)
            risk_ranks = zeros(Int, length(sorted_indices))
            for (rank, idx) in enumerate(sorted_indices)
                risk_ranks[idx] = rank
            end
            for scenario_idx in 1:num_joint_scenarios
                risk_ranks[scenario_idx] <= min_num_of_risk_scenarios ? risk_indicator_list[scenario_idx] = true : nothing
            end
        end
    end

    for scenario_idx in 1:num_joint_scenarios
        if risk_indicator_list[scenario_idx] == true
            push!(selected_joint_data, joint_key_list[scenario_idx])
        end
    end

    return selected_joint_data
end

function parse_ra_renewable_scenarios!(ALEAF_setting, RA_Info)
    raw_scenarios = ALEAF_setting["RA Scenarios"]

    parsed = Dict{String, Dict{String, Any}}()
    weights = Dict{String, Float64}()

    for row_data in values(raw_scenarios)
        if row_data["Enabled"] != true
            continue
        end

        scenario_id = string(row_data["Scenario_ID"])
        weight = Float64(row_data["Weight"])

        parsed[scenario_id] = Dict{String, Any}(
            "weight" => weight,
            "wind_ons_file_id" => string(row_data["Wind_Ons_File_ID"]),
            "pv_file_id" => string(row_data["PV_File_ID"]),
        )
        weights[scenario_id] = weight
    end

    total_weight = sum(values(weights))
    for scenario_id in keys(parsed)
        normalized_weight = weights[scenario_id] / total_weight
        parsed[scenario_id]["weight"] = normalized_weight
        weights[scenario_id] = normalized_weight
    end

    renewable_scenario_ids = sort(collect(keys(parsed)))
    base_scenario_id = findfirst(x -> uppercase(x) == "BASE", renewable_scenario_ids)
    if base_scenario_id !== nothing
        base_scenario = renewable_scenario_ids[base_scenario_id]
        deleteat!(renewable_scenario_ids, base_scenario_id)
        renewable_scenario_ids = vcat([base_scenario], renewable_scenario_ids)
    end

    RA_Info["setting"]["renewable_scenarios"] = parsed
    RA_Info["setting"]["renewable_scenario_ids"] = renewable_scenario_ids
    RA_Info["setting"]["renewable_scenario_weights"] = weights
    RA_Info["setting"]["num_renewable_scenarios"] = length(renewable_scenario_ids)

    return nothing
end

function _get_ra_wind_ons_scenario_file_path(ALEAF_setting, scenario_id, wind_ons_file_id)
    
    data_location = ALEAF_setting["data_location"] 

    if uppercase(scenario_id) == "BASE"
        return joinpath(data_location, ALEAF_setting["File Path"]["timeseries_data_wind_ons_path"])
    end

    data_path = joinpath(data_location, "timeseries_data_files", "0_additional_scenarios", "WIND")
    file_name = "timeseries_wind_ons_hourly_$(wind_ons_file_id).csv"
    return joinpath(data_path, file_name)
end

function _get_ra_pv_scenario_file_path(ALEAF_setting, scenario_id, pv_file_id)
    data_location = ALEAF_setting["data_location"]

    if uppercase(scenario_id) == "BASE"
        return joinpath(data_location, ALEAF_setting["File Path"]["timeseries_data_pv_path"])
    end

    data_path = joinpath(data_location, "timeseries_data_files", "0_additional_scenarios", "PV")
    file_name = "timeseries_pv_hourly_$(pv_file_id).csv"
    return joinpath(data_path, file_name)
end

function _read_ra_renewable_timeseries_by_tag(file_path)
    time_series_df = CSV.read(file_path, DataFrame)
    time_columns = Set(["Year", "Month", "Day", "Period"])
    time_series_by_tag = Dict{String, Vector{Float64}}()

    for col_name in names(time_series_df)
        if col_name in time_columns
            continue
        end
        time_series_by_tag[string(col_name)] = Float64.(time_series_df[!, col_name])
    end

    return time_series_by_tag
end

function load_ra_renewable_scenario_timeseries!(ALEAF_setting, RA_Info)
    renewable_scenarios = RA_Info["setting"]["renewable_scenarios"]
    renewable_scenario_ids = RA_Info["setting"]["renewable_scenario_ids"]

    renewable_scenario_data = Dict{String, Any}(
        "wind_ons" => Dict{String, Dict{String, Vector{Float64}}}(),
        "pv" => Dict{String, Dict{String, Vector{Float64}}}(),
    )

    for scenario_id in renewable_scenario_ids
        scenario_info = renewable_scenarios[scenario_id]

        wind_ons_path = _get_ra_wind_ons_scenario_file_path(
            ALEAF_setting,
            scenario_id,
            scenario_info["wind_ons_file_id"],
        )
        pv_path = _get_ra_pv_scenario_file_path(
            ALEAF_setting,
            scenario_id,
            scenario_info["pv_file_id"],
        )

        renewable_scenario_data["wind_ons"][scenario_id] = _read_ra_renewable_timeseries_by_tag(wind_ons_path)
        renewable_scenario_data["pv"][scenario_id] = _read_ra_renewable_timeseries_by_tag(pv_path)
    end

    RA_Info["renewable_scenario_data"] = renewable_scenario_data

    return nothing
end

function build_ra_renewable_daygroup_cache!(RA_Info, ALEAF_model_instance)
    renewable_scenario_data = RA_Info["renewable_scenario_data"]
    renewable_scenario_ids = RA_Info["setting"]["renewable_scenario_ids"]
    repday_groups = ALEAF_model_instance.ref[:nw][0][:repday_groups]
    repdays = ALEAF_model_instance.ref[:nw][0][:repdays]
    ids_h = collect(ALEAF_model_instance.setting["run_H"])

    renewable_scenario_daygroup_data = Dict{Int, Any}()

    for day_group_id in keys(repday_groups)
        ids_d = repday_groups[day_group_id]["Day_Idx_List"]
        renewable_scenario_daygroup_data[day_group_id] = Dict{String, Any}()

        for scenario_id in renewable_scenario_ids
            renewable_scenario_daygroup_data[day_group_id][scenario_id] = Dict{String, Any}()

            for tech in keys(renewable_scenario_data)
                daygroup_tech_data = Dict{String, Vector{Float64}}()

                for (timeseries_tag, annual_profile) in renewable_scenario_data[tech][scenario_id]
                    hourly_profile = Vector{Float64}(undef, length(ids_d) * length(ids_h))
                    idx = 1
                    for repday_id in ids_d
                        reference_day = repdays[repday_id]["Day"]
                        for h in ids_h
                            annual_hour_idx = (reference_day - 1) * 24 + h
                            hourly_profile[idx] = annual_profile[annual_hour_idx]
                            idx += 1
                        end
                    end
                    daygroup_tech_data[timeseries_tag] = hourly_profile
                end

                renewable_scenario_daygroup_data[day_group_id][scenario_id][tech] = daygroup_tech_data
            end
        end
    end

    RA_Info["renewable_scenario_daygroup_data"] = renewable_scenario_daygroup_data
    RA_Info["setting"]["renewable_scenario_daygroup_data"] = renewable_scenario_daygroup_data

    return nothing
end

function build_ra_joint_scenario_map!(RA_Info, simulation_list)
    renewable_scenario_ids = RA_Info["setting"]["renewable_scenario_ids"]
    num_risk_scenario = RA_Info["setting"]["num_risk_scenario"]
    joint_scenario_map = Dict{Int, Dict{Int, Dict{String, Any}}}()

    for day_group_id in simulation_list
        joint_scenario_map[day_group_id] = Dict{Int, Dict{String, Any}}()
        joint_id = 1

        for outage_risk_id in 1:num_risk_scenario
            for renewable_scenario_id in renewable_scenario_ids
                joint_scenario_map[day_group_id][joint_id] = Dict{String, Any}(
                    "outage_risk_id" => outage_risk_id,
                    "renewable_scenario_id" => renewable_scenario_id,
                )
                joint_id += 1
            end
        end
    end

    RA_Info["RA_reference_risk_data"]["joint_scenario_map"] = joint_scenario_map
    RA_Info["setting"]["joint_scenario_map"] = joint_scenario_map

    return nothing
end

function update_RA_setting(ALEAF_setting, case_id, RA_input, current_year)

    resolve_network_input_files!(ALEAF_setting, case_id)
    data_location = ALEAF_setting["data_location"]
    load_file_path_sheet!(ALEAF_setting, ALEAF_setting["network_data_file_location"])
    
    # RA assessment summary dictionary
    RA_Info = Dict{String, Any}()
    RA_Info["setting"] = Dict{String, Any}()
    
    # Propagate global logging mode into RA setting so RA internals honor detailed logging.
    global_logging_level = get(ALEAF_setting["Simulation Setting"], "logging_level_value", "simple")

    ra_setting_default = ALEAF_setting["RA Setting"]
    default_elcc_resource_list = String[]
    for tech_idx in keys(ALEAF_setting["Gen Technology"])
        if ALEAF_setting["Gen Technology"][tech_idx]["ELCC_Flag"] == true
            push!(default_elcc_resource_list, ALEAF_setting["Gen Technology"][tech_idx]["UNITGROUP"])
        end
    end

    if haskey(RA_input, "setting_type")
        RA_Info["setting"]["setting_type"] = RA_input["setting_type"]
        RA_Info["setting"]["num_risk_scenario"] = get(RA_input, "num_risk_scenario", ra_setting_default["num_risk_scenario"])
        RA_Info["setting"]["metrics"] = Dict{String, Any}()
        RA_Info["setting"]["RA_method"] = get(RA_input, "RA_method", ra_setting_default["RA_method"])
        RA_Info["setting"]["current_year"] = current_year
        RA_Info["setting"]["logging_level_value"] = get(RA_input, "logging_level_value", global_logging_level)

        RA_Info["setting"]["reference_temp"] = get(RA_input, "reference_temp", ra_setting_default["reference_temp"])
        RA_Info["setting"]["system_peak_scale"] = get(RA_input, "system_peak_scale", ra_setting_default["system_peak_scale"])
        RA_Info["setting"]["risk_tol"] = get(RA_input, "risk_tol", ra_setting_default["risk_tol_value"])
        RA_Info["setting"]["min_num_risk_in_each_day_value"] = get(RA_input, "min_num_risk_in_each_day_value", ra_setting_default["min_num_risk_in_each_day_value"])
        RA_Info["setting"]["repair_time_bound"] = get(RA_input, "repair_time_bound", ra_setting_default["repair_time_bound"])

        RA_Info["setting"]["ELCC_resource_list"] = get(RA_input, "ELCC_resource_list", default_elcc_resource_list)

        RA_Info["setting"]["distributed_run_flag"] = get(RA_input, "distributed_run_flag", ra_setting_default["distributed_run_flag"])
        RA_Info["setting"]["num_distributed_scenarios_per_worker_value"] = get(RA_input, "num_distributed_scenarios_per_worker_value", ra_setting_default["num_distributed_scenarios_per_worker_value"])
        RA_Info["setting"]["output_verbose_level"] = lowercase(string(get(RA_input, "output_verbose_level", "compact")))
        RA_Info["setting"]["export_dispatch_results_threshold_value"] = get(RA_input, "export_dispatch_results_threshold_value", ra_setting_default["export_dispatch_results_threshold_value"])
        RA_Info["setting"]["export_reference_dispatch_results_flag"] = get(RA_input, "export_reference_dispatch_results_flag", ra_setting_default["export_reference_dispatch_results_flag"])
        RA_Info["setting"]["export_baseline_dispatch_results_flag"] = get(RA_input, "export_baseline_dispatch_results_flag", ra_setting_default["export_baseline_dispatch_results_flag"])
        RA_Info["setting"]["export_ELCC_dispatch_results_flag"] = get(RA_input, "export_ELCC_dispatch_results_flag", ra_setting_default["export_ELCC_dispatch_results_flag"])

        RA_Info["setting"]["risk_filtering_flag"] = get(RA_input, "risk_filtering_flag", ra_setting_default["risk_filtering_flag"])
        if haskey(RA_input, "preselected_days_list")
            preselected_days = RA_input["preselected_days_list"]
            if preselected_days isa AbstractString
                if preselected_days == "[]"
                    RA_Info["setting"]["preselected_days_list"] = []
                else
                    RA_Info["setting"]["preselected_days_list"] = parse.(Int, split(strip(preselected_days, ['[', ']']), ","))
                end
            else
                RA_Info["setting"]["preselected_days_list"] = preselected_days
            end
        else
            if ra_setting_default["preselected_days_list"] == "[]"
                RA_Info["setting"]["preselected_days_list"] = []
            else
                RA_Info["setting"]["preselected_days_list"] = parse.(Int, split(strip(ra_setting_default["preselected_days_list"], ['[', ']']), ","))
            end
        end

        # ELCC setting
        RA_Info["setting"]["calculate_capacity_credit_flag"] = get(RA_input, "calculate_capacity_credit_flag", ra_setting_default["calculate_capacity_credit_flag"])
        RA_Info["setting"]["capacity_credit_type"] = get(RA_input, "capacity_credit_type", get(ra_setting_default, "capacity_credit_type", ""))
        if RA_Info["setting"]["calculate_capacity_credit_flag"] == true
            RA_Info["setting"]["capacity_credit_RA_simulation_method"] = get(RA_input, "capacity_credit_RA_simulation_method", ra_setting_default["capacity_credit_RA_simulation_method"])
            RA_Info["setting"]["capacity_credit_assessment_mode"] = get(RA_input, "capacity_credit_assessment_mode", "add_new")
            RA_Info["setting"]["capacity_credit_reference_RA_metric"] = get(RA_input, "capacity_credit_reference_RA_metric", ra_setting_default["capacity_credit_reference_RA_metric"])
            RA_Info["setting"]["capacity_credit_reference_RA_metric_spatial_resolution"] = get(RA_input, "capacity_credit_reference_RA_metric_spatial_resolution", ra_setting_default["capacity_credit_reference_RA_metric_spatial_resolution"])
            RA_Info["setting"]["ELCC_new_resource_location"] = get(RA_input, "ELCC_new_resource_location", 1)

            RA_Info["setting"]["capacity_credit_rel_tol_value"] = get(RA_input, "capacity_credit_rel_tol_value", ra_setting_default["capacity_credit_rel_tol_value"])
            RA_Info["setting"]["capacity_credit_abs_tol_value"] = get(RA_input, "capacity_credit_abs_tol_value", ra_setting_default["capacity_credit_abs_tol_value"])
            RA_Info["setting"]["capacity_credit_max_iteration_value"] = get(RA_input, "capacity_credit_max_iteration_value", ra_setting_default["capacity_credit_max_iteration_value"])
        end

        # Storage Setting
        RA_Info["setting"]["Post_Contingency_Discharge_method"] = get(RA_input, "Post_Contingency_Discharge_method", ra_setting_default["Post_Contingency_Discharge_method"])
        RA_Info["setting"]["Post_Contingency_Charge_method"] = get(RA_input, "Post_Contingency_Charge_method", ra_setting_default["Post_Contingency_Charge_method"])

        # Sequential ED rolling-horizon look-ahead length. Default 6.
        RA_Info["setting"]["sequential_horizon_hours_value"] = Int(get(RA_input, "sequential_horizon_hours_value", get(ra_setting_default, "sequential_horizon_hours_value", 6)))
        # Redispatch-cost marginal-cost floor.
        RA_Info["setting"]["min_redispatch_mc_value"] = Float64(get(RA_input, "min_redispatch_mc_value", get(ra_setting_default, "min_redispatch_mc_value", 0.0001)))

    else
        RA_Info["setting"]["setting_type"] = "from Setting File"
        RA_Info["setting"]["num_risk_scenario"] = ALEAF_setting["RA Setting"]["num_risk_scenario"]
        RA_Info["setting"]["metrics"] = Dict{String, Any}()
        RA_Info["setting"]["RA_method"] = ALEAF_setting["RA Setting"]["RA_method"]
        RA_Info["setting"]["current_year"] = current_year
        RA_Info["setting"]["logging_level_value"] = global_logging_level
        
        RA_Info["setting"]["reference_temp"] = ALEAF_setting["RA Setting"]["reference_temp"]
        RA_Info["setting"]["system_peak_scale"] = ALEAF_setting["RA Setting"]["system_peak_scale"]  
        RA_Info["setting"]["risk_tol"] = ALEAF_setting["RA Setting"]["risk_tol_value"]  
        RA_Info["setting"]["min_num_risk_in_each_day_value"] = ALEAF_setting["RA Setting"]["min_num_risk_in_each_day_value"]  
        RA_Info["setting"]["repair_time_bound"] = ALEAF_setting["RA Setting"]["repair_time_bound"]  

        RA_Info["setting"]["distributed_run_flag"] = ALEAF_setting["RA Setting"]["distributed_run_flag"]  
        RA_Info["setting"]["num_distributed_scenarios_per_worker_value"] = ALEAF_setting["RA Setting"]["num_distributed_scenarios_per_worker_value"] 
        RA_Info["setting"]["output_verbose_level"] = lowercase(string(get(ALEAF_setting["RA Setting"], "output_verbose_level", "compact")))

        RA_Info["setting"]["export_dispatch_results_threshold_value"] = get(ALEAF_setting["RA Setting"], "export_dispatch_results_threshold_value", Inf)
        RA_Info["setting"]["export_reference_dispatch_results_flag"] = get(ALEAF_setting["RA Setting"], "export_reference_dispatch_results_flag", false)
        RA_Info["setting"]["export_baseline_dispatch_results_flag"] = get(ALEAF_setting["RA Setting"], "export_baseline_dispatch_results_flag", false)
        RA_Info["setting"]["export_ELCC_dispatch_results_flag"] = get(ALEAF_setting["RA Setting"], "export_ELCC_dispatch_results_flag", false)

        RA_Info["setting"]["ELCC_resource_list"] = []
        for tech_idx in keys(ALEAF_setting["Gen Technology"])
            if ALEAF_setting["Gen Technology"][tech_idx]["ELCC_Flag"] == true
                push!(RA_Info["setting"]["ELCC_resource_list"], ALEAF_setting["Gen Technology"][tech_idx]["UNITGROUP"])
            end
        end

        RA_Info["setting"]["risk_filtering_flag"] = ALEAF_setting["RA Setting"]["risk_filtering_flag"]
        if ALEAF_setting["RA Setting"]["preselected_days_list"] == "[]"
            RA_Info["setting"]["preselected_days_list"] = []
        else 
            RA_Info["setting"]["preselected_days_list"] = parse.(Int, split(strip(ALEAF_setting["RA Setting"]["preselected_days_list"], ['[', ']']), ","))
        end

        # ELCC setting
        RA_Info["setting"]["calculate_capacity_credit_flag"] = ALEAF_setting["RA Setting"]["calculate_capacity_credit_flag"]
        RA_Info["setting"]["capacity_credit_type"] = get(ALEAF_setting["RA Setting"], "capacity_credit_type", "")
        if RA_Info["setting"]["calculate_capacity_credit_flag"] == true
            RA_Info["setting"]["capacity_credit_RA_simulation_method"] = ALEAF_setting["RA Setting"]["capacity_credit_RA_simulation_method"]
            RA_Info["setting"]["capacity_credit_assessment_mode"] = get(ALEAF_setting["RA Setting"], "capacity_credit_assessment_mode", "add_new")
            RA_Info["setting"]["capacity_credit_reference_RA_metric"] = ALEAF_setting["RA Setting"]["capacity_credit_reference_RA_metric"]
            RA_Info["setting"]["capacity_credit_reference_RA_metric_spatial_resolution"] = ALEAF_setting["RA Setting"]["capacity_credit_reference_RA_metric_spatial_resolution"]

            RA_Info["setting"]["capacity_credit_rel_tol_value"] = ALEAF_setting["RA Setting"]["capacity_credit_rel_tol_value"]
            RA_Info["setting"]["capacity_credit_abs_tol_value"] = ALEAF_setting["RA Setting"]["capacity_credit_abs_tol_value"]
            RA_Info["setting"]["capacity_credit_max_iteration_value"] = ALEAF_setting["RA Setting"]["capacity_credit_max_iteration_value"]
            
            # temporal
            RA_Info["setting"]["ELCC_new_resource_location"] = 1
        end

        # Storage Setting
        RA_Info["setting"]["Post_Contingency_Discharge_method"] = ALEAF_setting["RA Setting"]["Post_Contingency_Discharge_method"]
        RA_Info["setting"]["Post_Contingency_Charge_method"] = ALEAF_setting["RA Setting"]["Post_Contingency_Charge_method"]

        # Sequential ED rolling-horizon look-ahead length. Default 6.
        RA_Info["setting"]["sequential_horizon_hours_value"] = Int(get(ALEAF_setting["RA Setting"], "sequential_horizon_hours_value", 6))
        # Redispatch-cost marginal-cost floor.
        RA_Info["setting"]["min_redispatch_mc_value"] = Float64(get(ALEAF_setting["RA Setting"], "min_redispatch_mc_value", 0.0001))

    end
    
    RA_Info["RA_reference_risk_data"] = Dict{String, Any}()
    RA_Info["RA_solutions"] = Dict{String, Any}()
    RA_Info["RA_metrics"] = Dict{String, Any}()
    RA_Info["system_info"] = Dict{String, Any}()
    RA_Info["capacity_credit_result"] = Dict{String, Any}()

    # gen disaggregation flag
    if ALEAF_setting["Simulation Setting"]["generation_representation_option"] == "Aggregated_by_Tech"
        RA_Info["setting"]["gen_disaggregation_flag"] = true
    elseif ALEAF_setting["Simulation Setting"]["generation_representation_option"] == "Individual_Plants"
        RA_Info["setting"]["gen_disaggregation_flag"] = false
    end

    if haskey(RA_input, "prior_soc_adjust")
        RA_Info["setting"]["prior_soc_adjust"] = RA_input["prior_soc_adjust"]
    end
    
    # Parse enabled RA renewable scenarios into settings.
    parse_ra_renewable_scenarios!(ALEAF_setting, RA_Info)
    load_ra_renewable_scenario_timeseries!(ALEAF_setting, RA_Info)

    # Scenario name
    case = ALEAF_setting["Simulation Setting"]["test_system_name"]
    case_name = ALEAF_setting["Simulation Configuration"][string(case_id)]["Case_ID"]

    RA_Info["scenario"] = case_name

    return RA_Info

end


function get_risk_scenario_of_new_gen_for_RA!(am_reference_data, am_setting_data, RA_setting, gen_info, ELCC_Info)
    
    reference_temp = RA_setting["reference_temp"]
    simulation_list = [i for i in 1:am_setting_data["Simulation Configuration"]["NDAY_Groups_RA"]]

    ELCC_Info["new_gen_risk_data"] = generate_outage_scenarios_for_single_unit(am_reference_data[:nw][0], simulation_list, reference_temp, RA_setting["num_risk_scenario"], gen_info, RA_setting; output_option="Status")
    
    return ELCC_Info
end


function find_gen_tech_info(am::Abstract_ALEAF_Model, new_resource_unit_group)

    gen_tech_id = 0
    for idx in keys(am.ref[:nw][0][:gen_technology])
        if (am.ref[:nw][0][:gen_technology][idx]["UNITGROUP"] == new_resource_unit_group) 
            gen_tech_id = idx
        end
    end

    if gen_tech_id == 0
        gen_tech_id = nothing
        return gen_tech_id
    else
        update_regional_fuel_prices_for_new_gen_tech_with_pu!(am, am.ref[:nw][0][:gen_technology][gen_tech_id])
        return am.ref[:nw][0][:gen_technology][gen_tech_id]
    end

            
end


function update_regional_fuel_prices_for_new_gen_tech_with_pu!(am::Abstract_ALEAF_Model, gen_tech_info)
                                                                                                                                                    
    pu_power_base = am.setting["Simulation Setting"]["per_unit_base_value"]                                                                       
    pu_econ_base = am.setting["Simulation Setting"]["per_unit_econ_base_value"] / pu_power_base    # per unit economic base
    region_data_identifier = get(am.setting["Network Setting"], "regional_fuel_zone_resolution_type", "system")
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

    gen_ref = gen_tech_info

    # check region based on data_identifier
    fuel_region = "System-wide"
    if gen_tech_info["FC"] == "Fuel-Regional"                    
        fuel_region = am.ref[:nw][0][:bus][bus_id]["region_mapping_info"][region_data_identifier]
    end

    fuel_type = gen_tech_info["FUEL"]
    
    for year_id in keys(am.ref[:nw][0][:planning_stages])
        
        year_string = string(year_id)
        gen_tech_info["Annual_FC"][year_string] = Dict{Int64, Any}() # initialize
        gen_tech_info["Annual_MC"][year_string] = Dict{Int64, Any}() # initialize
        
        vom = gen_tech_info["Annual_VOM"][year_string] * pu_econ_base # $/MWh

        for day_id in keys(am.ref[:nw][0][:repdays])

            month = am.ref[:nw][0][:repdays][day_id]["Month"]

            fuel_cost = gen_tech_info["FC"]
            if fuel_cost in ["Fuel", "Fuel-Regional"]
                fuel_cost = find_monthly_fuel_cost(fuel_region, fuel_type, am.ref[:nw][0][:planning_stages][year_id]["year"], month) 
            end

            # add fuel cost
            gen_tech_info["Annual_FC"][year_string][day_id] = fuel_cost 
            gen_tech_info["Annual_MC"][year_string][day_id] = gen_tech_info["HR"] * fuel_cost + vom
                                                    
            # apply per unit
            gen_tech_info["Annual_FC"][year_string][day_id] /= pu_econ_base
            gen_tech_info["Annual_MC"][year_string][day_id] /= pu_econ_base

        end
        
    end    
                                                                                                                                             
end


function build_ALEAF_model_instance_for_RA(ALEAF_setting, case_id, current_year, network_data, recorded_investment_decisions)

    # build common ALEAF model instance structure
    ALEAF_model_instance = build_model_structure_GTEP(ALEAF_setting, network_data, Abstract_LC_GTEP_Model)
    add_ref_LC_RA_model!(ALEAF_model_instance, ALEAF_setting, recorded_investment_decisions)

    return ALEAF_model_instance
end


function add_ref_LC_RA_model!(am::Abstract_ALEAF_Model, ALEAF_setting::Dict{String, Any}, recorded_investment_decisions)
    
    # model setting
    add_model_setting_GTEP!(am, ALEAF_setting)
    
    # run period
    add_run_period_GTEP!(am)

    # add common reference
    add_common_ref_GTEP!(am)

    # generation aggregation
    if ALEAF_setting["Simulation Setting"]["generation_representation_option"] == "Aggregated_by_Tech"
        aggregate_generation_GTEP!(am)
    elseif ALEAF_setting["Simulation Setting"]["generation_representation_option"] == "Individual_Plants"

        if ALEAF_setting["Simulation Setting"]["generation_parameter_source_option"] == "Individual Plant"
            update_individual_plant_info_GTEP_using_plant_info!(am, ALEAF_setting; recorded_investment_decisions) 
        elseif ALEAF_setting["Simulation Setting"]["generation_parameter_source_option"] == "Gen Technology"
            update_individual_plant_info_GTEP_using_tech_info!(am; recorded_investment_decisions)
        end
    end

    # add gen index
    add_gen_index_ref_GTEP!(am)

    # add local gen index
    add_local_gen_idx_GTEP!(am)

    # adjust gen capacity for UC
    update_CAP_i_noEXP!(am)

    # calculate vre class average generation (per year, bus, tech)
    aggregate_vre_timeseries_per_bus_LC_GTEP!(am)

    # calculate vre class average generation (per year, zone, tech,)
    aggregate_vre_timeseries_per_policy_zone_LC_GTEP!(am)

    # calculate vre class average generation (per year, zone, day, time, tech)
    aggregate_vre_timeseries_zdt_GTEP!(am)

    # resource-bearing sub-region masks reused when RA re-aggregates vre shapes per renewable scenario
    add_ra_vre_region_masks!(am)

    # water level management data
    add_water_management_info_LC_GTEP!(am)

    # add per unit to the market parameters
    apply_per_unit_to_market_parameters!(am, ALEAF_setting)

    # update water budget values
    update_water_budget_values_with_pu!(am)

    # update regional fuel prices and marginal costs
    update_regional_fuel_prices_with_pu!(am, ALEAF_setting)

end


function _ra_hourly_vre_profile(hourly_data, tech::String)
    labels = Dict(
        "csp" => haskey(hourly_data, "csp_BA") ? "csp_BA" : "csp",
        "wind_ons" => haskey(hourly_data, "wind_ons_BA") ? "wind_ons_BA" : "wind_ons",
        "wind_ofs" => haskey(hourly_data, "wind_ofs_BA") ? "wind_ofs_BA" : "wind_ofs",
        "pv" => "pv",
        "rtpv" => "rtpv",
        "hydro" => "hydro",
    )
    label = get(labels, tech, tech)
    return haskey(hourly_data, label) ? hourly_data[label] : Dict{Any, Any}()
end


function add_ra_vre_region_masks!(am::Abstract_ALEAF_Model)

    num_stages = am.setting["Planning Design"]["num_stages_value"]

    if haskey(am.ref[:nw][0], :operation_year)
        ids_y = [am.ref[:nw][0][:operation_year]]
    else
        ids_y = [i for i in 1:num_stages]
    end

    ids_zone = [(z) for (z) in get_index(am, :bus, 0)]
    ids_d = [(d) for (d) in get_index(am, :repdays, 0)]
    vre_techs = String["hydro", "csp", "wind_ons", "wind_ofs", "rtpv", "pv"]

    # Only regions that carry the resource at some hour belong in the zone average; a region that is
    # zero (or absent) across every rep-day hour has no resource and must stay out of the denominator.
    region_masks = Dict{Int64, Dict{Int64, Dict{String, Vector{String}}}}()

    for year in ids_y

        repday_data = am.ref[:nw][0][:planning_stages][year]["repdays"]
        resource_bearing = Dict{String, Set{String}}(tech => Set{String}() for tech in vre_techs)

        for d in ids_d
            hourly_data = repday_data[string(d)]["data"]
            for h_key in keys(hourly_data)
                sub_period_data = hourly_data[h_key]["1"]
                for tech in vre_techs
                    bearing_regions = resource_bearing[tech]
                    for (region_id, value) in _ra_hourly_vre_profile(sub_period_data, tech)
                        if value != 0
                            push!(bearing_regions, string(region_id))
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
            for tech in vre_techs
                bearing_regions = resource_bearing[tech]
                # Mirror the GTEP builder: keep the finest region id as the mask entry, but test
                # resource-bearing at its DATA-REGION (hourly data is keyed by data-region, not finest id).
                zone_masks[tech] = String[string(region_id) for region_id in region_list if string(profile_data_region(region_map, tech, region_id)) in bearing_regions]
            end
        end

    end

    am.ref[:nw][0][:vre_region_masks] = region_masks

end


function calculate_max_units(am, bus_idx, genco_tech_id)
    """Calculate maximum number of individual plant units for a technology"""
    
    current_U_G_iy = get(am.ref[:nw][bus_idx][:gen_bus][genco_tech_id], "U_G_iy", am.ref[:nw][bus_idx][:gen_bus][genco_tech_id]["EXUNITS"])

    # Get original existing units before expansion if original was overwritten
    if haskey(am.ref[:nw][bus_idx][:gen_bus][genco_tech_id], "Original_EXUNITS")
        Original_EXUNITS = am.ref[:nw][bus_idx][:gen_bus][genco_tech_id]["Original_EXUNITS"]
    else
        Original_EXUNITS = am.ref[:nw][bus_idx][:gen_bus][genco_tech_id]["EXUNITS"]
    end
    
    # Respect RA option for including investment options in generator indexing.
    # If disabled, size only from existing units to avoid pre-allocating zero-cap slots.
    include_investment_options = false
    if haskey(am.setting, "RA Setting")
        include_investment_options = get(am.setting["RA Setting"], "generate_gen_index_with_investment_options_flag", false)
    end

    if include_investment_options
        # Per bus max investment (same as bounds from variable_u_newG_iy_integer_real)
        MAXINVEST = parameter(am, bus_idx, :gen_bus, genco_tech_id, "MAXINVEST")
        max_U_G_iy = max(current_U_G_iy, Original_EXUNITS + MAXINVEST)
    else
        max_U_G_iy = max(current_U_G_iy, Original_EXUNITS)
    end

    max_num_units = ceil(Int, max_U_G_iy)
    return max(1, max_num_units)
end


function prepare_RA_gen!(am)

    # create a new gen index dict
    am.ref[:nw][0][:RA_gen_index] = Dict{Int64, Any}()
    for bus_idx in keys(am.ref[:nw][0][:bus])
        am.ref[:nw][bus_idx][:RA_gen_bus] = Dict{Int64, Any}()
        am.ref[:nw][0][:bus][bus_idx]["new_disaggregated_local_gen_idx"] = []
    end
    RA_gen_idx = 1

    gen_disaggregation_flag = false

    for gen_idx in keys(am.ref[:nw][0][:gen_index])

        bus_idx = am.ref[:nw][0][:gen_index][gen_idx]["bus_idx"]
        genco_tech_id = am.ref[:nw][0][:gen_index][gen_idx]["genco_tech_id"]
        U_G_iy = parameter(am, bus_idx, :gen_bus, genco_tech_id, "U_G_iy")
        RA_FOR = parameter(am, bus_idx, :gen_bus, genco_tech_id, "RA_FOR")
        ICAP = parameter(am, bus_idx, :gen_bus, genco_tech_id, "CAP")
        CAPCRED = parameter(am, bus_idx, :gen_bus, genco_tech_id, "CAPCRED")
        if CAPCRED isa String
            data_identifier = am.ref[:nw][0][:zone]["capacity_credit"]["data_identifier"]
            data_region = am.ref[:nw][0][:bus][bus_idx]["region_mapping_info"][data_identifier]
            data_region = data_region isa AbstractString ? data_region : string(data_region)
            CAPCRED = am.ref[:nw][0][:zone]["capacity_credit"][data_region][CAPCRED]
        end
        max_units = calculate_max_units(am, bus_idx, genco_tech_id)

        # Calculate actual capacity distribution
        if U_G_iy > 1.0
            gen_disaggregation_flag = true
            EXCAP = U_G_iy * ICAP
            num_units = round(U_G_iy, RoundDown)
            remainder = mod(EXCAP, ICAP)
            capacity_to_add = remainder / num_units
            final_capacity = ICAP + capacity_to_add
        elseif U_G_iy > 0.0
            num_units = 1
            final_capacity = ICAP * U_G_iy
        else
            num_units = 0
            final_capacity = 0.0
        end

        # Pre-allocate all slots for this technology
        for slot_num in 1:max_units
            
            am.ref[:nw][0][:RA_gen_index][RA_gen_idx] = deepcopy(am.ref[:nw][0][:gen_index][gen_idx])
            # copy gen_bus data as well
            # define new_genco_tech_id - note this produces a new id even if max_units = 1
            new_genco_tech_id = parse(Int, string(genco_tech_id) * "00" * string(slot_num))
            am.ref[:nw][bus_idx][:RA_gen_bus][new_genco_tech_id] = deepcopy(am.ref[:nw][bus_idx][:gen_bus][genco_tech_id])
            am.ref[:nw][bus_idx][:RA_gen_bus][new_genco_tech_id]["gen_idx"] = RA_gen_idx

            # Update RA_gen_index with fields that are independent of U_G_iy
            am.ref[:nw][0][:RA_gen_index][RA_gen_idx]["bus_name"] = am.ref[:nw][0][:bus][bus_idx]["bus_i"]
            am.ref[:nw][0][:RA_gen_index][RA_gen_idx]["original_gen_index"] = gen_idx
            am.ref[:nw][0][:RA_gen_index][RA_gen_idx]["RA_FOR"] = RA_FOR
            am.ref[:nw][0][:RA_gen_index][RA_gen_idx]["genco_tech_id"] = new_genco_tech_id
            am.ref[:nw][0][:RA_gen_index][RA_gen_idx]["slot_number"] = slot_num
            
            # Fill first num_units slots, leave rest at zero; this is needed for the indexing consistency
            if slot_num <= num_units
                am.ref[:nw][0][:RA_gen_index][RA_gen_idx]["ICAP"] = final_capacity
                am.ref[:nw][0][:RA_gen_index][RA_gen_idx]["CAP"] = final_capacity
                am.ref[:nw][0][:RA_gen_index][RA_gen_idx]["UCAP"] = final_capacity * CAPCRED
                am.ref[:nw][0][:RA_gen_index][RA_gen_idx]["U_G_iy"] = 1.0
                # These last two assignments weren't made when U_G_iy > 1.0 previously
                am.ref[:nw][0][:RA_gen_index][RA_gen_idx]["EXUNITS"] = 1.0
                am.ref[:nw][0][:RA_gen_index][RA_gen_idx]["EXCAPS"] = final_capacity

                # Assign gen_bus capacity
                am.ref[:nw][bus_idx][:RA_gen_bus][new_genco_tech_id]["ICAP"] = final_capacity
                am.ref[:nw][bus_idx][:RA_gen_bus][new_genco_tech_id]["CAP"] = final_capacity
                am.ref[:nw][bus_idx][:RA_gen_bus][new_genco_tech_id]["U_G_iy"] = 1.0
                am.ref[:nw][bus_idx][:RA_gen_bus][new_genco_tech_id]["EXUNITS"] = 1.0
                am.ref[:nw][bus_idx][:RA_gen_bus][new_genco_tech_id]["EXCAPS"] = final_capacity

                # update storage charge cap
                if am.ref[:nw][bus_idx][:RA_gen_bus][new_genco_tech_id]["Charge_CAP"] > 0
                    am.ref[:nw][bus_idx][:RA_gen_bus][new_genco_tech_id]["Charge_CAP"] = final_capacity
                end
            else
                am.ref[:nw][0][:RA_gen_index][RA_gen_idx]["ICAP"] = 0.0
                am.ref[:nw][0][:RA_gen_index][RA_gen_idx]["CAP"] = 0.0
                am.ref[:nw][0][:RA_gen_index][RA_gen_idx]["UCAP"] = 0.0
                am.ref[:nw][0][:RA_gen_index][RA_gen_idx]["U_G_iy"] = 0.0
                am.ref[:nw][0][:RA_gen_index][RA_gen_idx]["EXUNITS"] = 0.0
                am.ref[:nw][0][:RA_gen_index][RA_gen_idx]["EXCAPS"] = 0.0

                am.ref[:nw][bus_idx][:RA_gen_bus][new_genco_tech_id]["ICAP"] = 0.0
                am.ref[:nw][bus_idx][:RA_gen_bus][new_genco_tech_id]["CAP"] = 0.0
                am.ref[:nw][bus_idx][:RA_gen_bus][new_genco_tech_id]["U_G_iy"] = 0.0
                am.ref[:nw][bus_idx][:RA_gen_bus][new_genco_tech_id]["EXUNITS"] = 0.0
                am.ref[:nw][bus_idx][:RA_gen_bus][new_genco_tech_id]["EXCAPS"] = 0.0

                if am.ref[:nw][bus_idx][:RA_gen_bus][new_genco_tech_id]["Charge_CAP"] > 0
                    am.ref[:nw][bus_idx][:RA_gen_bus][new_genco_tech_id]["Charge_CAP"] = 0.0
                end
            end

            # update storage MWh
            if am.ref[:nw][0][:gen_index][gen_idx]["UNIT_CATEGORY"] == "STORAGE"
                original_ES_MWh = parameter(am, bus_idx, :gen_bus, genco_tech_id, "ES_MWh")
                if slot_num <= num_units
                    # original_ES_MWh is total energy for the active fleet. Allocate it
                    # to each RA slot in proportion to the slot's power capacity.
                    total_capacity = U_G_iy * ICAP
                    scaled_ES_MWh = total_capacity > 0 ? original_ES_MWh * (final_capacity / total_capacity) : 0.0
                    am.ref[:nw][bus_idx][:RA_gen_bus][new_genco_tech_id]["new_u_ESE_iy"] = scaled_ES_MWh
                    am.ref[:nw][bus_idx][:RA_gen_bus][new_genco_tech_id]["ES_MWh"] = scaled_ES_MWh
                else
                    am.ref[:nw][bus_idx][:RA_gen_bus][new_genco_tech_id]["new_u_ESE_iy"] = 0.0
                    am.ref[:nw][bus_idx][:RA_gen_bus][new_genco_tech_id]["ES_MWh"] = 0.0
                end
            end

            # update new local gen index list
            push!(am.ref[:nw][0][:bus][bus_idx]["new_disaggregated_local_gen_idx"], RA_gen_idx)

            RA_gen_idx += 1
        end
    end

    # Rebuild hybrid linkage indices in RA index space.
    for ra_idx in keys(am.ref[:nw][0][:RA_gen_index])
        if haskey(am.ref[:nw][0][:RA_gen_index][ra_idx], "hybrid_main_gen_idx")
            delete!(am.ref[:nw][0][:RA_gen_index][ra_idx], "hybrid_main_gen_idx")
        end
        if haskey(am.ref[:nw][0][:RA_gen_index][ra_idx], "hybrid_ES_gen_idx")
            delete!(am.ref[:nw][0][:RA_gen_index][ra_idx], "hybrid_ES_gen_idx")
        end
    end

    for ra_idx in keys(am.ref[:nw][0][:RA_gen_index])
        hybrid_type = am.ref[:nw][0][:RA_gen_index][ra_idx]["hybrid_type"]
        if hybrid_type == "GEN"
            bus_idx = am.ref[:nw][0][:RA_gen_index][ra_idx]["bus_idx"]
            genco_tech_id = am.ref[:nw][0][:RA_gen_index][ra_idx]["genco_tech_id"]
            hybrid_id = am.ref[:nw][bus_idx][:RA_gen_bus][genco_tech_id]["PLANT_NAME"]

            for other_ra_idx in am.ref[:nw][0][:bus][bus_idx]["new_disaggregated_local_gen_idx"]
                if other_ra_idx != ra_idx
                    other_hybrid_type = am.ref[:nw][0][:RA_gen_index][other_ra_idx]["hybrid_type"]
                    if other_hybrid_type == "ES"
                        other_genco_tech_id = am.ref[:nw][0][:RA_gen_index][other_ra_idx]["genco_tech_id"]
                        other_hybrid_id = am.ref[:nw][bus_idx][:RA_gen_bus][other_genco_tech_id]["PLANT_NAME"]
                        if other_hybrid_id == hybrid_id
                            am.ref[:nw][0][:RA_gen_index][other_ra_idx]["hybrid_main_gen_idx"] = ra_idx
                            am.ref[:nw][0][:RA_gen_index][ra_idx]["hybrid_ES_gen_idx"] = other_ra_idx
                        end
                    end
                end
            end
        end
    end

    am.ref[:nw][0][:RA_gen_index] = sort(am.ref[:nw][0][:RA_gen_index])

    for bus_idx in keys(am.ref[:nw][0][:bus])
        am.ref[:nw][0][:bus][bus_idx]["aggregation_info"]["new_local_gen_idx"] = deepcopy(am.ref[:nw][0][:bus][bus_idx]["new_disaggregated_local_gen_idx"])
    end

    if gen_disaggregation_flag == true
        @aleaf_info "[ALEAF RA Model]: Prepare generator data: performed generation disaggregation"

    end

end

#= RA via ED: filter to critical risk scenarios, run reserve-free ED (respecting
   topology/power flow) on them, and derive metrics from the load-shedding results. =#

function execute_RA_model_with_ED!(RA_Info, ALEAF_model_instance, risk_scenario_dict; report_log=true) 
    
    # Run RA model in parallel
    RA_Info["RA_solutions"] = perform_ED_analysis_for_RA(ALEAF_model_instance.ref, ALEAF_model_instance.setting, RA_Info["setting"], risk_scenario_dict, RA_Info["RA_reference_risk_data"]["filtered_joint_scenario_map"]; report_log=report_log)
    if report_log @aleaf_info "[ALEAF RA Model]: ED Analysis: Completed." end

    # Calculate RA metrics (EUE, LOLE)
    RA_Info = calculate_risk_metrics_ED!(ALEAF_model_instance, RA_Info)
    if report_log @aleaf_info "[ALEAF RA Model]: RA Metrics Calculation: Completed." end

    return RA_Info
end


function perform_ELCC_analysis_ED(ALEAF_model_instance, RA_Info, simulation_setting, risk_scenario_dict)

    # copy RA info for ELCC calculations
    ELCC_RA_Info = Dict{String, Any}()
    ELCC_RA_Info["setting"] = RA_Info["setting"]
    ELCC_RA_Info["RA_metrics"] = deepcopy(RA_Info["RA_metrics"])
    ELCC_RA_Info["system_info"] = Dict{String, Any}()
    ELCC_RA_Info["RA_reference_risk_data"] = get(RA_Info, "RA_reference_risk_data", Dict{String, Any}())
    
    # RA method
    RA_method = ELCC_RA_Info["setting"]["RA_method"]
    new_resource_unit_group = simulation_setting["new_resource_unit_group"]
    reference_RA_metric = simulation_setting["reference_RA_metric"]
    reference_RA_metric_spatial_resolution = simulation_setting["reference_RA_metric_spatial_resolution"] 
    new_resource_location = simulation_setting["new_resource_location"]
    metric_scale_mw = ALEAF_model_instance.setting["Simulation Setting"]["per_unit_base_value"]
    constant_load_distribution = reference_RA_metric_spatial_resolution == "Systemwide" ? "load_weighted_systemwide" : "single_bus"
    constant_load_bus_idx = constant_load_distribution == "single_bus" ? new_resource_location : nothing
    constant_load_distribution_profile = Dict{Int, Float64}()

    @aleaf_info "[ALEAF RA Model]\t ELCC Process for ($new_resource_unit_group at bus $new_resource_location) using RA method: $RA_method"
    @aleaf_info "[ALEAF RA Model]\t ELCC load adjustment mode: $constant_load_distribution (profile_hours=$(length(constant_load_distribution_profile)))"
    
    # define Capacity Credit assessment summary dictionary
    ELCC_Info = Dict{String, Any}()
    ELCC_Info["reference_unit_group"] = new_resource_unit_group
    ELCC_Info["reference_unit_location"] = new_resource_location

    # reference RA metric of the reference system
    ref_RA_metric = 0.0
    new_resource_bus_name = ALEAF_model_instance.ref[:nw][0][:bus][new_resource_location]["bus_i"]
    if reference_RA_metric_spatial_resolution == "Systemwide"
        ref_RA_metric = ELCC_RA_Info["RA_metrics"][reference_RA_metric]["Annual"]["Systemwide"]
    elseif reference_RA_metric_spatial_resolution == "Regional"
        ref_RA_metric = ELCC_RA_Info["RA_metrics"][reference_RA_metric]["Annual"]["Regional"][string(new_resource_bus_name)]
    end

    # Add a resource to the system
    # 1) get ref gen info
    original_gen_index = deepcopy(ALEAF_model_instance.ref[:nw][0][:RA_gen_index])

    ref_tech_info = find_gen_tech_info(ALEAF_model_instance, new_resource_unit_group)
    if ref_tech_info == nothing
        ELCC_Info["ELCC"] = 0.0
        ELCC_Info["status"] = "tartet technology not found"
        @aleaf_warn "[ALEAF RA Model]\t ELCC Existing-Asset Process: target technology not found for unit_group=$new_resource_unit_group"
        return ELCC_Info
    end

    ELCC_Info["reference_unit_ICAP"] = ref_tech_info["CAP"]

    # 2) add new resource to the system
    new_gen_idx = maximum(get_index(ALEAF_model_instance, :RA_gen_index, 0)) + 1
    new_genco_tech_id = parse(Int, string("100", new_gen_idx))
    original_local_gen_idx = nothing
    try
        ALEAF_model_instance.ref[:nw][0][:RA_gen_index][new_gen_idx] = Dict{String, Any}(
            
            "bus_idx"              => new_resource_location,
            "UNIT_GROUP"           => ref_tech_info["UNITGROUP"],
            "ICAP"                 => ref_tech_info["CAP"],
            "CAP"                  => ref_tech_info["CAP"],
            "U_G_iy"               => 1.0,
            "genco_tech_UNIT_REPORT_LABEL_1" => ref_tech_info["UNIT_REPORT_LABEL_1"],
            "genco_tech_UNIT_REPORT_LABEL_2" => ref_tech_info["UNIT_REPORT_LABEL_2"],
            "Profile_Type"         => ref_tech_info["Profile_Type"],
            "bus_name"             => ALEAF_model_instance.ref[:nw][0][:bus][new_resource_location]["bus_i"],
            "genco_tech_id"        => new_genco_tech_id,
            "UNIT_CATEGORY"        => ref_tech_info["UNIT_CATEGORY"],
            "original_gen_index"   => 0,
            "RA_FOR"               => ref_tech_info["RA_FOR"],
            "new_unit"             => true,
            "PLANT_NAME"           => string("new_ELCC_", ref_tech_info["UNIT_CATEGORY"]),

            # Hybrid systems are not supported in ELCC analysis yet.
            "hybrid_type"          => "NA",
            "hybrid_main_gen_idx"  => 0,
            "hybrid_ES_gen_idx"    => 0
        )

        # update gen dict
        ALEAF_model_instance.ref[:nw][new_resource_location][:RA_gen_bus][new_genco_tech_id] = deepcopy(ref_tech_info)
        ALEAF_model_instance.ref[:nw][new_resource_location][:RA_gen_bus][new_genco_tech_id]["Timeseries_Tag"] = "LOCAL"
        ALEAF_model_instance.ref[:nw][new_resource_location][:RA_gen_bus][new_genco_tech_id]["EXUNITS"] = 1.0
        ALEAF_model_instance.ref[:nw][new_resource_location][:RA_gen_bus][new_genco_tech_id]["PLANT_NAME"] = string("new_ELCC_", ref_tech_info["UNIT_CATEGORY"])
        if ref_tech_info["UNIT_CATEGORY"] == "STORAGE"
            ALEAF_model_instance.ref[:nw][new_resource_location][:RA_gen_bus][new_genco_tech_id]["ES_MWh"] = ref_tech_info["STOHR_MAX"] * ref_tech_info["CAP"]
        end

        # create ELCC output path
        output_path = string(ALEAF_model_instance.ref[:nw][0][:output_path], "ELCC_of_", ref_tech_info["UNIT_REPORT_LABEL_2"], "_at_bus_", new_resource_location, "/")
        ALEAF.check_and_create_path(output_path)
        RA_Info["setting"]["ELCC_output_path"] = output_path

        # update local gen index
        original_local_gen_idx = deepcopy(ALEAF_model_instance.ref[:nw][0][:bus][new_resource_location]["aggregation_info"]["new_local_gen_idx"])
        push!(ALEAF_model_instance.ref[:nw][0][:bus][new_resource_location]["aggregation_info"]["new_local_gen_idx"], new_gen_idx)
        
        # Reference RA metrics
        @aleaf_info "[ALEAF RA Model]\t ELCC Process: iteration=0, target_ref_metric=$(round(ref_RA_metric, digits=3)), added_unit_icap=$(ELCC_Info["reference_unit_ICAP"]*metric_scale_mw) MW"
        
        # get risk scenario of the new resource
        ELCC_Info = get_risk_scenario_of_new_gen_for_RA!(ALEAF_model_instance.ref, ALEAF_model_instance.setting, RA_Info["setting"], ALEAF_model_instance.ref[:nw][0][:RA_gen_index][new_gen_idx], ELCC_Info)
        ELCC_Info["filtered_joint_scenario_map"] = RA_Info["RA_reference_risk_data"]["filtered_joint_scenario_map"]
        @aleaf_info "[ALEAF RA Model]\t ELCC Process:  Obtained risk scenario of the new resource"

        # Solve for ELCC via bracketed search on added load (see search block below).
        ELCC_RA_Info["RA_solutions"] = perform_ED_analysis_for_RA(ALEAF_model_instance.ref, ALEAF_model_instance.setting, RA_Info["setting"], risk_scenario_dict, ELCC_Info["filtered_joint_scenario_map"]; risk_scenario_dict_of_new_gen = ELCC_Info["new_gen_risk_data"], constant_load_bus_idx, constant_load_distribution, constant_load_distribution_profile, report_log=false)
        ELCC_RA_Info = calculate_risk_metrics_ED!(ALEAF_model_instance, ELCC_RA_Info)

        # get new RA metric
        initial_RA_metric = 0.0
        if reference_RA_metric_spatial_resolution == "Systemwide"
            initial_RA_metric = ELCC_RA_Info["RA_metrics"][reference_RA_metric]["Annual"]["Systemwide"]
        elseif reference_RA_metric_spatial_resolution == "Regional"
            initial_RA_metric = ELCC_RA_Info["RA_metrics"][reference_RA_metric]["Annual"]["Regional"][string(new_resource_bus_name)]
        end
            
        RA_improvement = ref_RA_metric - initial_RA_metric
        @aleaf_info "[ALEAF RA Model]\t ELCC Process: iteration=1, new_metric=$(round(initial_RA_metric, digits=3)), target_ref_metric=$(round(ref_RA_metric, digits=3)), improvement_vs_target=$(round(RA_improvement, digits=3))"
        if RA_improvement > 0

            # Bracketed search on added load to return to the reference RA metric;
            # avoids Newton failures when storage dispatch changes discontinuously.
            best_constant_load = 0.0
            min_relative_gap = Inf
            iter = 2
            tol = RA_Info["setting"]["capacity_credit_rel_tol_value"] * 0.01
            abs_tol = parse(Float64, string(RA_Info["setting"]["capacity_credit_abs_tol_value"]))
            max_iter = RA_Info["setting"]["capacity_credit_max_iteration_value"]
            RA_metric_gap = initial_RA_metric - ref_RA_metric
            relative_RA_metric_gap = abs(ref_RA_metric) > 1e-9 ? abs(RA_metric_gap) / abs(ref_RA_metric) : abs(RA_metric_gap)
            converged = false

            function evaluate_added_load(added_load, iter)
                constant_load = max(0.0, added_load)
                eval_RA_Info = ELCC_RA_Info
                eval_RA_Info["RA_solutions"] = perform_ED_analysis_for_RA(ALEAF_model_instance.ref, ALEAF_model_instance.setting, RA_Info["setting"], risk_scenario_dict, ELCC_Info["filtered_joint_scenario_map"]; constant_load, constant_load_bus_idx, constant_load_distribution, constant_load_distribution_profile, risk_scenario_dict_of_new_gen = ELCC_Info["new_gen_risk_data"], report_log=false, iter)
                eval_RA_Info = calculate_risk_metrics_ED!(ALEAF_model_instance, eval_RA_Info)

                new_RA_metric = 0.0
                if reference_RA_metric_spatial_resolution == "Systemwide"
                    new_RA_metric = eval_RA_Info["RA_metrics"][reference_RA_metric]["Annual"]["Systemwide"]
                elseif reference_RA_metric_spatial_resolution == "Regional"
                    new_RA_metric = eval_RA_Info["RA_metrics"][reference_RA_metric]["Annual"]["Regional"][string(new_resource_bus_name)]
                end

                metric_gap = new_RA_metric - ref_RA_metric
                if abs(ref_RA_metric) > 1e-9
                    rel_gap = abs(metric_gap) / abs(ref_RA_metric)
                else
                    rel_gap = abs(metric_gap)
                end

                @aleaf_info "[ALEAF RA Model]\t ELCC Process: iteration=$iter, added_load=$(round(added_load*metric_scale_mw, digits=3)) MW, new_metric=$(round(new_RA_metric, digits=3)), target_ref_metric=$(round(ref_RA_metric, digits=3)), abs_gap_to_target=$(round(metric_gap, digits=3)), rel_gap_to_target=$(round(rel_gap, digits=5))"
                return new_RA_metric, metric_gap, rel_gap
            end

            lower_load = 0.0
            lower_gap = initial_RA_metric - ref_RA_metric   # gap at added_load=0 (from improved-system RA at iter=1)
            upper_load = ELCC_Info["reference_unit_ICAP"]
            _, upper_gap, upper_rel_gap = evaluate_added_load(upper_load, iter)
            RA_metric_gap = upper_gap
            relative_RA_metric_gap = upper_rel_gap
            best_constant_load = upper_load
            min_relative_gap = upper_rel_gap

            if abs(upper_gap) <= abs_tol || upper_rel_gap < tol
                converged = true
            else
                iter += 1
            end

            while !converged && upper_gap < 0.0 && iter <= max_iter
                lower_load = upper_load
                lower_gap = upper_gap                       # track gap at new lower for regula falsi
                upper_load *= 2.0
                _, upper_gap, upper_rel_gap = evaluate_added_load(upper_load, iter)
                RA_metric_gap = upper_gap
                relative_RA_metric_gap = upper_rel_gap
                if upper_rel_gap < min_relative_gap
                    min_relative_gap = upper_rel_gap
                    best_constant_load = upper_load
                end
                if abs(upper_gap) <= abs_tol || upper_rel_gap < tol
                    converged = true
                    break
                end
                iter += 1
            end

            if !converged && upper_gap >= 0.0
                # Illinois-protected regula falsi: interpolate between bracket endpoints
                # by gap; faster than bisection, weights demote a stuck endpoint.
                w_lower = 1.0
                w_upper = 1.0
                last_side = :none
                while iter <= max_iter
                    denom = w_upper * lower_gap - w_lower * upper_gap
                    if abs(denom) < 1e-12
                        trial_load = 0.5 * (lower_load + upper_load)    # degenerate: fall back to bisection
                    else
                        trial_load = (w_upper * lower_gap * upper_load - w_lower * upper_gap * lower_load) / denom
                    end
                    bracket_eps = 1e-9 * max(1.0, upper_load - lower_load)
                    trial_load = clamp(trial_load, lower_load + bracket_eps, upper_load - bracket_eps)

                    _, trial_gap, trial_rel_gap = evaluate_added_load(trial_load, iter)
                    RA_metric_gap = trial_gap
                    relative_RA_metric_gap = trial_rel_gap
                    if trial_rel_gap < min_relative_gap
                        min_relative_gap = trial_rel_gap
                        best_constant_load = trial_load
                    end
                    if abs(trial_gap) <= abs_tol || trial_rel_gap < tol
                        converged = true
                        break
                    end
                    if trial_gap < 0.0
                        lower_load = trial_load
                        lower_gap = trial_gap
                        if last_side == :lower
                            w_upper *= 0.5      # Illinois: demote stuck upper endpoint
                        else
                            w_upper = 1.0
                        end
                        last_side = :lower
                    else
                        upper_load = trial_load
                        upper_gap = trial_gap
                        if last_side == :upper
                            w_lower *= 0.5
                        else
                            w_lower = 1.0
                        end
                        last_side = :upper
                    end
                    iter += 1
                end
            elseif !converged && upper_gap < 0.0
                @aleaf_warn "[ALEAF RA Model]\t ELCC Process: failed to bracket target metric within $max_iter iterations."
            end

            ELCC_Info["ELCC"] = best_constant_load / ELCC_Info["reference_unit_ICAP"]
            ELCC_Info["equivalent_load_MW"] = best_constant_load * metric_scale_mw
            ELCC_Info["status"] = converged ? "converged" : "best_effort"

            if ELCC_Info["status"] == "converged"
                if abs(RA_metric_gap) <= abs_tol
                    @aleaf_info "[ALEAF RA Model]\t ELCC Process: convergence achieved with absolute gap=$(round(abs(RA_metric_gap), digits=3)) <= abs_tol=$(abs_tol)."
                elseif relative_RA_metric_gap < tol
                    @aleaf_info "[ALEAF RA Model]\t ELCC Process: convergence achieved with relative gap=$(round(relative_RA_metric_gap, digits=5)) < rel_tol=$(tol)."
                end
                @aleaf_info "[ALEAF RA Model]\t ELCC Process Completed. ELCC = $(ELCC_Info["ELCC"] * 100)%"
            else
                @aleaf_info "[ALEAF RA Model]\t ELCC bracketed search did not converge in $max_iter iterations. Use best added load to calculate ELCC = $(ELCC_Info["ELCC"] * 100)%"
            end

        else

            ELCC_Info["ELCC"] = 0.0
            ELCC_Info["status"] = "no_improvement_after_addition"
            @aleaf_info "[ALEAF RA Model]\t No RA improvement after adding $(ELCC_Info["reference_unit_ICAP"]*metric_scale_mw) MW of $new_resource_unit_group at bus $new_resource_location,\t"

        end

        # remove unnecessary data from the RA_Info dict
        delete!(ELCC_Info, "new_gen_risk_data")
        return ELCC_Info
    finally
        # Ensure temporary ELCC resources are always removed, even on errors.
        if haskey(ALEAF_model_instance.ref[:nw][0][:RA_gen_index], new_gen_idx)
            delete!(ALEAF_model_instance.ref[:nw][0][:RA_gen_index], new_gen_idx)
        end
        if haskey(ALEAF_model_instance.ref[:nw][new_resource_location][:RA_gen_bus], new_genco_tech_id)
            delete!(ALEAF_model_instance.ref[:nw][new_resource_location][:RA_gen_bus], new_genco_tech_id)
        end
        if original_local_gen_idx !== nothing
            ALEAF_model_instance.ref[:nw][0][:bus][new_resource_location]["aggregation_info"]["new_local_gen_idx"] = original_local_gen_idx
        end
        # restore gen index
        ALEAF_model_instance.ref[:nw][0][:gen_index] = original_gen_index
    end
end


function perform_ELCC_analysis_existing_asset_ED(ALEAF_model_instance, RA_Info, simulation_setting, risk_scenario_dict)

    ELCC_RA_Info = Dict{String, Any}()
    ELCC_RA_Info["setting"] = RA_Info["setting"]
    ELCC_RA_Info["RA_metrics"] = deepcopy(RA_Info["RA_metrics"])
    ELCC_RA_Info["system_info"] = Dict{String, Any}()
    ELCC_RA_Info["RA_reference_risk_data"] = get(RA_Info, "RA_reference_risk_data", Dict{String, Any}())

    RA_method = ELCC_RA_Info["setting"]["RA_method"]
    reference_RA_metric = simulation_setting["reference_RA_metric"]
    reference_RA_metric_spatial_resolution = simulation_setting["reference_RA_metric_spatial_resolution"]
    target_plant_name = string(get(simulation_setting, "deactivated_asset_plant_name", ""))
    target_unit_group = lowercase(string(get(simulation_setting, "deactivated_asset_unit_group", "")))
    metric_scale_mw = ALEAF_model_instance.setting["Simulation Setting"]["per_unit_base_value"]
    @aleaf_info "[ALEAF RA Model]\t ELCC Existing-Asset Process for (plant=$target_plant_name, unit_group=$target_unit_group) using RA method: $RA_method"
    
    ELCC_Info = Dict{String, Any}()
    ELCC_Info["assessment_type"] = "deactivate_existing"
    ELCC_Info["target_PLANT_NAME"] = target_plant_name
    ELCC_Info["target_UNIT_GROUP"] = target_unit_group
    ELCC_Info["deactivated_assets"] = Int[]
    ELCC_Info["removed_ICAP_MW"] = 0.0

    # identify existing assets
    target_gen_indices = Int[]
    for gen_idx in keys(ALEAF_model_instance.ref[:nw][0][:RA_gen_index])
        gen_data = ALEAF_model_instance.ref[:nw][0][:RA_gen_index][gen_idx]
        
        bus_idx = gen_data["bus_idx"]
        tech_idx = gen_data["genco_tech_id"]
        
        gen_bus_data = ALEAF_model_instance.ref[:nw][bus_idx][:RA_gen_bus][tech_idx]
        plant_name = get(gen_bus_data, "PLANT_NAME", get(gen_data, "PLANT_NAME", ""))
        unit_group = gen_data["UNIT_GROUP"]
        
        if (lowercase(string(unit_group)) == target_unit_group) && (string(plant_name) == target_plant_name)
            push!(target_gen_indices, gen_idx)
        end
    end
    
    if isempty(target_gen_indices)
        ELCC_Info["ELCC"] = 0.0
        ELCC_Info["status"] = "target_asset_not_found"
        @aleaf_warn "[ALEAF RA Model]\t ELCC Existing-Asset Process: no matching assets found for plant=$target_plant_name unit_group=$target_unit_group"
        return ELCC_Info
    end

    # identify bus names of the target assets for regional metric lookup
    target_bus_idx_list = sort!(unique([
        ALEAF_model_instance.ref[:nw][0][:RA_gen_index][idx]["bus_idx"] for idx in target_gen_indices
    ]))
    target_bus_names = sort!(unique([
        string(ALEAF_model_instance.ref[:nw][0][:RA_gen_index][idx]["bus_name"]) for idx in target_gen_indices
    ]))
    ELCC_Info["target_bus_indices"] = target_bus_names

    # reference RA metric
    ref_RA_metric = 0.0
    regional_metric_bus_idx = nothing
    if reference_RA_metric_spatial_resolution == "Systemwide"
        ref_RA_metric = ELCC_RA_Info["RA_metrics"][reference_RA_metric]["Annual"]["Systemwide"]

    elseif reference_RA_metric_spatial_resolution == "Regional"
        
        if length(target_bus_names) == 1
            regional_metric_bus_idx = first(target_bus_names)
            ref_RA_metric = ELCC_RA_Info["RA_metrics"][reference_RA_metric]["Annual"]["Regional"][string(regional_metric_bus_idx)]
        else
            ELCC_Info["ELCC"] = 0.0
            ELCC_Info["status"] = "regional_metric_bus_ambiguous"
            @aleaf_warn "[ALEAF RA Model]\t ELCC Existing-Asset Process: Regional metric requires one bus, but matched assets span buses=$(target_bus_names). Use Systemwide metric or narrow the target."
            return ELCC_Info
        end
    end

    # backup and deactivate selected assets (no index deletion)
    backup_gen_index = Dict{Int, Dict{String, Any}}()
    backup_gen_bus = Dict{Tuple{Int, Int}, Dict{String, Any}}()
    removed_icap = 0.0
    bus_tag = length(target_bus_names) == 1 ? string("bus_", first(target_bus_names)) : "matched_buses"
    output_path = string(ALEAF_model_instance.ref[:nw][0][:output_path], "ELCC_existing_asset_", replace(target_plant_name, " " => "_"), "_", bus_tag, "/")
    ALEAF.check_and_create_path(output_path)
    RA_Info["setting"]["ELCC_output_path"] = output_path

    try
        for gen_idx in target_gen_indices
            gen_data = ALEAF_model_instance.ref[:nw][0][:RA_gen_index][gen_idx]
            bus_idx = gen_data["bus_idx"]
            tech_idx = gen_data["genco_tech_id"]
            gen_bus_data = ALEAF_model_instance.ref[:nw][bus_idx][:RA_gen_bus][tech_idx]

            index_fields_to_zero = [
                "CAP", "ICAP", "UCAP", "U_G_iy", "EXUNITS", "EXCAPS"
            ]
            bus_fields_to_zero = [
                "CAP", "ICAP", "UCAP", "U_G_iy", "EXUNITS", "EXCAPS",
                "Charge_CAP", "ES_MWh", "new_u_ESE_iy", "Existing_ES_MWh"
            ]
            backup_gen_index[gen_idx] = Dict{String, Any}(
                field => gen_data[field] for field in index_fields_to_zero if haskey(gen_data, field)
            )
            backup_gen_bus[(bus_idx, tech_idx)] = Dict{String, Any}(
                field => gen_bus_data[field] for field in bus_fields_to_zero if haskey(gen_bus_data, field)
            )

            removed_icap += Float64(get(gen_data, "ICAP", get(gen_data, "CAP", 0.0)))
            for field in index_fields_to_zero
                if haskey(gen_data, field)
                    gen_data[field] = 0.0
                end
            end
            for field in bus_fields_to_zero
                if haskey(gen_bus_data, field)
                    gen_bus_data[field] = 0.0
                end
            end

            push!(ELCC_Info["deactivated_assets"], gen_idx)
        end

        ELCC_Info["removed_ICAP_MW"] = removed_icap * metric_scale_mw
        @aleaf_info "[ALEAF RA Model]\t ELCC Process: iteration=0, target_ref_metric=$(round(ref_RA_metric, digits=3)), removed_unit_icap=$(round(ELCC_Info["removed_ICAP_MW"], digits=3)) MW"

        # run RA after deactivation
        filtered_list = RA_Info["RA_reference_risk_data"]["filtered_joint_scenario_map"]
        constant_load_distribution = reference_RA_metric_spatial_resolution == "Systemwide" ? "load_weighted_systemwide" : "single_bus"
        constant_load_distribution_profile = Dict{Int, Float64}()
        load_relief_bus_idx = constant_load_distribution == "single_bus" && !isempty(target_bus_idx_list) ? first(target_bus_idx_list) : nothing
        @aleaf_info "[ALEAF RA Model]\t ELCC load adjustment mode: $constant_load_distribution"

        ELCC_RA_Info["RA_solutions"] = perform_ED_analysis_for_RA(
            ALEAF_model_instance.ref, ALEAF_model_instance.setting, RA_Info["setting"],
            risk_scenario_dict, filtered_list; constant_load_bus_idx=load_relief_bus_idx, constant_load_distribution, constant_load_distribution_profile, report_log=false, force_ELCC_flag=true
        )
        ELCC_RA_Info = calculate_risk_metrics_ED!(ALEAF_model_instance, ELCC_RA_Info)

        deactivated_RA_metric = 0.0
        if reference_RA_metric_spatial_resolution == "Systemwide"
            deactivated_RA_metric = ELCC_RA_Info["RA_metrics"][reference_RA_metric]["Annual"]["Systemwide"]
        else
            deactivated_RA_metric = ELCC_RA_Info["RA_metrics"][reference_RA_metric]["Annual"]["Regional"][string(regional_metric_bus_idx)]
        end

        RA_degradation = deactivated_RA_metric - ref_RA_metric
        @aleaf_info "[ALEAF RA Model]\t ELCC Process: iteration=1, new_metric=$(round(deactivated_RA_metric, digits=3)), target_ref_metric=$(round(ref_RA_metric, digits=3)), degradation_vs_target=$(round(RA_degradation, digits=3))"
        if RA_degradation <= 0 || removed_icap <= 1e-9
            ELCC_Info["ELCC"] = 0.0
            ELCC_Info["equivalent_load_relief_MW"] = 0.0
            ELCC_Info["status"] = "no_degradation_after_deactivation"
            @aleaf_info "[ALEAF RA Model]\t No RA degradation after deactivation; ELCC = 0.0%"
            return ELCC_Info
        end

        # Bracketed search on load reduction (negative constant_load) to recover the
        # reference RA metric; ED/storage response is too non-smooth for Newton updates.
        best_relief = 0.0
        min_relative_gap = Inf
        iter = 2
        tol = RA_Info["setting"]["capacity_credit_rel_tol_value"] * 0.01
        abs_tol = parse(Float64, string(RA_Info["setting"]["capacity_credit_abs_tol_value"]))
        max_iter = RA_Info["setting"]["capacity_credit_max_iteration_value"]
        RA_metric_gap = RA_degradation
        relative_RA_metric_gap = abs(ref_RA_metric) > 1e-9 ? abs(RA_metric_gap) / abs(ref_RA_metric) : abs(RA_metric_gap)
        converged = false

        function evaluate_load_relief(relief, iter)
            constant_load = -max(0.0, relief)
            eval_RA_Info = ELCC_RA_Info
            eval_RA_Info["RA_solutions"] = perform_ED_analysis_for_RA(
                ALEAF_model_instance.ref, ALEAF_model_instance.setting, RA_Info["setting"],
                risk_scenario_dict, filtered_list; constant_load, constant_load_bus_idx=load_relief_bus_idx, constant_load_distribution, constant_load_distribution_profile, report_log=false, iter, force_ELCC_flag=true
            )
            eval_RA_Info = calculate_risk_metrics_ED!(ALEAF_model_instance, eval_RA_Info)

            new_RA_metric = 0.0
            if reference_RA_metric_spatial_resolution == "Systemwide"
                new_RA_metric = eval_RA_Info["RA_metrics"][reference_RA_metric]["Annual"]["Systemwide"]
            else
                new_RA_metric = eval_RA_Info["RA_metrics"][reference_RA_metric]["Annual"]["Regional"][string(regional_metric_bus_idx)]
            end

            metric_gap = new_RA_metric - ref_RA_metric
            if abs(ref_RA_metric) > 1e-9
                rel_gap = abs(metric_gap) / abs(ref_RA_metric)
            else
                rel_gap = abs(metric_gap)
            end

            @aleaf_info "[ALEAF RA Model]\t ELCC Process: iteration=$iter, added_load_relief=$(round(relief*metric_scale_mw, digits=3)) MW, new_metric=$(round(new_RA_metric, digits=3)), target_ref_metric=$(round(ref_RA_metric, digits=3)), abs_gap_to_target=$(round(metric_gap, digits=3)), rel_gap_to_target=$(round(rel_gap, digits=5))"
            return new_RA_metric, metric_gap, rel_gap
        end

        lower_relief = 0.0
        lower_gap = RA_degradation                   # gap at relief=0 (from deactivation RA at iter=1)
        upper_relief = removed_icap
        _, upper_gap, upper_rel_gap = evaluate_load_relief(upper_relief, iter)
        RA_metric_gap = upper_gap
        relative_RA_metric_gap = upper_rel_gap
        best_relief = upper_relief
        min_relative_gap = upper_rel_gap

        if abs(upper_gap) <= abs_tol || upper_rel_gap < tol
            converged = true
        else
            iter += 1
        end

        while !converged && upper_gap > 0.0 && iter <= max_iter
            lower_relief = upper_relief
            lower_gap = upper_gap                    # track gap at new lower for regula falsi
            upper_relief *= 2.0
            _, upper_gap, upper_rel_gap = evaluate_load_relief(upper_relief, iter)
            RA_metric_gap = upper_gap
            relative_RA_metric_gap = upper_rel_gap
            if upper_rel_gap < min_relative_gap
                min_relative_gap = upper_rel_gap
                best_relief = upper_relief
            end
            if abs(upper_gap) <= abs_tol || upper_rel_gap < tol
                converged = true
                break
            end
            iter += 1
        end

        if !converged && upper_gap <= 0.0
            # Illinois-protected regula falsi: interpolate between bracket endpoints by
            # gap (near-linear response); faster than bisection, weights demote a stuck endpoint.
            w_lower = 1.0
            w_upper = 1.0
            last_side = :none
            while iter <= max_iter
                denom = w_upper * lower_gap - w_lower * upper_gap
                if abs(denom) < 1e-12
                    trial_relief = 0.5 * (lower_relief + upper_relief)   # degenerate: fall back to bisection
                else
                    trial_relief = (w_upper * lower_gap * upper_relief - w_lower * upper_gap * lower_relief) / denom
                end
                # Keep trial strictly inside the bracket (defensive against numerical drift)
                bracket_eps = 1e-9 * max(1.0, upper_relief - lower_relief)
                trial_relief = clamp(trial_relief, lower_relief + bracket_eps, upper_relief - bracket_eps)

                _, trial_gap, trial_rel_gap = evaluate_load_relief(trial_relief, iter)
                RA_metric_gap = trial_gap
                relative_RA_metric_gap = trial_rel_gap
                if trial_rel_gap < min_relative_gap
                    min_relative_gap = trial_rel_gap
                    best_relief = trial_relief
                end
                if abs(trial_gap) <= abs_tol || trial_rel_gap < tol
                    converged = true
                    break
                end
                if trial_gap > 0.0
                    lower_relief = trial_relief
                    lower_gap = trial_gap
                    if last_side == :lower
                        w_upper *= 0.5      # Illinois: demote stuck upper endpoint
                    else
                        w_upper = 1.0       # reset weight on side switch
                    end
                    last_side = :lower
                else
                    upper_relief = trial_relief
                    upper_gap = trial_gap
                    if last_side == :upper
                        w_lower *= 0.5
                    else
                        w_lower = 1.0
                    end
                    last_side = :upper
                end
                iter += 1
            end
        elseif !converged && upper_gap > 0.0
            @aleaf_warn "[ALEAF RA Model]\t ELCC Process: failed to bracket target metric within $max_iter iterations."
        end

        equivalent_load_relief = best_relief
        ELCC_Info["equivalent_load_relief_MW"] = equivalent_load_relief * metric_scale_mw
        ELCC_Info["ELCC"] = removed_icap > 0 ? equivalent_load_relief / removed_icap : 0.0
        ELCC_Info["status"] = converged ? "converged" : "best_effort"

        if ELCC_Info["status"] == "converged"
            @aleaf_info "[ALEAF RA Model]\t ELCC Process Completed. ELCC = $(ELCC_Info["ELCC"] * 100)%"
        else
            @aleaf_info "[ALEAF RA Model]\t ELCC bracketed search did not converge in $max_iter iterations. Use best load relief to calculate ELCC = $(ELCC_Info["ELCC"] * 100)%"
        end
        return ELCC_Info
    finally
        # restore all deactivated assets
        for gen_idx in target_gen_indices
            if haskey(backup_gen_index, gen_idx) && haskey(ALEAF_model_instance.ref[:nw][0][:RA_gen_index], gen_idx)
                for (field, value) in backup_gen_index[gen_idx]
                    ALEAF_model_instance.ref[:nw][0][:RA_gen_index][gen_idx][field] = value
                end
            end
            if haskey(ALEAF_model_instance.ref[:nw][0][:RA_gen_index], gen_idx)
                bus_idx = ALEAF_model_instance.ref[:nw][0][:RA_gen_index][gen_idx]["bus_idx"]
                tech_idx = ALEAF_model_instance.ref[:nw][0][:RA_gen_index][gen_idx]["genco_tech_id"]
                key = (bus_idx, tech_idx)
                if haskey(backup_gen_bus, key) && haskey(ALEAF_model_instance.ref[:nw][bus_idx][:RA_gen_bus], tech_idx)
                    for (field, value) in backup_gen_bus[key]
                        ALEAF_model_instance.ref[:nw][bus_idx][:RA_gen_bus][tech_idx][field] = value
                    end
                end
            end
        end
    end
end


const _RA_WORKER_CACHE = Dict{Symbol, Any}()

function _compute_ra_static_cache_tag(am_reference_data, am_setting_data, RA_setting)
    output_path = ""
    if haskey(am_reference_data, :nw) && haskey(am_reference_data[:nw], 0) && haskey(am_reference_data[:nw][0], :output_path)
        output_path = string(am_reference_data[:nw][0][:output_path])
    end
    current_year = get(RA_setting, "current_year", "NA")
    ra_method = get(RA_setting, "RA_method", "NA")
    nday_groups = get(am_setting_data["Simulation Configuration"], "NDAY_Groups_RA", "NA")
    capacity_credit_type = get(RA_setting, "capacity_credit_type", "NA")
    capacity_credit_assessment_mode = get(RA_setting, "capacity_credit_assessment_mode", "NA")
    elcc_output_path = get(RA_setting, "ELCC_output_path", "NA")
    export_elcc_dispatch = get(RA_setting, "export_ELCC_dispatch_results_flag", "NA")
    export_baseline_dispatch = get(RA_setting, "export_baseline_dispatch_results_flag", "NA")
    return string(
        output_path, "|", current_year, "|", ra_method, "|", nday_groups, "|",
        capacity_credit_type, "|", capacity_credit_assessment_mode, "|",
        elcc_output_path, "|", export_elcc_dispatch, "|", export_baseline_dispatch
    )
end

function _emit_worker_log!(worker_log_ch, worker_id::Int, message::String)
    if worker_log_ch !== nothing
        try
            put!(worker_log_ch, (worker_id, message))
        catch
            # ignore logging channel failures
        end
    end
    return nothing
end

function _set_ra_worker_cache_path!(outage_cache_path::String)
    prior_path = get(_RA_WORKER_CACHE, :outage_cache_path, nothing)
    _RA_WORKER_CACHE[:outage_cache_path] = outage_cache_path
    if !haskey(_RA_WORKER_CACHE, :risk_day_index) || prior_path != outage_cache_path
        # Index maps are keyed by day_group only; they must be rebuilt when
        # switching to a different case/cache path.
        _RA_WORKER_CACHE[:risk_day_index] = Dict{Int, Dict{Int, Int}}()
    end
    return nothing
end

function _set_ra_reference_cache_path!(reference_cache_path::String)
    _RA_WORKER_CACHE[:reference_cache_path] = reference_cache_path
    _RA_WORKER_CACHE[:reference_day_cache] = Dict{Tuple{Int, String}, Any}()
    _RA_WORKER_CACHE[:reference_day_cache_order] = Tuple{Int, String}[]
    return nothing
end

function _atomic_jld2_write!(writer_fn::Function, file_path::String)
    dir_path = dirname(file_path)
    mkpath(dir_path)
    tmp_path = joinpath(dir_path, "." * basename(file_path) * ".tmp." * string(myid()) * "." * randstring(8))
    try
        jldopen(tmp_path, "w") do f
            writer_fn(f)
        end
        mv(tmp_path, file_path; force=true)
    catch err
        if isfile(tmp_path)
            rm(tmp_path; force=true)
        end
        rethrow(err)
    end
    return nothing
end

function _write_risk_scenario_cache!(outage_cache_path::String, risk_scenario_dict::Dict{Int64, Any}; prefix::String="day")
    mkpath(outage_cache_path)
    # Materialize iteration order for Threads.@threads; each iteration writes a unique file
    # and owns its `scenarios` BitArray, sharing no mutable state, so it is safe to parallelize.
    day_group_ids = collect(keys(risk_scenario_dict))
    Threads.@threads for idx in eachindex(day_group_ids)
        day_group_id = day_group_ids[idx]
        day_dict = risk_scenario_dict[day_group_id]
        risk_ids = sort(collect(keys(day_dict)))

        # Keep only scenarios that carry an outage matrix; filtered-out scenarios are often stored as [].
        valid_risk_ids = [rid for rid in risk_ids if day_dict[rid] isa AbstractMatrix]
        risk_id_to_idx = Dict{Int, Int}()

        if isempty(valid_risk_ids)
            file_path = joinpath(outage_cache_path, "$(prefix)_$(day_group_id).jld2")
            _atomic_jld2_write!(file_path) do f
                f["scenarios"] = falses(0, 0, 0)
                f["risk_id_to_idx"] = risk_id_to_idx
            end
            continue
        end

        sample = day_dict[valid_risk_ids[1]]
        num_hours, num_gens = size(sample)

        # Enforce consistent matrix size within the same day-group cache file.
        valid_risk_ids = [rid for rid in valid_risk_ids if size(day_dict[rid]) == (num_hours, num_gens)]
        num_scenarios = length(valid_risk_ids)
        scenarios = falses(num_hours, num_gens, num_scenarios)
        for (i, rid) in enumerate(valid_risk_ids)
            scenarios[:, :, i] = day_dict[rid]
            risk_id_to_idx[rid] = i
        end

        file_path = joinpath(outage_cache_path, "$(prefix)_$(day_group_id).jld2")
        _atomic_jld2_write!(file_path) do f
            f["scenarios"] = scenarios
            f["risk_id_to_idx"] = risk_id_to_idx
        end
    end
    return nothing
end

function _load_risk_scenario_from_cache(day_group_id::Int, risk_id::Int; newgen::Bool=false)
    outage_cache_path = _RA_WORKER_CACHE[:outage_cache_path]
    prefix = newgen ? "day_newgen" : "day"
    file_path = joinpath(outage_cache_path, "$(prefix)_$(day_group_id).jld2")
    if !isfile(file_path)
        error("Outage cache file not found: $file_path")
    end
    index_cache = _RA_WORKER_CACHE[:risk_day_index]
    idx_map = get(index_cache, day_group_id, nothing)
    if idx_map === nothing
        try
            jldopen(file_path, "r") do f
                idx_map = f["risk_id_to_idx"]
            end
        catch err
            error("Failed to read outage index cache: $file_path | $(sprint(showerror, err))")
        end
        index_cache[day_group_id] = idx_map
    end
    if !haskey(idx_map, risk_id)
        error("Outage scenario cache index missing risk_id=$risk_id for day_group_id=$day_group_id (file=$file_path).")
    end
    idx = idx_map[risk_id]
    t1 = time()
    try
        jldopen(file_path, "r") do f
            return f["scenarios"][:, :, idx]
        end
    catch err
        error("Failed to read outage scenario cache: $file_path (day_group_id=$day_group_id, risk_id=$risk_id) | $(sprint(showerror, err))")
    end
end

function _load_risk_scenarios_from_cache_batch(risk_index_list; newgen::Bool=false)
    risk_scenario = Dict{Tuple{Int64, Int64}, Any}()
    outage_cache_path = _RA_WORKER_CACHE[:outage_cache_path]
    prefix = newgen ? "day_newgen" : "day"
    index_cache = _RA_WORKER_CACHE[:risk_day_index]

    risk_ids_by_day = Dict{Int, Vector{Int}}()
    for risk_id in risk_index_list
        day_group_id, rid = risk_id
        if !haskey(risk_ids_by_day, day_group_id)
            risk_ids_by_day[day_group_id] = Int[]
        end
        push!(risk_ids_by_day[day_group_id], rid)
    end

    for (day_group_id, risk_ids) in risk_ids_by_day
        file_path = joinpath(outage_cache_path, "$(prefix)_$(day_group_id).jld2")
        if !isfile(file_path)
            error("Outage cache file not found: $file_path")
        end

        idx_map = get(index_cache, day_group_id, nothing)
        local scenarios
        if idx_map === nothing
            try
                jldopen(file_path, "r") do f
                    idx_map = f["risk_id_to_idx"]
                    scenarios = f["scenarios"]
                end
            catch err
                error("Failed to read outage scenario cache: $file_path (day_group_id=$day_group_id) | $(sprint(showerror, err))")
            end
            index_cache[day_group_id] = idx_map
        else
            try
                jldopen(file_path, "r") do f
                    scenarios = f["scenarios"]
                end
            catch err
                error("Failed to read outage scenario cache: $file_path (day_group_id=$day_group_id) | $(sprint(showerror, err))")
            end
        end

        for rid in risk_ids
            if !haskey(idx_map, rid)
                error("Outage scenario cache index missing risk_id=$rid for day_group_id=$day_group_id (file=$file_path).")
            end
            idx = idx_map[rid]
            risk_scenario[(day_group_id, rid)] = scenarios[:, :, idx]
        end
    end
    return risk_scenario
end

function _get_filtered_joint_scenario(filtered_joint_scenario_map, day_group_id::Int, joint_id::Int)
    if !haskey(filtered_joint_scenario_map, day_group_id) || !haskey(filtered_joint_scenario_map[day_group_id], joint_id)
        error("Filtered joint scenario is missing day_group_id=$day_group_id, joint_id=$joint_id.")
    end
    return filtered_joint_scenario_map[day_group_id][joint_id]
end

function _write_reference_ed_cache!(reference_cache_path::String, daily_reference_ED_solutions::Dict{Int64, Any})
    mkpath(reference_cache_path)
    # Flatten the nested dict so Threads.@threads can balance all (day, scenario) writes;
    # unique target paths plus per-call random tmp paths mean no races, safe to parallelize.
    write_tasks = Tuple{String, Any}[]
    for (day_group_id, scenario_solution_dict) in daily_reference_ED_solutions
        for (renewable_scenario_id, day_solution) in scenario_solution_dict
            file_path = joinpath(reference_cache_path, "day_$(day_group_id)__$(renewable_scenario_id).jld2")
            push!(write_tasks, (file_path, day_solution))
        end
    end
    Threads.@threads for i in eachindex(write_tasks)
        file_path, day_solution = write_tasks[i]
        _atomic_jld2_write!(file_path) do f
            f["reference_solution"] = day_solution
        end
    end
    return nothing
end

function _load_reference_ed_solution_from_cache(day_group_id::Int, renewable_scenario_id::String)
    reference_cache_path = _RA_WORKER_CACHE[:reference_cache_path]
    day_cache = _RA_WORKER_CACHE[:reference_day_cache]
    day_cache_order = get!(_RA_WORKER_CACHE, :reference_day_cache_order, Tuple{Int, String}[])
    cache_limit = max(1, Int(get(_RA_WORKER_CACHE, :reference_day_cache_size, 8)))
    cache_key = (day_group_id, renewable_scenario_id)
    if haskey(day_cache, cache_key)
        idx = findfirst(==(cache_key), day_cache_order)
        if idx !== nothing
            deleteat!(day_cache_order, idx)
        end
        push!(day_cache_order, cache_key)
        return day_cache[cache_key]
    end
    file_path = joinpath(reference_cache_path, "day_$(day_group_id)__$(renewable_scenario_id).jld2")
    if !isfile(file_path)
        error("Reference ED cache file not found: $file_path")
    end
    t0 = time()
    solution = try
        jldopen(file_path, "r") do f
            f["reference_solution"]
        end
    catch err
        error("Failed to read reference ED cache: $file_path (day_group_id=$day_group_id, renewable_scenario_id=$renewable_scenario_id) | $(sprint(showerror, err))")
    end
    day_cache[cache_key] = solution
    idx = findfirst(==(cache_key), day_cache_order)
    if idx !== nothing
        deleteat!(day_cache_order, idx)
    end
    push!(day_cache_order, cache_key)
    while length(day_cache_order) > cache_limit
        evict_key = popfirst!(day_cache_order)
        if haskey(day_cache, evict_key)
            delete!(day_cache, evict_key)
        end
    end
    # silent; caller emits compact cache-load summaries
    return solution
end

function _start_ra_worker_progress_logger(report_log::Bool, detailed_log::Bool, last_simulation_id::Int, total_num_risk_scenarios::Int)
    if !report_log
        return nothing, nothing
    end

    worker_log_ch = RemoteChannel(()->Channel{Tuple{Int, String}}(10000), myid())
    worker_log_task = @async begin
        total_batches = max(1, last_simulation_id)
        total_scenarios = max(1, total_num_risk_scenarios)
        completed_batches = 0
        completed_scenarios = 0
        seen_batch_done = Set{Int}()
        t_progress_start = time()
        next_simple_pct = 5.0
        last_simple_heartbeat = t_progress_start
        heartbeat_sec = 600.0
        while true
            wait_status = timedwait(() -> isready(worker_log_ch) || !isopen(worker_log_ch), heartbeat_sec; pollint=1.0)
            if wait_status == :timed_out
                if !detailed_log
                    now_t = time()
                    elapsed = now_t - t_progress_start
                    batch_rate = completed_batches / max(elapsed, 1e-6)
                    remaining_batches = max(0, total_batches - completed_batches)
                    eta_sec = remaining_batches / max(batch_rate, 1e-9)
                    global_pct = 100.0 * completed_batches / total_batches
                    @aleaf_info "[ALEAF RA ED Analysis]: progress: batches=$(completed_batches)/$(total_batches) ($(round(global_pct, digits=1))%), risk scenarios completed=$(completed_scenarios), eta=$(round(eta_sec / 60.0, digits=1))m"
                    last_simple_heartbeat = now_t
                end
                continue
            end

            if !isopen(worker_log_ch) && !isready(worker_log_ch)
                break
            end

            try
                worker_id, msg = take!(worker_log_ch)
                m_done = match(r"master_done sim=(\d+) scenarios=(\d+)", msg)
                m = match(r"sim=(\d+) progress=(\d+)% \((\d+)/(\d+)\) t=([0-9.]+)s", msg)
                if m_done !== nothing
                    sim_id = parse(Int, m_done.captures[1])
                    done = parse(Int, m_done.captures[2])
                    if !(sim_id in seen_batch_done)
                        push!(seen_batch_done, sim_id)
                        completed_batches += 1
                        completed_scenarios += done
                        now_t = time()
                        elapsed = now_t - t_progress_start
                        batch_rate = completed_batches / max(elapsed, 1e-6)
                        remaining_batches = max(0, total_batches - completed_batches)
                        eta_sec = remaining_batches / max(batch_rate, 1e-9)
                        global_pct = 100.0 * completed_batches / total_batches
                        if detailed_log || global_pct >= next_simple_pct || completed_batches == total_batches
                            @aleaf_info "[ALEAF RA ED Analysis]: progress: batches=$(completed_batches)/$(total_batches) ($(round(global_pct, digits=1))%), risk scenarios completed=$(completed_scenarios), eta=$(round(eta_sec / 60.0, digits=1))m"
                            if !detailed_log
                                while next_simple_pct <= global_pct
                                    next_simple_pct += 5.0
                                end
                                last_simple_heartbeat = now_t
                            end
                        end
                    end
                elseif m !== nothing
                    sim_id = parse(Int, m.captures[1])
                    pct = parse(Int, m.captures[2])
                    done = parse(Int, m.captures[3])
                    total = parse(Int, m.captures[4])

                    if detailed_log && (pct == 25 || pct == 50 || pct == 75 || pct == 100)
                        @aleaf_info "[ALEAF RA ED Analysis][worker $(worker_id)]: sim=$(sim_id) progress=$(pct)% ($(done)/$(total))"
                    end
                elseif detailed_log
                    @aleaf_info "[ALEAF RA ED Analysis][worker $worker_id]: $msg"
                end
            catch
                break
            end
        end
    end

    return worker_log_ch, worker_log_task
end

function _set_ra_worker_cache!(am_reference_data, am_setting_data, RA_setting, daily_reference_ED_solutions, risk_scenario_dict, risk_scenario_dict_of_new_gen)
    _RA_WORKER_CACHE[:am_reference_data] = am_reference_data
    _RA_WORKER_CACHE[:am_setting_data] = am_setting_data
    _RA_WORKER_CACHE[:RA_setting] = RA_setting
    _RA_WORKER_CACHE[:static_cache_tag] = _compute_ra_static_cache_tag(am_reference_data, am_setting_data, RA_setting)
    _RA_WORKER_CACHE[:reference_day_cache_size] = max(1, Int(get(RA_setting, "reference_day_cache_size_value", 8)))
    if daily_reference_ED_solutions !== nothing
        _RA_WORKER_CACHE[:daily_reference_ED_solutions] = daily_reference_ED_solutions
    elseif haskey(_RA_WORKER_CACHE, :daily_reference_ED_solutions)
        delete!(_RA_WORKER_CACHE, :daily_reference_ED_solutions)
    end
    if risk_scenario_dict !== nothing
        _RA_WORKER_CACHE[:risk_scenario_dict] = risk_scenario_dict
    end
    if risk_scenario_dict_of_new_gen !== nothing
        _RA_WORKER_CACHE[:risk_scenario_dict_of_new_gen] = risk_scenario_dict_of_new_gen
    end
    return nothing
end

function _get_required_ra_worker_cache!(keys_needed::Vector{Symbol})
    missing = Symbol[]
    for k in keys_needed
        if !haskey(_RA_WORKER_CACHE, k)
            push!(missing, k)
        end
    end
    if !isempty(missing)
        error("Missing worker cache keys: $(missing)")
    end
    return nothing
end

function perform_ED_analysis_for_RA(am_reference_data, am_setting_data, RA_setting, risk_scenario_dict, filtered_joint_scenario_map; risk_scenario_dict_of_new_gen = [], constant_load = 0.0, constant_load_bus_idx=nothing, constant_load_distribution="single_bus", constant_load_distribution_profile=Dict{Tuple{Int, Int, Int}, Dict{Int, Float64}}(), report_log=true, iter=1, force_ELCC_flag=false)
    
    # ALEAF_model_instance build-------------------------------------------------
    detailed_log = lowercase(string(get(RA_setting, "logging_level_value", "simple"))) == "detailed"

    # Check ELCC flag
    ELCC_flag = force_ELCC_flag || (length(risk_scenario_dict_of_new_gen) > 0) || (abs(constant_load) > 1e-12)
    run_DLOL_flag = get(RA_setting, "calculate_capacity_credit_flag", false) && (string(get(RA_setting, "capacity_credit_type", "")) == "DLOL")

    # During ELCC iterations use the dedicated simulation method so the user can run
    # a cheaper Sequential ED for RA but a full PF method for capacity-credit evaluation.
    effective_RA_method = (ELCC_flag && haskey(RA_setting, "capacity_credit_RA_simulation_method")) ?
        RA_setting["capacity_credit_RA_simulation_method"] : RA_setting["RA_method"]

    # initialize
    RA_solution = Dict{String, Any}()
    for day_id in [d for d in 1:am_setting_data["Simulation Configuration"]["NDAY_Groups_RA"]]
        RA_solution[string(day_id)] = Dict{Int64, Any}()
    end
   
    # simulation set
    # Keep only joint IDs whose outage scenario maps to an outage matrix in each day-group.
    function valid_joint_ids_for_day(day_group_id::Int)
        valid_ids = Int[]
        if !haskey(filtered_joint_scenario_map, day_group_id) || !haskey(risk_scenario_dict, day_group_id)
            return valid_ids
        end
        for joint_id in keys(filtered_joint_scenario_map[day_group_id])
            outage_risk_id = filtered_joint_scenario_map[day_group_id][joint_id]["outage_risk_id"]
            if haskey(risk_scenario_dict[day_group_id], outage_risk_id) && (risk_scenario_dict[day_group_id][outage_risk_id] isa AbstractMatrix)
                push!(valid_ids, joint_id)
            end
        end
        return valid_ids
    end

    valid_joint_scenario_list = Dict{Int, Vector{Int}}()
    for day_group_id in keys(filtered_joint_scenario_map)
        valid_joint_scenario_list[day_group_id] = valid_joint_ids_for_day(day_group_id)
    end
    
    total_num_risk_scenarios = sum(length(risks) for risks in values(valid_joint_scenario_list))
    if total_num_risk_scenarios == 0
        return RA_solution
    end
    RA_setting["filtered_joint_scenario_map"] = filtered_joint_scenario_map

    simulation_set = Dict{Int64, Any}()
    simulation_set_solution = Dict{Int64, Any}()
    DLOL_unit_accumulator_global = Dict{Int, Dict{Symbol, Any}}()
            
    # identify day groups that have risk scenarios
    day_group_list = Int[]
    for day_group_id in keys(valid_joint_scenario_list)
        if isempty(valid_joint_scenario_list[day_group_id]) == false
            push!(day_group_list, day_group_id)
        end
    end

    # Run base case operations for the identified day groups
    daily_reference_ED_solutions = perform_reference_ED_simulation_for_RA(filtered_joint_scenario_map, am_reference_data, am_setting_data, RA_setting; constant_load, constant_load_bus_idx, constant_load_distribution, constant_load_distribution_profile, ELCC_flag)
    if report_log @aleaf_info "[ALEAF RA ED Analysis]: Obtained reference ED solutions" end
    if haskey(am_reference_data[:nw][0], :output_path)
        reference_cache_path = joinpath(am_reference_data[:nw][0][:output_path], "RA_reference_cache")
        _write_reference_ed_cache!(reference_cache_path, daily_reference_ED_solutions)
        RA_setting["reference_cache_path"] = reference_cache_path
        if report_log
            @aleaf_info "[ALEAF RA ED Analysis]: Reference ED cache written to $reference_cache_path"
        end
    end

    # cache large static data on workers to reduce communication (static only)
    worker_list = workers()
    if isempty(worker_list)
        error("No distributed workers are available for RA ED analysis.")
    end
    
    # Build one global pair list; dynamic queue handles balancing at runtime.
    all_pairs = Tuple{Int, Int}[]
    for day_group_id in sort(collect(day_group_list))
        for joint_id in sort(valid_joint_scenario_list[day_group_id])
            push!(all_pairs, (day_group_id, joint_id))
        end
    end

    # push static cache to workers (reference ED solutions are loaded from disk cache on-demand)
    static_cache_tag = _compute_ra_static_cache_tag(am_reference_data, am_setting_data, RA_setting)
    _set_ra_worker_cache!(am_reference_data, am_setting_data, RA_setting, nothing, nothing, nothing)
    if haskey(RA_setting, "outage_cache_path")
        _set_ra_worker_cache_path!(RA_setting["outage_cache_path"])
    end
    if haskey(RA_setting, "reference_cache_path")
        _set_ra_reference_cache_path!(RA_setting["reference_cache_path"])
    end
    @sync for p in worker_list
        if p != myid()
            @async begin
                needs_refresh = !remotecall_fetch(p, static_cache_tag) do expected_tag
                    get(_RA_WORKER_CACHE, :static_cache_tag, nothing) == expected_tag
                end
                if needs_refresh
                    remotecall_fetch(p, am_reference_data, am_setting_data, RA_setting) do am_ref, am_set, ra_set
                        _set_ra_worker_cache!(am_ref, am_set, ra_set, nothing, nothing, nothing)
                        nothing
                    end
                end
                remotecall_fetch(p, RA_setting) do ra_set
                    if haskey(_RA_WORKER_CACHE, :daily_reference_ED_solutions)
                        delete!(_RA_WORKER_CACHE, :daily_reference_ED_solutions)
                    end
                    if haskey(ra_set, "outage_cache_path")
                        _set_ra_worker_cache_path!(ra_set["outage_cache_path"])
                    end
                    if haskey(ra_set, "reference_cache_path")
                        _set_ra_reference_cache_path!(ra_set["reference_cache_path"])
                    end
                    _RA_WORKER_CACHE[:reference_day_cache_size] = max(1, Int(get(ra_set, "reference_day_cache_size_value", 8)))
                    nothing
                end
            end
        end
    end

    # Adaptive dynamic memory-probing batch sizing.
    configured_batch_size = max(1, Int(get(RA_setting, "num_distributed_scenarios_per_worker_value", 200)))
    worker_batch_size = Dict{Int, Int}()
    sample_scenario_bytes = 1
    if !isempty(day_group_list)
        sample_day_group = first(day_group_list)
        sample_joint_id = first(valid_joint_scenario_list[sample_day_group])
        sample_risk_id = filtered_joint_scenario_map[sample_day_group][sample_joint_id]["outage_risk_id"]
        try
            if haskey(RA_setting, "outage_cache_path")
                sample_scenario = _load_risk_scenario_from_cache(sample_day_group, sample_risk_id; newgen=false)
                sample_scenario_bytes = max(1, Base.summarysize(sample_scenario))
            end
        catch
            # fallback below
        end
        if sample_scenario_bytes == 1 && haskey(risk_scenario_dict, sample_day_group) && haskey(risk_scenario_dict[sample_day_group], sample_risk_id)
            sample_scenario_bytes = max(1, Base.summarysize(risk_scenario_dict[sample_day_group][sample_risk_id]))
        end
    end

    # Check representative reference day size for memory budgeting.
    representative_ref_day_bytes = 0
    if !isempty(day_group_list)
        rep_day = first(day_group_list)
        if haskey(daily_reference_ED_solutions, rep_day)
            representative_ref_day_bytes = Base.summarysize(daily_reference_ED_solutions[rep_day])
        end
    end

    unique_day_count_all_pairs = length(unique(p[1] for p in all_pairs))
    unique_day_count_per_worker_est = max(1, ceil(Int, unique_day_count_all_pairs / max(1, length(worker_list))))
    mem_headroom_ratio = 0.30
    for p in worker_list
        free_mem_p = try
            p == myid() ? Sys.free_memory() : remotecall_fetch(() -> Sys.free_memory(), p)
        catch
            Sys.free_memory()
        end
        available_budget = Int(floor(free_mem_p * (1.0 - mem_headroom_ratio)))
        ref_cache_bytes = representative_ref_day_bytes * unique_day_count_per_worker_est
        worker_budget = max(1, available_budget - ref_cache_bytes)
        adaptive_batch = max(1, Int(floor(worker_budget / sample_scenario_bytes)))
        worker_batch_size[p] = min(configured_batch_size, adaptive_batch)
    end
    avg_batch_size = max(1, Int(floor(sum(values(worker_batch_size)) / length(worker_batch_size))))
    
    # Keep enough chunks to utilize workers; otherwise a high configured limit can collapse
    # almost all work into 1-2 batches and starve parallelism.
    target_chunks = max(1, min(length(worker_list), length(all_pairs)))
    worker_balanced_limit = max(1, ceil(Int, length(all_pairs) / target_chunks))
    effective_batch_size = min(configured_batch_size, avg_batch_size, worker_balanced_limit)
    if report_log && detailed_log
        free_mem_gb = round(Sys.free_memory() / 1024^3, digits=1)
        @aleaf_info "[ALEAF RA ED Analysis]: Adaptive batch sizing (configured_limit=$(configured_batch_size), adaptive_avg=$(avg_batch_size), worker_balanced_limit=$(worker_balanced_limit), effective_batch_size=$(effective_batch_size), sample_scenario_bytes=$(sample_scenario_bytes), free_mem_gb=$(free_mem_gb), headroom=$(mem_headroom_ratio))"
    end

    simulation_idx = 1
    pair_idx = 1
    while pair_idx <= length(all_pairs)
        simulation_set[simulation_idx] = Tuple{Int, Int}[]
        while pair_idx <= length(all_pairs) && length(simulation_set[simulation_idx]) < effective_batch_size
            day_group_id, joint_id = all_pairs[pair_idx]
            push!(simulation_set[simulation_idx], (day_group_id, joint_id))
            RA_solution[string(day_group_id)][joint_id] = Dict()
            pair_idx += 1
        end
        simulation_idx += 1
    end

    last_simulation_id = simulation_idx - 1
    if report_log
        @aleaf_info "[ALEAF RA ED Analysis]: Total number of simulation batches: $last_simulation_id, Total number of workers: $(length(worker_list))"
    end

    # write new-gen outage cache if provided (so workers can load locally)
    if length(risk_scenario_dict_of_new_gen) > 0 && haskey(RA_setting, "outage_cache_path")
        _write_risk_scenario_cache!(RA_setting["outage_cache_path"], risk_scenario_dict_of_new_gen; prefix="day_newgen")
    end

    # central worker log stream for near real-time worker messages
    worker_log_ch, worker_log_task = _start_ra_worker_progress_logger(report_log, detailed_log, last_simulation_id, total_num_risk_scenarios)

    # submit batches through a dynamic queue to reduce tail latency
    dispatch_order = sort(collect(keys(simulation_set)), by = sid -> length(simulation_set[sid]), rev = true)
    simulation_id_ch = Channel{Int}(max(1, length(dispatch_order)))
    @async begin
        for simulation_id in dispatch_order
            put!(simulation_id_ch, simulation_id)
        end
        close(simulation_id_ch)
    end
    try
        @sync for p in worker_list
            @async begin
                for simulation_id in simulation_id_ch
                    try
                        if ELCC_flag == true
                            if effective_RA_method == "Economic Dispatch"
                                simulation_set_solution[simulation_id] = remotecall_fetch(run_ra_perfect_foresight_batch!, p, simulation_set[simulation_id], p, simulation_id; constant_load, constant_load_bus_idx, constant_load_distribution, constant_load_distribution_profile, new_gen_risk_scenario = nothing, report_log, iter, worker_log_ch, force_ELCC_flag = ELCC_flag)
                            elseif effective_RA_method == "Sequential Economic Dispatch"
                                simulation_set_solution[simulation_id] = remotecall_fetch(run_ra_imperfect_foresight_sequential_snapshot_batch!, p, simulation_set[simulation_id], p, simulation_id; constant_load, constant_load_bus_idx, constant_load_distribution, constant_load_distribution_profile, new_gen_risk_scenario = nothing, report_log, iter, worker_log_ch, force_ELCC_flag = ELCC_flag)
                            end
                        else
                            if effective_RA_method == "Economic Dispatch"
                                simulation_set_solution[simulation_id] = remotecall_fetch(run_ra_perfect_foresight_batch!, p, simulation_set[simulation_id], p, simulation_id; report_log, worker_log_ch)
                            elseif effective_RA_method == "Sequential Economic Dispatch"
                                simulation_set_solution[simulation_id] = remotecall_fetch(run_ra_imperfect_foresight_sequential_snapshot_batch!, p, simulation_set[simulation_id], p, simulation_id; report_log, worker_log_ch)
                            end
                        end
                        _emit_worker_log!(worker_log_ch, p, "master_done sim=$simulation_id scenarios=$(length(simulation_set[simulation_id]))")
                    catch err
                        err_msg = sprint(showerror, err)
                        @aleaf_error "ED analysis failed for simulation id $simulation_id on worker $p: $err_msg"
                        rethrow(err)
                    end
                end
            end
        end
    finally
        if worker_log_ch !== nothing
            close(worker_log_ch)
        end
        if worker_log_task !== nothing
            wait(worker_log_task)
        end
    end

    # Accept both payload formats: legacy scenario result dict, or DLOL wrapped dict
    # containing scenario_results + DLOL partials.
    for simulation_idx in keys(simulation_set_solution)
        worker_payload = simulation_set_solution[simulation_idx]
        scenario_results = worker_payload
        if (worker_payload isa AbstractDict) && haskey(worker_payload, "scenario_results")
            scenario_results = worker_payload["scenario_results"]
        end

        # Merge worker-level DLOL partial accumulator when available.
        if run_DLOL_flag
            _merge_DLOL_worker_partials!(DLOL_unit_accumulator_global, worker_payload)
        end

        for risk_index in keys(scenario_results)
            day_group_id = risk_index[1]
            joint_id = risk_index[2]
            RA_solution[string(day_group_id)][joint_id] = scenario_results[risk_index]
        end
    end

    # DLOL eligibility: unit group needs ELCC_Flag in Gen Technology; bus needs
    # RA_ELCC_Calculation_Flag for bus-level aggregation.
    eligible_unit_groups = Set{String}()
    for tech_idx in keys(am_reference_data[:nw][0][:gen_technology])
        tech_data = am_reference_data[:nw][0][:gen_technology][tech_idx]
        if get(tech_data, "ELCC_Flag", false) == true
            push!(eligible_unit_groups, string(get(tech_data, "UNITGROUP", "")))
        end
    end

    eligible_buses = Set{Int}()
    for (bus_idx, bus_data) in am_reference_data[:nw][0][:bus]
        if get(bus_data, "RA_ELCC_Calculation_Flag", false) == true
            push!(eligible_buses, Int(bus_idx))
        end
    end

    # Finalize DLOL outputs on master once all worker payloads are merged.
    if run_DLOL_flag
        RA_solution["DLOL"] = _finalize_DLOL_outputs!(DLOL_unit_accumulator_global, am_reference_data, RA_setting, ELCC_flag, eligible_unit_groups, eligible_buses; report_log)
    end

    return RA_solution
end


function _merge_DLOL_worker_partials!(DLOL_unit_accumulator_global::Dict{Int, Dict{Symbol, Any}}, worker_payload)
    if !((worker_payload isa AbstractDict) && haskey(worker_payload, "DLOL_unit_accumulator_partial"))
        return
    end

    partial = worker_payload["DLOL_unit_accumulator_partial"]
    for (unit_id, unit_data) in partial
        if !haskey(DLOL_unit_accumulator_global, unit_id)
            DLOL_unit_accumulator_global[unit_id] = Dict{Symbol, Any}(
                :plant_name => unit_data[:plant_name],
                :bus_idx => unit_data[:bus_idx],
                :unit_group => unit_data[:unit_group],
                :unit_category => unit_data[:unit_category],
                :ICAP_MW => unit_data[:ICAP_MW],
                :DLOL_sum_output_stress_weighted => 0.0,
                :DLOL_sum_stress_hour_weights => 0.0
            )
        end
        DLOL_unit_accumulator_global[unit_id][:DLOL_sum_output_stress_weighted] += unit_data[:DLOL_sum_output_stress_weighted]
        DLOL_unit_accumulator_global[unit_id][:DLOL_sum_stress_hour_weights] += unit_data[:DLOL_sum_stress_hour_weights]
    end
end


function _finalize_DLOL_outputs!(DLOL_unit_accumulator_global::Dict{Int, Dict{Symbol, Any}}, am_reference_data, RA_setting, ELCC_flag::Bool, eligible_unit_groups::Set{String}, eligible_buses::Set{Int}; report_log::Bool=true)
    unit_rows = Tuple{Int, String, Int, String, String, String, Float64, Float64, Float64, Float64}[]
    unit_json = Dict{String, Any}()
    bus_name_map = Dict{Int, String}()
    for (bus_idx, bus_data) in am_reference_data[:nw][0][:bus]
        bus_name = get(bus_data, "bus_i", string(bus_idx))
        bus_name_map[Int(bus_idx)] = isa(bus_name, String) ? bus_name : string(bus_name)
    end

    for unit_id in sort(collect(keys(DLOL_unit_accumulator_global)))
        u = DLOL_unit_accumulator_global[unit_id]
        unit_group = string(get(u, :unit_group, ""))
        if !(unit_group in eligible_unit_groups)
            continue
        end

        icap = Float64(get(u, :ICAP_MW, 0.0))
        sum_out = Float64(get(u, :DLOL_sum_output_stress_weighted, 0.0))
        sum_w = Float64(get(u, :DLOL_sum_stress_hour_weights, 0.0))
        avg_output_stress = sum_w > 0.0 ? (sum_out / sum_w) : 0.0
        dlol = icap > 0.0 ? (avg_output_stress / icap) : 0.0
        bus_idx = Int(get(u, :bus_idx, 0))
        bus_name = get(bus_name_map, bus_idx, string(bus_idx))
        plant_name = string(get(u, :plant_name, ""))
        unit_category = string(get(u, :unit_category, ""))

        push!(unit_rows, (unit_id, plant_name, bus_idx, bus_name, unit_group, unit_category, icap, avg_output_stress, dlol, sum_w))

        unit_json[string(unit_id)] = Dict{String, Any}(
            "plant_name" => plant_name,
            "bus_idx" => bus_idx,
            "bus_name" => bus_name,
            "UNIT_GROUP" => unit_group,
            "UNIT_CATEGORY" => unit_category,
            "ICAP_MW" => icap,
            "AvgOutputStress_MW" => avg_output_stress,
            "DLOL" => dlol,
            "stress_hour_weights" => sum_w
        )
    end

    unit_df = DataFrame(
        unit_id = Int[],
        plant_name = String[],
        bus_idx = Int[],
        bus_name = String[],
        UNIT_GROUP = String[],
        UNIT_CATEGORY = String[],
        ICAP_MW = Float64[],
        AvgOutputStress_MW = Float64[],
        DLOL = Float64[],
        stress_hour_weights = Float64[]
    )
    for row in unit_rows
        push!(unit_df, row)
    end

    # Bus-level summary by UNIT_GROUP (simple and ICAP-weighted averages)
    bus_group_stats = Dict{Tuple{String, String}, Dict{Symbol, Float64}}()
    for row in eachrow(unit_df)
        if !(Int(row.bus_idx) in eligible_buses)
            continue
        end
        key = (String(row.bus_name), String(row.UNIT_GROUP))
        if !haskey(bus_group_stats, key)
            bus_group_stats[key] = Dict{Symbol, Float64}(
                :unit_count => 0.0,
                :sum_dlol => 0.0,
                :sum_icap => 0.0,
                :sum_dlol_icap => 0.0
            )
        end
        bus_group_stats[key][:unit_count] += 1.0
        bus_group_stats[key][:sum_dlol] += Float64(row.DLOL)
        bus_group_stats[key][:sum_icap] += Float64(row.ICAP_MW)
        bus_group_stats[key][:sum_dlol_icap] += Float64(row.DLOL) * Float64(row.ICAP_MW)
    end

    bus_group_rows = Tuple{String, String, Int, Float64, Float64}[]
    bus_group_json = Dict{String, Any}()
    for ((bus_name, unit_group), st) in sort(collect(bus_group_stats), by = x -> x[1])
        unit_count = Int(round(st[:unit_count]))
        avg_dlol = st[:unit_count] > 0.0 ? (st[:sum_dlol] / st[:unit_count]) : 0.0
        weighted_avg_dlol = st[:sum_icap] > 0.0 ? (st[:sum_dlol_icap] / st[:sum_icap]) : 0.0
        push!(bus_group_rows, (bus_name, unit_group, unit_count, avg_dlol, weighted_avg_dlol))

        if !haskey(bus_group_json, bus_name)
            bus_group_json[bus_name] = Dict{String, Any}()
        end
        bus_group_json[bus_name][unit_group] = Dict{String, Any}(
            "unit_count" => unit_count,
            "Avg_DLOL" => avg_dlol,
            "ICAP_Weighted_Avg_DLOL" => weighted_avg_dlol
        )
    end
    bus_group_df = DataFrame(
        bus_name = String[],
        UNIT_GROUP = String[],
        unit_count = Int[],
        Avg_DLOL = Float64[],
        ICAP_Weighted_Avg_DLOL = Float64[]
    )
    for row in bus_group_rows
        push!(bus_group_df, row)
    end

    # System-level summary by UNIT_GROUP (simple and ICAP-weighted averages)
    system_group_stats = Dict{String, Dict{Symbol, Float64}}()
    for row in eachrow(unit_df)
        unit_group = String(row.UNIT_GROUP)
        if !haskey(system_group_stats, unit_group)
            system_group_stats[unit_group] = Dict{Symbol, Float64}(
                :unit_count => 0.0,
                :sum_dlol => 0.0,
                :sum_icap => 0.0,
                :sum_dlol_icap => 0.0
            )
        end
        system_group_stats[unit_group][:unit_count] += 1.0
        system_group_stats[unit_group][:sum_dlol] += Float64(row.DLOL)
        system_group_stats[unit_group][:sum_icap] += Float64(row.ICAP_MW)
        system_group_stats[unit_group][:sum_dlol_icap] += Float64(row.DLOL) * Float64(row.ICAP_MW)
    end

    system_group_rows = Tuple{String, Int, Float64, Float64}[]
    system_group_json = Dict{String, Any}()
    for (unit_group, st) in sort(collect(system_group_stats), by = x -> x[1])
        unit_count = Int(round(st[:unit_count]))
        avg_dlol = st[:unit_count] > 0.0 ? (st[:sum_dlol] / st[:unit_count]) : 0.0
        weighted_avg_dlol = st[:sum_icap] > 0.0 ? (st[:sum_dlol_icap] / st[:sum_icap]) : 0.0
        push!(system_group_rows, (unit_group, unit_count, avg_dlol, weighted_avg_dlol))
        system_group_json[unit_group] = Dict{String, Any}(
            "unit_count" => unit_count,
            "Avg_DLOL" => avg_dlol,
            "ICAP_Weighted_Avg_DLOL" => weighted_avg_dlol
        )
    end
    system_group_df = DataFrame(
        UNIT_GROUP = String[],
        unit_count = Int[],
        Avg_DLOL = Float64[],
        ICAP_Weighted_Avg_DLOL = Float64[]
    )
    for row in system_group_rows
        push!(system_group_df, row)
    end

    if report_log
        stressed_units = count(row -> row.stress_hour_weights > 0.0, eachrow(unit_df))
        @aleaf_info "[ALEAF RA DLOL]: Finalized DLOL (units_total=$(length(DLOL_unit_accumulator_global)), units_filtered=$(nrow(unit_df)), stressed_units=$stressed_units, eligible_buses=$(length(eligible_buses)))"
        if nrow(unit_df) == 0
            @aleaf_warn "[ALEAF RA DLOL]: No units passed ELCC_Flag filtering; DLOL outputs are empty."
        end
    end

    return Dict{String, Any}(
        "unit_level" => unit_json,
        "bus_level_UNIT_GROUP" => bus_group_json,
        "system_level_UNIT_GROUP" => system_group_json
    )
end



function _build_reference_ed_pair_list(filtered_joint_scenario_map)
    pair_set = Set{Tuple{Int, String}}()
    for day_group_id in sort(collect(keys(filtered_joint_scenario_map)))
        for joint_id in sort(collect(keys(filtered_joint_scenario_map[day_group_id])))
            renewable_scenario_id = filtered_joint_scenario_map[day_group_id][joint_id]["renewable_scenario_id"]
            push!(pair_set, (day_group_id, renewable_scenario_id))
        end
    end
    return sort!(collect(pair_set); by = x -> (x[1], x[2]))
end

function _build_ra_local_reference_for_scenario(am_reference_data, RA_setting, run_H, day_group_id::Int, renewable_scenario_id::String, y::Int, bus_keys)
    local_ref = copy(am_reference_data)
    local_ref[:nw] = copy(am_reference_data[:nw])
    local_ref[:nw][0] = copy(am_reference_data[:nw][0])
    for bus_idx in bus_keys
        local_ref[:nw][bus_idx] = copy(am_reference_data[:nw][bus_idx])
    end
    local_ref[:nw][0][:gen_index] = copy(am_reference_data[:nw][0][:RA_gen_index])
    for bus_idx in bus_keys
        local_ref[:nw][bus_idx][:gen_bus] = copy(am_reference_data[:nw][bus_idx][:RA_gen_bus])
    end

    local_ref[:nw][0][:planning_stages] = copy(am_reference_data[:nw][0][:planning_stages])
    local_ref[:nw][0][:planning_stages][y] = copy(am_reference_data[:nw][0][:planning_stages][y])
    local_ref[:nw][0][:planning_stages][y]["repdays"] = copy(am_reference_data[:nw][0][:planning_stages][y]["repdays"])
    repday_ids = local_ref[:nw][0][:repday_groups][day_group_id]["Day_Idx_List"]

    # Fresh dense zdt for exactly the (year, bus, repday, hour) this local build fills and later reads.
    d_pos = Dict{Int, Int}(rid => pos for (pos, rid) in enumerate(sort(collect(repday_ids))))
    h_pos = Dict{Int, Int}(h => pos for (pos, h) in enumerate(run_H))
    nd = length(d_pos)
    nh = length(h_pos)
    local_ref[:nw][0][:vre_zdt_index] = (d = d_pos, h = h_pos)
    local_ref[:nw][0][:vre_aggregated_data_zdt] = Dict{Int, Dict{Int, Array{Float64,3}}}()
    local_ref[:nw][0][:vre_aggregated_data_zdt][y] = Dict{Int, Array{Float64,3}}()
    for repday_id in repday_ids
        repday_key = string(repday_id)
        local_ref[:nw][0][:planning_stages][y]["repdays"][repday_key] = deepcopy(am_reference_data[:nw][0][:planning_stages][y]["repdays"][repday_key])
    end

    scenario_cache = RA_setting["renewable_scenario_daygroup_data"][day_group_id][renewable_scenario_id]
    row_idx = 1
    for repday_id in repday_ids
        repday_key = string(repday_id)
        for h in run_H
            hourly_data = local_ref[:nw][0][:planning_stages][y]["repdays"][repday_key]["data"][string(h)]["1"]

            if haskey(hourly_data, "wind_ons")
                for (timeseries_tag, hourly_profile) in scenario_cache["wind_ons"]
                    if haskey(hourly_data["wind_ons"], timeseries_tag)
                        hourly_data["wind_ons"][timeseries_tag] = hourly_profile[row_idx]
                    end
                end
            end

            if haskey(hourly_data, "pv")
                for (timeseries_tag, hourly_profile) in scenario_cache["pv"]
                    if haskey(hourly_data["pv"], timeseries_tag)
                        hourly_data["pv"][timeseries_tag] = hourly_profile[row_idx]
                    end
                end
            end

            row_idx += 1
        end
    end

    year_region_masks = am_reference_data[:nw][0][:vre_region_masks][y]
    region_map = get(am_reference_data[:nw][0], :profile_data_region_map, nothing)

    local_nw0 = local_ref[:nw][0]

    # Zero-init dense store for every bus (empty-drs cells stay 0.0, matching an empty mask).
    for bus_idx in bus_keys
        local_nw0[:vre_aggregated_data_zdt][y][bus_idx] = zeros(Float64, nd, nh, 6)
    end

    # tech -> zdt shape key
    tech_specs = (
        ("hydro",    "hydro_shape"),
        ("csp",      "csp_shape"),
        ("wind_ons", "wind_ons_shape"),
        ("wind_ofs", "wind_ofs_shape"),
        ("rtpv",     "rtpv_shape"),
        ("pv",       "pv_shape"),
    )

    for (tech, shape_key) in tech_specs
        # Per-bus ordered DATA-REGION id list: map each finest mask id to its data-region (mirror GTEP).
        # Then group buses by identical drs so the masked average is computed once per distinct group.
        group_drs = Vector{Vector{String}}()
        group_buses = Vector{Vector{Int}}()
        sig_to_group = Dict{Vector{String}, Int}()
        for bus_idx in bus_keys
            mask = year_region_masks[bus_idx][tech]
            drs = String[string(profile_data_region(region_map, tech, region_id)) for region_id in mask]
            isempty(drs) && continue   # dense store already 0.0; matches an empty mask -> 0.0
            gi = get(sig_to_group, drs, 0)
            if gi == 0
                push!(group_drs, drs)
                push!(group_buses, Int[bus_idx])
                sig_to_group[drs] = length(group_drs)
            else
                push!(group_buses[gi], bus_idx)
            end
        end

        for repday_id in repday_ids
            repday_key = string(repday_id)
            for h in run_H
                hourly_data = local_nw0[:planning_stages][y]["repdays"][repday_key]["data"][string(h)]["1"]
                profile_dict = _ra_hourly_vre_profile(hourly_data, tech)   # bus-independent: fetch once
                for gi in eachindex(group_drs)
                    drs = group_drs[gi]
                    n = length(drs)
                    acc = 0.0
                    for dr in drs
                        acc += Float64(get(profile_dict, dr, 0.0))
                    end
                    val = acc / n
                    for bus_idx in group_buses[gi]
                        set_vre_zdt_shape!(local_nw0, y, bus_idx, repday_id, h, shape_key, val)
                    end
                end
            end
        end
    end

    return local_ref
end

function perform_reference_ED_simulation_for_RA(filtered_joint_scenario_map, am_reference_data, am_setting_data, RA_setting; constant_load = 0.0, constant_load_bus_idx=nothing, constant_load_distribution="single_bus", constant_load_distribution_profile=Dict{Tuple{Int, Int, Int}, Dict{Int, Float64}}(), ELCC_flag=false)
    
    # ALEAF_model_instance build-------------------------------------------------
 
    # initialize
    daily_reference_ED_solutions = Dict{Int64, Any}()
    dispatch_report = false
    if ELCC_flag && RA_setting["export_ELCC_dispatch_results_flag"] == true
        dispatch_report = true
    elseif (!ELCC_flag) && RA_setting["export_reference_dispatch_results_flag"] == true
        dispatch_report = true
    end
       
    reference_pair_list = _build_reference_ed_pair_list(filtered_joint_scenario_map)
    num_reference_pairs = length(reference_pair_list)

    # update num_sim_limit_per_call
    configured_limit = max(1, RA_setting["num_distributed_scenarios_per_worker_value"])
    num_sim_limit_per_call = configured_limit

    # Cache large static data on workers; tag-based skip avoids re-sending the multi-MB
    # reference data a worker already holds, and the push is parallelized via @sync @async.
    static_cache_tag = _compute_ra_static_cache_tag(am_reference_data, am_setting_data, RA_setting)
    _set_ra_worker_cache!(am_reference_data, am_setting_data, RA_setting, Dict{Int64, Any}(), nothing, nothing)

    futures = Vector{Tuple{Future, Int}}()
    worker_list = workers()
    if isempty(worker_list)
        error("No distributed workers are available for reference ED analysis.")
    end
    @sync for p in worker_list
        if p != myid()
            @async begin
                needs_refresh = !remotecall_fetch(p, static_cache_tag) do expected_tag
                    get(_RA_WORKER_CACHE, :static_cache_tag, nothing) == expected_tag
                end
                if needs_refresh
                    remotecall_fetch(p, am_reference_data, am_setting_data, RA_setting) do am_ref, am_set, ra_set
                        _set_ra_worker_cache!(am_ref, am_set, ra_set, Dict{Int64, Any}(), nothing, nothing)
                        nothing
                    end
                end
            end
        end
    end

    # Ensure reference ED utilizes all available workers by targeting at least
    # one simulation chunk per worker when daygroups are sufficient.
    target_chunks = max(1, min(length(worker_list), num_reference_pairs))
    worker_balanced_limit = max(1, ceil(Int, num_reference_pairs / target_chunks))
    num_sim_limit_per_call = min(configured_limit, worker_balanced_limit)
    total_num_required_simulations = ceil(Int, num_reference_pairs / num_sim_limit_per_call)

    # rebuild simulation sets based on finalized chunk size
    simulation_set = [Tuple{Int, String}[] for _ in 1:total_num_required_simulations]
    simulation_idx = 1
    check_count = 1
    for reference_pair in reference_pair_list
        push!(simulation_set[simulation_idx], reference_pair)
        if check_count == num_sim_limit_per_call
            simulation_idx += 1
            check_count = 1
        else
            check_count += 1
        end
    end

    for simulation_id in 1:total_num_required_simulations
        p = worker_list[(simulation_id - 1) % length(worker_list) + 1]
        fut = remotecall(run_reference_ED_analysis, p, simulation_set[simulation_id]; constant_load, constant_load_bus_idx, constant_load_distribution, constant_load_distribution_profile, dispatch_report, ELCC_flag)
        push!(futures, (fut, p))
    end
    simulation_set_solution = Vector{Any}(undef, total_num_required_simulations)
    for simulation_id in 1:total_num_required_simulations
        fut, worker_id = futures[simulation_id]
        try
            simulation_set_solution[simulation_id] = fetch(fut)
        catch err
            err_msg = sprint(showerror, err)
            @aleaf_error "Reference ED failed for simulation id $simulation_id (worker $worker_id): $err_msg"
            rethrow(err)
        end
    end

    for simulation_idx in eachindex(simulation_set)
        worker_result = simulation_set_solution[simulation_idx]
        has_valid_worker_result = worker_result isa AbstractDict
        if !has_valid_worker_result
            @aleaf_warn "Reference ED worker result for simulation $simulation_idx is $(typeof(worker_result)); falling back to local build_and_run_reference_ED."
        end

        for (day_group_id, renewable_scenario_id) in simulation_set[simulation_idx]
            if !haskey(daily_reference_ED_solutions, day_group_id)
                daily_reference_ED_solutions[day_group_id] = Dict{String, Any}()
            end

            if has_valid_worker_result && haskey(worker_result, (day_group_id, renewable_scenario_id))
                daily_reference_ED_solutions[day_group_id][renewable_scenario_id] = worker_result[(day_group_id, renewable_scenario_id)]
            else
                daily_reference_ED_solutions[day_group_id][renewable_scenario_id] = build_and_run_reference_ED(day_group_id, renewable_scenario_id, am_reference_data, am_setting_data, RA_setting["current_year"], RA_setting["system_peak_scale"], RA_setting; constant_load, constant_load_bus_idx, constant_load_distribution, constant_load_distribution_profile, dispatch_report)
                if dispatch_report == true
                    report_result_dispatch_reference_RA(am_reference_data, am_setting_data, day_group_id, renewable_scenario_id, RA_setting["current_year"], daily_reference_ED_solutions[day_group_id][renewable_scenario_id], RA_setting; ELCC_flag)
                end
            end
                
        end
    end

    return daily_reference_ED_solutions
end


function run_reference_ED_analysis(reference_pair_list, am_reference_data, am_setting_data, RA_setting; constant_load=0.0, constant_load_bus_idx=nothing, constant_load_distribution="single_bus", constant_load_distribution_profile=Dict{Tuple{Int, Int, Int}, Dict{Int, Float64}}(), dispatch_report=false, ELCC_flag=false)

    reference_ED_solutions = Dict{Tuple{Int, String}, Any}()

    # Run reference ED for each day group and renewable scenario.
    for (day_group_id, renewable_scenario_id) in reference_pair_list
        reference_ED_solutions[(day_group_id, renewable_scenario_id)] = build_and_run_reference_ED(day_group_id, renewable_scenario_id, am_reference_data, am_setting_data, RA_setting["current_year"], RA_setting["system_peak_scale"], RA_setting; constant_load, constant_load_bus_idx, constant_load_distribution, constant_load_distribution_profile, dispatch_report)
    end

    # report reference results
    for (day_group_id, renewable_scenario_id) in reference_pair_list 
        if dispatch_report == true
            report_result_dispatch_reference_RA(am_reference_data, am_setting_data, day_group_id, renewable_scenario_id, RA_setting["current_year"], reference_ED_solutions[(day_group_id, renewable_scenario_id)], RA_setting; ELCC_flag)
        end
    end

    return reference_ED_solutions
end

# lightweight worker overload: pull cached data on worker
function run_reference_ED_analysis(reference_pair_list; constant_load=0.0, constant_load_bus_idx=nothing, constant_load_distribution="single_bus", constant_load_distribution_profile=Dict{Tuple{Int, Int, Int}, Dict{Int, Float64}}(), dispatch_report=false, ELCC_flag=false)
    am_reference_data = _RA_WORKER_CACHE[:am_reference_data]
    am_setting_data = _RA_WORKER_CACHE[:am_setting_data]
    RA_setting = _RA_WORKER_CACHE[:RA_setting]
    return run_reference_ED_analysis(reference_pair_list, am_reference_data, am_setting_data, RA_setting; constant_load, constant_load_bus_idx, constant_load_distribution, constant_load_distribution_profile, dispatch_report, ELCC_flag)
end


function build_and_run_reference_ED(day_group_id, renewable_scenario_id, am_reference_data, am_setting_data, current_year, system_peak_scale, RA_setting; constant_load=0.0, constant_load_bus_idx=nothing, constant_load_distribution="single_bus", constant_load_distribution_profile=Dict{Tuple{Int, Int, Int}, Dict{Int, Float64}}(), dispatch_report=false) 
    
    # build common ALEAF model instance structure
    RA_ALEAF_model_instance = initialize_model_instance_RA(Abstract_LC_GTEP_Model, [1])
    RA_ALEAF_model_instance.setting = am_setting_data
    bus_keys = collect(keys(am_reference_data[:nw][0][:bus]))
    RA_ALEAF_model_instance.ref = _build_ra_local_reference_for_scenario(am_reference_data, RA_setting, am_setting_data["run_H"], day_group_id, renewable_scenario_id, current_year, bus_keys)

    # build ED instance
    build_LCO_GTEP_reference_ED_for_RA_instance!(RA_ALEAF_model_instance, 1, day_group_id, current_year; constant_load, constant_load_bus_idx, constant_load_distribution, constant_load_distribution_profile, system_peak_scale, dispatch_report)

    # Define JuMP_model, solution_list, and solver setting
    JuMP_model = RA_ALEAF_model_instance.model[:nw][1][1]
    solution_list = RA_ALEAF_model_instance.sol[:nw][1][1]
    solver_setting = RA_ALEAF_model_instance.setting["Solver Setting"]
    iteration = 1

    # Solve
    ED_solution = solve_model_GTEP!(JuMP_model, day_group_id, solution_list, solver_setting; iteration, PH_flag=false)

    RA_ALEAF_model_instance = 0.0

    return ED_solution
end


function run_ra_perfect_foresight_batch!(risk_index_list, worker_id, simulation_id; constant_load=0.0, constant_load_bus_idx=nothing, constant_load_distribution="single_bus", constant_load_distribution_profile=Dict{Tuple{Int, Int, Int}, Dict{Int, Float64}}(), new_gen_risk_scenario=nothing, report_log=true, iter=1, worker_log_ch=nothing, force_ELCC_flag=false)

    RA_solution = Dict()
    if report_log
        if worker_log_ch === nothing
            @aleaf_info "[ALEAF RA ED Analysis]: Worker $worker_id start simulation id $simulation_id (risk_count=$(length(risk_index_list)), threads=$(Threads.nthreads()))"
        end
        _emit_worker_log!(worker_log_ch, worker_id, "start simulation_id=$simulation_id risk_count=$(length(risk_index_list))")
    end

    _get_required_ra_worker_cache!([:am_reference_data, :am_setting_data, :RA_setting, :outage_cache_path])
    am_reference_data = _RA_WORKER_CACHE[:am_reference_data]
    am_setting_data = _RA_WORKER_CACHE[:am_setting_data]
    RA_setting = _RA_WORKER_CACHE[:RA_setting]
    filtered_joint_scenario_map = RA_setting["filtered_joint_scenario_map"]

    if haskey(_RA_WORKER_CACHE, :daily_reference_ED_solutions)
        daily_reference_ED_solutions = _RA_WORKER_CACHE[:daily_reference_ED_solutions]
    elseif haskey(_RA_WORKER_CACHE, :reference_cache_path)
        daily_reference_ED_solutions = Dict{Int64, Dict{String, Any}}()
        for joint_risk_id in risk_index_list
            day_group_id, joint_id = joint_risk_id
            joint_scenario = _get_filtered_joint_scenario(filtered_joint_scenario_map, day_group_id, joint_id)
            renewable_scenario_id = joint_scenario["renewable_scenario_id"]
            if !haskey(daily_reference_ED_solutions, day_group_id)
                daily_reference_ED_solutions[day_group_id] = Dict{String, Any}()
            end
            if !haskey(daily_reference_ED_solutions[day_group_id], renewable_scenario_id)
                daily_reference_ED_solutions[day_group_id][renewable_scenario_id] = _load_reference_ed_solution_from_cache(day_group_id, renewable_scenario_id)
            end
        end
    else
        error("Missing worker cache key: :reference_cache_path (and no in-memory :daily_reference_ED_solutions)")
    end

    outage_risk_index_list = Tuple{Int, Int}[]
    joint_scenario_lookup = Dict{Tuple{Int, Int}, Dict{String, Any}}()
    for joint_risk_id in risk_index_list
        day_group_id, joint_id = joint_risk_id
        joint_scenario = _get_filtered_joint_scenario(filtered_joint_scenario_map, day_group_id, joint_id)
        joint_scenario_lookup[joint_risk_id] = joint_scenario
        push!(outage_risk_index_list, (day_group_id, joint_scenario["outage_risk_id"]))
    end
    unique_outage_risk_index_list = unique(outage_risk_index_list)
    risk_scenario = _load_risk_scenarios_from_cache_batch(unique_outage_risk_index_list; newgen=false)
    if new_gen_risk_scenario === nothing && force_ELCC_flag
        first_day = risk_index_list[1][1]
        newgen_file = joinpath(_RA_WORKER_CACHE[:outage_cache_path], "day_newgen_$(first_day).jld2")
        if isfile(newgen_file)
            new_gen_risk_scenario = _load_risk_scenarios_from_cache_batch(unique_outage_risk_index_list; newgen=true)
        else
            new_gen_risk_scenario = Dict{Tuple{Int64, Int64}, Any}()
        end
    elseif new_gen_risk_scenario === nothing
        new_gen_risk_scenario = Dict{Tuple{Int64, Int64}, Any}()
    end
    if length(new_gen_risk_scenario) > 0
        RA_solution = run_ra_perfect_foresight_cases!(risk_index_list, joint_scenario_lookup, risk_scenario, am_reference_data, am_setting_data, daily_reference_ED_solutions, RA_setting; constant_load, constant_load_bus_idx, constant_load_distribution, constant_load_distribution_profile, new_gen_risk_matrix = new_gen_risk_scenario, report_log, iter, worker_id, simulation_id, worker_log_ch, force_ELCC_flag)
    else
        RA_solution = run_ra_perfect_foresight_cases!(risk_index_list, joint_scenario_lookup, risk_scenario, am_reference_data, am_setting_data, daily_reference_ED_solutions, RA_setting; constant_load, constant_load_bus_idx, constant_load_distribution, constant_load_distribution_profile, report_log, worker_id, simulation_id, worker_log_ch, force_ELCC_flag)
    end
    _emit_worker_log!(worker_log_ch, worker_id, "finish simulation_id=$simulation_id")


    return RA_solution
end


function run_ra_perfect_foresight_cases!(risk_index_list, joint_scenario_lookup, risk_scenario, am_reference_data, am_setting_data, daily_reference_ED_solutions, RA_setting; constant_load=0.0, constant_load_bus_idx=nothing, constant_load_distribution="single_bus", constant_load_distribution_profile=Dict{Tuple{Int, Int, Int}, Dict{Int, Float64}}(), new_gen_risk_matrix=[], report_log=true, iter=1, worker_id=0, simulation_id=0, worker_log_ch=nothing, force_ELCC_flag=false) 
    
    current_year = RA_setting["current_year"]
    system_peak_scale = RA_setting["system_peak_scale"]
    detailed_log = lowercase(string(get(RA_setting, "logging_level_value", "simple"))) == "detailed"
    
    ELCC_flag = force_ELCC_flag || (length(new_gen_risk_matrix) > 0) || (abs(constant_load) > 1e-12)
    run_DLOL_flag = get(RA_setting, "calculate_capacity_credit_flag", false) && (string(get(RA_setting, "capacity_credit_type", "")) == "DLOL")

    dispatch_report = false
    if ELCC_flag && (RA_setting["export_ELCC_dispatch_results_flag"])
        dispatch_report = true        
    elseif (!ELCC_flag) && (RA_setting["export_baseline_dispatch_results_flag"])
        dispatch_report = true
    end
    
    # Compute outage scenario list
    outage_scenario_list = eachindex(risk_index_list)
    bus_keys = collect(keys(am_reference_data[:nw][0][:bus]))
    num_units = length(am_reference_data[:nw][0][:RA_gen_index])
    results = Vector{Dict{String, Any}}(undef, length(risk_index_list))
    for i in eachindex(results)
        results[i] = Dict{String, Any}(
            "solution" => Dict{String, Any}(),
            "total ENS" => 0.0,
            "max ENS" => 0.0
        )
    end

    # Add DLOL accumulator initialization if needed (only on workers running ELCC cases with DLOL)
    DLOL_unit_accumulator = Dict{Int, Dict{Symbol, Any}}()
    # Key by actual thread_id to avoid BoundsError when thread_id space exceeds nthreads().
    DLOL_thread_accumulators = Dict{Int, Dict{Int, Dict{Symbol, Any}}}()
    if run_DLOL_flag
        am_ref = am_reference_data[:nw][0]
        pu_power_base = am_setting_data["Simulation Setting"]["per_unit_base_value"]

        for unit_id in sort(collect(keys(am_ref[:RA_gen_index])))
            bus_idx = am_ref[:RA_gen_index][unit_id]["bus_idx"]
            tech_idx = am_ref[:RA_gen_index][unit_id]["genco_tech_id"]
            gen_bus_data = am_reference_data[:nw][bus_idx][:RA_gen_bus][tech_idx]

            plant_name = string(get(gen_bus_data, "PLANT_NAME", get(gen_bus_data, "Tech_ID", string(unit_id))))
            unit_group = string(get(gen_bus_data, "UNITGROUP", ""))
            unit_category = string(get(gen_bus_data, "UNIT_CATEGORY", ""))
            icap_mw = Float64(get(gen_bus_data, "CAP", 0.0) * get(gen_bus_data, "EXUNITS", 0.0) * pu_power_base)

            DLOL_unit_accumulator[unit_id] = Dict{Symbol, Any}(
                :plant_name => plant_name,
                :bus_idx => Int(bus_idx),
                :unit_group => unit_group,
                :unit_category => unit_category,
                :ICAP_MW => icap_mw,
                :DLOL_sum_output_stress_weighted => 0.0,
                :DLOL_sum_stress_hour_weights => 0.0
            )
        end

        if report_log && detailed_log
            _emit_worker_log!(worker_log_ch, worker_id, "DLOL setup: simulation_id=$simulation_id, initialized worker-local unit metadata for $(length(DLOL_unit_accumulator)) units")
        end
    end

    # Check missing reference solutions
    missing_ref_pairs = Set{Tuple{Int, String}}()
    for risk_idx in outage_scenario_list
        joint_risk_id = risk_index_list[risk_idx]
        day_group_in_risk_idx = joint_risk_id[1]
        renewable_scenario_id = joint_scenario_lookup[joint_risk_id]["renewable_scenario_id"]
        if !haskey(daily_reference_ED_solutions, day_group_in_risk_idx) || !haskey(daily_reference_ED_solutions[day_group_in_risk_idx], renewable_scenario_id)
            push!(missing_ref_pairs, (day_group_in_risk_idx, renewable_scenario_id))
        end
    end
    if !isempty(missing_ref_pairs)
        missing_pairs_sorted = sort(collect(missing_ref_pairs))
        error("Missing reference ED cache for simulation_id=$simulation_id, day_group/scenario pairs=$(missing_pairs_sorted)")
    end

    # build and solve ED instance
    function run_ra_perfect_foresight_case!(risk_id, joint_scenario, risk_scenario_of_day, day_group_id, current_year, dispatch_report, system_peak_scale, new_gen_risk_matrix, RA_setting, daily_reference_ED_solutions, local_result, DLOL_thread_accumulator=nothing)
        joint_id = risk_id[2]
        outage_risk_id = joint_scenario["outage_risk_id"]
        renewable_scenario_id = joint_scenario["renewable_scenario_id"]
        model_instance = initialize_model_instance_RA(Abstract_LC_GTEP_Model, [joint_id])
        model_instance.setting = copy(am_setting_data)

        model_instance.ref = _build_ra_local_reference_for_scenario(am_reference_data, RA_setting, am_setting_data["run_H"], day_group_id, renewable_scenario_id, current_year, bus_keys)
    
        ###################################
        #------ Merge Risk Data
        ###################################
        if length(new_gen_risk_matrix) > 0
            risk_scenario_of_day = hcat(risk_scenario_of_day, new_gen_risk_matrix[(day_group_id, outage_risk_id)])
        end
    
        # Build and solve ED instance if there is a generator outage
        num_of_unit_outages_per_hour = vec(num_units .- sum(risk_scenario_of_day, dims=2))
        gen_outage_flag = any(>(0), num_of_unit_outages_per_hour)

        if gen_outage_flag == true
            # build ED instance
            build_ra_perfect_foresight_ed_model!(model_instance, joint_id, day_group_id, current_year, risk_scenario_of_day, daily_reference_ED_solutions[day_group_id][renewable_scenario_id]["solution"], RA_setting; constant_load, constant_load_bus_idx, constant_load_distribution, constant_load_distribution_profile, system_peak_scale, dispatch_report)
            # Define JuMP_model, solution_list, and solver setting
            JuMP_model = model_instance.model[:nw][1][joint_id]
            solution_list = model_instance.sol[:nw][1][joint_id]
            solver_setting = model_instance.setting["Solver Setting"]
            iteration = 1
    
            # Solve model and extract ENS
            ED_solution = solve_model_GTEP!(JuMP_model, day_group_id, solution_list, solver_setting; iteration, PH_flag=false)
    
            scenario_total_ens = sum(ED_solution["solution"]["scarcity"][time_idx]["ens_ndhty"] for time_idx in keys(ED_solution["solution"]["scarcity"]))
            max_ens = maximum(ED_solution["solution"]["scarcity"][time_idx]["ens_ndhty"] for time_idx in keys(ED_solution["solution"]["scarcity"]))
    
            local_result["total ENS"] = scenario_total_ens
            local_result["max ENS"] = max_ens
            if scenario_total_ens > 0
                local_result["solution"] = ED_solution["solution"]
            end

            # DLOL accumulation for this solved scenario:
            # use stressed hours only (system ENS > 0) and apply representative-day weights.
            if run_DLOL_flag && (DLOL_thread_accumulator !== nothing) && (scenario_total_ens > 0)
                ids_d = model_instance.ref[:nw][0][:repday_groups][day_group_id]["Day_Idx_List"]
                ids_h = model_instance.setting["run_H"]
                ids_t = model_instance.setting["run_T"]
                ids_i_local = collect(keys(model_instance.ref[:nw][0][:gen_index]))
                pu_power_base = model_instance.setting["Simulation Setting"]["per_unit_base_value"]

                for d in ids_d, h in ids_h, t in ids_t
                    hourly_system_ens = 0.0
                    for n in bus_keys
                        scarcity_idx = string("(", n, ", ", d, ", ", h, ", ", t, ", ", current_year, ")")
                        if haskey(ED_solution["solution"]["scarcity"], scarcity_idx)
                            hourly_system_ens += ED_solution["solution"]["scarcity"][scarcity_idx]["ens_ndhty"]
                        end
                    end

                    if hourly_system_ens > 0.0
                        num_days_weight = Float64(get(model_instance.ref[:nw][0][:repdays][d], "NumDays", 1.0))
                        for unit_id in ids_i_local
                            dispatch_idx = string("(", unit_id, ", ", d, ", ", h, ", ", t, ", ", current_year, ")")
                            if haskey(ED_solution["solution"]["dispatch"], dispatch_idx)
                                if !haskey(DLOL_thread_accumulator, unit_id)
                                    if haskey(DLOL_unit_accumulator, unit_id)
                                        base = DLOL_unit_accumulator[unit_id]
                                        DLOL_thread_accumulator[unit_id] = Dict{Symbol, Any}(
                                            :plant_name => base[:plant_name],
                                            :bus_idx => base[:bus_idx],
                                            :unit_group => base[:unit_group],
                                            :unit_category => base[:unit_category],
                                            :ICAP_MW => base[:ICAP_MW],
                                            :DLOL_sum_output_stress_weighted => 0.0,
                                            :DLOL_sum_stress_hour_weights => 0.0
                                        )
                                    else
                                        unit_ref = model_instance.ref[:nw][0][:gen_index][unit_id]
                                        bus_idx = unit_ref["bus_idx"]
                                        tech_idx = unit_ref["genco_tech_id"]
                                        gen_bus_data = model_instance.ref[:nw][bus_idx][:gen_bus][tech_idx]
                                        DLOL_thread_accumulator[unit_id] = Dict{Symbol, Any}(
                                            :plant_name => string(get(gen_bus_data, "PLANT_NAME", get(gen_bus_data, "Tech_ID", string(unit_id)))),
                                            :bus_idx => Int(bus_idx),
                                            :unit_group => string(get(gen_bus_data, "UNITGROUP", "")),
                                            :unit_category => string(get(gen_bus_data, "UNIT_CATEGORY", "")),
                                            :ICAP_MW => Float64(get(gen_bus_data, "CAP", 0.0) * get(gen_bus_data, "EXUNITS", 0.0) * pu_power_base),
                                            :DLOL_sum_output_stress_weighted => 0.0,
                                            :DLOL_sum_stress_hour_weights => 0.0
                                        )
                                    end
                                end

                                g_mw = ED_solution["solution"]["dispatch"][dispatch_idx]["g_idhty"] * pu_power_base
                                DLOL_thread_accumulator[unit_id][:DLOL_sum_output_stress_weighted] += g_mw * num_days_weight
                                DLOL_thread_accumulator[unit_id][:DLOL_sum_stress_hour_weights] += num_days_weight
                            end
                        end
                    end
                end
            end
            
            # report dispatch
            if dispatch_report == true
                if max_ens > RA_setting["export_dispatch_results_threshold_value"]    # MW threshold
                    report_result_dispatch_RA(local_result["solution"], model_instance.ref[:nw], model_instance.setting, day_group_id, joint_id, current_year, risk_scenario_of_day, daily_reference_ED_solutions[day_group_id][renewable_scenario_id], RA_setting; ELCC_flag, iter)
                end
            end

            # Drop per-scenario dispatch/storage before returning to master: it is
            # already in CSV and consumed by DLOL, and keeping it OOMs the master.
            if haskey(local_result["solution"], "dispatch")
                delete!(local_result["solution"], "dispatch")
            end
            if haskey(local_result["solution"], "storage")
                delete!(local_result["solution"], "storage")
            end
        end
    end

    # Main parallelized loop (dynamic scheduling)
    risk_index_ch = Channel{Int}(min(length(outage_scenario_list), max(32, Threads.nthreads() * 4)))
    @async begin
        for risk_index in outage_scenario_list
            put!(risk_index_ch, risk_index)
        end
        close(risk_index_ch)
    end

    _dlol_acc_lock = ReentrantLock()
    _task_counter = Threads.Atomic{Int}(0)
    worker_tasks = Task[]
    for _ in 1:Threads.nthreads()
        push!(worker_tasks, Threads.@spawn begin
            local _task_id = Threads.atomic_add!(_task_counter, 1)
            local _dlol_thread_acc = if run_DLOL_flag
                lock(_dlol_acc_lock) do
                    get!(DLOL_thread_accumulators, _task_id, Dict{Int, Dict{Symbol, Any}}())
                end
            else
                nothing
            end
            for risk_index in risk_index_ch
                try
                    risk_id = risk_index_list[risk_index]
                    joint_scenario = joint_scenario_lookup[risk_id]
                    day_group_id = risk_id[1]
                    outage_risk_id = joint_scenario["outage_risk_id"]
                    risk_scenario_of_day = risk_scenario[(day_group_id, outage_risk_id)]

                    local_result = Dict{String, Any}()
                    local_result["solution"] = Dict{String, Any}()
                    local_result["total ENS"] = 0.0
                    local_result["max ENS"] = 0.0

                    # Route DLOL updates to a task-local accumulator for thread-safe writes.
                    if run_DLOL_flag
                        run_ra_perfect_foresight_case!(risk_id, joint_scenario, risk_scenario_of_day, day_group_id, current_year, dispatch_report, system_peak_scale, new_gen_risk_matrix, RA_setting, daily_reference_ED_solutions, local_result, _dlol_thread_acc)
                    else
                        run_ra_perfect_foresight_case!(risk_id, joint_scenario, risk_scenario_of_day, day_group_id, current_year, dispatch_report, system_peak_scale, new_gen_risk_matrix, RA_setting, daily_reference_ED_solutions, local_result)
                    end
                    results[risk_index] = local_result
                catch e
                    local rid = risk_index_list[risk_index]
                    local day_group_id_log = (rid isa Tuple && length(rid) >= 1) ? rid[1] : "NA"
                    local joint_id_log = (rid isa Tuple && length(rid) >= 2) ? rid[2] : "NA"
                    local js = get(joint_scenario_lookup, rid, nothing)
                    local outage_risk_id_log = js === nothing ? "NA" : get(js, "outage_risk_id", "NA")
                    local renewable_scenario_id_log = js === nothing ? "NA" : get(js, "renewable_scenario_id", "NA")
                    local worker_id_log = try myid() catch; 0 end
                    local host_log = try gethostname() catch; "" end
                    err_msg = sprint(showerror, e)
                    @error "Perfect ED failure: worker=$worker_id_log host=$host_log simulation_id=$simulation_id day_group_id=$day_group_id_log joint_id=$joint_id_log outage_risk_id=$outage_risk_id_log renewable_scenario_id=$renewable_scenario_id_log err=$err_msg"
                    flush(stdout); flush(stderr)
                    throw(ErrorException("Perfect ED failed for day_group_id $(day_group_id_log) risk $(rid) (outage=$outage_risk_id_log renew=$renewable_scenario_id_log worker=$worker_id_log host=$host_log): $err_msg"))
                end
            end
        end)
    end
    foreach(wait, worker_tasks)

    # Merge per-thread DLOL partials into one worker-local partial accumulator.
    if run_DLOL_flag
        for DLOL_thread_acc in values(DLOL_thread_accumulators)
            for (unit_id, unit_data) in DLOL_thread_acc
                if !haskey(DLOL_unit_accumulator, unit_id)
                    DLOL_unit_accumulator[unit_id] = Dict{Symbol, Any}(
                        :plant_name => unit_data[:plant_name],
                        :bus_idx => unit_data[:bus_idx],
                        :unit_group => unit_data[:unit_group],
                        :unit_category => unit_data[:unit_category],
                        :ICAP_MW => unit_data[:ICAP_MW],
                        :DLOL_sum_output_stress_weighted => 0.0,
                        :DLOL_sum_stress_hour_weights => 0.0
                    )
                end
                DLOL_unit_accumulator[unit_id][:DLOL_sum_output_stress_weighted] += unit_data[:DLOL_sum_output_stress_weighted]
                DLOL_unit_accumulator[unit_id][:DLOL_sum_stress_hour_weights] += unit_data[:DLOL_sum_stress_hour_weights]
            end
        end

        if report_log && detailed_log
            units_with_stressed_hours = count(v -> v[:DLOL_sum_stress_hour_weights] > 0.0, values(DLOL_unit_accumulator))
            _emit_worker_log!(worker_log_ch, worker_id, "DLOL partial merge: simulation_id=$simulation_id, units_with_stress_hours=$units_with_stressed_hours")
        end
    end

    result_EA_analysis_of_day = Dict{Any, Any}()
    for i in eachindex(risk_index_list)
        result_EA_analysis_of_day[risk_index_list[i]] = results[i]
    end
    # DLOL mode returns wrapped payload so master can consume both scenario results
    # and worker-local DLOL partials.
    if run_DLOL_flag
        return Dict{String, Any}(
            "scenario_results" => result_EA_analysis_of_day,
            "DLOL_unit_accumulator_partial" => DLOL_unit_accumulator
        )
    end
    return result_EA_analysis_of_day
end


function run_ra_imperfect_foresight_sequential_snapshot_batch!(risk_index_list, worker_id, simulation_id; constant_load=0.0, constant_load_bus_idx=nothing, constant_load_distribution="single_bus", constant_load_distribution_profile=Dict{Tuple{Int, Int, Int}, Dict{Int, Float64}}(), new_gen_risk_scenario=nothing, report_log=true, iter=1, worker_log_ch=nothing, force_ELCC_flag=false)

    RA_solution = Dict()
    if report_log
        if worker_log_ch === nothing
            @aleaf_info "[ALEAF RA ED Analysis]: Worker $worker_id start sequential simulation id $simulation_id (risk_count=$(length(risk_index_list)), threads=$(Threads.nthreads()))"
        end
        _emit_worker_log!(worker_log_ch, worker_id, "start sequential simulation_id=$simulation_id risk_count=$(length(risk_index_list))")
    end

    _get_required_ra_worker_cache!([:am_reference_data, :am_setting_data, :RA_setting, :outage_cache_path])
    am_reference_data = _RA_WORKER_CACHE[:am_reference_data]
    am_setting_data = _RA_WORKER_CACHE[:am_setting_data]
    RA_setting = _RA_WORKER_CACHE[:RA_setting]
    filtered_joint_scenario_map = RA_setting["filtered_joint_scenario_map"]
    if haskey(_RA_WORKER_CACHE, :daily_reference_ED_solutions)
        daily_reference_ED_solutions = _RA_WORKER_CACHE[:daily_reference_ED_solutions]
    elseif haskey(_RA_WORKER_CACHE, :reference_cache_path)
        daily_reference_ED_solutions = Dict{Int64, Dict{String, Any}}()
        for joint_risk_id in risk_index_list
            day_group_id, joint_id = joint_risk_id
            joint_scenario = _get_filtered_joint_scenario(filtered_joint_scenario_map, day_group_id, joint_id)
            renewable_scenario_id = joint_scenario["renewable_scenario_id"]
            if !haskey(daily_reference_ED_solutions, day_group_id)
                daily_reference_ED_solutions[day_group_id] = Dict{String, Any}()
            end
            if !haskey(daily_reference_ED_solutions[day_group_id], renewable_scenario_id)
                daily_reference_ED_solutions[day_group_id][renewable_scenario_id] = _load_reference_ed_solution_from_cache(day_group_id, renewable_scenario_id)
            end
        end
    else
        error("Missing worker cache key: :reference_cache_path (and no in-memory :daily_reference_ED_solutions)")
    end
    outage_risk_index_list = Tuple{Int, Int}[]
    joint_scenario_lookup = Dict{Tuple{Int, Int}, Dict{String, Any}}()
    for joint_risk_id in risk_index_list
        day_group_id, joint_id = joint_risk_id
        joint_scenario = _get_filtered_joint_scenario(filtered_joint_scenario_map, day_group_id, joint_id)
        joint_scenario_lookup[joint_risk_id] = joint_scenario
        push!(outage_risk_index_list, (day_group_id, joint_scenario["outage_risk_id"]))
    end
    unique_outage_risk_index_list = unique(outage_risk_index_list)
    risk_scenario = _load_risk_scenarios_from_cache_batch(unique_outage_risk_index_list; newgen=false)
    if new_gen_risk_scenario === nothing && force_ELCC_flag
        first_day = risk_index_list[1][1]
        newgen_file = joinpath(_RA_WORKER_CACHE[:outage_cache_path], "day_newgen_$(first_day).jld2")
        if isfile(newgen_file)
            new_gen_risk_scenario = _load_risk_scenarios_from_cache_batch(unique_outage_risk_index_list; newgen=true)
        else
            new_gen_risk_scenario = Dict{Tuple{Int64, Int64}, Any}()
        end
    elseif new_gen_risk_scenario === nothing
        new_gen_risk_scenario = Dict{Tuple{Int64, Int64}, Any}()
    end
    if length(new_gen_risk_scenario) > 0
        RA_solution = run_ra_imperfect_foresight_sequential_snapshot_cases!(risk_index_list, joint_scenario_lookup, risk_scenario, am_reference_data, am_setting_data, daily_reference_ED_solutions, RA_setting, simulation_id; constant_load, constant_load_bus_idx, constant_load_distribution, constant_load_distribution_profile, new_gen_risk_matrix = new_gen_risk_scenario, report_log, iter, worker_id, worker_log_ch, force_ELCC_flag)
    else
        RA_solution = run_ra_imperfect_foresight_sequential_snapshot_cases!(risk_index_list, joint_scenario_lookup, risk_scenario, am_reference_data, am_setting_data, daily_reference_ED_solutions, RA_setting, simulation_id; constant_load, constant_load_bus_idx, constant_load_distribution, constant_load_distribution_profile, report_log, worker_id, worker_log_ch, force_ELCC_flag)

    end
    _emit_worker_log!(worker_log_ch, worker_id, "finish sequential simulation_id=$simulation_id")


    return RA_solution
end


function run_ra_imperfect_foresight_sequential_snapshot_cases!(risk_index_list, joint_scenario_lookup, risk_scenario, am_reference_data, am_setting_data, daily_reference_ED_solutions, RA_setting, simulation_id; constant_load=0.0, constant_load_bus_idx=nothing, constant_load_distribution="single_bus", constant_load_distribution_profile=Dict{Tuple{Int, Int, Int}, Dict{Int, Float64}}(), new_gen_risk_matrix=[], report_log=true, iter=1, worker_id=0, worker_log_ch=nothing, force_ELCC_flag=false) 
    
    current_year = RA_setting["current_year"]
    system_peak_scale = RA_setting["system_peak_scale"]
    run_DLOL_flag = get(RA_setting, "calculate_capacity_credit_flag", false) && (string(get(RA_setting, "capacity_credit_type", "")) == "DLOL")

    ELCC_flag = force_ELCC_flag || (length(new_gen_risk_matrix) > 0) || (abs(constant_load) > 1e-12)
    dispatch_report = false
    if ELCC_flag && (RA_setting["export_ELCC_dispatch_results_flag"] == true)
        dispatch_report = true
    elseif (!ELCC_flag) && (RA_setting["export_baseline_dispatch_results_flag"] == true)
        dispatch_report = true
    end
    
    outage_scenario_list = eachindex(risk_index_list)
    
    results = Vector{Dict{String, Any}}(undef, length(risk_index_list))
    for i in eachindex(results)
        # initialize to avoid unassigned slots when individual scenarios fail
        results[i] = Dict{String, Any}(
            "solution" => Dict{String, Any}(),
            "total ENS" => 0.0,
            "max ENS" => 0.0
        )
    end

    missing_ref_pairs = Set{Tuple{Int, String}}()
    for risk_idx in outage_scenario_list
        joint_risk_id = risk_index_list[risk_idx]
        day_group_in_risk_idx = joint_risk_id[1]
        renewable_scenario_id = joint_scenario_lookup[joint_risk_id]["renewable_scenario_id"]
        if !haskey(daily_reference_ED_solutions, day_group_in_risk_idx) || !haskey(daily_reference_ED_solutions[day_group_in_risk_idx], renewable_scenario_id)
            push!(missing_ref_pairs, (day_group_in_risk_idx, renewable_scenario_id))
        end
    end
    if !isempty(missing_ref_pairs)
        missing_pairs_sorted = sort(collect(missing_ref_pairs))
        error("Missing reference ED cache for simulation_id=$simulation_id, day_group/scenario pairs=$(missing_pairs_sorted)")
    end

    bus_keys = collect(keys(am_reference_data[:nw][0][:bus]))
    run_h = collect(am_setting_data["run_H"])
    run_t = collect(am_setting_data["run_T"])
    num_units = length(am_reference_data[:nw][0][:RA_gen_index])

    # DLOL worker-local accumulators (thread-safe merge pattern, consistent with perfect-foresight path)
    DLOL_unit_accumulator = Dict{Int, Dict{Symbol, Any}}()
    DLOL_thread_accumulators = Dict{Int, Dict{Int, Dict{Symbol, Any}}}()
    if run_DLOL_flag
        am_ref = am_reference_data[:nw][0]
        pu_power_base = am_setting_data["Simulation Setting"]["per_unit_base_value"]
        for unit_id in sort(collect(keys(am_ref[:RA_gen_index])))
            bus_idx = am_ref[:RA_gen_index][unit_id]["bus_idx"]
            tech_idx = am_ref[:RA_gen_index][unit_id]["genco_tech_id"]
            gen_bus_data = am_reference_data[:nw][bus_idx][:RA_gen_bus][tech_idx]
            plant_name = string(get(gen_bus_data, "PLANT_NAME", get(gen_bus_data, "Tech_ID", string(unit_id))))
            unit_group = string(get(gen_bus_data, "UNITGROUP", ""))
            unit_category = string(get(gen_bus_data, "UNIT_CATEGORY", ""))
            icap_mw = Float64(get(gen_bus_data, "CAP", 0.0) * get(gen_bus_data, "EXUNITS", 0.0) * pu_power_base)
            DLOL_unit_accumulator[unit_id] = Dict{Symbol, Any}(
                :plant_name => plant_name,
                :bus_idx => Int(bus_idx),
                :unit_group => unit_group,
                :unit_category => unit_category,
                :ICAP_MW => icap_mw,
                :DLOL_sum_output_stress_weighted => 0.0,
                :DLOL_sum_stress_hour_weights => 0.0
            )
        end
    end

    # Rolling-horizon sequential ED: event hours solve a k-hour MPC LP (commit only h),
    # non-event hours synthesize SOC from ref. k = sequential_horizon_hours_value (default 6).
    function run_ra_imperfect_foresight_sequential_snapshot_case!(risk_id, joint_scenario, risk_scenario_of_day, day_group_id, current_year, dispatch_report, system_peak_scale, new_gen_risk_matrix, RA_setting, daily_reference_ED_solutions, local_result, DLOL_thread_accumulator=nothing)
        joint_id = risk_id[2]
        outage_risk_id = joint_scenario["outage_risk_id"]
        renewable_scenario_id = joint_scenario["renewable_scenario_id"]

        # build common ALEAF model instance structure (local to avoid races)
        model_instance = initialize_model_instance_RA(Abstract_LC_GTEP_Model, [joint_id])
        model_instance.setting = copy(am_setting_data)

        model_instance.ref = _build_ra_local_reference_for_scenario(am_reference_data, RA_setting, am_setting_data["run_H"], day_group_id, renewable_scenario_id, current_year, bus_keys)

        ###################################
        #------ Merge Risk Data
        ###################################
        if length(new_gen_risk_matrix) > 0
            risk_scenario_of_day = hcat(risk_scenario_of_day, new_gen_risk_matrix[(day_group_id, outage_risk_id)])
        end

        ###################################
        #------ Check Risk Events
        ###################################
        ids_d = [(d) for (d) in model_instance.ref[:nw][0][:repday_groups][day_group_id]["Day_Idx_List"]]
        ids_h = run_h
        ids_t = run_t
        ids_dht = [(d,h,t) for d in ids_d for h in ids_h for t in ids_t]

        num_of_unit_outages_per_hour = vec(num_units .- sum(risk_scenario_of_day, dims=2))
        indices_of_outage_hours = findall(x -> x > 0, num_of_unit_outages_per_hour)

        if length(indices_of_outage_hours) > 0

            h_first_event_idx = indices_of_outage_hours[1]
            simulated_event_hours = []         # accumulates across the whole window (not reset per event)
            failed_this_window = false

            # Reference ED solution and BATEFF cache for this scenario
            local_ref_solution = daily_reference_ED_solutions[day_group_id][renewable_scenario_id]["solution"]
            card_T_local = length(am_setting_data["run_T"])
            storage_meta = Dict{Int, Tuple{Int, Int, Float64}}()   # i => (bus_idx, tech_idx, sqrt(BATEFF))
            for i in 1:num_units
                if am_reference_data[:nw][0][:RA_gen_index][i]["UNIT_CATEGORY"] == "STORAGE"
                    b = am_reference_data[:nw][0][:RA_gen_index][i]["bus_idx"]
                    t = am_reference_data[:nw][0][:RA_gen_index][i]["genco_tech_id"]
                    bateff = am_reference_data[:nw][b][:RA_gen_bus][t]["BATEFF"]
                    storage_meta[i] = (b, t, sqrt(bateff))
                end
            end

            # Rolling-horizon length k (default 6): k=1 is single-hour snapshot, larger k gives
            # charging look-ahead without true foresight (later hours use a projected outage mask).
            horizon_k = Int(get(RA_setting, "sequential_horizon_hours_value", 6))
            if horizon_k < 1
                horizon_k = 1
            end

            # When neither CSV export nor DLOL need full per-hour history, drop
            # synthesized non-event-hour keys once the next event hour's LP consumes them.
            incremental_drop = !dispatch_report && (DLOL_thread_accumulator === nothing)
            synthesized_keys = String[]

            for window_offset in 0:(length(ids_dht) - h_first_event_idx)
                row_idx = h_first_event_idx + window_offset
                event_id_set = ids_dht[row_idx]
                event_day = event_id_set[1]
                event_hour = event_id_set[2]
                event_time = event_id_set[3]
                is_event_hour = num_of_unit_outages_per_hour[row_idx] > 0
                daygroup_first_hour_flag = (row_idx == 1)

                if is_event_hour
                    # Build the rolling-horizon window [h, h+k-1], truncated to end-of-day-group.
                    actual_horizon = min(horizon_k, length(ids_dht) - row_idx + 1)
                    horizon_dht_list = Vector{Tuple{Int,Int,Int}}(undef, actual_horizon)
                    for r in 1:actual_horizon
                        horizon_dht_list[r] = ids_dht[row_idx + r - 1]
                    end

                    # Outage projection: units out at h use their real recovery profile, units
                    # available at h stay available across the horizon (standard SCED look-ahead).
                    horizon_gen_risk_matrix = BitMatrix(undef, actual_horizon, num_units)
                    for i in 1:num_units
                        unit_avail_at_h = Bool(risk_scenario_of_day[row_idx, i])
                        @inbounds horizon_gen_risk_matrix[1, i] = unit_avail_at_h
                        for r in 2:actual_horizon
                            if unit_avail_at_h
                                @inbounds horizon_gen_risk_matrix[r, i] = true
                            else
                                actual_row = row_idx + r - 1
                                @inbounds horizon_gen_risk_matrix[r, i] = Bool(risk_scenario_of_day[actual_row, i])
                            end
                        end
                    end

                    # Build/solve the rolling-horizon LP via PF-ED; prior_state_solution
                    # supplies committed prior-hour SoC to the first-hour SoC balance helper.
                    build_ra_perfect_foresight_ed_model!(model_instance, joint_id, day_group_id, current_year, horizon_gen_risk_matrix, local_ref_solution, RA_setting; constant_load, constant_load_bus_idx, constant_load_distribution, constant_load_distribution_profile, system_peak_scale, dispatch_report, horizon_dht_list, prior_state_solution=local_result["solution"], horizon_daygroup_first_hour_flag=daygroup_first_hour_flag)

                    JuMP_model = model_instance.model[:nw][1][joint_id]
                    solution_list = model_instance.sol[:nw][1][joint_id]
                    solver_setting = model_instance.setting["Solver Setting"]
                    iteration = 1

                    ED_solution = solve_model_GTEP!(JuMP_model, joint_id, solution_list, solver_setting; iteration, PH_flag=false, terminate_if_error=false)

                    if ED_solution == "No Solution"
                        lp_file_name = string("RA_day_id_", risk_id[1], "_joint_id_", joint_id, "_dispatch.lp")
                        JuMP.write_to_file(JuMP_model, lp_file_name)
                        @aleaf_warn "[ALEAF RA ED Analysis]: Failed to Find a Solution: Day Group id: $day_group_id, Joint ID: $risk_id, Event Hour: $event_hour; aborting redispatch window; LP file generated: $lp_file_name"
                        failed_this_window = true
                        break
                    end

                    # MPC commit: only hour h is kept; the look-ahead at h+1..h+k-1 is
                    # discarded and resolved next iteration with the real outage mask.
                    commit_key_suffix = string(", ", event_day, ", ", event_hour, ", ", event_time, ", ", current_year, ")")

                    ens = 0.0
                    max_ens = 0.0
                    if haskey(ED_solution["solution"], "scarcity")
                        for time_idx in keys(ED_solution["solution"]["scarcity"])
                            if !endswith(time_idx, commit_key_suffix)
                                continue
                            end
                            if ED_solution["solution"]["scarcity"][time_idx]["ens_ndhty"] > 0
                                ens += ED_solution["solution"]["scarcity"][time_idx]["ens_ndhty"]
                                if ED_solution["solution"]["scarcity"][time_idx]["ens_ndhty"] > max_ens
                                    max_ens = ED_solution["solution"]["scarcity"][time_idx]["ens_ndhty"]
                                end
                            end
                        end
                    end

                    local_result["total ENS"] += ens
                    if max_ens > local_result["max ENS"]
                        local_result["max ENS"] = max_ens
                    end

                    # Merge only commit-hour entries; look-ahead entries (hours > h) are
                    # dropped so downstream prior-state lookups never see tentative values.
                    for solutionID in keys(ED_solution["solution"])
                        src = ED_solution["solution"][solutionID]
                        if !(src isa AbstractDict)
                            continue
                        end
                        if !haskey(local_result["solution"], solutionID)
                            local_result["solution"][solutionID] = Dict{String, Any}()
                        end
                        dst = local_result["solution"][solutionID]
                        for (k_str, v) in src
                            if k_str isa AbstractString && endswith(k_str, commit_key_suffix)
                                # For scarcity, also drop zero-ENS entries (existing
                                # convention; avoids inflating local_result with zeros).
                                if solutionID == "scarcity"
                                    if v isa AbstractDict && get(v, "ens_ndhty", 0.0) > 0
                                        dst[k_str] = v
                                    end
                                else
                                    dst[k_str] = v
                                end
                            end
                        end
                    end

                    push!(simulated_event_hours, event_id_set)

                    # Drop synthesized non-event-hour entries from the gap before this event hour;
                    # the LP already consumed them and later hours chain from solved values.
                    if incremental_drop && !isempty(synthesized_keys)
                        dispatch_d = get(local_result["solution"], "dispatch", nothing)
                        storage_d  = get(local_result["solution"], "storage",  nothing)
                        for k in synthesized_keys
                            dispatch_d === nothing || delete!(dispatch_d, k)
                            storage_d  === nothing || delete!(storage_d,  k)
                        end
                        empty!(synthesized_keys)
                    end

                else
                    # Non-event hour: no LP. Synthesize prior-state entries for downstream
                    # event hours (g, chg = ref; SoC chained from prior via SoC balance).
                    if !haskey(local_result["solution"], "dispatch")
                        local_result["solution"]["dispatch"] = Dict{String, Any}()
                    end
                    if !haskey(local_result["solution"], "storage")
                        local_result["solution"]["storage"] = Dict{String, Any}()
                    end

                    # Prior idx for SoC chaining (same convention as constraint_ES_SOC_Balance_Inter_Hour_NoReserve_Cont_Snapshot_idhty)
                    prior_idx_template = if event_hour == 1 && event_time == 1
                        (event_day - 1, 24, 1)
                    elseif event_hour != 1 && event_time == 1
                        (event_day, event_hour - 1, event_time)
                    else
                        (event_day, event_hour, event_time)   # sub-hourly fallback; same hour, will pick up ref via lookup
                    end

                    for i in 1:num_units
                        idx_str = string("(", i, ", ", event_day, ", ", event_hour, ", ", event_time, ", ", current_year, ")")
                        if !haskey(local_ref_solution, "dispatch") || !haskey(local_ref_solution["dispatch"], idx_str)
                            continue
                        end
                        ref_g = local_ref_solution["dispatch"][idx_str]["g_idhty"]
                        local_result["solution"]["dispatch"][idx_str] = Dict{String, Any}("g_idhty" => ref_g)
                        # Track synthesized key for incremental drop after the next event hour
                        # solves (no-op when incremental_drop is false).
                        incremental_drop && push!(synthesized_keys, idx_str)

                        if haskey(storage_meta, i)
                            (_, _, sqrt_bateff) = storage_meta[i]
                            ref_chg = local_ref_solution["storage"][idx_str]["chg_idhty"]

                            prior_idx_str = string("(", i, ", ", prior_idx_template[1], ", ", prior_idx_template[2], ", ", prior_idx_template[3], ", ", current_year, ")")
                            prior_soc = if haskey(local_result["solution"]["storage"], prior_idx_str)
                                local_result["solution"]["storage"][prior_idx_str]["soc_idhty"]
                            elseif haskey(local_ref_solution, "storage") && haskey(local_ref_solution["storage"], prior_idx_str)
                                local_ref_solution["storage"][prior_idx_str]["soc_idhty"]
                            else
                                # Last-resort fallback (e.g., day-group first hour): use current-hour ref soc
                                local_ref_solution["storage"][idx_str]["soc_idhty"]
                            end

                            new_soc = prior_soc + (sqrt_bateff * ref_chg - ref_g / sqrt_bateff) / card_T_local
                            local_result["solution"]["storage"][idx_str] = Dict{String, Any}("chg_idhty" => ref_chg, "soc_idhty" => new_soc)
                        end
                    end
                end
            end

            if failed_this_window
                @aleaf_warn "[ALEAF RA ED Analysis]: Aborted continuous redispatch window for day_group_id=$day_group_id, joint_id=$joint_id"
            end
    
            if dispatch_report == true
                if local_result["max ENS"] > RA_setting["export_dispatch_results_threshold_value"]     
                    report_result_dispatch_RA(local_result["solution"], model_instance.ref[:nw], model_instance.setting, day_group_id, joint_id, current_year, risk_scenario_of_day, daily_reference_ED_solutions[day_group_id][renewable_scenario_id], RA_setting; ELCC_flag, iter)
                end
            end

            # DLOL accumulation for this solved sequential scenario:
            # use stressed hours only (system ENS > 0) and representative-day weights.
            if run_DLOL_flag && (DLOL_thread_accumulator !== nothing) && (local_result["total ENS"] > 0)
                ids_d = model_instance.ref[:nw][0][:repday_groups][day_group_id]["Day_Idx_List"]
                ids_h = run_h
                ids_t = run_t
                ids_i_local = collect(keys(model_instance.ref[:nw][0][:gen_index]))
                pu_power_base = model_instance.setting["Simulation Setting"]["per_unit_base_value"]

                for d in ids_d, h in ids_h, t in ids_t
                    hourly_system_ens = 0.0
                    for n in bus_keys
                        scarcity_idx = string("(", n, ", ", d, ", ", h, ", ", t, ", ", current_year, ")")
                        if haskey(local_result["solution"]["scarcity"], scarcity_idx)
                            hourly_system_ens += local_result["solution"]["scarcity"][scarcity_idx]["ens_ndhty"]
                        end
                    end

                    if hourly_system_ens > 0.0
                        num_days_weight = Float64(get(model_instance.ref[:nw][0][:repdays][d], "NumDays", 1.0))
                        for unit_id in ids_i_local
                            dispatch_idx = string("(", unit_id, ", ", d, ", ", h, ", ", t, ", ", current_year, ")")
                            if haskey(local_result["solution"]["dispatch"], dispatch_idx)
                                if !haskey(DLOL_thread_accumulator, unit_id)
                                    if haskey(DLOL_unit_accumulator, unit_id)
                                        base = DLOL_unit_accumulator[unit_id]
                                        DLOL_thread_accumulator[unit_id] = Dict{Symbol, Any}(
                                            :plant_name => base[:plant_name],
                                            :bus_idx => base[:bus_idx],
                                            :unit_group => base[:unit_group],
                                            :unit_category => base[:unit_category],
                                            :ICAP_MW => base[:ICAP_MW],
                                            :DLOL_sum_output_stress_weighted => 0.0,
                                            :DLOL_sum_stress_hour_weights => 0.0
                                        )
                                    else
                                        unit_ref = model_instance.ref[:nw][0][:gen_index][unit_id]
                                        bus_idx = unit_ref["bus_idx"]
                                        tech_idx = unit_ref["genco_tech_id"]
                                        gen_bus_data = model_instance.ref[:nw][bus_idx][:gen_bus][tech_idx]
                                        DLOL_thread_accumulator[unit_id] = Dict{Symbol, Any}(
                                            :plant_name => string(get(gen_bus_data, "PLANT_NAME", get(gen_bus_data, "Tech_ID", string(unit_id)))),
                                            :bus_idx => Int(bus_idx),
                                            :unit_group => string(get(gen_bus_data, "UNITGROUP", "")),
                                            :unit_category => string(get(gen_bus_data, "UNIT_CATEGORY", "")),
                                            :ICAP_MW => Float64(get(gen_bus_data, "CAP", 0.0) * get(gen_bus_data, "EXUNITS", 0.0) * pu_power_base),
                                            :DLOL_sum_output_stress_weighted => 0.0,
                                            :DLOL_sum_stress_hour_weights => 0.0
                                        )
                                    end
                                end

                                g_idhty_pu = local_result["solution"]["dispatch"][dispatch_idx]["g_idhty"]
                                g_idhty_mw = g_idhty_pu * pu_power_base
                                DLOL_thread_accumulator[unit_id][:DLOL_sum_output_stress_weighted] += g_idhty_mw * num_days_weight
                                DLOL_thread_accumulator[unit_id][:DLOL_sum_stress_hour_weights] += num_days_weight
                            end
                        end
                    end
                end
            end

            # Drop per-scenario dispatch/storage unconditionally before returning to
            # master: already in CSV and consumed by DLOL, and keeping it OOMs the master.
            if haskey(local_result["solution"], "dispatch")
                delete!(local_result["solution"], "dispatch")
            end
            if haskey(local_result["solution"], "storage")
                delete!(local_result["solution"], "storage")
            end
        end

    end

    # Initialize progress tracking
    total_work = length(outage_scenario_list)
    pct_step    = 25              # log every x%
    t0 = time()
    
    # Channel/logger only when reporting is enabled.
    progress_ch = report_log ? Channel{Bool}(Inf) : nothing
    logger_task = nothing
    if report_log
        logger_task = @async begin
            progress   = 0
            last_pct   = -1

            while progress < total_work
                # Exit cleanly if the channel closes during shutdown/error paths.
                try
                    take!(progress_ch)                 # wait for a tick
                catch e
                    if e isa InvalidStateException
                        break
                    else
                        rethrow(e)
                    end
                end
                progress += 1

                pct = Int(floor(100 * progress / max(total_work, 1)))

                if pct >= last_pct + pct_step
                    elapsed = time() - t0
                    if worker_log_ch === nothing
                        @aleaf_info "worker progress=$pct% ($progress/$total_work), elapsed= $(round(elapsed, digits=1))s"
                    end
                    _emit_worker_log!(worker_log_ch, worker_id, "sim=$simulation_id progress=$pct% ($progress/$total_work) t=$(round(elapsed, digits=1))s")
                    last_pct   = pct - (pct % pct_step)
                end
            end
        end
    end

    # Dynamic case assignment across local threads.
    # Buffer is sized to keep workers fed without over-buffering very large runs.
    risk_index_ch = Channel{Int}(min(length(outage_scenario_list), max(32, Threads.nthreads() * 4)))
    @async begin
        # Producer: enqueue all case indices once, then close to signal completion.
        for risk_index in outage_scenario_list
            put!(risk_index_ch, risk_index)
        end
        close(risk_index_ch)
    end

    _dlol_acc_lock = ReentrantLock()
    _task_counter = Threads.Atomic{Int}(0)
    worker_tasks = Task[]
    for _ in 1:Threads.nthreads()
        push!(worker_tasks, Threads.@spawn begin
            local _task_id = Threads.atomic_add!(_task_counter, 1)
            local _dlol_thread_acc = if run_DLOL_flag
                lock(_dlol_acc_lock) do
                    get!(DLOL_thread_accumulators, _task_id, Dict{Int, Dict{Symbol, Any}}())
                end
            else
                nothing
            end
            # Each thread pulls next available case until channel closes.
            for risk_index in risk_index_ch
                local ticked = false
                try
                    risk_id = risk_index_list[risk_index]
                    joint_scenario = joint_scenario_lookup[risk_id]
                    day_group_id = risk_id[1]
                    outage_risk_id = joint_scenario["outage_risk_id"]
                    risk_scenario_of_day = risk_scenario[(day_group_id, outage_risk_id)]

                    local_result = Dict{String, Any}()
                    local_result["solution"] = Dict{String, Any}()
                    local_result["total ENS"] = 0.0
                    local_result["max ENS"] = 0.0
                    if run_DLOL_flag
                        run_ra_imperfect_foresight_sequential_snapshot_case!(risk_id, joint_scenario, risk_scenario_of_day, day_group_id,
                            current_year, dispatch_report, system_peak_scale, new_gen_risk_matrix,
                            RA_setting, daily_reference_ED_solutions, local_result, _dlol_thread_acc)
                    else
                        run_ra_imperfect_foresight_sequential_snapshot_case!(risk_id, joint_scenario, risk_scenario_of_day, day_group_id,
                            current_year, dispatch_report, system_peak_scale, new_gen_risk_matrix,
                            RA_setting, daily_reference_ED_solutions, local_result)
                    end

                    results[risk_index] = local_result
                    if report_log
                        # Exactly one progress tick per attempted case.
                        if isopen(progress_ch)
                            try
                                put!(progress_ch, true)
                                ticked = true
                            catch e
                                if !(e isa InvalidStateException)
                                    rethrow(e)
                                end
                            end
                        end
                    end

                catch e
                    local rid = risk_index_list[risk_index]
                    local day_group_id_log = (rid isa Tuple && length(rid) >= 1) ? rid[1] : "NA"
                    err_msg = sprint(showerror, e)
                    @aleaf_error "Sequential ED failed for day_group_id $(day_group_id_log) risk $(rid): $err_msg"
                    rethrow()
                finally
                    if report_log && !ticked
                        # Ensure logger does not stall if an exception occurs before tick.
                        if isopen(progress_ch)
                            try
                                put!(progress_ch, true)
                            catch e
                                if !(e isa InvalidStateException)
                                    rethrow(e)
                                end
                            end
                        end
                    end
                end
            end
        end)
    end
    
    worker_error = nothing
    try
        # Wait all workers so progress channel is not closed while workers are still reporting.
        for t in worker_tasks
            try
                wait(t)
            catch e
                if worker_error === nothing
                    if e isa TaskFailedException && e.task.exception !== nothing
                        worker_error = e.task.exception
                    else
                        worker_error = ErrorException(sprint(showerror, e))
                    end
                end
            end
        end
    finally
        if report_log
            # Close progress channel after all workers complete, then drain logger task.
            close(progress_ch)
            wait(logger_task)
        end
    end
    if worker_error !== nothing
        throw(worker_error)
    end

    result_EA_analysis_of_day = Dict()
    for i in eachindex(risk_index_list)
        result_EA_analysis_of_day[risk_index_list[i]] = results[i]
    end

    # Merge per-thread DLOL partials into one worker-local partial accumulator.
    if run_DLOL_flag
        for DLOL_thread_acc in values(DLOL_thread_accumulators)
            for (unit_id, unit_data) in DLOL_thread_acc
                if !haskey(DLOL_unit_accumulator, unit_id)
                    DLOL_unit_accumulator[unit_id] = Dict{Symbol, Any}(
                        :plant_name => unit_data[:plant_name],
                        :bus_idx => unit_data[:bus_idx],
                        :unit_group => unit_data[:unit_group],
                        :unit_category => unit_data[:unit_category],
                        :ICAP_MW => unit_data[:ICAP_MW],
                        :DLOL_sum_output_stress_weighted => 0.0,
                        :DLOL_sum_stress_hour_weights => 0.0
                    )
                end
                DLOL_unit_accumulator[unit_id][:DLOL_sum_output_stress_weighted] += unit_data[:DLOL_sum_output_stress_weighted]
                DLOL_unit_accumulator[unit_id][:DLOL_sum_stress_hour_weights] += unit_data[:DLOL_sum_stress_hour_weights]
            end
        end
        return Dict{String, Any}(
            "scenario_results" => result_EA_analysis_of_day,
            "DLOL_unit_accumulator_partial" => DLOL_unit_accumulator
        )
    end

    return result_EA_analysis_of_day
end


function report_result_dispatch_RA(result::Dict{String, Any}, am_ref_data, am_setting_data, day_group_id, risk_id, current_year, risk_scenario_of_day, daily_reference_ED_solutions, RA_setting; ELCC_flag=false, iter=1)
    
    output_path = am_ref_data[0][:output_path]
    system_peak_scale = RA_setting["system_peak_scale"]

    pu_power_base = am_setting_data["Simulation Setting"]["per_unit_base_value"]
    pu_econ_base = am_setting_data["Simulation Setting"]["per_unit_econ_base_value"] / pu_power_base

    case_name = get(am_setting_data["Simulation Configuration"], "Case_ID", "RA")

    # Indices and variables
    ids_h = [(h) for (h) in am_setting_data["run_H"]]
    ids_t = [(t) for (t) in am_setting_data["run_T"]]
    ids_i = [(i) for (i) in keys(am_ref_data[0][:RA_gen_index])]
    ids_k = [(k) for (k) in keys(am_ref_data[0][:branch])]
    ids_n = [(k) for (k) in keys(am_ref_data[0][:bus])]
    y = current_year
    total_scarcity = 0.0

    # define days group index
    ids_d = am_ref_data[0][:repday_groups][day_group_id]["Day_Idx_List"]

    # -----------------
    # Dispatch
    # -----------------
    dispatch_data = Tuple{
        Int,      # day
        Int,      # hour
        Int,      # time
        Int,      # unit_id
        String,   # PLANT_NAME
        Int,      # bus_id
        String,   # bus_name
        String,   # Tech_ID
        String,   # UnitGroup
        String,   # Unit_Category
        String,   # Unit_Report_Label_1
        String,   # Unit_Report_Label_2
        Float64,  # u_G_iy
        Float64,  # u_ESE_iy
        Float64,  # ICAP
        Bool,     # status
        Float64,  # g_idhty
        Float64,  # chg_idhty
        Float64,  # soc_idhty
        Float64,  # g_ref_idht
        Float64,  # chg_ref_idhty
        Float64,  # soc_ref_idhty
        Bool      # risk_status
    }[]

    total_ICAP_loss_per_bus = zeros(length(ids_n), length(ids_d)*length(ids_h)*length(ids_t))

    for i in ids_i

        check_count = 1

        unit_id = i
        bus_idx = am_ref_data[0][:RA_gen_index][i]["bus_idx"]
        tech_idx = am_ref_data[0][:RA_gen_index][i]["genco_tech_id"]
        bus_name = am_ref_data[0][:bus][bus_idx]["bus_i"]
        bus_name = isa(bus_name, String) ? bus_name : string(bus_name)  # make sure that bus_name is a string
        
        Tech_ID = am_ref_data[bus_idx][:RA_gen_bus][tech_idx]["Tech_ID"]
        UnitGroup = am_ref_data[bus_idx][:RA_gen_bus][tech_idx]["UNITGROUP"]
        Unit_Category = am_ref_data[bus_idx][:RA_gen_bus][tech_idx]["UNIT_CATEGORY"]
        Unit_Report_Label_1 = am_ref_data[bus_idx][:RA_gen_bus][tech_idx]["UNIT_REPORT_LABEL_1"]
        Unit_Report_Label_2 = am_ref_data[bus_idx][:RA_gen_bus][tech_idx]["UNIT_REPORT_LABEL_2"]
        PLANT_NAME = haskey(am_ref_data[bus_idx][:RA_gen_bus][tech_idx], "PLANT_NAME") ? string(am_ref_data[bus_idx][:RA_gen_bus][tech_idx]["PLANT_NAME"]) : ""

        u_G_iy = am_ref_data[bus_idx][:RA_gen_bus][tech_idx]["EXUNITS"]
        ICAP = am_ref_data[bus_idx][:RA_gen_bus][tech_idx]["CAP"] * u_G_iy * pu_power_base

        for d in ids_d
            for h in ids_h
                for t in ids_t

                    idx_idhty = string("(", i, ", ", d, ", ", h, ", ", t, ", ", y, ")")

                    if haskey(result["dispatch"], idx_idhty)

                        day = d
                        hour = h
                        time = t

                        g_idhty = result["dispatch"][idx_idhty]["g_idhty"] * pu_power_base
                        g_ref_idht = daily_reference_ED_solutions["solution"]["dispatch"][idx_idhty]["g_idhty"] * pu_power_base
                        
                        if Unit_Category == "STORAGE"
                            u_ESE_iy = am_ref_data[bus_idx][:RA_gen_bus][tech_idx]["ES_MWh"] * pu_power_base
                            chg_idhty = result["storage"][idx_idhty]["chg_idhty"] * pu_power_base
                            chg_ref_idhty = daily_reference_ED_solutions["solution"]["storage"][idx_idhty]["chg_idhty"] * pu_power_base
                            soc_idhty = result["storage"][idx_idhty]["soc_idhty"] * pu_power_base
                            soc_ref_idhty = daily_reference_ED_solutions["solution"]["storage"][idx_idhty]["soc_idhty"] * pu_power_base
                            
                        else
                            u_ESE_iy = 0.0
                            chg_idhty = 0.0
                            chg_ref_idhty = 0.0
                            soc_idhty = 0.0
                            soc_ref_idhty = 0.0
                        end

                        status = risk_scenario_of_day[check_count, i]

                        status == false ? total_ICAP_loss_per_bus[bus_idx, check_count] += ICAP : nothing

                        risk_status = false
                        if length(ids_i) - sum(risk_scenario_of_day[check_count, :]) > 0
                            risk_status = true
                        end

                        push!(dispatch_data, (day, hour, time, unit_id, PLANT_NAME, bus_idx, bus_name, Tech_ID, UnitGroup, Unit_Category, Unit_Report_Label_1, Unit_Report_Label_2,
                                 u_G_iy, u_ESE_iy, ICAP, status, g_idhty, chg_idhty, soc_idhty,
                                 g_ref_idht, chg_ref_idhty, soc_ref_idhty, risk_status))
                        
                    else    # Default values when no data is found for the given idx_idhty

                        day = d
                        hour = h
                        time = t

                        g_idhty = 0.0
                        g_ref_idht = 0.0
                        u_ESE_iy = 0.0
                        chg_idhty = 0.0
                        chg_ref_idhty = 0.0
                        soc_idhty = 0.0
                        soc_ref_idhty = 0.0
                        
                        if Unit_Category == "STORAGE"
                            if haskey(daily_reference_ED_solutions["solution"], "expansion")
                                u_ESE_iy = daily_reference_ED_solutions["solution"]["expansion"][string("(", i, ", ", y, ")")]["u_ESE_iy"] * pu_power_base
                            else
                                u_ESE_iy = am_ref_data[bus_idx][:RA_gen_bus][tech_idx]["ES_MWh"] * pu_power_base
                            end
                        end
                        
                        status = risk_scenario_of_day[check_count, i]

                        risk_status = false
                        if length(ids_i) - sum(risk_scenario_of_day[check_count, :]) > 0
                            risk_status = true
                        end

                        push!(dispatch_data, (day, hour, time, unit_id, PLANT_NAME, bus_idx, bus_name, Tech_ID, UnitGroup, Unit_Category, Unit_Report_Label_1, Unit_Report_Label_2,
                                 u_G_iy, u_ESE_iy, ICAP, status, g_idhty, chg_idhty, soc_idhty,
                                 g_ref_idht, chg_ref_idhty, soc_ref_idhty, risk_status))
                    end

                    check_count += 1
                end
            end
        end

    end

    dispatch_df = DataFrame(dispatch_data, ra_report_names([:day, :hour, :time, :unit_id, :PLANT_NAME, :bus_id, :bus_name, :Tech_ID, :UnitGroup, :Unit_Category,
                            :Unit_Report_Label_1, :Unit_Report_Label_2, :u_G_iy, :u_ESE_iy, :ICAP, :status, :g_idhty, :chg_idhty,
                            :soc_idhty, :g_ref_idht, :chg_ref_idhty, :soc_ref_idhty, :risk_status]))
    dispatch_data = nothing

    if ELCC_flag == true
        file_name = string(case_name, "__RA_stage_", current_year, "_ELCC_iteration_", iter, "_day_id_", day_group_id, "_risk_id_", risk_id, "_dispatch.csv")
        output_path = get(RA_setting, "ELCC_output_path", output_path)
        CSV.write(joinpath(output_path, file_name), dispatch_df)
        dispatch_df = nothing
    else
        file_name = string(case_name, "__RA_stage_", current_year, "_day_id_", day_group_id, "_risk_id_", risk_id, "_dispatch.csv")
        CSV.write(joinpath(am_ref_data[0][:output_path], file_name), dispatch_df)
        dispatch_df = nothing 
    end


    # -----------------
    # System
    # -----------------
    system_data = Tuple{
        Int, # day
        Int, # hour
        Int, # time
        Int, # bus_id
        String, # bus_name
        Float64, # Scarcity
        Float64, # demand
        Float64  # total_gen_loss_ICAP
        }[]

    for n in ids_n
        
        check_count = 1
        bus_id = n 
        bus_name = am_ref_data[0][:bus][bus_id]["bus_i"]
        bus_name = isa(bus_name, String) ? bus_name : string(bus_name)  # make sure that bus_name is a string
        
        for d in ids_d
            for h in ids_h
                for t in ids_t

                    Demand = get_bus_demand_with_growth(am_ref_data, bus_id, d, h, t, y; nw=0, system_peak_scale) * pu_power_base

                    idx_string_ndhty = string("(", bus_id, ", ", d, ", ", h, ", ", t, ", ", y,")")
                    scarcity_E = 0.0
                    if haskey(result["scarcity"], idx_string_ndhty)
                        scarcity_E = result["scarcity"][idx_string_ndhty]["ens_ndhty"] * pu_power_base
                        total_scarcity += scarcity_E
                    end

                    total_gen_loss_ICAP = total_ICAP_loss_per_bus[bus_id, check_count]

                    push!(system_data, (d, h, t, bus_id, bus_name, scarcity_E, Demand, total_gen_loss_ICAP))

                    check_count += 1
                end
            end
        end
    end

    system_data_df = DataFrame(system_data, ra_report_names([:day, :hour, :time, :bus_id, :bus_name, :ens, :demand, :total_ICAP_loss]))
    system_data = nothing

    if ELCC_flag == true
        file_name = string(case_name, "__RA_stage_", current_year, "_ELCC_iteration_", iter, "_day_id_", day_group_id, "_risk_id_", risk_id, "_system.csv")
        CSV.write(joinpath(output_path, file_name), system_data_df)
        system_data_df = nothing         
    else
        file_name = string(case_name, "__RA_stage_", current_year, "_day_id_", day_group_id, "_risk_id_", risk_id, "_system.csv")
        CSV.write(joinpath(am_ref_data[0][:output_path], file_name), system_data_df)
        system_data_df = nothing 
    end
    

end


function report_result_dispatch_reference_RA(am_reference_data, am_setting_data, day_group_id, renewable_scenario_id, current_year, daily_reference_ED_solutions, RA_setting; ELCC_flag=false)

    case_name = get(am_setting_data["Simulation Configuration"], "Case_ID", "RA")

    am_ref_data = am_reference_data[:nw][0]
    system_peak_scale = RA_setting["system_peak_scale"]

    pu_power_base = am_setting_data["Simulation Setting"]["per_unit_base_value"]
    pu_econ_base  = am_setting_data["Simulation Setting"]["per_unit_econ_base_value"] / pu_power_base

    ids_h = [(h) for (h) in am_setting_data["run_H"]]
    ids_t = [(t) for (t) in am_setting_data["run_T"]]
    ids_n = [(k) for (k) in get_index(am_reference_data[:nw], :bus, 0)]
    ids_i = [(i) for (i) in get_index(am_reference_data[:nw], :RA_gen_index, 0)]
    ids_k = [(k) for (k) in get_index(am_reference_data[:nw], :branch, 0) if parameter(am_reference_data[:nw], 0, :branch, "model_flag", k) == true]
    y = current_year

    ids_d = am_ref_data[:repday_groups][day_group_id]["Day_Idx_List"]


    # remove power flow data
    delete!(daily_reference_ED_solutions["solution"], "powerflow")
    
    # -----------------
    # Dispatch
    # -----------------
    dispatch_data = Tuple{
        Int,      # day
        Int,      # hour
        Int,      # time
        Int,      # unit_id
        String,   # PLANT_NAME
        Int,      # bus_id
        String,   # bus_name
        String,   # Tech_ID
        String,   # UnitGroup
        String,   # Unit_Category
        String,   # Unit_Report_Label_1
        String,   # Unit_Report_Label_2
        Float64,  # u_G_iy
        Float64,  # u_ESE_iy
        Float64,  # ICAP
        Float64,  # MC
        Float64,  # g_idhty
        Float64,  # chg_idhty
        Float64,  # soc_idhty
    }[]

    for i in ids_i

        check_count = 1

        unit_id = i
        bus_idx = am_reference_data[:nw][0][:RA_gen_index][i]["bus_idx"]
        tech_idx = am_reference_data[:nw][0][:RA_gen_index][i]["genco_tech_id"]
        bus_name = am_reference_data[:nw][0][:bus][bus_idx]["bus_i"]
        bus_name = isa(bus_name, String) ? bus_name : string(bus_name)  # make sure that bus_name is a string
        
        Tech_ID = am_reference_data[:nw][bus_idx][:RA_gen_bus][tech_idx]["Tech_ID"]
        UnitGroup = am_reference_data[:nw][bus_idx][:RA_gen_bus][tech_idx]["UNITGROUP"]
        Unit_Category = am_reference_data[:nw][bus_idx][:RA_gen_bus][tech_idx]["UNIT_CATEGORY"]
        Unit_Report_Label_1 = am_reference_data[:nw][bus_idx][:RA_gen_bus][tech_idx]["UNIT_REPORT_LABEL_1"]
        Unit_Report_Label_2 = am_reference_data[:nw][bus_idx][:RA_gen_bus][tech_idx]["UNIT_REPORT_LABEL_2"]
        PLANT_NAME = haskey(am_reference_data[:nw][bus_idx][:RA_gen_bus][tech_idx], "PLANT_NAME") ? string(am_reference_data[:nw][bus_idx][:RA_gen_bus][tech_idx]["PLANT_NAME"]) : ""

        u_G_iy = am_reference_data[:nw][bus_idx][:RA_gen_bus][tech_idx]["EXUNITS"]
        ICAP = am_reference_data[:nw][bus_idx][:RA_gen_bus][tech_idx]["CAP"] * u_G_iy * pu_power_base
        
        for d in ids_d
            for h in ids_h
                for t in ids_t

                    idx_idhty = string("(", i, ", ", d, ", ", h, ", ", t, ", ", y, ")")
                    idx_idh = string("(", i, ", ", d, ", ", h, ")")
                    idx_string_dht = string("(", d, ", ", h, ", ", t, ")")
                    
                    day = d
                    hour = h
                    time = t

                    MC = parameter(am_reference_data[:nw], bus_idx, :RA_gen_bus, tech_idx, "Annual_MC")[string(y)][d] * pu_econ_base
                    g_ref_idht = daily_reference_ED_solutions["solution"]["dispatch"][idx_idhty]["g_idhty"] * pu_power_base
                    
                    if Unit_Category == "STORAGE"
                        u_ESE_iy = am_reference_data[:nw][bus_idx][:RA_gen_bus][tech_idx]["ES_MWh"] * pu_power_base
                        chg_idhty = daily_reference_ED_solutions["solution"]["storage"][idx_idhty]["chg_idhty"] * pu_power_base
                        soc_idhty = daily_reference_ED_solutions["solution"]["storage"][idx_idhty]["soc_idhty"] * pu_power_base
                        
                    else
                        u_ESE_iy = 0.0
                        chg_idhty = 0.0
                        soc_idhty = 0.0
                    end

                    push!(dispatch_data, (
                                        day,
                                        hour,
                                        time,

                                        unit_id,
                                        PLANT_NAME,
                                        bus_idx,
                                        bus_name,

                                        Tech_ID,
                                        UnitGroup,
                                        Unit_Category,
                                        Unit_Report_Label_1,
                                        Unit_Report_Label_2,

                                        u_G_iy,
                                        u_ESE_iy,
                                        ICAP,
                                        MC,

                                        g_ref_idht,
                                        chg_idhty,
                                        soc_idhty
                        ))
                end
            end
        end

    end


    dispatch_df = DataFrame(dispatch_data, ra_report_names([:day, :hour, :time, :unit_id, :PLANT_NAME, :bus_id, :bus_name, :Tech_ID, :UnitGroup, :Unit_Category,
                            :Unit_Report_Label_1, :Unit_Report_Label_2, :u_G_iy, :u_ESE_iy, :ICAP, :MC, :g_idhty, :chg_idhty, :soc_idhty]))
    dispatch_data = nothing 

    if ELCC_flag == true
        file_name = string(case_name, "__Reference_RA_stage_", current_year, "_day_id_", day_group_id, "_scenario_", renewable_scenario_id, "_dispatch.csv")
        output_path = get(RA_setting, "ELCC_output_path", am_reference_data[:nw][0][:output_path])
        CSV.write(joinpath(output_path, file_name), dispatch_df)
        dispatch_df = nothing

    else

        file_name = string(case_name, "__Reference_RA_stage_", current_year, "_day_id_", day_group_id, "_scenario_", renewable_scenario_id, "_dispatch.csv")
        CSV.write(joinpath(am_reference_data[:nw][0][:output_path], file_name), dispatch_df)
        dispatch_df = nothing
    end

    
    # -----------------
    # System
    # -----------------
    system_data = Tuple{
        Int, # day
        Int, # hour
        Int, # time
        Int, # bus_id
        String, # bus_name
        Float64, # Scarcity
        Float64 # demand
        }[]

    for n in ids_n
        
        bus_id = n 
        bus_name = am_reference_data[:nw][0][:bus][bus_id]["bus_i"]
        bus_name = isa(bus_name, String) ? bus_name : string(bus_name)  # make sure that bus_name is a string
        
        for d in ids_d
            for h in ids_h
                for t in ids_t

                    Demand = get_bus_demand_with_growth(am_reference_data[:nw], bus_id, d, h, t, y; nw=0, system_peak_scale) * pu_power_base

                    idx_string_ndhty = string("(", bus_id, ", ", d, ", ", h, ", ", t, ", ", y,")")
                    scarcity_E = daily_reference_ED_solutions["solution"]["scarcity"][idx_string_ndhty]["ens_ndhty"] * pu_power_base
                    
                    push!(system_data, (d, h, t, bus_id, bus_name, scarcity_E, Demand))

                end
            end
        end
    end

    system_data_df = DataFrame(system_data, ra_report_names([:day, :hour, :time, :bus_id, :bus_name, :ens, :demand]))
    system_data = nothing

    if ELCC_flag == true
        file_name = string(case_name, "__Reference_RA_stage_", current_year, "_day_id_", day_group_id, "_scenario_", renewable_scenario_id, "_system.csv")
        output_path = get(RA_setting, "ELCC_output_path", am_reference_data[:nw][0][:output_path])
        CSV.write(joinpath(output_path, file_name), system_data_df)
        system_data_df = nothing

    else

        file_name = string(case_name, "__Reference_RA_stage_", current_year, "_day_id_", day_group_id, "_scenario_", renewable_scenario_id, "_system.csv")
        CSV.write(joinpath(am_reference_data[:nw][0][:output_path], file_name), system_data_df)
        system_data_df = nothing
    end

    # remove scarcity data
    delete!(daily_reference_ED_solutions["solution"], "scarcity")
end


function build_LCO_GTEP_reference_ED_for_RA_instance!(am::Abstract_ALEAF_Model, risk_id::Int, day_group_id::Int, y::Int; nw::Int=am.cnw, constant_load=0.0, constant_load_bus_idx=nothing, constant_load_distribution="single_bus", constant_load_distribution_profile=Dict{Tuple{Int, Int, Int}, Dict{Int, Float64}}(), system_peak_scale=1.0, dispatch_report=false)

    # Reset a JuMP model
    if (am.setting["Solver Setting"]["solver_name"] == "CPLEX") & (am.setting["Solver Setting"]["1"]["Value"] == true)
        
        JuMP_model = cplex_direct_model(Val(:CPLEX))
    else
        JuMP_model = JuMP.Model()
    end

    const_name_flag = am.setting["Simulation Setting"]["const_name_flag"]
    if const_name_flag == false
        JuMP.set_string_names_on_creation(JuMP_model, false)
    end

    #------ run for each day group
    ids_d = [(d) for (d) in am.ref[:nw][0][:repday_groups][day_group_id]["Day_Idx_List"]]
    ids_y = [(y)]

    ###################################
    #------ Pre-processing
    ###################################
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

    ids_h = [(h) for (h) in am.setting["run_H"]]
    ids_t = [(h) for (h) in am.setting["run_T"]]
    ids_k = [(k) for (k) in get_index(am, :branch, 0) if parameter(am, 0, :branch, "model_flag", k) == true]
    ids_n = [(k) for (k) in get_index(am, :bus, 0)]
    ids_z = [(z) for (z) in get_index(am, :zone, 0, "reserve")]
    ids_p = [(z) for (z) in get_index(am, :zone, 0, "policy")]
    ids_lfl = [(lfl) for (lfl) in get_index(am, :demand, 0)]    # Large Flexible Load

    ###################################
    #------ DEFINE DECISION VARIABLES
    ###################################
    # Power flow variables
    if am.setting["Simulation Setting"]["power_flow_mode_flag"] in ["PTDF", "B-theta"]
        variable_f_kdhty_real(JuMP_model, am, :powerflow, "f_kdhty", risk_id, ids_k, ids_d, ids_h, ids_t, ids_y; bounded=false, report=dispatch_report)    # Power flow (unbounded for DC power flow)
    else # power_flow_mode_flag == "Network_Flow"
        variable_f_kdhty_real(JuMP_model, am, :powerflow, "f_kdhty", risk_id, ids_k, ids_d, ids_h, ids_t, ids_y; report=dispatch_report)                   # Power flow (bounded >=0 for network flow)
    end

    # Dispatch Variables
    variable_g_idhty_real(JuMP_model, am, :dispatch, "g_idhty", risk_id, ids_i, ids_d, ids_h, ids_t, ids_y; report=true)                  # Unit Generation (MW)
    define_variable_idhty_real(JuMP_model, am, :dispatch, "redispatch_abs_dev_idhty", risk_id, ids_i, ids_d, ids_h, ids_t, ids_y; bounded_lower=true, lower_bound=0.0, report=dispatch_report)
    variable_curt_idhty_real(JuMP_model, am, :curtailment, "curt_idhty", risk_id, ids_i, ids_d, ids_h, ids_t, ids_y; report=false)         # Curtailment (MW)

    # Scarcity variables
    variable_ens_ndhty_real(JuMP_model, am, :scarcity, "ens_ndhty", risk_id, ids_n, ids_d, ids_h, ids_t, ids_y; report=dispatch_report)                              # Energy Not Served (MW).
    
    # Commitment Variables
    if am.setting["Simulation Configuration"]["Dispatch_Mode_in_OP"] == "Unit Commitment"
        variable_c_idhty_integer(JuMP_model, am, :commitment, "c_idhty", risk_id, ids_i_commit, ids_d, ids_h, ids_t, ids_y; report=dispatch_report)                     # Number of Units Committed
        variable_su_idhty_real(JuMP_model, am, :commitment, "su_idhty", risk_id, ids_i_commit, ids_d, ids_h, ids_t, ids_y; report=dispatch_report)               # Number of Units Started Up
    end
    
    # Storage
    variable_soc_idhty_real(JuMP_model, am, :storage, "soc_idhty", risk_id, ids_i_sto, ids_d, ids_h, ids_t, ids_y; report=true)            # Storage charge level at end of the period
    variable_chg_idhty_real(JuMP_model, am, :storage, "chg_idhty", risk_id, ids_i_sto, ids_d, ids_h, ids_t, ids_y; report=true)            # Storage charging MW
    variable_sto_c_idhty_integer(JuMP_model, am, :storage_commitment, "sto_c_idhty", risk_id, ids_i, ids_d, ids_h, ids_t, ids_y; report=dispatch_report)               

    # Large Flexible Load
    define_variable_idhty_real(JuMP_model, am, :demand, "lfl_ldhty", risk_id, ids_lfl, ids_d, ids_h, ids_t, ids_y; bounded_lower=true, lower_bound=0, report=true)  
    define_variable_idhty_real(JuMP_model, am, :demand, "lfl_DR_ldhty", risk_id, ids_lfl, ids_d, ids_h, ids_t, ids_y; bounded_lower=true, lower_bound=0, report=dispatch_report)  
    variable_lfl_seg_lsdhty_real(JuMP_model, am, :demand, "lfl_seg_lsdhty", risk_id, ids_lfl, ids_d, ids_h, ids_t, ids_y; bounded_lower=true, lower_bound=0, report=dispatch_report)  
    variable_lfl_ind_lsdhty_binary(JuMP_model, am, :demand, "lfl_ind_lsdhty", risk_id, ids_lfl, ids_d, ids_h, ids_t, ids_y, report=dispatch_report)  

    # Define LFL onsite hybrid generation variables if applicable
    for lfl in ids_lfl 
        
        if am.ref[:nw][0][:demand][lfl]["Hybrid_Gen"] != "NA"
            define_variable_idhty_real(JuMP_model, am, :demand, "lfl_g_G_LFL_ldhty", risk_id, ids_lfl, ids_d, ids_h, ids_t, ids_y; bounded_lower=true, lower_bound=0, report=dispatch_report)  
            define_variable_idhty_real(JuMP_model, am, :demand, "lfl_g_G_Grid_ldhty", risk_id, ids_lfl, ids_d, ids_h, ids_t, ids_y; bounded_lower=true, lower_bound=0, report=dispatch_report)  
            define_variable_idhty_real(JuMP_model, am, :demand, "lfl_g_G_ES_ldhty", risk_id, ids_lfl, ids_d, ids_h, ids_t, ids_y; bounded_lower=true, lower_bound=0, report=dispatch_report)  
        end

        if am.ref[:nw][0][:demand][lfl]["Hybrid_ES"] != "NA"
            define_variable_idhty_real(JuMP_model, am, :demand, "lfl_g_ES_LFL_ldhty", risk_id, ids_lfl, ids_d, ids_h, ids_t, ids_y; bounded_lower=true, lower_bound=0, report=dispatch_report)  
            define_variable_idhty_real(JuMP_model, am, :demand, "lfl_g_ES_Grid_ldhty", risk_id, ids_lfl, ids_d, ids_h, ids_t, ids_y; bounded_lower=true, lower_bound=0, report=dispatch_report)  
            define_variable_idhty_real(JuMP_model, am, :demand, "lfl_chg_Grid_ES_ldhty", risk_id, ids_lfl, ids_d, ids_h, ids_t, ids_y; bounded_lower=true, lower_bound=0, report=dispatch_report)   
            define_variable_idhty_real(JuMP_model, am, :demand, "lfl_soc_ldhty", risk_id, ids_lfl, ids_d, ids_h, ids_t, ids_y; bounded_lower=true, lower_bound=0, report=dispatch_report)   
        end
    end

    # Hybrid Generation to ES variable (report=true: post-contingency SOC balance reads it from ref_ED).
    define_variable_idhty_real(JuMP_model, am, :hybrid, "g_G_ES_idhty", risk_id, ids_i_hybrid_gen, ids_d, ids_h, ids_t, ids_y; bounded_lower=true, lower_bound=0, report=true)
    

    ###################################
    #------ DEFINE Constraints
    ###################################

    # System Balance & power flow
    if am.setting["Simulation Setting"]["power_flow_mode_flag"] == "PTDF"
        
        variable_p_inj_ndhty_real(JuMP_model, am, :powerflow, "p_inj_ndhty", risk_id, ids_n, ids_d, ids_h, ids_t, ids_y)                              # power injection at node n
        variable_demand_ndhty_real(JuMP_model, am, :powerflow, "demand_ndhty", risk_id, ids_n, ids_d, ids_h, ids_t, ids_y)                              # delivered demand

        for d in ids_d, h in ids_h, t in ids_t

            constraint_sum_p_injdhty_RA(JuMP_model, am, "constraint_sum_p_injdhty_RA", risk_id, d, h, t, y)

            for n in ids_n
                constraint_PTDF_Power_injection_ndhty_real(JuMP_model, am, "constraint_PTDF_Power_injection_ndhty_real", risk_id, n, d, h, t, y)
                constraint_PTDF_LoadBalance_ndhty_real_RA(JuMP_model, am, "constraint_PTDF_LoadBalance_ndhty_real_RA", risk_id, n, d, h, t, y; constant_load, constant_load_bus_idx, constant_load_distribution, constant_load_distribution_profile, system_peak_scale)
            end

            for k in ids_k
                constraint_dc_power_flow_kdhty_RA(JuMP_model, am, "constraint_dc_power_flow_kdhty_RA", risk_id, k, d, h, t, y; constant_load, system_peak_scale)
                constraint_dc_power_flow_max_kdhty_RA(JuMP_model, am, "constraint_dc_power_flow_max_kdhty_RA", risk_id, k, d, h, t, y)
                constraint_dc_power_flow_min_kdhty_RA(JuMP_model, am, "constraint_dc_power_flow_min_kdhty_RA", risk_id, k, d, h, t, y)
            end
        end
    elseif am.setting["Simulation Setting"]["power_flow_mode_flag"] == "Network_Flow"
        for d in ids_d, h in ids_h, t in ids_t
            for n in ids_n
                constraint_LoadBalance_ndhty_real_RA(JuMP_model, am, "constraint_LoadBalance_ndhty_real_RA", risk_id, n, d, h, t, y; constant_load, constant_load_bus_idx, constant_load_distribution, constant_load_distribution_profile, system_peak_scale)
            end
                    
            for k in ids_k
                constraint_power_flow_max_kdhty_OP(JuMP_model, am, "constraint_power_flow_max_kdhty_OP", risk_id, k, d, h, t, y; const_name_flag)
                constraint_power_flow_min_kdhty_OP(JuMP_model, am, "constraint_power_flow_min_kdhty_OP", risk_id, k, d, h, t, y; const_name_flag)
            end
        end
    
    elseif am.setting["Simulation Setting"]["power_flow_mode_flag"] == "B-theta"

        # One slack per AC island so DC ties don't tie separate synchronous areas together.
        ac_ref_bus_ids = get_ac_reference_buses(am; nw=0)
        variable_bus_angle_ndhty_real(JuMP_model, am, :powerflow, "bus_angle_ndhty", risk_id, ids_n, ids_d, ids_h, ids_t, ids_y; ref_bus_ids = ac_ref_bus_ids)                              # power injection at node n

        for d in ids_d, h in ids_h, t in ids_t
            for n in ids_n
                constraint_LoadBalance_ndhty_real_RA(JuMP_model, am, "constraint_LoadBalance_ndhty_real_RA", risk_id, n, d, h, t, y; constant_load, constant_load_bus_idx, constant_load_distribution, constant_load_distribution_profile, system_peak_scale)
            end

            for k in ids_k
                # DC ties carry no B-theta angle coupling; flow limits still apply to all branches.
                if get(am.ref[:nw][0][:branch][k], "dc_line", false) != true
                    constraint_b_theta_power_flow_kdhty_RA(JuMP_model, am, "constraint_b_theta_power_flow_kdhty_RA", risk_id, k, d, h, t, y)
                end
                constraint_dc_power_flow_max_kdhty_RA(JuMP_model, am, "constraint_dc_power_flow_max_kdhty_RA", risk_id, k, d, h, t, y)
                constraint_dc_power_flow_min_kdhty_RA(JuMP_model, am, "constraint_dc_power_flow_min_kdhty_RA", risk_id, k, d, h, t, y)
            end
        end

    end

    # Unit Dispatch and commitment
    for i in ids_i

        unit_category = parameter(am, 0, :gen_index, "UNIT_CATEGORY", i)

        for d in ids_d, h in ids_h, t in ids_t
            
            if (i in ids_i_commit) && (am.setting["Simulation Configuration"]["Dispatch_Mode_in_OP"] == "Unit Commitment") # if unit i requires a commitment decision
            
                # Unit commitment constraints
                constraint_Commit_Limit_idhty(JuMP_model, am, "constraint_CommitLimit_OP_idhty", risk_id, i, d, h, t, y; const_name_flag)
                constraint_Start_Up_Limit_idhty(JuMP_model, am, "constraint_Start_Up_Limit_OP_idhty", risk_id, i, d, h, t, y; const_name_flag)
                constraint_Start_Up_Status_Dn_idhty(JuMP_model, am, "constraint_Start_Up_Status_Dn_idhty", risk_id, i, d, h, t, y; const_name_flag)
                constraint_Start_Up_Status_Up_idhty(JuMP_model, am, "constraint_Start_Up_Status_Up_idhty", risk_id, i, d, h, t, y; const_name_flag)

                # Unit dispatch constraints
                constraint_TherMin_Dispatch_UC_NoReserve_idhty(JuMP_model, am, "constraint_TherMin_Dispatch_UC_NoReserve_idhty", risk_id, i, d, h, t, y; const_name_flag)
                constraint_TherMax_Dispatch_UC_NoReserve_idhty(JuMP_model, am, "constraint_TherMax_Dispatch_UC_NoReserve_idhty", risk_id, i, d, h, t, y; const_name_flag)
                
                # Inter-temporal Ramping Constraints (with commitment)
                if unit_category in ["THERMAL", "NUCLEAR", "OTHER"]
                    constraint_RampUp_InterTemporal_Hour_UC_idhy(JuMP_model, am, "constraint_RampUp_InterTemporal_Hour_UC_idhy", risk_id, i, d, h, t, y; const_name_flag, apply_PMAX=false)
                    constraint_RampDn_InterTemporal_Hour_UC_idhy(JuMP_model, am, "constraint_RampDn_InterTemporal_Hour_UC_idhy", risk_id, i, d, h, t, y; const_name_flag, apply_PMAX=false)
                end

            else
                
                # Unit dispatch constraint for thermal resources
                constraint_TherMaxDispatch_NoReserve_RA_idhty(JuMP_model, am, "constraint_TherMaxDispatch_NoReserve_RA_idhty", risk_id, i, d, h, t, y)
                    
                # Inter-temporal Ramping Constraints for thermal resources (w/o commitment)
                if unit_category in ["THERMAL", "NUCLEAR", "OTHER"]
                    constraint_RampUp_InterTemporal_Hour_idhy_RA(JuMP_model, am, "constraint_RampUp_InterTemporal_Hour_idhy_RA", risk_id, i, d, h, t, y, day_group_id; const_name_flag, apply_PMAX=false)
                    constraint_RampDn_InterTemporal_Hour_idhy_RA(JuMP_model, am, "constraint_RampDn_InterTemporal_Hour_idhy_RA", risk_id, i, d, h, t, y, day_group_id; const_name_flag, apply_PMAX=false)
                end

            end 
        end
    end

    # VRE Balance and Budget
    for i in ids_i

        bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
        tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)

        if parameter(am, bus_idx, :gen_bus, tech_idx, "VRE_Flag") == true

            unit_group = parameter(am, 0, :gen_index, "UNIT_GROUP", i)
            fuel_limit = parameter(am, bus_idx, :gen_bus, tech_idx, "FUEL_LIMIT")

            if (fuel_limit == "Fixed Profile")
                for d in ids_d, h in ids_h, t in ids_t
                    constraint_VREBalance_Fixed_Profile_NoReserve_idhty_RA(JuMP_model, am, "constraint_VREBalance_Fixed_Profile_NoReserve_idhty_RA", risk_id, i, d, h, t, y, unit_group; apply_PMAX=false)
                end

            elseif fuel_limit != "NA"

                if am.setting["Simulation Configuration"]["Hydro_Flexibility_Flag"] == true
                    for d in ids_d, h in ids_h, t in ids_t
                        constraint_VREBalance_Fixed_Profile_wifh_Flexibilty_NoReserve_idhty_RA(JuMP_model, am, "constraint_VREBalance_Fixed_Profile_wifh_Flexibilty_NoReserve_idhty_RA", risk_id, i, d, h, t, y, unit_group; const_name_flag, apply_PMAX=false)
                    end
                end

                if am.setting["Simulation Configuration"]["Hydro_Budget_Flag"] == true
                    constraint_VREBalance_Budget_NonReserve_iy_RA(JuMP_model, am, "constraint_VREBalance_Budget_NonReserve_iy_RA", risk_id, i, day_group_id, y, unit_group; const_name_flag)
                end
            end
        end

    end


    # Storage Balance and Power Constraints
    for i in ids_i_sto
        hybrid_type = am.ref[:nw][0][:gen_index][i]["hybrid_type"] 

        if hybrid_type == "ES"      # limit hybrid ES grid charge
            constraint_hybrid_ES_Charge_from_Grid_idhty(JuMP_model, am, "constraint_hybrid_ES_Charge_from_Grid_idhty", risk_id, i, ids_d, ids_h, ids_t, ids_y; const_name_flag)
        end

        for d in ids_d, h in ids_h, t in ids_t
            
            # charge limit
            if hybrid_type == "ES"
                constraint_hybrid_ES_Charge_Max_NonReserve_idhty(JuMP_model, am, "constraint_hybrid_ES_Charge_Max_NonReserve_idhty", risk_id, i, d, h, t, y; const_name_flag)
            else
                constraint_ES_Charge_Max_NonReserve_idhty_RA(JuMP_model, am, "constraint_ES_Charge_Max_NonReserve_idhty_RA", risk_id, i, d, h, t, y)
            end
                
            # SOC bounds
            constraint_ES_SOC_Max_NonReserve_idhty_RA(JuMP_model, am, "constraint_ES_SOC_Max_NonReserve_idhty_RA", risk_id, i, d, h, t, y)
            constraint_ES_SOC_Min_NonReserve_idhty_RA(JuMP_model, am, "constraint_ES_SOC_Min_NonReserve_idhty_RA", risk_id, i, d, h, t, y)
            
            # Additional SOC bound for SATA
            constraint_ES_SOC_Min_NonReserve_SATA_idhty(JuMP_model, am, "constraint_ES_SOC_Min_NonReserve_SATA_idhty", risk_id, i, d, h, t, y; const_name_flag)

            # Storage commitment constraints
            constraint_ES_Charge_Max_Sto_UC_NonReserve_idhty(JuMP_model, am, "constraint_ES_Charge_Max_Sto_UC_NonReserve_idhty", risk_id, i, d, h, t, y; const_name_flag)
            constraint_ES_DisCharge_Max_Sto_UC_NonReserve_idhty(JuMP_model, am, "constraint_ES_DisCharge_Max_Sto_UC_NonReserve_idhty", risk_id, i, d, h, t, y; const_name_flag)
            
            if hybrid_type == "ES"
                constraint_hybrid_ES_Charge_Max_Sto_UC_NoReserve_idhty(JuMP_model, am, "constraint_hybrid_ES_Charge_Max_Sto_UC_NoReserve_idhty", risk_id, i, d, h, t, y; const_name_flag)
            else
                constraint_ES_Sto_UC_Limit_idhty_RA(JuMP_model, am, "constraint_ES_Sto_UC_Limit_idhty_RA", risk_id, i, d, h, t, y; const_name_flag)
            end
        end

    end


    # Storage SOC Balance
    for d in ids_d
        for i in ids_i_sto
            hybrid_type = am.ref[:nw][0][:gen_index][i]["hybrid_type"] 

            if hybrid_type == "ES"
                for h in ids_h, t in ids_t
                    constraint_hybrid_ES_SOC_Balance_Inter_Hour_idhty(JuMP_model, am, "constraint_hybrid_ES_SOC_Balance_Inter_Hour_idhty", risk_id, i, d, h, t, y, day_group_id; const_name_flag)
                end
            else
                for h in ids_h, t in ids_t
                    constraint_ES_SOC_Balance_Inter_Hour_NoReserve_with_SATA_idhty(JuMP_model, am, "constraint_ES_SOC_Balance_Inter_Hour_NoReserve_with_SATA_idhty", risk_id, i, d, h, t, y, day_group_id)
                end 
            end
        end
    end


    # Storage SOC Neutral
    for d in ids_d
        for i in ids_i_sto
            hybrid_type = am.ref[:nw][0][:gen_index][i]["hybrid_type"] 
            end_day = am.ref[:nw][0][:repday_groups][day_group_id]["End_Day_Id"]
            start_day_idx = ids_d[1]
            if am.ref[:nw][0][:repdays][d]["Day"] == end_day
                if hybrid_type == "ES"
                    constraint_hybrid_ES_SOC_Neutral_idy(JuMP_model, am, "constraint_hybrid_ES_SOC_Neutral_idy", risk_id, i, d, y, start_day_idx; const_name_flag)
                else
                    constraint_ES_SOC_Neutral_idy(JuMP_model, am, "constraint_ES_SOC_Neutral_idy", risk_id, i, d, y, start_day_idx; const_name_flag)
                end
            end
        end
    end

    # Large Flexible Load
    for d in ids_d, h in ids_h, t in ids_t
        for lfl in ids_lfl

            constraint_LFL_power_balance_ldhty(JuMP_model, am, "constraint_LFL_power_balance_ldhty", risk_id, lfl, d, h, t, y; const_name_flag)
            constraint_LFL_inter_connection_limit_injection_ldhty(JuMP_model, am, "constraint_LFL_inter_connection_limit_injection_ldhty", risk_id, lfl, d, h, t, y; const_name_flag)
            constraint_LFL_inter_connection_limit_withdraw_ldhty(JuMP_model, am, "constraint_LFL_inter_connection_limit_withdraw_ldhty", risk_id, lfl, d, h, t, y; const_name_flag)
            
            constraint_LFL_segment_bound_lsdhty(JuMP_model, am, "constraint_LFL_segment_bound_lsdhty", risk_id, lfl, d, h, t, y; const_name_flag)
            constraint_LFL_segment_relation_lsdhty(JuMP_model, am, "constraint_LFL_segment_relation_lsdhty", risk_id, lfl, d, h, t, y; const_name_flag)

            constraint_lfl_DR_balance_ldhty(JuMP_model, am, "constraint_lfl_DR_balance_ldhty", risk_id, lfl, d, h, t, y; const_name_flag)

            if am.ref[:nw][0][:demand][lfl]["Hybrid_Gen"] != "NA"
                constraint_lfl_onsite_gen_thermal_cap_ldhty(JuMP_model, am, "constraint_lfl_onsite_gen_thermal_cap_ldhty", risk_id, lfl, d, h, t, y; const_name_flag)
            end

            if am.ref[:nw][0][:demand][lfl]["Hybrid_ES"] != "NA"
                constraint_lfl_onsite_ES_charge_cap_ldhty(JuMP_model, am, "constraint_lfl_onsite_ES_charge_cap_ldhty", risk_id, lfl, d, h, t, y; const_name_flag)
                constraint_lfl_onsite_ES_discharge_cap_ldhty(JuMP_model, am, "constraint_lfl_onsite_ES_discharge_cap_ldhty", risk_id, lfl, d, h, t, y; const_name_flag)
                constraint_lfl_onsite_ES_SOC_cap_ldhty(JuMP_model, am, "constraint_lfl_onsite_ES_SOC_cap_ldhty", risk_id, lfl, d, h, t, y; const_name_flag)
            end
        end
    end

    # Large Flexible Load DR Daily Limit
    for d in ids_d
        for lfl in ids_lfl
            constraint_lfl_DR_daily_limit_ldhty(JuMP_model, am, "constraint_lfl_DR_daily_limit_ldhty", risk_id, lfl, d, y; const_name_flag)

        end
    end

    # Large Flexible Load Onsite Storage SOC Balance
    for lfl in ids_lfl
        if am.ref[:nw][0][:demand][lfl]["Hybrid_ES"] != "NA"
            for d in am.ref[:nw][0][:repday_groups][day_group_id]["Day_Idx_List"]
            
                for h in ids_h, t in ids_t
                    constraint_lfl_onsite_ES_SOC_balance_ldhty(JuMP_model, am, "constraint_lfl_onsite_ES_SOC_balance_ldhty", risk_id, lfl, d, h, t, y, day_group_id; const_name_flag)
                end

                if am.ref[:nw][0][:repdays][d]["Day"] == am.ref[:nw][0][:repday_groups][day_group_id]["End_Day_Id"]
                    start_day_idx = am.ref[:nw][0][:repday_groups][day_group_id]["Day_Idx_List"][1]
                    constraint_lfl_onsite_ES_SOC_neutral_ldhty(JuMP_model, am, "constraint_lfl_onsite_ES_SOC_neutral_ldhty", risk_id, lfl, d, y, start_day_idx; const_name_flag)
                end 

            end
        end
    end

    # Hybrid Plant Generation Constraint
    for i in ids_i_hybrid_gen
        for d in ids_d, h in ids_h, t in ids_t
            constraint_hybrid_onsite_gen_thermal_cap_ldhty(JuMP_model, am, "constraint_hybrid_onsite_gen_thermal_cap_ldhty", risk_id, i, d, h, t, y; const_name_flag, reserve_flag=false)
            constraint_hybrid_inter_connection_limit_injection_ldhty(JuMP_model, am, "constraint_hybrid_inter_connection_limit_injection_ldhty", risk_id, i, d, h, t, y; const_name_flag, reserve_flag=false)
            constraint_hybrid_inter_connection_limit_withdraw_ldhty(JuMP_model, am, "constraint_hybrid_inter_connection_limit_withdraw_ldhty", risk_id, i, d, h, t, y; const_name_flag, reserve_flag=false)
        end
    end

    # Policy Constraints

    # Storage annual energy throughput limits
    if am.setting["Simulation Configuration"]["Energy_Storage_AET_Limit_Flag"] == true
        for i in ids_i_sto
            constraint_ES_AET_OP_RA_y(JuMP_model, am, "constraint_ES_AET_y", risk_id, i, y, ids_d; reserve_flag=false)
        end
    end

    # Clean Energy Generation Constraint
    if am.setting["Simulation Configuration"]["Clean_Energy_Generation_Target_OP_Flag"] == true
        # slack_CEG_ndy must exist before the constraint references it
        variable_slack_CEG_ndy_real(JuMP_model, am, :slack, "slack_CEG_ndy", risk_id, ids_p, ids_d, ids_y)
        if am.ref[:nw][0][:planning_stages][y]["year"] >= am.setting["Simulation Configuration"]["Clean_Energy_Generation_Target_Start_Year"]
            for n in ids_p
                for d in ids_d
                    constraint_clean_energy_generation_ndy(JuMP_model, am, "constraint_clean_energy_generation_ndy", risk_id, n, d, y)
                end
            end
        end
    end

    ###################################
    #------ DEFINE Objective
    ###################################
    ED_objective_function_operation_for_RA(JuMP_model, am, day_group_id, ids_i, ids_d, ids_y, ids_i_commit, ids_i_hybrid_gen, ids_h, ids_t, ids_k, ids_n, ids_z, ids_p, ids_lfl, risk_id)

    am.model[:nw][nw][risk_id] = JuMP_model
    
end


function ED_objective_function_operation_for_RA(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, day_group_id, ids_i, ids_d, ids_y, ids_i_commit, ids_i_hybrid_gen, ids_h, ids_t, ids_k, ids_n, ids_z, ids_p, ids_lfl, risk_id; nw::Int=am.cnw, report::Bool=true)
    
    objective = JuMP.AffExpr(0.0)    

    y = ids_y[1]
    future_year = am.ref[:nw][0][:planning_stages][y]["year"]

    # Random adder
    Random.seed!(12345)
    epsilon = 1e-7

    ids_order = sort(collect(ids_i))
    jitter = Dict(id => rand() * epsilon for id in ids_order)

    # Hourly Costs, divide by number of timesteps in an hour (card(T))
    inv_card_T = 1 / length(am.setting["run_T"])

    # LFL demand response cost
    for lfl in ids_lfl
        local_lfl = am.ref[:nw][0][:demand][lfl]
        num_seg = local_lfl["Num_DR_Segments"]

        # skip s==1 by construction
        for s in 2:num_seg
            dr_cost = local_lfl["Price_$(s)"]  # computed once per s

            for d in ids_d 
                num_days = parameter(am, 0, :repdays, "NumDays", d)
                coef = -(num_days * dr_cost * inv_card_T)
                for h in ids_h, t in ids_t
                    JuMP.add_to_expression!(objective, coef, variable(am, nw, risk_id, :lfl_seg_lsdhty, (lfl,s,d,h,t,y)))  
                end
            end
        end

        # Ancillary penalties for hybrid ES grid charging
        if am.ref[:nw][0][:demand][lfl]["Hybrid_ES"] != "NA"
            for d in ids_d, h in ids_h, t in ids_t
                JuMP.add_to_expression!(objective, 0.0000001, variable(am, nw, risk_id, :lfl_chg_Grid_ES_ldhty, (lfl,d,h,t,y)))  
            end
        end
    end

    for i in ids_i, d in ids_d, h in ids_h, t in ids_t
        bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
        tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)
        CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")
        unit_category = parameter(am, 0, :gen_index, "UNIT_CATEGORY", i)

        # Marginal cost
        MC = parameter(am, bus_idx, :gen_bus, tech_idx, "Annual_MC")[string(y)][d] + jitter[i]

        if haskey(am.ref[:nw][:0][:gen_index][i], "new_unit")   # ensure that the added unit is dispatched first
            if am.ref[:nw][:0][:gen_index][i]["new_unit"] == true
                if MC > 0 
                    MC = MC * 0.9
                else
                    MC = -0.0001 
                end
            end
        end

        # Generation Cost
        JuMP.add_to_expression!(objective, MC, variable(am, nw, risk_id, :g_idhty, (i,d,h,t,y)))  
        
    end

    # Start up and no load cost
    if am.setting["Simulation Configuration"]["Dispatch_Mode_in_OP"] == "Unit Commitment"
        for i in ids_i, d in ids_d, h in ids_h
            if (i in ids_i_commit)
                bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
                tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)
                CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")
                
                # Startup Cost ($/MW)
                JuMP.add_to_expression!(objective, CAP * parameter(am, bus_idx, :gen_bus, tech_idx, "SUC"), variable(am, nw, risk_id, :su_idhty, (i,d,h,1,y)))  

                # No Load Cost ($/MW)
                JuMP.add_to_expression!(objective, CAP * parameter(am, bus_idx, :gen_bus, tech_idx, "NLC") * inv_card_T, variable(am, nw, risk_id, :c_idhty, (i,d,h,1,y)))  
            end
        end
    end

    ids_n_order = sort(collect(ids_n))
    n_jitter = Dict(id => rand() * epsilon for id in ids_n_order)
    
    # Load shedding
    for n in ids_n, d in ids_d, h in ids_h, t in ids_t
        JuMP.add_to_expression!(objective, (am.setting["Simulation Configuration"]["VOLL"] + n_jitter[n]) * inv_card_T, variable(am, nw, risk_id, :ens_ndhty, (n,d,h,t,y)))
    end

    # Emissions
    for i in ids_i, d in ids_d, h in ids_h, t in ids_t
        bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
        tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)
        Emission_rate = parameter(am, bus_idx, :gen_bus, tech_idx, "Emission_CO2")
        # Emissions Cost ($/MW) = $/tonne * kg/MMbtu * MMbtu/MWh * MWh * 1 tonne/1000kg = $/tonne * metric tonnes CO2 per timestep

        JuMP.add_to_expression!(objective, am.setting["Simulation Configuration"]["CTAX"] * Emission_rate * inv_card_T, variable(am, nw, risk_id, :g_idhty, (i,d,h,t,y)))  
    end

     # Ancillary penalties for storage commitment and hybrid ES grid charging
    for i in ids_i
        bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
        tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)
        hybrid_type = am.ref[:nw][0][:gen_index][i]["hybrid_type"] 

        if parameter(am, 0, :gen_index, "UNIT_CATEGORY", i) == "STORAGE"
            if parameter(am, bus_idx, :gen_bus, tech_idx, "Storage Commitment") == true
                for d in ids_d, h in ids_h, t in ids_t
                    JuMP.add_to_expression!(objective, 0.0000001, variable(am, nw, risk_id, :sto_c_idhty, (i,d,h,t,y)))    
                end 
            end

            if hybrid_type == "ES"
                for d in ids_d, h in ids_h, t in ids_t
                    JuMP.add_to_expression!(objective, 0.0000001, variable(am, nw, risk_id, :chg_idhty, (i,d,h,t,y)))    
                end 
            end
        end
    end


    # DEFINE OBJECTIVE FUNCTION
    return JuMP.@objective(JuMP_model, Min, objective     
    )

end


function build_ra_imperfect_foresight_sequential_snapshot_ed_model!(am::Abstract_ALEAF_Model, risk_id::Int, event_idx::Int, event_day::Int, event_hour::Int, event_time::Int, day_group_id::Int, y::Int, gen_risk_matrix, ref_ED_solution, risk_ED_solution, RA_setting, event_first_hour_flag, daygroup_first_hour_flag, simulated_event_hours; nw::Int=am.cnw, constant_load=0.0, constant_load_bus_idx=nothing, constant_load_distribution="single_bus", constant_load_distribution_profile=Dict{Tuple{Int, Int, Int}, Dict{Int, Float64}}(), system_peak_scale=1.0, dispatch_report=false)

    # Reset a JuMP model
    if (am.setting["Solver Setting"]["solver_name"] == "CPLEX") & (am.setting["Solver Setting"]["1"]["Value"] == true)
        
        JuMP_model = cplex_direct_model(Val(:CPLEX))
    else
        JuMP_model = JuMP.Model()
    end

    const_name_flag = am.setting["Simulation Setting"]["const_name_flag"]
    if const_name_flag == false
        JuMP.set_string_names_on_creation(JuMP_model, false)
    end

    am.var = Dict{Symbol,Any}(:nw => Dict{Int,Any}())
    am.con = Dict{Symbol,Any}(:nw => Dict{Int,Any}())
    am.sol = Dict{Symbol,Any}(:nw => Dict{Int,Any}())

    for (nw_id, nw) in am.ref[:nw]
        am.var[:nw][nw_id] = Dict{Int,Any}()
        am.con[:nw][nw_id] = Dict{Int,Any}()
        am.sol[:nw][nw_id] = Dict{Int,Any}()
        
        am.var[:nw][nw_id][risk_id] = Dict{Symbol,Any}()
        am.con[:nw][nw_id][risk_id] = Dict{Symbol,Any}()
        am.sol[:nw][nw_id][risk_id] = Dict{Symbol,Any}()
        
    end

    #------ run for each day group
    ids_day_group = [(d) for (d) in am.ref[:nw][0][:repday_groups][day_group_id]["Day_Idx_List"]]
    ids_d = [event_day]
    ids_h = [event_hour]
    ids_t = [event_time]
    ids_y = [(y)]

    ###################################
    #------ Pre-processing
    ###################################
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

    ids_k = [(k) for (k) in get_index(am, :branch, 0) if parameter(am, 0, :branch, "model_flag", k) == true]
    ids_n = [(k) for (k) in get_index(am, :bus, 0)]
    ids_z = [(z) for (z) in get_index(am, :zone, 0, "reserve")]
    ids_p = [(z) for (z) in get_index(am, :zone, 0, "policy")]
    ids_lfl = [(lfl) for (lfl) in get_index(am, :demand, 0)]    # Large Flexible Load

    ids_dht = [(d,h,t) for d in ids_d for h in ids_h for t in ids_t]

    ###################################
    #------ DEFINE DECISION VARIABLES
    ###################################
    # Power flow variables
    if am.setting["Simulation Setting"]["power_flow_mode_flag"] in ["PTDF", "B-theta"]
        variable_f_kdhty_real(JuMP_model, am, :powerflow, "f_kdhty", risk_id, ids_k, ids_d, ids_h, ids_t, ids_y; bounded=false, report=dispatch_report)    # Power flow (unbounded for DC power flow)
    else # power_flow_mode_flag == "Network_Flow"
        variable_f_kdhty_real(JuMP_model, am, :powerflow, "f_kdhty", risk_id, ids_k, ids_d, ids_h, ids_t, ids_y; report=dispatch_report)                   # Power flow (bounded >=0 for network flow)
    end

    # Dispatch Variables
    variable_g_idhty_real(JuMP_model, am, :dispatch, "g_idhty", risk_id, ids_i, ids_d, ids_h, ids_t, ids_y; report=true)                  # Unit Generation (MW)
    define_variable_idhty_real(JuMP_model, am, :dispatch, "redispatch_abs_dev_idhty", risk_id, ids_i, ids_d, ids_h, ids_t, ids_y; bounded_lower=true, lower_bound=0.0, report=dispatch_report)

    # Scarcity variables
    variable_ens_ndhty_real(JuMP_model, am, :scarcity, "ens_ndhty", risk_id, ids_n, ids_d, ids_h, ids_t, ids_y)                              # Energy Not Served (MW).
    
    # Storage
    variable_soc_idhty_real(JuMP_model, am, :storage, "soc_idhty", risk_id, ids_i_sto, ids_d, ids_h, ids_t, ids_y; report=true)            # Storage charge level at end of the period
    variable_chg_idhty_real(JuMP_model, am, :storage, "chg_idhty", risk_id, ids_i_sto, ids_d, ids_h, ids_t, ids_y; report=true)            # Storage charging MW
    variable_sto_c_idhty_integer(JuMP_model, am, :storage_commitment, "sto_c_idhty", risk_id, ids_i, ids_d, ids_h, ids_t, ids_y; report=dispatch_report)       
    
    # Large Flexible Load
    define_variable_idhty_real(JuMP_model, am, :demand, "lfl_ldhty", risk_id, ids_lfl, ids_d, ids_h, ids_t, ids_y; bounded_lower=true, lower_bound=0, report=dispatch_report)    # Large Flexible Load demand (MW)

    ###################################
    #------ DEFINE Constraints
    ###################################

    # System Balance & power flow
    if am.setting["Simulation Setting"]["power_flow_mode_flag"] == "PTDF"
        
        variable_p_inj_ndhty_real(JuMP_model, am, :powerflow, "p_inj_ndhty", risk_id, ids_n, ids_d, ids_h, ids_t, ids_y)                              # power injection at node n
        variable_demand_ndhty_real(JuMP_model, am, :powerflow, "demand_ndhty", risk_id, ids_n, ids_d, ids_h, ids_t, ids_y)                              # delivered demand

        for d in ids_d, h in ids_h, t in ids_t
            constraint_sum_p_injdhty_RA(JuMP_model, am, "constraint_sum_p_injdhty_RA", risk_id, d, h, t, y)

            for n in ids_n
                constraint_PTDF_Power_injection_ndhty_real(JuMP_model, am, "constraint_PTDF_Power_injection_ndhty_real", risk_id, n, d, h, t, y)
                constraint_PTDF_LoadBalance_ndhty_real_RA(JuMP_model, am, "constraint_PTDF_LoadBalance_ndhty_real_RA", risk_id, n, d, h, t, y; constant_load, constant_load_bus_idx, constant_load_distribution, constant_load_distribution_profile, system_peak_scale)
            end

            for k in ids_k
                constraint_dc_power_flow_kdhty_RA(JuMP_model, am, "constraint_dc_power_flow_kdhty_RA", risk_id, k, d, h, t, y; constant_load, system_peak_scale)
                constraint_dc_power_flow_max_kdhty_RA(JuMP_model, am, "constraint_dc_power_flow_max_kdhty", risk_id, k, d, h, t, y)
                constraint_dc_power_flow_min_kdhty_RA(JuMP_model, am, "constraint_dc_power_flow_min_kdhty", risk_id, k, d, h, t, y)
            end
        end


    elseif am.setting["Simulation Setting"]["power_flow_mode_flag"] == "Network_Flow"
        for d in ids_d, h in ids_h, t in ids_t
            for n in ids_n
                constraint_LoadBalance_ndhty_real_RA(JuMP_model, am, "constraint_LoadBalance_ndhty_real_RA", risk_id, n, d, h, t, y; constant_load, constant_load_bus_idx, constant_load_distribution, constant_load_distribution_profile, system_peak_scale)
            end
                    
            for k in ids_k
                constraint_power_flow_max_kdhty_OP(JuMP_model, am, "constraint_power_flow_max_kdhty_OP", risk_id, k, d, h, t, y)
                constraint_power_flow_min_kdhty_OP(JuMP_model, am, "constraint_power_flow_min_kdhty_OP", risk_id, k, d, h, t, y)
            end
        end

    elseif am.setting["Simulation Setting"]["power_flow_mode_flag"] == "B-theta"

        # One slack per AC island so DC ties don't tie separate synchronous areas together.
        ac_ref_bus_ids = get_ac_reference_buses(am; nw=0)
        variable_bus_angle_ndhty_real(JuMP_model, am, :powerflow, "bus_angle_ndhty", risk_id, ids_n, ids_d, ids_h, ids_t, ids_y; ref_bus_ids = ac_ref_bus_ids)                              # power injection at node n

        for d in ids_d, h in ids_h, t in ids_t
            for n in ids_n
                constraint_LoadBalance_ndhty_real_RA(JuMP_model, am, "constraint_LoadBalance_ndhty_real_RA", risk_id, n, d, h, t, y; constant_load, constant_load_bus_idx, constant_load_distribution, constant_load_distribution_profile, system_peak_scale)
            end

            for k in ids_k
                # DC ties carry no B-theta angle coupling; flow limits still apply to all branches.
                if get(am.ref[:nw][0][:branch][k], "dc_line", false) != true
                    constraint_b_theta_power_flow_kdhty_RA(JuMP_model, am, "constraint_b_theta_power_flow_kdhty_RA", risk_id, k, d, h, t, y)
                end
                constraint_dc_power_flow_max_kdhty_RA(JuMP_model, am, "constraint_dc_power_flow_max_kdhty_RA", risk_id, k, d, h, t, y)
                constraint_dc_power_flow_min_kdhty_RA(JuMP_model, am, "constraint_dc_power_flow_min_kdhty_RA", risk_id, k, d, h, t, y)
            end
        end

    end

    # VRE Budget
    for i in ids_i
        bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
        tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)
        if parameter(am, bus_idx, :gen_bus, tech_idx, "Hydro_Flag") in ("ROR", "IMPOUNDMENT")

            unit_group = parameter(am, 0, :gen_index, "UNIT_GROUP", i)
            fuel_limit = parameter(am, bus_idx, :gen_bus, tech_idx, "FUEL_LIMIT")

            if fuel_limit != "NA" && fuel_limit != "Fixed Profile"
                if am.setting["Simulation Configuration"]["Hydro_Budget_Flag"] == true
                    constraint_VREBalance_Budget_NonReserve_RA_snapshot_iy(JuMP_model, am, "constraint_VREBalance_Budget_NonReserve_RA_snapshot_iy", risk_id, i, day_group_id, y, ids_d, ids_h, ids_t, simulated_event_hours, risk_ED_solution, unit_group; const_name_flag)
                end
            end
        end
    end

    # Unit Dispatch : Fix the dispatch and related variables to 0 when unit is on outage based on the gen_risk_matrix for this scenario
    for i in ids_i
        if gen_risk_matrix[event_idx, i] == 0

            for (d, h, t) in ids_dht
                g_idhty = variable(am, nw, risk_id, :g_idhty, (i,d,h,t,y))
                JuMP.fix(g_idhty, 0.0; force = true)

                redispatch_abs_dev_idhty = variable(am, nw, risk_id, :redispatch_abs_dev_idhty, (i,d,h,t,y))
                JuMP.fix(redispatch_abs_dev_idhty, 0.0; force = true)

                if i in ids_i_sto
                    chg_idhty = variable(am, nw, risk_id, :chg_idhty, (i,d,h,t,y))
                    JuMP.fix(chg_idhty, 0.0; force = true)  # we don't allow charging during outage hours

                    # Pin SOC during outage to ref ED to prevent fabricated energy at return-from-outage.
                    soc_idhty = variable(am, nw, risk_id, :soc_idhty, (i,d,h,t,y))
                    idx_string = string("(", i, ", ", d, ", ", h, ", ", t, ", ", y, ")")
                    ref_soc = ref_ED_solution["storage"][idx_string]["soc_idhty"]
                    JuMP.fix(soc_idhty, ref_soc; force = true)
                end
            end
        end
    end

    # Unit Re-dispatch
    for (d, h, t) in ids_dht

        # System-wide cap: available units cannot collectively exceed reference net output
        # plus the outage shortfall; mirrors constraint_cont_deploy_dhty in the PF builder.
        constraint_total_cont_deploy_snapshot_dhty(JuMP_model, am, "constraint_total_cont_deploy_snapshot_dhty", risk_id, ids_n, ids_i, ids_i_sto, d, h, t, y, gen_risk_matrix, event_idx, ref_ED_solution)

        for i in ids_i

            risk_flag = gen_risk_matrix[event_idx, i]

            # SOC bounds apply for ALL storage regardless of outage (LP would otherwise fabricate energy).
            if i in ids_i_sto
                constraint_ES_SOC_Max_NonReserve_idhty_RA(JuMP_model, am, "constraint_ES_SOC_Max_NonReserve_idhty_RA", risk_id, i, d, h, t, y)
                constraint_ES_SOC_Min_NonReserve_idhty_RA(JuMP_model, am, "constraint_ES_SOC_Min_NonReserve_idhty_RA", risk_id, i, d, h, t, y)
            end

            if risk_flag == 1   # unit is available, we allow redispatch
                # apply redispatch power balance, ramping, and inter-temporal constraints for non-outaged units
                constraint_cont_deploy_snapshot_idhty(JuMP_model, am, "constraint_cont_deploy_snapshot_idhty", risk_id, i, d, h, t, y, ref_ED_solution, risk_ED_solution, event_first_hour_flag)
                constraint_redispatch_abs_dev_idhty(JuMP_model, am, "constraint_redispatch_abs_dev_idhty", risk_id, i, d, h, t, y, ref_ED_solution)

                # Storage Balance and Power Constraints
                if i in ids_i_sto

                    # Post-contingency ES charging limits; no ramp rates needed since ES
                    # responds fast and model resolution is hourly.
                    constraint_cont_deploy_ES_snapshot_idhty(JuMP_model, am, "constraint_cont_deploy_ES_snapshot_idhty", risk_id, i, d, h, t, y, ref_ED_solution, risk_ED_solution, RA_setting, event_first_hour_flag)

                    # ES discharge/charge pinned to reference when no ENS exists at this hour.
                    constraint_cont_deploy_ES_limit_snapshot_idhty(JuMP_model, am, "constraint_cont_deploy_ES_limit_snapshot_idhty", risk_id, i, d, h, t, y, risk_flag, ref_ED_solution)

                    # Storage commitment constraints
                    if i in ids_i_commit
                        constraint_ES_Charge_Max_Sto_UC_NonReserve_idhty(JuMP_model, am, "constraint_ES_Charge_Max_Sto_UC_NonReserve_idhty", risk_id, i, d, h, t, y)
                        constraint_ES_DisCharge_Max_Sto_UC_NonReserve_idhty(JuMP_model, am, "constraint_ES_DisCharge_Max_Sto_UC_NonReserve_idhty", risk_id, i, d, h, t, y)
                        constraint_ES_Sto_UC_Limit_idhty_RA(JuMP_model, am, "constraint_ES_Sto_UC_Limit_idhty_RA", risk_id, i, d, h, t, y)
                    end

                    # SOC balance: outage-event first hour starts from the reference prior-hour
                    # SOC; subsequent hours chain SOC within the RA model.
                    constraint_ES_SOC_Balance_Inter_Hour_NoReserve_Cont_Snapshot_idhty(JuMP_model, am, "constraint_ES_SOC_Balance_Inter_Hour_NoReserve_Cont_Snapshot_idhty", risk_id, i, d, h, t, y, ref_ED_solution, risk_ED_solution, daygroup_first_hour_flag, RA_setting)
                end
            end
        end
    end

    # LFL and Hybrid Constraints
    for (d, h, t) in ids_dht
        # We assume LFL can be curtailed during post-contingency hours.
        for l in ids_lfl
            constraint_LFL_Limit_ldhty(JuMP_model, am, "constraint_LFL_Limit_ldhty", risk_id, l, d, h, t, y)
        end

        # we still enforce inter-connection limit constraints for hybrid plants during post-contingency hours
        for i in ids_i_hybrid_gen
            constraint_hybrid_inter_connection_limit_injection_ldhty(JuMP_model, am, "constraint_hybrid_inter_connection_limit_injection_ldhty", risk_id, i, d, h, t, y; const_name_flag, reserve_flag=false)
            constraint_hybrid_inter_connection_limit_withdraw_ldhty(JuMP_model, am, "constraint_hybrid_inter_connection_limit_withdraw_ldhty", risk_id, i, d, h, t, y; const_name_flag, reserve_flag=false)
        end
    end

    ###################################
    #------ DEFINE Objective
    ###################################
    ED_objective_function_sequential_ED_Imperfect_Foresight_for_RA(JuMP_model, am, ids_i, ids_n, ids_dht, ids_y, ref_ED_solution, risk_id, event_idx, gen_risk_matrix)
    
    am.model[:nw][nw][risk_id] = JuMP_model

end


function ED_objective_function_ED_Perfect_Foresight_for_RA(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, ids_i, ids_n, ids_dht, ids_y, risk_id, gen_risk_matrix, redispatch_window_dht_set; nw::Int=am.cnw, report::Bool=true)

    objective = JuMP.AffExpr(0.0)
    y = ids_y[1]

    # Random adder
    Random.seed!(12345)
    epsilon = 1e-7

    ids_order = sort(collect(ids_i))
    jitter = Dict(id => rand() * epsilon for id in ids_order)

    # Redispatch-cost MC floor: deters arbitrary cycling of zero-MC units (storage, must-take VRE).
    # Workbook value is in pu cost; default 0.0001 ≈ $1/MWh real, above the jitter floor.
    ra_setting_lookup = get(am.setting, "RA Setting", Dict{String,Any}())
    min_redispatch_mc = Float64(get(ra_setting_lookup, "min_redispatch_mc_value", 0.0001))

    # Redispatch deviation cost applies at every continuous-window hour (pre-window dev is
    # fixed to 0, so the gate just keeps the LP expression compact).
    for i in ids_i

        bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
        tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)

        for (row_idx, (d,h,t)) in enumerate(ids_dht)

            if (d,h,t) in redispatch_window_dht_set

                MC_raw = parameter(am, bus_idx, :gen_bus, tech_idx, "Annual_MC")[string(y)][d]
                MC = max(MC_raw, min_redispatch_mc) + jitter[i]

                if haskey(am.ref[:nw][:0][:gen_index][i], "new_unit")   # ensure that the added unit is dispatched first
                    if am.ref[:nw][:0][:gen_index][i]["new_unit"] == true
                        MC = MC * 0.9
                    end
                end

                if (haskey(am.ref[:nw][bus_idx][:gen_bus][tech_idx], "asset_type")) && (parameter(am, bus_idx, :gen_bus, tech_idx, "asset_type") == "Transmission")   # this is a new battery added as SATA
                    JuMP.add_to_expression!(objective, - 1.0, variable(am, nw, risk_id, :redispatch_abs_dev_idhty, (i,d,h,t,y)))
                else
                    JuMP.add_to_expression!(objective, MC, variable(am, nw, risk_id, :redispatch_abs_dev_idhty, (i,d,h,t,y)))
                end
            end
        end
    end

    ids_n_order = sort(collect(ids_n))
    n_jitter = Dict(id => rand() * epsilon for id in ids_n_order)

    # Load shedding
    for (row_idx, (d,h,t)) in enumerate(ids_dht)
        for n in ids_n
            JuMP.add_to_expression!(objective, (am.setting["Simulation Configuration"]["VOLL"] + n_jitter[n]), variable(am, nw, risk_id, :ens_ndhty, (n,d,h,t,y)))
        end
    end

    # DEFINE OBJECTIVE FUNCTION
    return JuMP.@objective(JuMP_model, Min, objective     
    )


end


function ED_objective_function_sequential_ED_Imperfect_Foresight_for_RA(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, ids_i, ids_n, ids_dht, ids_y, ref_ED_solution, risk_id, event_idx, gen_risk_matrix; nw::Int=am.cnw, report::Bool=true)

    objective = JuMP.AffExpr(0.0)
    y = ids_y[1]

    # Random adder
    Random.seed!(12345)
    epsilon = 1e-7

    ids_order = sort(collect(ids_i))
    jitter = Dict(id => rand() * epsilon for id in ids_order)

    # Hourly Costs, divide by number of timesteps in an hour (card(T))
    card_T = length(am.setting["run_T"])

    # Redispatch-cost MC floor (mirror of perfect-foresight twin).
    ra_setting_lookup = get(am.setting, "RA Setting", Dict{String,Any}())
    min_redispatch_mc = Float64(get(ra_setting_lookup, "min_redispatch_mc_value", 0.0001))

    # we already have outages in these hours
    for i in ids_i
        bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
        tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)

        for (d, h, t) in ids_dht

            MC_raw = parameter(am, bus_idx, :gen_bus, tech_idx, "Annual_MC")[string(y)][d]
            MC = max(MC_raw, min_redispatch_mc) + jitter[i]

            if haskey(am.ref[:nw][:0][:gen_index][i], "new_unit")   # ensure that the added unit is dispatched first
                if am.ref[:nw][:0][:gen_index][i]["new_unit"] == true
                    MC = MC * 0.9
                end
            end

            if (haskey(am.ref[:nw][bus_idx][:gen_bus][tech_idx], "asset_type")) && (parameter(am, bus_idx, :gen_bus, tech_idx, "asset_type") == "Transmission")   # this is a new battery added as SATA
                JuMP.add_to_expression!(objective, - 1.0, variable(am, nw, risk_id, :redispatch_abs_dev_idhty, (i,d,h,t,y)))
            else
                JuMP.add_to_expression!(objective, MC, variable(am, nw, risk_id, :redispatch_abs_dev_idhty, (i,d,h,t,y)))
            end
        end

    end
    
    ids_n_order = sort(collect(ids_n))
    n_jitter = Dict(id => rand() * epsilon for id in ids_n_order)

    # Load shedding
    for (d, h, t) in ids_dht
        for n in ids_n
            JuMP.add_to_expression!(objective, (am.setting["Simulation Configuration"]["VOLL"] + n_jitter[n]), variable(am, nw, risk_id, :ens_ndhty, (n,d,h,t,y)))
        end
    end
    
    # DEFINE OBJECTIVE FUNCTION
    return JuMP.@objective(JuMP_model, Min, objective     
    )


end


function convert_LOLH_to_LOLE(lolh::Float64, num_hours::Int64, total_duration::Float64)
    # LOLE is reported as expected outage days per year.
    lole = lolh / 24.0
    return lole
end


function initialize_RA_dict(am, ids_n; day_group_flag = false, description = "")

    # index
    ids_ht = [(h,t) for h in am.setting["run_H"] for t in am.setting["run_T"]]

    bus_name_dict = Dict()
    for n in ids_n
        bus_name = am.ref[:nw][0][:bus][n]["bus_i"]
        bus_name = isa(bus_name, String) ? bus_name : string(bus_name)  # make sure that bus_name is a string

        bus_name_dict[n] = bus_name
    end

    # Nested by Annual/Day-Group x Systemwide/Regional.
    RA_metric_Info = Dict{String, Any}()

    RA_metric_Info["Annual"] = Dict{String, Any}()
    RA_metric_Info["Annual"]["Systemwide"] = 0.0
    RA_metric_Info["Annual"]["Regional"] = Dict{String, Any}()
    for n in ids_n
        RA_metric_Info["Annual"]["Regional"][bus_name_dict[n]] = 0.0
    end

    if day_group_flag == true
    
        RA_metric_Info["Day Group"] = Dict{String, Any}()
        for d in [(d) for (d) in get_index(am, :repday_groups, 0)]
            RA_metric_Info["Day Group"][string(d)] = Dict{String, Any}()
            RA_metric_Info["Day Group"][string(d)]["Systemwide"] = 0.0
            RA_metric_Info["Day Group"][string(d)]["Regional"] = Dict{String, Any}()
            RA_metric_Info["Day Group"][string(d)]["Regional_by_Hour"] = Dict{String, Any}()
            
            for n in ids_n

                bus_name = bus_name_dict[n]
                
                RA_metric_Info["Day Group"][string(d)]["Regional"][bus_name] = 0.0
                RA_metric_Info["Day Group"][string(d)]["Regional_by_Hour"][bus_name] = zeros(length(ids_ht))
            
            end
        end

    end

    RA_metric_Info["Description"] = description

    return RA_metric_Info
end


function initialize_RA_dict(am; day_group_flag = false, description = "")

    # Nested by Annual/Day-Group x Systemwide/Regional.
    RA_metric_Info = Dict{String, Any}()

    RA_metric_Info["Annual"] = Dict{String, Any}()
    RA_metric_Info["Annual"]["Systemwide"] = 0.0

    if day_group_flag == true
    
        RA_metric_Info["Day Group"] = Dict{String, Any}()
        for d in [(d) for (d) in get_index(am, :repday_groups, 0)]
            RA_metric_Info["Day Group"][string(d)] = Dict{String, Any}()
            RA_metric_Info["Day Group"][string(d)]["Systemwide"] = 0.0
        end
    end

    RA_metric_Info["Description"] = description

    return RA_metric_Info
end


function calculate_risk_metrics_ED!(am, RA_Info)

    if !haskey(RA_Info, "system_info")
        RA_Info["system_info"] = Dict{String, Any}()
    end

    pu_power_base = am.setting["Simulation Setting"]["per_unit_base_value"]
    
    num_risk_scenario = RA_Info["setting"]["num_risk_scenario"]
    renewable_scenario_weights = get(RA_Info["setting"], "renewable_scenario_weights", Dict{String, Float64}())
    filtered_joint_scenario_map = get(RA_Info["RA_reference_risk_data"], "filtered_joint_scenario_map", Dict{Int, Dict{Int, Dict{String, Any}}}())
    ids_n = [(k) for (k) in get_index(am, :bus, 0)]
    ids_ht = [(h,t) for h in am.setting["run_H"] for t in am.setting["run_T"]]
    ids_daygroup = [d for d in get_index(am, :repday_groups, 0)]
    current_year = RA_Info["setting"]["current_year"]
    
    # Define and initialize RA metrics (see description strings for units/definitions).

    description = "Expected Unserved Energy: Total expected energy (in MWh) that could not be supplied to demand, averaged over risk scenarios."
    RA_Info["RA_metrics"]["EUE"] = initialize_RA_dict(am, ids_n; day_group_flag = true, description)

    description = "Normalized Expected Unserved Energy: Expected Unserved Energy normalized by regional demand, expressed in parts per million (ppm)."
    RA_Info["RA_metrics"]["NEUE"] = initialize_RA_dict(am, ids_n; day_group_flag = true, description)

    description = "Loss of Load Hours: Total hours in a year where demand exceeds supply, averaged over risk scenarios."
    RA_Info["RA_metrics"]["LOLH"] = initialize_RA_dict(am, ids_n; description)

    description = "Loss of Load Expectation (days/year): Average expected number of days in a year when load exceeds generation capacity."
    RA_Info["RA_metrics"]["LOLE"] = initialize_RA_dict(am, ids_n; description)

    description = "Maximum number of consecutive hours with unserved energy in a given region or system."
    RA_Info["RA_metrics"]["Max_Consecutive_Outage_Hours"] = initialize_RA_dict(am, ids_n; description)   

    description = "Maximum amount of energy (in MWh) that could not be served during a specific period."
    RA_Info["RA_metrics"]["Max_MWh_Loss"] = initialize_RA_dict(am, ids_n; description)   

    description = "Maximum power loss (in MW) experienced during a specific period."
    RA_Info["RA_metrics"]["Max_MW_Loss"] = initialize_RA_dict(am, ids_n; description)   

    # Get Demand
    demand_list = Dict{Tuple{Int, Int, Int, Int}, Float64}()
    total_system_demand = 0.0
    regional_demand_list = Dict{Int, Float64}()
    daygroup_demand_list = Dict{Int, Float64}(d => 0.0 for d in ids_daygroup)
    system_peak_scale = RA_Info["setting"]["system_peak_scale"]
    
    for n in ids_n

        bus_data = am.ref[:nw][0][:bus][n]
        original_load = bus_data["aggregation_info"]["original_load_(bus_i, MW)"]
        aggregated_ba_ids = bus_data["aggregation_info"]["aggregated_regions_bus_i"]
        regional_demand = 0.0

        for day_group_id in ids_daygroup

            day_idx_list = am.ref[:nw][0][:repday_groups][day_group_id]["Day_Idx_List"]

            for d in day_idx_list
                num_days = parameter(am, 0, :repdays, "NumDays", d)
                for (h, t) in ids_ht
                    
                    # Keep demand units consistent with ENS (both scaled by pu_power_base),
                    # and annualize representative days by NumDays.
                    demand = get_bus_demand_with_growth(am, n, d, h, t, current_year; nw=0, system_peak_scale) * pu_power_base * num_days

                    demand_list[(n, d, h, t)] = demand
                    total_system_demand += demand
                    regional_demand += demand
                    daygroup_demand_list[day_group_id] += demand
                    
                end
            end 
        end
        
        regional_demand_list[n] = regional_demand
    end


    # calculate metrics
    regional_daygroup_demand = Dict{Tuple{Int, Int}, Float64}()
    for n in ids_n
        for day_group_id in ids_daygroup
            regional_daygroup_demand[(n, day_group_id)] = 0.0
            day_idx_list = am.ref[:nw][0][:repday_groups][day_group_id]["Day_Idx_List"]
            for d in day_idx_list
                for (h, t) in ids_ht
                    regional_daygroup_demand[(n, day_group_id)] += demand_list[(n, d, h, t)]
                end
            end
        end
    end
    
    # Initialize weighted expectation counters across joint scenarios.
    annual_unserved_energy_event_hours_count = 0.0
    annual_unserved_energy_MWh_count = 0.0

    # Track systemwide maxima and where they occur.
    system_max_consecutive_outage_hours = 0.0
    system_max_mwh_loss = 0.0
    system_max_mw_loss = 0.0

    # Initialize arrays to store the day_group_id for max losses
    day_group_id_max_consecutive_outages = 0
    day_group_id_max_mwh_loss = 0
    day_group_id_max_mw_loss = 0

    # Prevent double counting system LOLH across buses within each solved joint scenario.
    system_hour_counted = Dict{Tuple{Int, Int}, Set{Tuple{Int, Int, Int}}}()
        
    # bus_name dict
    bus_name_dict = Dict{Int, String}()
    for n in ids_n
        bus_name = am.ref[:nw][0][:bus][n]["bus_i"]
        bus_name = isa(bus_name, String) ? bus_name : string(bus_name)  # make sure that bus_name is a string

        bus_name_dict[n] = bus_name 
    end
    
    
    for n in ids_n

        regional_annual_unserved_energy_event_hours_count = 0.0
        regional_unserved_energy_MWh_count = 0.0

        # Track maximum consecutive outage hours, max MW loss, and max MWh loss
        regional_max_consecutive_outage_hours = 0
        regional_max_mwh_loss = 0.0  # Track maximum cumulative MWh loss
        regional_max_mw_loss = 0.0   # Track maximum instantaneous MW loss

        bus_name = bus_name_dict[n]
        
        for day_group_id in ids_daygroup
            day_group_key = string(day_group_id)
            if !haskey(RA_Info["RA_solutions"], day_group_key)
                continue
            end

            day_group_solution = RA_Info["RA_solutions"][day_group_key]
            day_group_eue = RA_Info["RA_metrics"]["EUE"]["Day Group"][day_group_key]
            day_group_neue = RA_Info["RA_metrics"]["NEUE"]["Day Group"][day_group_key]
            day_idx_list = am.ref[:nw][0][:repday_groups][day_group_id]["Day_Idx_List"]
            daygroup_expected_unserved_energy = 0.0

            for joint_id in keys(day_group_solution)
                if !haskey(filtered_joint_scenario_map, day_group_id) || !haskey(filtered_joint_scenario_map[day_group_id], joint_id)
                    @aleaf_warn "Skipped RA metric accumulation for unexpected joint_id=$joint_id in day_group=$day_group_id"
                    continue
                end

                joint_scenario = filtered_joint_scenario_map[day_group_id][joint_id]
                renewable_scenario_id = joint_scenario["renewable_scenario_id"]
                renewable_weight = get(renewable_scenario_weights, renewable_scenario_id, 0.0)
                joint_weight = renewable_weight / num_risk_scenario
                if joint_weight <= 0.0
                    continue
                end

                scenario_solution = day_group_solution[joint_id]
                if !(haskey(scenario_solution, "solution") && haskey(scenario_solution["solution"], "scarcity"))
                    continue
                end

                scarcity = scenario_solution["solution"]["scarcity"]

                # Temporary counters for consecutive outages in this (bus, day-group, scenario).
                consecutive_outage_hours = 0
                cumulative_mwh_loss = 0.0

                for d in day_idx_list
                    num_days = parameter(am, 0, :repdays, "NumDays", d)
                    check_count = 1

                    for (h, t) in ids_ht
                        string_key = string("(", n, ", ", d, ", ", h, ", ", t, ", ", current_year, ")")
                        unserved_energy = haskey(scarcity, string_key) ? scarcity[string_key]["ens_ndhty"] * pu_power_base : 0.0
                        
                        if unserved_energy > 0.0001 
                            weighted_unserved_energy = unserved_energy * num_days
                            # Accumulate annual and day-group expected unserved energy.
                            expected_unserved_energy = weighted_unserved_energy * joint_weight
                            annual_unserved_energy_MWh_count += expected_unserved_energy
                            regional_unserved_energy_MWh_count += expected_unserved_energy
                            daygroup_expected_unserved_energy += expected_unserved_energy
                            day_group_eue["Regional"][bus_name] += expected_unserved_energy
                            day_group_eue["Regional_by_Hour"][bus_name][check_count] += expected_unserved_energy

                            demand = demand_list[(n, d, h, t)]
                            if demand > 0
                                # Keep hourly profile as ratio at each hour; day-group NEUE is set from ratio-of-totals later.
                                expected_neue_hourly = (expected_unserved_energy / demand) * 1000000
                                day_group_neue["Regional_by_Hour"][bus_name][check_count] += expected_neue_hourly
                            end

                            # Track outage streak and severity.
                            consecutive_outage_hours += 1
                            cumulative_mwh_loss += unserved_energy
                            regional_max_mw_loss = max(regional_max_mw_loss, unserved_energy)

                            # Count system LOLH once per scenario/day/hour across all buses.
                            hour_key = (d, h, t)
                            scenario_key = (day_group_id, joint_id)
                            if !haskey(system_hour_counted, scenario_key)
                                system_hour_counted[scenario_key] = Set{Tuple{Int, Int, Int}}()
                            end
                            if !(hour_key in system_hour_counted[scenario_key])
                                push!(system_hour_counted[scenario_key], hour_key)
                                annual_unserved_energy_event_hours_count += num_days * joint_weight
                            end
                            regional_annual_unserved_energy_event_hours_count += num_days * joint_weight

                        else
                            # No outage in this hour: close and reset the current streak.
                            regional_max_consecutive_outage_hours = max(regional_max_consecutive_outage_hours, consecutive_outage_hours)
                            regional_max_mwh_loss = max(regional_max_mwh_loss, cumulative_mwh_loss)
                            consecutive_outage_hours = 0
                            cumulative_mwh_loss = 0.0
                        end

                        check_count += 1
                    end
                end

                # Final flush if the last simulated hour ended inside an outage streak.
                regional_max_consecutive_outage_hours = max(regional_max_consecutive_outage_hours, consecutive_outage_hours)
                regional_max_mwh_loss = max(regional_max_mwh_loss, cumulative_mwh_loss)

                # Update systemwide maxima and record where each max occurred.
                if regional_max_consecutive_outage_hours > system_max_consecutive_outage_hours
                    system_max_consecutive_outage_hours = regional_max_consecutive_outage_hours
                    day_group_id_max_consecutive_outages = day_group_id
                end
                if regional_max_mwh_loss > system_max_mwh_loss
                    system_max_mwh_loss = regional_max_mwh_loss
                    day_group_id_max_mwh_loss = day_group_id
                end
                if regional_max_mw_loss > system_max_mw_loss
                    system_max_mw_loss = regional_max_mw_loss
                    day_group_id_max_mw_loss = day_group_id
                end
            end

            # Day-group regional NEUE as ratio-of-totals for consistency with systemwide NEUE.
            daygroup_demand = regional_daygroup_demand[(n, day_group_id)]
            if daygroup_demand > 0
                day_group_neue["Regional"][bus_name] = (daygroup_expected_unserved_energy / daygroup_demand) * 1000000
            else
                day_group_neue["Regional"][bus_name] = 0.0
            end
        end

        # Save the regional LOLH, EUE, NEUE, Max MWh, and Max MW loss metrics
        RA_Info["RA_metrics"]["LOLH"]["Annual"]["Regional"][bus_name] = regional_annual_unserved_energy_event_hours_count
        RA_Info["RA_metrics"]["EUE"]["Annual"]["Regional"][bus_name] = regional_unserved_energy_MWh_count
        RA_Info["RA_metrics"]["Max_Consecutive_Outage_Hours"]["Annual"]["Regional"][bus_name] = regional_max_consecutive_outage_hours
        RA_Info["RA_metrics"]["Max_MWh_Loss"]["Annual"]["Regional"][bus_name] = regional_max_mwh_loss
        RA_Info["RA_metrics"]["Max_MW_Loss"]["Annual"]["Regional"][bus_name] = regional_max_mw_loss

        if regional_demand_list[n] > 0
            RA_Info["RA_metrics"]["NEUE"]["Annual"]["Regional"][bus_name] = (regional_unserved_energy_MWh_count / regional_demand_list[n]) * 1000000
        end

    end

    # Update system-wide reliability metrics

    # System wide LOLH
    RA_Info["RA_metrics"]["LOLH"]["Annual"]["Systemwide"] = annual_unserved_energy_event_hours_count

    # EUE
    RA_Info["RA_metrics"]["EUE"]["Annual"]["Systemwide"] = annual_unserved_energy_MWh_count

    for day_group_id in ids_daygroup
        for n in ids_n
            RA_Info["RA_metrics"]["EUE"]["Day Group"][string(day_group_id)]["Systemwide"] += RA_Info["RA_metrics"]["EUE"]["Day Group"][string(day_group_id)]["Regional"][bus_name_dict[n]]
        end
    end

    # NEUE
    if total_system_demand > 0
        RA_Info["RA_metrics"]["NEUE"]["Annual"]["Systemwide"] = (annual_unserved_energy_MWh_count / total_system_demand) * 1000000
    else
        RA_Info["RA_metrics"]["NEUE"]["Annual"]["Systemwide"] = 0.0
    end

    for day_group_id in ids_daygroup
        daygroup_demand = daygroup_demand_list[day_group_id]
        if daygroup_demand > 0
            RA_Info["RA_metrics"]["NEUE"]["Day Group"][string(day_group_id)]["Systemwide"] = (RA_Info["RA_metrics"]["EUE"]["Day Group"][string(day_group_id)]["Systemwide"] / daygroup_demand) * 1000000
        else
            RA_Info["RA_metrics"]["NEUE"]["Day Group"][string(day_group_id)]["Systemwide"] = 0.0
        end
    end

    # Calculate System wide Max Consecutive Outage Hours
    RA_Info["RA_metrics"]["Max_Consecutive_Outage_Hours"]["Annual"]["Systemwide"] = system_max_consecutive_outage_hours
    RA_Info["RA_metrics"]["Max_Consecutive_Outage_Hours"]["Annual"]["day_group_id"] = day_group_id_max_consecutive_outages

    # Calculate System wide Max MWh Loss
    RA_Info["RA_metrics"]["Max_MWh_Loss"]["Annual"]["Systemwide"] = system_max_mwh_loss
    RA_Info["RA_metrics"]["Max_MWh_Loss"]["Annual"]["day_group_id"] = day_group_id_max_mwh_loss

    # Calculate System wide Max MW Loss
    RA_Info["RA_metrics"]["Max_MW_Loss"]["Annual"]["Systemwide"] = system_max_mw_loss
    RA_Info["RA_metrics"]["Max_MW_Loss"]["Annual"]["day_group_id"] = day_group_id_max_mw_loss

    # Convert LOLH (hours/year) to LOLE (days/year)
    total_duration = 1.0 # years
    num_hours = 8760 # total number of hours in a year
    RA_Info["RA_metrics"]["LOLE"]["Annual"]["Systemwide"] = convert_LOLH_to_LOLE(RA_Info["RA_metrics"]["LOLH"]["Annual"]["Systemwide"], num_hours, total_duration)  

    for n in ids_n
        RA_Info["RA_metrics"]["LOLE"]["Annual"]["Regional"][bus_name_dict[n]] = convert_LOLH_to_LOLE(RA_Info["RA_metrics"]["LOLH"]["Annual"]["Regional"][bus_name_dict[n]], num_hours, total_duration)  
    end

    # Always export compact demand diagnostics for EUE/NEUE verification.
    regional_demand_by_bus = Dict{String, Float64}()
    regional_daygroup_demand_by_bus = Dict{String, Dict{String, Float64}}()
    for n in ids_n
        bus_name = bus_name_dict[n]
        regional_demand_by_bus[bus_name] = get(regional_demand_list, n, 0.0)
        by_day = Dict{String, Float64}()
        for day_group_id in ids_daygroup
            by_day[string(day_group_id)] = get(regional_daygroup_demand, (n, day_group_id), 0.0)
        end
        regional_daygroup_demand_by_bus[bus_name] = by_day
    end
    daygroup_demand = Dict{String, Float64}()
    solved_scenarios_by_daygroup = Dict{String, Int}()
    for day_group_id in ids_daygroup
        day_key = string(day_group_id)
        daygroup_demand[day_key] = get(daygroup_demand_list, day_group_id, 0.0)
        solved_scenarios_by_daygroup[day_key] = haskey(RA_Info["RA_solutions"], day_key) ? length(keys(RA_Info["RA_solutions"][day_key])) : 0
    end
    RA_Info["system_info"]["demand"] = Dict{String, Any}(
        "num_risk_scenario_setting" => num_risk_scenario,
        "total_system_demand_mwh" => total_system_demand,
        "regional_demand_mwh" => regional_demand_by_bus,
        "daygroup_demand_mwh" => daygroup_demand,
        "regional_daygroup_demand_mwh" => regional_daygroup_demand_by_bus,
        "solved_scenarios_by_daygroup" => solved_scenarios_by_daygroup
    )

    return RA_Info
end


"""
    build_ra_perfect_foresight_ed_model!(am, risk_id, day_group_id, y, gen_risk_matrix, ref_ED_solution, RA_setting; ...)

Build the perfect-foresight post-contingency Economic Dispatch JuMP model for a single
joint risk scenario over a day-group. The LP partitions day-group hours into:

  - **Pre-window** (hours strictly before the first contingency): g, chg, soc, lfl,
    redispatch_abs_dev are pinned to the reference dispatch (steady-state — the
    system hasn't been disturbed yet).
  - **Continuous redispatch window** (first contingency hour through end-of-day-group):
    joint LP redispatch with full physics. SOC chains continuously across event AND
    non-event hours (the unit's own outage hours hold SOC constant via the standard
    balance with g=chg=0).

This design (2026-06-12, commit f1b0958) replaces a prior per-event windowed
construction where each event's first hour read prior_soc from the reference dispatch
— letting the LP "free-ride" on reference SOC refills between events. With the
continuous-window design, the LP-tracked SOC is honest end-to-end through the day-group
and storage energy budget is binding in the natural way.

See also: `run_ra_imperfect_foresight_sequential_snapshot_case!` (the sequential
counterpart, applying the same continuous-window pattern hour-by-hour with virtualised
SOC bookkeeping at non-event hours).
"""
function build_ra_perfect_foresight_ed_model!(am::Abstract_ALEAF_Model, risk_id::Int, day_group_id::Int, y::Int, gen_risk_matrix, ref_ED_solution, RA_setting; nw::Int=am.cnw, constant_load=0.0, constant_load_bus_idx=nothing, constant_load_distribution="single_bus", constant_load_distribution_profile=Dict{Tuple{Int, Int, Int}, Dict{Int, Float64}}(), system_peak_scale=1.0, dispatch_report=false, horizon_dht_list::Union{Nothing, Vector}=nothing, prior_state_solution=nothing, horizon_daygroup_first_hour_flag::Bool=false)

    # Rolling-horizon mode: with `horizon_dht_list` the LP covers only that (d,h,t) subset,
    # skips pre-window pin-to-ref, and seeds the first hour from `prior_state_solution`.

    # Reset a JuMP model
    if (am.setting["Solver Setting"]["solver_name"] == "CPLEX") & (am.setting["Solver Setting"]["1"]["Value"] == true)

        JuMP_model = cplex_direct_model(Val(:CPLEX))
    else
        JuMP_model = JuMP.Model()
    end

    const_name_flag = am.setting["Simulation Setting"]["const_name_flag"]
    if const_name_flag == false
        JuMP.set_string_names_on_creation(JuMP_model, false)
    end

    horizon_mode = (horizon_dht_list !== nothing)

    # Reset registries when reusing a model_instance across LP builds (rolling-horizon
    # SeqED reuses it per event hour); add_sol_component asserts !haskey otherwise.
    if horizon_mode
        am.var = Dict{Symbol,Any}(:nw => Dict{Int,Any}())
        am.con = Dict{Symbol,Any}(:nw => Dict{Int,Any}())
        am.sol = Dict{Symbol,Any}(:nw => Dict{Int,Any}())
        for (nw_id, _) in am.ref[:nw]
            am.var[:nw][nw_id] = Dict{Int,Any}()
            am.con[:nw][nw_id] = Dict{Int,Any}()
            am.sol[:nw][nw_id] = Dict{Int,Any}()
            am.var[:nw][nw_id][risk_id] = Dict{Symbol,Any}()
            am.con[:nw][nw_id][risk_id] = Dict{Symbol,Any}()
            am.sol[:nw][nw_id][risk_id] = Dict{Symbol,Any}()
        end
    end

    #------ run for each day group
    ids_d_full = [(d) for (d) in am.ref[:nw][0][:repday_groups][day_group_id]["Day_Idx_List"]]
    ids_y = [(y)]

    ###################################
    #------ Pre-processing
    ###################################
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

    ids_h_full = [(h) for (h) in am.setting["run_H"]]
    ids_t_full = [(h) for (h) in am.setting["run_T"]]
    ids_k = [(k) for (k) in get_index(am, :branch, 0) if parameter(am, 0, :branch, "model_flag", k) == true]
    ids_n = [(k) for (k) in get_index(am, :bus, 0)]
    ids_z = [(z) for (z) in get_index(am, :zone, 0, "reserve")]
    ids_p = [(z) for (z) in get_index(am, :zone, 0, "policy")]
    ids_lfl = [(lfl) for (lfl) in get_index(am, :demand, 0)]    # Large Flexible Load

    if horizon_mode
        # Horizon-mode: derive per-axis ids_d/ids_h/ids_t (builders expect them) as the
        # union of axes touched; variables outside the horizon tuples are never referenced.
        ids_dht = collect(horizon_dht_list)
        ids_d = sort!(unique([dht[1] for dht in ids_dht]))
        ids_h = sort!(unique([dht[2] for dht in ids_dht]))
        ids_t = sort!(unique([dht[3] for dht in ids_dht]))
    else
        ids_d = ids_d_full
        ids_h = ids_h_full
        ids_t = ids_t_full
        ids_dht = [(d,h,t) for d in ids_d for h in ids_h for t in ids_t]
    end

    ###################################
    #------ Identify the continuous redispatch window
    ###################################
    # Continuous redispatch window: pre-window hours pin g/chg/soc/lfl to ref, window hours
    # (first contingency onward) redispatch with SOC chaining. Horizon mode: horizon = window.
    if horizon_mode
        h_first_event_idx = 1
        redispatch_window_dht_set = Set{Tuple{Int,Int,Int}}(ids_dht)
    else
        num_of_unit_outages_per_hour = vec(length(ids_i) .- sum(gen_risk_matrix, dims=2))
        indices_of_outage_hours = findall(x -> x > 0, num_of_unit_outages_per_hour)

        if !isempty(indices_of_outage_hours)
            h_first_event_idx = indices_of_outage_hours[1]
            redispatch_window_dht_set = Set{Tuple{Int,Int,Int}}(ids_dht[h_first_event_idx:end])
        else
            # No contingency in this day-group: entire day-group is pre-window (no LP
            # deviation from reference). Sentinel h_first_event_idx beyond ids_dht.
            h_first_event_idx = length(ids_dht) + 1
            redispatch_window_dht_set = Set{Tuple{Int,Int,Int}}()
        end
    end

    ###################################
    #------ DEFINE DECISION VARIABLES
    ###################################
    # Power flow variables
    if am.setting["Simulation Setting"]["power_flow_mode_flag"] in ["PTDF", "B-theta"]
        variable_f_kdhty_real(JuMP_model, am, :powerflow, "f_kdhty", risk_id, ids_k, ids_d, ids_h, ids_t, ids_y; bounded=false, report=dispatch_report)    # Power flow (unbounded for DC power flow)
    else # power_flow_mode_flag == "Network_Flow"
        variable_f_kdhty_real(JuMP_model, am, :powerflow, "f_kdhty", risk_id, ids_k, ids_d, ids_h, ids_t, ids_y; report=dispatch_report)                   # Power flow (bounded >=0 for network flow)
    end

    # Dispatch Variables
    variable_g_idhty_real(JuMP_model, am, :dispatch, "g_idhty", risk_id, ids_i, ids_d, ids_h, ids_t, ids_y; report=true)                  # Unit Generation (MW)
    define_variable_idhty_real(JuMP_model, am, :dispatch, "redispatch_abs_dev_idhty", risk_id, ids_i, ids_d, ids_h, ids_t, ids_y; bounded_lower=true, lower_bound=0.0, report=dispatch_report)

    # Scarcity variables
    variable_ens_ndhty_real(JuMP_model, am, :scarcity, "ens_ndhty", risk_id, ids_n, ids_d, ids_h, ids_t, ids_y)                              # Energy Not Served (MW).
    
    # Storage
    variable_soc_idhty_real(JuMP_model, am, :storage, "soc_idhty", risk_id, ids_i_sto, ids_d, ids_h, ids_t, ids_y; report=dispatch_report)            # Storage charge level at end of the period
    variable_chg_idhty_real(JuMP_model, am, :storage, "chg_idhty", risk_id, ids_i_sto, ids_d, ids_h, ids_t, ids_y; report=dispatch_report)            # Storage charging MW
    variable_sto_c_idhty_integer(JuMP_model, am, :storage_commitment, "sto_c_idhty", risk_id, ids_i, ids_d, ids_h, ids_t, ids_y; report=dispatch_report)               

    # Large Flexible Load
    define_variable_idhty_real(JuMP_model, am, :demand, "lfl_ldhty", risk_id, ids_lfl, ids_d, ids_h, ids_t, ids_y; bounded_lower=true, lower_bound=0, report=dispatch_report)    # Large Flexible Load demand (MW)

    ###################################
    #------ DEFINE Constraints
    ###################################

    # System Balance & power flow
    if am.setting["Simulation Setting"]["power_flow_mode_flag"] == "PTDF"
        
        variable_p_inj_ndhty_real(JuMP_model, am, :powerflow, "p_inj_ndhty", risk_id, ids_n, ids_d, ids_h, ids_t, ids_y)                              # power injection at node n
        variable_demand_ndhty_real(JuMP_model, am, :powerflow, "demand_ndhty", risk_id, ids_n, ids_d, ids_h, ids_t, ids_y)                              # delivered demand

        for (d, h, t) in ids_dht
            constraint_sum_p_injdhty_RA(JuMP_model, am, "constraint_sum_p_injdhty_RA", risk_id, d, h, t, y)

            for n in ids_n
                constraint_PTDF_Power_injection_ndhty_real(JuMP_model, am, "constraint_PTDF_Power_injection_ndhty_real", risk_id, n, d, h, t, y)
                constraint_PTDF_LoadBalance_ndhty_real_RA(JuMP_model, am, "constraint_PTDF_LoadBalance_ndhty_real_RA", risk_id, n, d, h, t, y; constant_load, constant_load_bus_idx, constant_load_distribution, constant_load_distribution_profile, system_peak_scale)
            end

            for k in ids_k
                constraint_dc_power_flow_kdhty_RA(JuMP_model, am, "constraint_dc_power_flow_kdhty_RA", risk_id, k, d, h, t, y; constant_load, system_peak_scale)
                constraint_dc_power_flow_max_kdhty_RA(JuMP_model, am, "constraint_dc_power_flow_max_kdhty", risk_id, k, d, h, t, y)
                constraint_dc_power_flow_min_kdhty_RA(JuMP_model, am, "constraint_dc_power_flow_min_kdhty", risk_id, k, d, h, t, y)
            end
        end


    elseif am.setting["Simulation Setting"]["power_flow_mode_flag"] == "Network_Flow"
        for (d, h, t) in ids_dht
            for n in ids_n
                constraint_LoadBalance_ndhty_real_RA(JuMP_model, am, "constraint_LoadBalance_ndhty_real_RA", risk_id, n, d, h, t, y; constant_load, constant_load_bus_idx, constant_load_distribution, constant_load_distribution_profile, system_peak_scale)
            end

            for k in ids_k
                constraint_power_flow_max_kdhty_OP(JuMP_model, am, "constraint_power_flow_max_kdhty_OP", risk_id, k, d, h, t, y)
                constraint_power_flow_min_kdhty_OP(JuMP_model, am, "constraint_power_flow_min_kdhty_OP", risk_id, k, d, h, t, y)
            end
        end

    elseif am.setting["Simulation Setting"]["power_flow_mode_flag"] == "B-theta"

        # One slack per AC island so DC ties don't tie separate synchronous areas together.
        ac_ref_bus_ids = get_ac_reference_buses(am; nw=0)
        variable_bus_angle_ndhty_real(JuMP_model, am, :powerflow, "bus_angle_ndhty", risk_id, ids_n, ids_d, ids_h, ids_t, ids_y; ref_bus_ids = ac_ref_bus_ids)                              # power injection at node n

        for (d, h, t) in ids_dht
            for n in ids_n
                constraint_LoadBalance_ndhty_real_RA(JuMP_model, am, "constraint_LoadBalance_ndhty_real_RA", risk_id, n, d, h, t, y; constant_load, constant_load_bus_idx, constant_load_distribution, constant_load_distribution_profile, system_peak_scale)
            end

            for k in ids_k
                # DC ties carry no B-theta angle coupling; flow limits still apply to all branches.
                if get(am.ref[:nw][0][:branch][k], "dc_line", false) != true
                    constraint_b_theta_power_flow_kdhty_RA(JuMP_model, am, "constraint_b_theta_power_flow_kdhty_RA", risk_id, k, d, h, t, y)
                end
                constraint_dc_power_flow_max_kdhty_RA(JuMP_model, am, "constraint_dc_power_flow_max_kdhty_RA", risk_id, k, d, h, t, y)
                constraint_dc_power_flow_min_kdhty_RA(JuMP_model, am, "constraint_dc_power_flow_min_kdhty_RA", risk_id, k, d, h, t, y)
            end
        end

    end

    # Hydro reservoir budget not enforced in horizon mode: an event spans a few hours so
    # daily reservoir balance is unchanged and per-hour max_g binds (perfect-foresight enforces it).
    if !horizon_mode
        for i in ids_i
            bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
            tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)
            if parameter(am, bus_idx, :gen_bus, tech_idx, "Hydro_Flag") in ("ROR", "IMPOUNDMENT")

                unit_group = parameter(am, 0, :gen_index, "UNIT_GROUP", i)
                fuel_limit = parameter(am, bus_idx, :gen_bus, tech_idx, "FUEL_LIMIT")

                if fuel_limit != "NA" && fuel_limit != "Fixed Profile"
                    if am.setting["Simulation Configuration"]["Hydro_Budget_Flag"] == true
                        constraint_VREBalance_Budget_NonReserve_iy_RA(JuMP_model, am, "constraint_VREBalance_Budget_NonReserve_iy_RA", risk_id, i, day_group_id, y, unit_group; const_name_flag)
                    end
                end
            end
        end
    end

    # Unit Dispatch : Fix the dispatch and related variables to 0 when unit is on outage based on the gen_risk_matrix for this scenario
    for i in ids_i
        risk_of_i = @view gen_risk_matrix[:, i]  # BitVector view for unit i
        for (row_idx, (d,h,t)) in enumerate(ids_dht)
            if !risk_of_i[row_idx]

                g_idhty = variable(am, nw, risk_id, :g_idhty, (i,d,h,t,y))
                JuMP.fix(g_idhty, 0.0; force = true)

                redispatch_abs_dev_idhty = variable(am, nw, risk_id, :redispatch_abs_dev_idhty, (i,d,h,t,y))
                JuMP.fix(redispatch_abs_dev_idhty, 0.0; force = true)

                if i in ids_i_sto
                    chg_idhty = variable(am, nw, risk_id, :chg_idhty, (i,d,h,t,y))
                    JuMP.fix(chg_idhty, 0.0; force = true)  # we don't allow charging during outage hours

                    # SOC is intentionally NOT pinned to ref during outage — with g=chg=0 fixed above
                    # and the hybrid onsite-charge injection gated to zero (es_avail), SOC holds constant.
                end
            end
        end
    end

    for row_idx in 1:(h_first_event_idx - 1)
        (d, h, t) = ids_dht[row_idx]
        for i in ids_i
            idx_string = string("(", i, ", ", d, ", ", h, ", ", t, ", ", y, ")")
            ref_g_idht = ref_ED_solution["dispatch"][idx_string]["g_idhty"]
            g_idhty = variable(am, nw, risk_id, :g_idhty, (i,d,h,t,y))
            JuMP.fix(g_idhty, ref_g_idht; force = true)

            redispatch_abs_dev_idhty = variable(am, nw, risk_id, :redispatch_abs_dev_idhty, (i,d,h,t,y))
            JuMP.fix(redispatch_abs_dev_idhty, 0.0; force = true)

            if i in ids_i_sto
                ref_chg_idht = ref_ED_solution["storage"][idx_string]["chg_idhty"]
                ref_soc_idht = ref_ED_solution["storage"][idx_string]["soc_idhty"]
                chg_idhty = variable(am, nw, risk_id, :chg_idhty, (i,d,h,t,y))
                soc_idhty = variable(am, nw, risk_id, :soc_idhty, (i,d,h,t,y))
                JuMP.fix(chg_idhty, ref_chg_idht; force = true)
                JuMP.fix(soc_idhty, ref_soc_idht; force = true)
            end
        end
        for lfl in ids_lfl
            idx_string = string("(", lfl, ", ", d, ", ", h, ", ", t, ", ", y, ")")
            ref_lfl_ldhty = ref_ED_solution["demand"][idx_string]["lfl_ldhty"]
            lfl_ldhty = variable(am, nw, risk_id, :lfl_ldhty, (lfl,d,h,t,y))
            JuMP.fix(lfl_ldhty, ref_lfl_ldhty; force = true)
        end
    end

    for window_offset in 0:(length(ids_dht) - h_first_event_idx)
        row_idx = h_first_event_idx + window_offset
        (d, h, t) = ids_dht[row_idx]

        # First redispatch-window hour: prior g and soc come from the reference dispatch at
        # the preceding hour (pre-window pinned to ref); later hours chain LP variables.
        window_first_hour_flag = (window_offset == 0)

        # First hour of the day-group: SOC balance first-hour helper reverse-calcs
        # the initial SOC state from ref at hour 1.
        # Horizon mode: the LP's first hour may be anywhere in the day-group; the
        # caller sets daygroup_first_hour_flag for the window's first hour, false after.
        daygroup_first_hour_flag = horizon_mode ? (window_first_hour_flag && horizon_daygroup_first_hour_flag) : (row_idx == 1)

        # System-wide redispatch limit at this hour
        constraint_cont_deploy_dhty(JuMP_model, am, "constraint_cont_deploy_dhty", risk_id, ids_n, ids_i, ids_i_sto, d, h, t, y, gen_risk_matrix, row_idx, ref_ED_solution)

        for i in ids_i

            risk_flag = gen_risk_matrix[row_idx, i]

            # SOC max/min bounds apply to ALL storage at every window hour, including
            # the storage unit's own outage hours (chained SOC stays constant when g=chg=0).
            if i in ids_i_sto
                constraint_ES_SOC_Max_NonReserve_idhty_RA(JuMP_model, am, "constraint_ES_SOC_Max_NonReserve_idhty_RA", risk_id, i, d, h, t, y)
                constraint_ES_SOC_Min_NonReserve_idhty_RA(JuMP_model, am, "constraint_ES_SOC_Min_NonReserve_idhty_RA", risk_id, i, d, h, t, y)

                # SOC balance at every window hour: only the first hour reads prior_soc from
                # ref, later hours chain from the LP at h-1 (horizon mode reads committed prior_soc).
                if window_first_hour_flag
                    constraint_post_cont_es_soc_balance_perfect_foresignt_first_hour_idhty(JuMP_model, am, "constraint_post_cont_es_soc_balance_perfect_foresignt_first_hour_idhty", risk_id, i, d, h, t, y, ref_ED_solution, daygroup_first_hour_flag, RA_setting; prior_state_solution=prior_state_solution, es_avail = Float64(risk_flag))
                else
                    constraint_post_cont_es_soc_balance_perfect_foresignt_idhty(JuMP_model, am, "constraint_post_cont_es_soc_balance_perfect_foresignt_idhty", risk_id, i, d, h, t, y, ref_ED_solution, daygroup_first_hour_flag; es_avail = Float64(risk_flag))
                end
            end

            # Dispatch + ramp + storage-charge constraints only for AVAILABLE units
            # at this hour (outaged units have g=chg=0 fixed in the outage block above).
            if risk_flag == 1
                # Ramp linkage: first window hour links to ref at h-1, later hours to
                # the LP at h-1; horizon mode anchors first hour on committed dispatch.
                constraint_cont_deploy_idhty(JuMP_model, am, "constraint_cont_deploy_idhty", risk_id, i, d, h, t, y, ref_ED_solution, window_first_hour_flag; prior_state_solution=prior_state_solution)
                constraint_redispatch_abs_dev_idhty(JuMP_model, am, "constraint_redispatch_abs_dev_idhty", risk_id, i, d, h, t, y, ref_ED_solution)

                if i in ids_i_sto
                    # Post-contingency ES charging limit (no ramp needed; storage responds fast).
                    constraint_cont_deploy_ES_idhty(JuMP_model, am, "constraint_cont_deploy_ES_idhty", risk_id, i, d, h, t, y, ref_ED_solution, window_first_hour_flag, RA_setting)

                    if i in ids_i_commit
                        constraint_ES_Charge_Max_Sto_UC_NonReserve_idhty(JuMP_model, am, "constraint_ES_Charge_Max_Sto_UC_NonReserve_idhty", risk_id, i, d, h, t, y)
                        constraint_ES_DisCharge_Max_Sto_UC_NonReserve_idhty(JuMP_model, am, "constraint_ES_DisCharge_Max_Sto_UC_NonReserve_idhty", risk_id, i, d, h, t, y)
                        constraint_ES_Sto_UC_Limit_idhty_RA(JuMP_model, am, "constraint_ES_Sto_UC_Limit_idhty_RA", risk_id, i, d, h, t, y)
                    end
                end
            end
        end
    end

    # LFL and hybrid constraints apply at every window hour (LFL curtailable, hybrids respect
    # inter-connection limits); pre-window LFL/hybrid dispatch is pinned to ref above.
    for row_idx in h_first_event_idx:length(ids_dht)
        (d, h, t) = ids_dht[row_idx]
        for l in ids_lfl
            constraint_LFL_Limit_ldhty(JuMP_model, am, "constraint_LFL_Limit_ldhty", risk_id, l, d, h, t, y)
        end
        for i in ids_i_hybrid_gen
            constraint_hybrid_inter_connection_limit_injection_ldhty(JuMP_model, am, "constraint_hybrid_inter_connection_limit_injection_ldhty", risk_id, i, d, h, t, y; const_name_flag, reserve_flag=false)
            constraint_hybrid_inter_connection_limit_withdraw_ldhty(JuMP_model, am, "constraint_hybrid_inter_connection_limit_withdraw_ldhty", risk_id, i, d, h, t, y; const_name_flag, reserve_flag=false)
        end
    end

    ###################################
    #------ DEFINE Objective
    ###################################
    # redispatch_window_dht_set: the LP may deviate from ref at every continuous-window
    # hour, so MC*|g - g_ref| covers them all; pre-window dev is pinned to 0.
    ED_objective_function_ED_Perfect_Foresight_for_RA(JuMP_model, am, ids_i, ids_n, ids_dht, ids_y, risk_id, gen_risk_matrix, redispatch_window_dht_set)


    am.model[:nw][nw][risk_id] = JuMP_model

end


function get_gen_outage_data(ALEAF_setting, case_id; RA_input::Dict{String,<:Any} = Dict{String,Any}(), current_year=1, recorded_investment_decisions::Dict{String,<:Any} = Dict{String,Any}(), external_GTEP_multi_round_info::Dict{String,<:Any} = Dict{String,Any}(), save_output::Bool=false)
    
    @aleaf_info "[ALEAF RA Model]: Generate outage sample data"

    # get recorded GTEP_multi_round_info if available
    output_path = define_output_path(ALEAF_setting, case_id)

    # Keep round/year selection consistent with execute_ALEAF_RA_model_using_external_data
    round_id = get(ALEAF_setting["RA Setting"], "round_id_to_run_RA", 1)
    if length(external_GTEP_multi_round_info) != 0

        @aleaf_info "[ALEAF RA Model]\tRunning RA using given generation mix"
        info = is_single_round_expansion_record(external_GTEP_multi_round_info) ?
            convert_single_round_record_to_multi_round_info(external_GTEP_multi_round_info, round_id) :
            external_GTEP_multi_round_info
        recorded_investment_decisions = info[string(round_id)]["updated_investment_decisions"]
        current_year = info[string(round_id)]["round_ids_y_decision"][1]

    else
        if ispath(joinpath(output_path, "GTEP_multi_round_info.json"))
            prior_GTEP_multi_round_info = JSON.parse(open(joinpath(output_path, "GTEP_multi_round_info.json")))
            recorded_investment_decisions = prior_GTEP_multi_round_info[string(round_id)]["updated_investment_decisions"]
            current_year = prior_GTEP_multi_round_info[string(round_id)]["round_ids_y_decision"][1]
            @aleaf_info "[ALEAF RA Model]\tRunning RA using saved generation mix"
        end
    end

    RA_ALEAF_setting = deepcopy(ALEAF_setting)

    update_solver_setting!(RA_ALEAF_setting)

    RA_Info = update_RA_setting(RA_ALEAF_setting, case_id, RA_input, current_year)

    network_data = generate_networkdata_LC_GTEP(RA_ALEAF_setting, case_id, "RA"; print_output_flag=true, year_id=current_year)
    @aleaf_info "[ALEAF RA Model]: Network Generation"

    network_data["output_path"] = define_output_path(ALEAF_setting, case_id)

    ALEAF_model_instance = ALEAF.build_ALEAF_model_instance_for_RA(RA_ALEAF_setting, case_id, RA_Info["setting"]["current_year"], network_data, recorded_investment_decisions);
    @aleaf_info "[ALEAF RA Model]: Build Reference System"

    if length(recorded_investment_decisions) > 0
        update_expansion_results!(ALEAF_model_instance, recorded_investment_decisions)
    end

    # Disaggregation happens here.
    prepare_RA_gen!(ALEAF_model_instance)
    @aleaf_info "[ALEAF RA Model]: Update Network Data for RA"

    risk_scenario_dict = Dict{Int, Any}()
    risk_scenario_dict, RA_Info["RA_reference_risk_data"]["filtered_joint_scenario_map"] = generate_and_filter_risk_scenarios(RA_ALEAF_setting, RA_Info, ALEAF_model_instance; RA_input)
    @aleaf_info "[ALEAF RA Model]: Risk Sampling Completed"

    # Save generated outage data so the same RA_input can be reused later.
    if save_output == true
        output_path = ALEAF_model_instance.ref[:nw][0][:output_path]
        outage_cache_path = joinpath(output_path, "RA_outage_cache")
        risk_scenario_dict_int64 = Dict{Int64, Any}(Int64(k) => v for (k, v) in risk_scenario_dict)
        _write_risk_scenario_cache!(outage_cache_path, risk_scenario_dict_int64; prefix="day")
        ra_input_output_path = joinpath(output_path, "RA_input_fixed_risk_sample.jld2")
        filtered_joint_scenario_map = RA_Info["RA_reference_risk_data"]["filtered_joint_scenario_map"]
        @save ra_input_output_path risk_scenario_dict filtered_joint_scenario_map
        @aleaf_info "[ALEAF RA Model]: Saved fixed-risk RA_input to $ra_input_output_path"
    end

    total_sce = 0
    total_num_days = 0
    for day_idx in keys(RA_Info["RA_reference_risk_data"]["filtered_joint_scenario_map"])
        num_of_risk = length(RA_Info["RA_reference_risk_data"]["filtered_joint_scenario_map"][day_idx])
        if num_of_risk > 0
            total_sce += length(RA_Info["RA_reference_risk_data"]["filtered_joint_scenario_map"][day_idx])
            total_num_days += 1
        end
    end
    total_num_risk = RA_ALEAF_setting["Simulation Configuration"][string(case_id)]["NDAY_Groups_RA"] * RA_Info["setting"]["num_risk_scenario"] * RA_Info["setting"]["num_renewable_scenarios"]
    identified_risk_ratio = total_num_risk > 0 ? (total_sce / total_num_risk) * 100 : 0.0
    @aleaf_info "[ALEAF RA Model]: Number of joint risk scenarios: $total_sce ($(round(identified_risk_ratio, digits=3)) %, in $total_num_days days)"

    return risk_scenario_dict, RA_Info["RA_reference_risk_data"]["filtered_joint_scenario_map"], network_data
end


function update_solver_setting!(RA_ALEAF_setting)

    if haskey(RA_ALEAF_setting, "CPLEX Setting")
        for setting_id in keys(RA_ALEAF_setting["CPLEX Setting"])
            if !(setting_id in ["optimizer", "solver_name", "1", "solution_type"])
                param = RA_ALEAF_setting["CPLEX Setting"][setting_id]["Parameter"]
                if param == "CPX_PARAM_THREADS"
                    RA_ALEAF_setting["CPLEX Setting"][setting_id]["Value"] = 1
                end
                if param == "CPX_PARAM_SCRIND"
                    RA_ALEAF_setting["CPLEX Setting"][setting_id]["Value"] = 0
                end
                if param == "CPXPARAM_SolutionType"
                    RA_ALEAF_setting["CPLEX Setting"][setting_id]["Flag"] = false
                end
            end
        end

    elseif haskey(RA_ALEAF_setting,"HiGHS Setting")
        for setting_id in keys(RA_ALEAF_setting["HiGHS Setting"])
            if !(setting_id in ["optimizer", "solver_name", "1", "solution_type"])
                param = RA_ALEAF_setting["HiGHS Setting"][setting_id]["Parameter"]
                # `threads` is deliberately not overridden: HiGHS fixes its global thread pool at the first solve in a
                # process (the expansion stage), and later solves asking for a different count fail with "error calling optimize!".
                if param == "output_flag"
                    RA_ALEAF_setting["HiGHS Setting"][setting_id]["Value"] = false
                end
            end
        end
    end

end
