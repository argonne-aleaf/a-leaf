# Copyright © 2026, UChicago Argonne, LLC. All Rights Reserved.
# SPDX-License-Identifier: BSD-3-Clause  (see LICENSE)

# Builds the network dict and time series for planning and RA workflows.


# Memoizes the year-invariant topology from get_network_data! so multi-round reruns skip
# the ~6s rebuild; deepcopy on every hit keeps each round isolated (sequential access).
const _TOPOLOGY_CACHE = Dict{String,Any}()

# The 6 topology keys produced by get_network_data! inside network_data.
const _TOPOLOGY_CACHE_NETWORK_KEYS = ["bus", "plant", "branch", "demand", "hybrid", "supplement_network_data"]

# ATB Technology names carrying a "Variable O&M" row; both the 2024 and 2025 spellings are listed
# because a Tech missing here reads as zero VOM rather than raising.
const ATB_TECHS_WITH_VOM = [
    "Biopower", "Coal", "CSP", "Natural Gas", "Nuclear", "Pumped Storage Hydropower",
    "Coal_FE", "NaturalGas_FE", "Geothermal", "Hydropower", "LandbasedWind",
    "OffShoreWind", "UtilityPV", "Utility-Scale Battery Storage",
]

# Cache key = config identity: network/config file paths + generation method/option.
# (Aggregation params live in the config file, so its path already disambiguates them.)
function _topology_cache_key(ALEAF_setting::Dict{String,<:Any})
    sim = get(ALEAF_setting, "Simulation Setting", Dict{String,Any}())
    return join(
        [
            string(get(ALEAF_setting, "network_data_file_location", "")),
            string(get(ALEAF_setting, "network_config_file_location", "")),
            string(get(sim, "power_flow_mode_flag", "")),   # PTDF mode attaches per-branch "ptdf" rows
        ],
        "||",
    )
end

# MISS runs get_network_data! and stores deepcopies; HIT restores fresh deepcopies and
# re-establishes the supplement alias + ALEAF_setting keys.
function get_network_data_cached!(network_data::Dict{String,<:Any}, ALEAF_setting::Dict{String,<:Any})
    # A selection-only build uses a coarsened aggregation resolution; never touch the shared topology
    # cache with it (would poison the full-run nodal topology). Build fresh and don't store.
    if get(network_data, "__selection_agg_override", nothing) !== nothing
        get_network_data!(network_data, ALEAF_setting)
        return nothing
    end
    key = _topology_cache_key(ALEAF_setting)
    cached = get(_TOPOLOGY_CACHE, key, nothing)
    if isnothing(cached)
        get_network_data!(network_data, ALEAF_setting)
        # Keep only the current config's topology (evict on config change) so the cache
        # can't accumulate across cases.
        empty!(_TOPOLOGY_CACHE)
        entry = Dict{String,Any}()
        for k in _TOPOLOGY_CACHE_NETWORK_KEYS
            entry[k] = deepcopy(network_data[k])
        end
        entry["network_resolution_level"] = deepcopy(ALEAF_setting["network_resolution_level"])
        entry["sub_area_mapping"] = deepcopy(ALEAF_setting["sub_area_mapping"])
        _TOPOLOGY_CACHE[key] = entry
    else
        for k in _TOPOLOGY_CACHE_NETWORK_KEYS
            network_data[k] = deepcopy(cached[k])
        end
        # additional_network_data is an alias (same object) of supplement_network_data.
        network_data["additional_network_data"] = network_data["supplement_network_data"]
        ALEAF_setting["network_resolution_level"] = deepcopy(cached["network_resolution_level"])
        ALEAF_setting["sub_area_mapping"] = deepcopy(cached["sub_area_mapping"])
    end
    return nothing
end


function get_timeseries_data!(ALEAF_setting::Dict{String,<:Any}, network_data::Dict{String,<:Any}, case_id)
    
    data_location = ALEAF_setting["data_location"]

    network_data["time_series_data"] = Dict{String, Any}()
    
    # Load
    file_path = joinpath(data_location, ALEAF_setting["File Path"]["timeseries_data_load_path"])
    if ALEAF_setting["Simulation Configuration"][string(case_id)]["Load_File_ID"] != "Base" 
        data_path = joinpath(data_location, "timeseries_data_files", "0_additional_scenarios", "Load")
        file_name = string("timeseries_load_hourly_", ALEAF_setting["Simulation Configuration"][string(case_id)]["Load_File_ID"], ".csv")
        file_path = joinpath(data_path, file_name)
    end
    timeSeries_df = CSV.read(file_path, DataFrame)
    network_data["time_series_data"]["load"] = timeSeries_df

    # Wind ons
    file_path = joinpath(data_location, ALEAF_setting["File Path"]["timeseries_data_wind_ons_path"])
    if ALEAF_setting["Simulation Configuration"][string(case_id)]["Wind_Ons_File_ID"] != "Base" 
        data_path = joinpath(data_location, "timeseries_data_files", "0_additional_scenarios", "WIND")
        file_name = string("timeseries_wind_ons_hourly_", ALEAF_setting["Simulation Configuration"][string(case_id)]["Wind_Ons_File_ID"], ".csv")
        file_path = joinpath(data_path, file_name)
    end
    timeSeries_df = CSV.read(file_path, DataFrame)
    network_data["time_series_data"]["wind_ons"] = timeSeries_df

    # Wind ofs
    timeSeries_df = CSV.read(joinpath(data_location, ALEAF_setting["File Path"]["timeseries_data_wind_ofs_path"]), DataFrame)
    network_data["time_series_data"]["wind_ofs"] = timeSeries_df

    # Fuel
    file_path = joinpath(data_location, ALEAF_setting["File Path"]["timeseries_data_fuel_price_path"])
    if ALEAF_setting["Simulation Configuration"][string(case_id)]["Fuel_ID"] != "Base" 
        data_path = joinpath(data_location, "timeseries_data_files", "0_additional_scenarios", "Fuel")
        file_name = string("timeseries_fuel_price_", ALEAF_setting["Simulation Configuration"][string(case_id)]["Fuel_ID"], ".csv")
        file_path = joinpath(data_path, file_name)
    end
    timeSeries_df = CSV.read(file_path, DataFrame)
    network_data["time_series_data"]["fuel_price"] = timeSeries_df

    # PV
    file_path = joinpath(data_location, ALEAF_setting["File Path"]["timeseries_data_pv_path"])
    if ALEAF_setting["Simulation Configuration"][string(case_id)]["PV_File_ID"] != "Base" 
        data_path = joinpath(data_location, "timeseries_data_files", "0_additional_scenarios", "PV")
        file_name = string("timeseries_pv_hourly_", ALEAF_setting["Simulation Configuration"][string(case_id)]["PV_File_ID"], ".csv")
        file_path = joinpath(data_path, file_name)
    end
    timeSeries_df = CSV.read(file_path, DataFrame)
    network_data["time_series_data"]["pv"] = timeSeries_df

    # RTPV
    timeSeries_df = CSV.read(joinpath(data_location, ALEAF_setting["File Path"]["timeseries_data_rtpv_path"]), DataFrame)
    network_data["time_series_data"]["rtpv"] = timeSeries_df

    # Hydro
    timeSeries_df = CSV.read(joinpath(data_location, ALEAF_setting["File Path"]["timeseries_data_hydro_path"]), DataFrame)
    network_data["time_series_data"]["hydro"] = timeSeries_df

    # CSP
    timeSeries_df = CSV.read(joinpath(data_location, ALEAF_setting["File Path"]["timeseries_data_csp_path"]), DataFrame)
    network_data["time_series_data"]["csp"] = timeSeries_df

    # DC (data-center large-flexible-load hourly profile) — optional; skip cleanly if unconfigured/absent
    dc_rel_path = get(ALEAF_setting["File Path"], "timeseries_data_dc_path", nothing)
    if dc_rel_path !== nothing
        dc_file_path = joinpath(data_location, dc_rel_path)
        if isfile(dc_file_path)
            network_data["time_series_data"]["dc"] = CSV.read(dc_file_path, DataFrame)
        end
    end

end


function load_timeseries_data!(data_location::String, ALEAF_setting::Dict{String,<:Any}, network_data::Dict{String,<:Any}, time_series_data_list, case_id)
    
    # Load
    if "load" in time_series_data_list
        timeSeries_df = network_data["time_series_data"]["load"]
        select_and_allocate_timeseries_data!(timeSeries_df, ALEAF_setting, network_data, "load")
    end

    # Wind
    if "wind" in time_series_data_list
        timeSeries_df = network_data["time_series_data"]["wind"]
        select_and_allocate_timeseries_data!(timeSeries_df, ALEAF_setting, network_data, "wind")
    end

    # Wind Ons
    if "wind_ons" in time_series_data_list
        timeSeries_df = network_data["time_series_data"]["wind_ons"]
        select_and_allocate_timeseries_data!(timeSeries_df, ALEAF_setting, network_data, "wind_ons")
    end

    # Wind Ofs
    if "wind_ofs" in time_series_data_list
        timeSeries_df = network_data["time_series_data"]["wind_ofs"]
        select_and_allocate_timeseries_data!(timeSeries_df, ALEAF_setting, network_data, "wind_ofs")
    end


    # Wind_BA
    if "wind_ons_BA" in time_series_data_list

        file_path = joinpath(data_location, ALEAF_setting["File Path"]["timeseries_data_wind_ons_path"])
        if ALEAF_setting["Simulation Configuration"][string(case_id)]["Wind_Ons_File_ID"] != "Base" 
            data_path = joinpath(data_location, "timeseries_data_files", "0_additional_scenarios", "WIND")
            file_name = string("timeseries_wind_ons_hourly_", ALEAF_setting["Simulation Configuration"][string(case_id)]["Wind_Ons_File_ID"], ".csv")
            file_path = joinpath(data_path, file_name)
        end

        timeSeries_df = CSV.read(file_path, DataFrame)
        select_and_allocate_timeseries_data!(timeSeries_df, ALEAF_setting, network_data, "wind_ons_BA")
    end

    # Wind_BA
    if "wind_ofs_BA" in time_series_data_list
        timeSeries_df = CSV.read(joinpath(data_location, ALEAF_setting["File Path"]["timeseries_data_wind_ofs_path"]), DataFrame)
        select_and_allocate_timeseries_data!(timeSeries_df, ALEAF_setting, network_data, "wind_ofs_BA")
    end

    # PV
    if "pv" in time_series_data_list
        timeSeries_df = network_data["time_series_data"]["pv"]
        select_and_allocate_timeseries_data!(timeSeries_df, ALEAF_setting, network_data, "pv")
    end

    # RTPV
    if "rtpv" in time_series_data_list
        timeSeries_df = network_data["time_series_data"]["rtpv"]
        select_and_allocate_timeseries_data!(timeSeries_df, ALEAF_setting, network_data, "rtpv")
    end

    # Hydro
    if "hydro" in time_series_data_list
        timeSeries_df = network_data["time_series_data"]["hydro"]
        select_and_allocate_timeseries_data!(timeSeries_df, ALEAF_setting, network_data, "hydro")
    end

    # CSP
    if "csp" in time_series_data_list
        timeSeries_df = network_data["time_series_data"]["csp"]
        select_and_allocate_timeseries_data!(timeSeries_df, ALEAF_setting, network_data, "csp")
    end

    # DC (data-center LFL hourly profile) — optional; only allocate if the file was loaded
    if "dc" in time_series_data_list && haskey(network_data["time_series_data"], "dc")
        timeSeries_df = network_data["time_series_data"]["dc"]
        select_and_allocate_timeseries_data!(timeSeries_df, ALEAF_setting, network_data, "dc")
    end

    # reg up
    if "reg_up" in time_series_data_list
        timeSeries_df = CSV.read(joinpath(data_location, ALEAF_setting["File Path"]["reserve_requirement_data_reg_up_path"]), DataFrame)
        select_and_allocate_timeseries_data!(timeSeries_df, ALEAF_setting, network_data, "reg_up")
    end

    # reg dn
    if "reg_down" in time_series_data_list
        timeSeries_df = CSV.read(joinpath(data_location, ALEAF_setting["File Path"]["reserve_requirement_data_reg_down_path"]), DataFrame)
        select_and_allocate_timeseries_data!(timeSeries_df, ALEAF_setting, network_data, "reg_down")
    end

    # reg
    if "reg" in time_series_data_list
        timeSeries_df = CSV.read(joinpath(data_location, ALEAF_setting["File Path"]["reserve_requirement_data_reg_path"]), DataFrame)
        select_and_allocate_timeseries_data!(timeSeries_df, ALEAF_setting, network_data, "reg")
    end

    # spin
    if "spin" in time_series_data_list
        timeSeries_df = CSV.read(joinpath(data_location, ALEAF_setting["File Path"]["reserve_requirement_data_spin_path"]), DataFrame)
        select_and_allocate_timeseries_data!(timeSeries_df, ALEAF_setting, network_data, "spin")
    end

    # non-spin
    if "nspin" in time_series_data_list
        timeSeries_df = CSV.read(joinpath(data_location, ALEAF_setting["File Path"]["reserve_requirement_data_nspin_path"]), DataFrame)
        select_and_allocate_timeseries_data!(timeSeries_df, ALEAF_setting, network_data, "nspin")
    end

    # flex up
    if "flex_up" in time_series_data_list
        timeSeries_df = CSV.read(joinpath(data_location, ALEAF_setting["File Path"]["reserve_requirement_data_flex_up_path"]), DataFrame)
        select_and_allocate_timeseries_data!(timeSeries_df, ALEAF_setting, network_data, "flex_up")
    end

    # flex dn
    if "flex_down" in time_series_data_list
        timeSeries_df = CSV.read(joinpath(data_location, ALEAF_setting["File Path"]["reserve_requirement_data_flex_down_path"]), DataFrame)
        select_and_allocate_timeseries_data!(timeSeries_df, ALEAF_setting, network_data, "flex_down")
    end


    
    ### 5 mins

    # Load 5mins
    if "load_5mins" in time_series_data_list
        timeSeries_df = CSV.read(joinpath(data_location, ALEAF_setting["File Path"]["timeseries_data_load_5mins_path"]), DataFrame)
        select_and_allocate_timeseries_data!(timeSeries_df, ALEAF_setting, network_data, "load")
    end

    # Wind
    if "wind_5mins" in time_series_data_list
        timeSeries_df = CSV.read(joinpath(data_location, ALEAF_setting["File Path"]["timeseries_data_wind_5mins_path"]), DataFrame)
        select_and_allocate_timeseries_data!(timeSeries_df, ALEAF_setting, network_data, "wind")
    end

    # PV
    if "pv_5mins" in time_series_data_list
        timeSeries_df = CSV.read(joinpath(data_location, ALEAF_setting["File Path"]["timeseries_data_pv_5mins_path"]), DataFrame)
        select_and_allocate_timeseries_data!(timeSeries_df, ALEAF_setting, network_data, "pv")
    end

    # reg
    if "reg_5mins" in time_series_data_list
        timeSeries_df = CSV.read(joinpath(data_location, ALEAF_setting["File Path"]["reserve_requirement_data_reg_5mins_path"]), DataFrame)
        select_and_allocate_timeseries_data!(timeSeries_df, ALEAF_setting, network_data, "reg")
    end

    # spin
    if "spin_5mins" in time_series_data_list
        timeSeries_df = CSV.read(joinpath(data_location, ALEAF_setting["File Path"]["reserve_requirement_data_spin_5mins_path"]), DataFrame)
        select_and_allocate_timeseries_data!(timeSeries_df, ALEAF_setting, network_data, "spin")
    end

    # non-spin
    if "nspin_5mins" in time_series_data_list
        timeSeries_df = CSV.read(joinpath(data_location, ALEAF_setting["File Path"]["reserve_requirement_data_nspin_5mins_path"]), DataFrame)
        select_and_allocate_timeseries_data!(timeSeries_df, ALEAF_setting, network_data, "nspin")
    end


    ### Reg up/down signal

    # flex dn
    if "reg_up_signal" in time_series_data_list
        # check if the file exists
        if haskey(ALEAF_setting["File Path"], "reserve_requirement_data_reg_up_signal_path")
            timeSeries_df = CSV.read(joinpath(data_location, ALEAF_setting["File Path"]["reserve_requirement_data_reg_up_signal_path"]), DataFrame)
            select_and_allocate_timeseries_data!(timeSeries_df, ALEAF_setting, network_data, "reg_up_signal")
        end
    end

    if "reg_down_signal" in time_series_data_list
        # check if the file exists
        if haskey(ALEAF_setting["File Path"], "reserve_requirement_data_reg_down_signal_path")
            timeSeries_df = CSV.read(joinpath(data_location, ALEAF_setting["File Path"]["reserve_requirement_data_reg_down_signal_path"]), DataFrame)
            select_and_allocate_timeseries_data!(timeSeries_df, ALEAF_setting, network_data, "reg_down_signal")
        end
    end

    ### Generation outages
    if "gen_outages" in time_series_data_list
        timeSeries_df = CSV.read(joinpath(data_location, ALEAF_setting["File Path"]["timeseries_data_gen_outages_path"]), DataFrame)
        select_and_allocate_timeseries_data!(timeSeries_df, ALEAF_setting, network_data, "gen_outages")
    end

    if "ngcc_outages" in time_series_data_list
        timeSeries_df = CSV.read(joinpath(data_location, ALEAF_setting["File Path"]["timeseries_data_gen_outages_path"]), DataFrame)
        select_and_allocate_timeseries_data!(timeSeries_df, ALEAF_setting, network_data, "NGCC")
    end

    if "ngct_outages" in time_series_data_list
        timeSeries_df = CSV.read(joinpath(data_location, ALEAF_setting["File Path"]["timeseries_data_gen_outages_path"]), DataFrame)
        select_and_allocate_timeseries_data!(timeSeries_df, ALEAF_setting, network_data, "NGCT")
    end

    if "nuclear_outages" in time_series_data_list
        timeSeries_df = CSV.read(joinpath(data_location, ALEAF_setting["File Path"]["timeseries_data_gen_outages_path"]), DataFrame)
        select_and_allocate_timeseries_data!(timeSeries_df, ALEAF_setting, network_data, "Nuclear")
    end


end



# Number of days in a time-series table, from its row count (365 or 366 days depending on the data year).
function infer_timeseries_num_days(timeSeries, steps_per_day::Int, time_series_type::AbstractString)
    n_rows = nrow(timeSeries)
    n_days, remainder = divrem(n_rows, steps_per_day)
    if remainder != 0 || !(n_days in (365, 366))
        error("[ALEAF]: time series '$time_series_type' has $n_rows rows; expected 365 or 366 days x $steps_per_day steps per day " *
              "(num_hours_per_day_value x num_sub_period_value).")
    end
    return n_days
end

function select_and_allocate_timeseries_data!(timeSeries, ALEAF_setting::Dict{String,<:Any}, network_data::Dict{String,<:Any}, time_series_type::String)
    # Add day id into the Dataframe (day starts at 1)

    day_id = []
    num_hours_per_day = ALEAF_setting["Planning Design"]["num_hours_per_day_value"]
    num_sub_period = ALEAF_setting["Planning Design"]["num_sub_period_value"]

    num_days = infer_timeseries_num_days(timeSeries, num_hours_per_day * num_sub_period, time_series_type)

    day = 1
    for d = 1:num_days
        for i = 1:(num_hours_per_day*num_sub_period)
            append!(day_id, day)
        end
        day += 1
    end
    insertcols!(timeSeries, 1, :day_id=>day_id)

    # Select and allocate timeseries data
    for stage in keys(network_data["planning_stages"])
        for id in keys(network_data["planning_stages"][stage]["repdays"])
            
            # Select timeseries data
            sub_dataframe = timeSeries[timeSeries.day_id.==network_data["planning_stages"][stage]["repdays"][id]["Day"],:]
            sub_dic = convert_dataFrame_to_dict_string_any(sub_dataframe)
    
            for sub_dic_id in keys(sub_dic)
                if haskey(sub_dic[sub_dic_id], "day_id") delete!(sub_dic[sub_dic_id], "day_id") end
                if haskey(sub_dic[sub_dic_id], "Year") delete!(sub_dic[sub_dic_id], "Year") end
                if haskey(sub_dic[sub_dic_id], "Month") delete!(sub_dic[sub_dic_id], "Month") end
                if haskey(sub_dic[sub_dic_id], "Day") delete!(sub_dic[sub_dic_id], "Day") end
                if haskey(sub_dic[sub_dic_id], "Period") delete!(sub_dic[sub_dic_id], "Period") end
                if haskey(sub_dic[sub_dic_id], "Period_Beginning") delete!(sub_dic[sub_dic_id], "Period_Beginning") end
    
                if haskey(sub_dic[sub_dic_id], "SolarShape") delete!(sub_dic[sub_dic_id], "SolarShape") end
                if haskey(sub_dic[sub_dic_id], "LoadShape") delete!(sub_dic[sub_dic_id], "LoadShape") end
                if haskey(sub_dic[sub_dic_id], "Load") delete!(sub_dic[sub_dic_id], "Load") end
                if haskey(sub_dic[sub_dic_id], "HydroShape") delete!(sub_dic[sub_dic_id], "HydroShape") end
                if haskey(sub_dic[sub_dic_id], "rtpvShape") delete!(sub_dic[sub_dic_id], "rtpvShape") end
                if haskey(sub_dic[sub_dic_id], "WindShape") delete!(sub_dic[sub_dic_id], "WindShape") end
            end
    
            # allocate/update timeseries data (hour/sub-hour)
            id_dic = 1
            if !haskey(network_data["planning_stages"][stage]["repdays"][id], "data") network_data["planning_stages"][stage]["repdays"][id]["data"] = Dict{String, Any}() end
            for h = 1:num_hours_per_day
                if !haskey(network_data["planning_stages"][stage]["repdays"][id]["data"], string(h)) network_data["planning_stages"][stage]["repdays"][id]["data"][string(h)] = Dict{String, Any}() end
                for t = 1:num_sub_period
                    if !haskey(network_data["planning_stages"][stage]["repdays"][id]["data"][string(h)] , string(t)) network_data["planning_stages"][stage]["repdays"][id]["data"][string(h)][string(t)] = Dict{String, Any}() end
                    network_data["planning_stages"][stage]["repdays"][id]["data"][string(h)][string(t)][time_series_type] = sub_dic[string(id_dic)]
                    id_dic += 1
                end
            end
        end
    end
    
end


function select_and_update_timeseries_data!(timeSeries, ALEAF_setting::Dict{String,<:Any}, network_data::Dict{String,<:Any}, time_series_type::String)
    # Add day id into the Dataframe (day starts at 1)

    day_id = []
    num_hours_per_day = ALEAF_setting["Planning Design"]["num_hours_per_day_value"]
    num_sub_period = ALEAF_setting["Planning Design"]["num_sub_period_value"]

    num_days = infer_timeseries_num_days(timeSeries, num_hours_per_day * num_sub_period, time_series_type)

    day = 1
    for d = 1:num_days
        for i = 1:(num_hours_per_day*num_sub_period)
            append!(day_id, day)
        end
        day += 1
    end
    insertcols!(timeSeries, 1, :day_id=>day_id)

    # Select and allocate timeseries data
    for stage in keys(network_data["planning_stages"])
        # Use the first present repday key, not a hardcoded "1": a distributed OP worker materializes
        # only its own day-group's repday, so "1" may be absent.
        if !isempty(network_data["planning_stages"][stage]["repdays"]) && haskey(network_data["planning_stages"][stage]["repdays"][first(keys(network_data["planning_stages"][stage]["repdays"]))], "data")
            for id in keys(network_data["planning_stages"][stage]["repdays"])
                
                # Select timeseries data
                sub_dataframe = timeSeries[timeSeries.day_id.==network_data["planning_stages"][stage]["repdays"][id]["Day"],:]
                sub_dic = convert_dataFrame_to_dict_string_any(sub_dataframe)

                for sub_dic_id in keys(sub_dic)
                    if haskey(sub_dic[sub_dic_id], "day_id") delete!(sub_dic[sub_dic_id], "day_id") end
                    if haskey(sub_dic[sub_dic_id], "Year") delete!(sub_dic[sub_dic_id], "Year") end
                    if haskey(sub_dic[sub_dic_id], "Month") delete!(sub_dic[sub_dic_id], "Month") end
                    if haskey(sub_dic[sub_dic_id], "Day") delete!(sub_dic[sub_dic_id], "Day") end
                    if haskey(sub_dic[sub_dic_id], "Period") delete!(sub_dic[sub_dic_id], "Period") end
                    if haskey(sub_dic[sub_dic_id], "Period_Beginning") delete!(sub_dic[sub_dic_id], "Period_Beginning") end

                    if haskey(sub_dic[sub_dic_id], "SolarShape") delete!(sub_dic[sub_dic_id], "SolarShape") end
                    if haskey(sub_dic[sub_dic_id], "LoadShape") delete!(sub_dic[sub_dic_id], "LoadShape") end
                    if haskey(sub_dic[sub_dic_id], "Load") delete!(sub_dic[sub_dic_id], "Load") end
                    if haskey(sub_dic[sub_dic_id], "HydroShape") delete!(sub_dic[sub_dic_id], "HydroShape") end
                    if haskey(sub_dic[sub_dic_id], "rtpvShape") delete!(sub_dic[sub_dic_id], "rtpvShape") end
                    if haskey(sub_dic[sub_dic_id], "WindShape") delete!(sub_dic[sub_dic_id], "WindShape") end
                end

                # allocate/update timeseries data (hour/sub-hour)
                id_dic = 1
                for h = 1:num_hours_per_day
                    for t = 1:num_sub_period
                        network_data["planning_stages"][stage]["repdays"][id]["data"][string(h)][string(t)][time_series_type] = sub_dic[string(id_dic)]
                        id_dic += 1
                    end
                end
            end
        end
    end
    
end


function select_and_update_stochastic_timeseries_data!(timeSeries, ALEAF_setting::Dict{String,<:Any}, network_data::Dict{String,<:Any}, time_series_type::String, scenario_id)
    # Add day id into the Dataframe (day starts at 1)

    day_id = []
    num_hours_per_day = ALEAF_setting["Planning Design"]["num_hours_per_day_value"]
    num_sub_period = ALEAF_setting["Planning Design"]["num_sub_period_value"]

    num_days = infer_timeseries_num_days(timeSeries, num_hours_per_day * num_sub_period, time_series_type)

    day = 1
    for d = 1:num_days
        for i = 1:(num_hours_per_day*num_sub_period)
            append!(day_id, day)
        end
        day += 1
    end
    insertcols!(timeSeries, 1, :day_id=>day_id)

    # Select and allocate timeseries data
    for stage in keys(network_data["planning_stages"])
        # Use the first present repday key, not a hardcoded "1": a distributed OP worker materializes
        # only its own day-group's repday, so "1" may be absent.
        if !isempty(network_data["planning_stages"][stage]["repdays"]) && haskey(network_data["planning_stages"][stage]["repdays"][first(keys(network_data["planning_stages"][stage]["repdays"]))], "data")
            for id in keys(network_data["planning_stages"][stage]["repdays"])

                if network_data["planning_stages"][stage]["repdays"][id]["Scenario_ID"] == scenario_id
                
                    # Select timeseries data
                    sub_dataframe = timeSeries[timeSeries.day_id.==network_data["planning_stages"][stage]["repdays"][id]["Day"],:]
                    sub_dic = convert_dataFrame_to_dict_string_any(sub_dataframe)

                    for sub_dic_id in keys(sub_dic)
                        if haskey(sub_dic[sub_dic_id], "day_id") delete!(sub_dic[sub_dic_id], "day_id") end
                        if haskey(sub_dic[sub_dic_id], "Year") delete!(sub_dic[sub_dic_id], "Year") end
                        if haskey(sub_dic[sub_dic_id], "Month") delete!(sub_dic[sub_dic_id], "Month") end
                        if haskey(sub_dic[sub_dic_id], "Day") delete!(sub_dic[sub_dic_id], "Day") end
                        if haskey(sub_dic[sub_dic_id], "Period") delete!(sub_dic[sub_dic_id], "Period") end
                        if haskey(sub_dic[sub_dic_id], "Period_Beginning") delete!(sub_dic[sub_dic_id], "Period_Beginning") end

                        if haskey(sub_dic[sub_dic_id], "SolarShape") delete!(sub_dic[sub_dic_id], "SolarShape") end
                        if haskey(sub_dic[sub_dic_id], "LoadShape") delete!(sub_dic[sub_dic_id], "LoadShape") end
                        if haskey(sub_dic[sub_dic_id], "Load") delete!(sub_dic[sub_dic_id], "Load") end
                        if haskey(sub_dic[sub_dic_id], "HydroShape") delete!(sub_dic[sub_dic_id], "HydroShape") end
                        if haskey(sub_dic[sub_dic_id], "rtpvShape") delete!(sub_dic[sub_dic_id], "rtpvShape") end
                        if haskey(sub_dic[sub_dic_id], "WindShape") delete!(sub_dic[sub_dic_id], "WindShape") end
                    end

                    # allocate/update timeseries data (hour/sub-hour)
                    id_dic = 1
                    for h = 1:num_hours_per_day
                        for t = 1:num_sub_period
                            network_data["planning_stages"][stage]["repdays"][id]["data"][string(h)][string(t)][time_series_type] = sub_dic[string(id_dic)]
                            id_dic += 1
                        end
                    end

                end
            end
        end
    end
    
end


function update_storage_technology_data!(ALEAF_setting::Dict{String,<:Any}, network_data::Dict{String,<:Any}, case_id)

    ESGC_ID = ALEAF_setting["Simulation Configuration"][string(case_id)]["ESGC_Setting_ID"]
    
    delete!(ALEAF_setting["Storage Cost and Performance"], "1") # delete: description
    delete!(ALEAF_setting["Storage Cost and Performance"], "2") # delete: units
    
    function find_ESGC_data(unit_group_id::String, data_label::String)
        for idx in keys(ALEAF_setting["Storage Cost and Performance"])
            if ALEAF_setting["Storage Cost and Performance"][idx]["ESGC_Setting_ID"] == ESGC_ID
                if ALEAF_setting["Storage Cost and Performance"][idx]["UNITGROUP"] == unit_group_id
                    return ALEAF_setting["Storage Cost and Performance"][idx][data_label]
                end
            end
        end
    end
    
    for id in keys(network_data["gen_technology"])
        if network_data["gen_technology"][id]["UNIT_CATEGORY"] == "STORAGE"

            if network_data["gen_technology"][id]["BATEFF"] == "ESGC"
                network_data["gen_technology"][id]["BATEFF_reference"] = "ESGC"
                network_data["gen_technology"][id]["BATEFF"] = find_ESGC_data(network_data["gen_technology"][id]["UNITGROUP"], "RTE")
            end

            if network_data["gen_technology"][id]["AET"] == "ESGC"
                network_data["gen_technology"][id]["AET_reference"] = "ESGC"
                network_data["gen_technology"][id]["AET"] = find_ESGC_data(network_data["gen_technology"][id]["UNITGROUP"], "AET")
            end
        end
    end
end


function update_branch_annual_cost_data!(ALEAF_setting::Dict{String,<:Any}, network_data)

    # Per-case so transmission cost can vary across cases; runs after the topology cache.
    case_config = ALEAF_setting["Simulation Configuration"][network_data["case_id"]]
    ac_present_value = case_config["transmission_cost_dollar_per_MW_mile_value"]
    dc_present_value = case_config["dc_tie_expansion_cost_dollar_per_MW_value"]
    WACC = ALEAF_setting["Planning Design"]["WACC_value"]
    crpyears = ALEAF_setting["Planning Design"]["transmission_investment_CRP_value"]
    # Bus centroids give straight-line distance, but real routes detour around terrain and land use.
    route_adder = get(ALEAF_setting["Planning Design"], "transmission_route_length_adder_value", 1.0)
    # Annual fixed O&M as a fraction of overnight capital (percent in the workbook, x0.01 here).
    fom_fraction = get(ALEAF_setting["Planning Design"], "transmission_FOM_percent_value", 0.0) * 0.01

    capital_recovery_factor = (WACC * ((1+WACC)^crpyears)) / ((1+WACC)^crpyears - 1)

    for k in keys(network_data["branch"])
        if get(network_data["branch"][k], "dc_line", false) == true
            # DC tie: 2 converter stations (one per side) + AC approach-line cost over branch length.
            branch_length = get(network_data["branch"][k], "length", get(network_data["branch"][k], "Length", 0.0))
            present_value = 2 * dc_present_value + ac_present_value * branch_length * route_adder
        else
            # AC length enters at objective time via len_factor, so the route adder folds in here.
            present_value = ac_present_value * route_adder
        end
        network_data["branch"][k]["transmission_expansion_cost"] = present_value * capital_recovery_factor
        # FOM shares the same overnight basis as capex (same len_factor at use), charged per operating year.
        network_data["branch"][k]["transmission_fom_cost"] = present_value * fom_fraction
    end
end


function update_hybrid_plant_technology_data!(network_data::Dict{String,<:Any})

    function find_hybrid_component_tech_idx(data::Dict{String,<:Any}, hybrid_id::String, LFL_id::String)
        
        hybrid_tech_idx = ""
        for idx in keys(data)
            if (data[idx]["UNITGROUP"] == hybrid_id) 
                hybrid_tech_idx = idx
            end
        end

        if (hybrid_tech_idx == "") 
            terminate_with_error(; msg="Failed to find hybrid component data for LFL id $LFL_id")
        else
            return hybrid_tech_idx
        end
    end
        
    for id in keys(network_data["hybrid"])
        if network_data["hybrid"][id]["Hybrid_Gen"] != "NA"
            # get components idx for generator
            hybrid_tech_idx = find_hybrid_component_tech_idx(network_data["gen_technology"], network_data["hybrid"][id]["Hybrid_Gen"], id)
            network_data["hybrid"][id]["Component_gen_tech_idx"] = hybrid_tech_idx
            network_data["hybrid"][id]["Component_gen_tech_data"] = deepcopy(network_data["gen_technology"][hybrid_tech_idx])
            network_data["hybrid"][id]["Component_gen_tech_data"]["CAP"] =  deepcopy(network_data["hybrid"][id]["Hybrid_Gen_CAP"])
        end 
        
        if network_data["hybrid"][id]["Hybrid_ES"] != "NA"
            # get components idx for storage
            hybrid_tech_idx = find_hybrid_component_tech_idx(network_data["gen_technology"], network_data["hybrid"][id]["Hybrid_ES"], id)
            network_data["hybrid"][id]["Component_ES_tech_idx"] = hybrid_tech_idx
            network_data["hybrid"][id]["Component_ES_tech_data"] = deepcopy(network_data["gen_technology"][hybrid_tech_idx])
            network_data["hybrid"][id]["Component_ES_tech_data"]["CAP"] =  deepcopy(network_data["hybrid"][id]["Hybrid_ES_CAP"])
            network_data["hybrid"][id]["Component_ES_tech_data"]["Charge_CAP"] =  deepcopy(network_data["hybrid"][id]["Hybrid_ES_CAP"])
        end
    end

    # Add hybrid plant to plant dict (separately)
    plant_id = length(network_data["plant"]) + 1
    plant_id_str = string(plant_id)

    # hybrid plant common info -> stay in hybrid dict 
    # add VRE and ES to plant for dispatch 
    common_info_list = ["PLANT_NAME", "bus_name", "bus_ID", "RetireYear", "Online_Year", "PLANT_ORIS_ID"]
    for id in keys(network_data["hybrid"])

        if network_data["hybrid"][id]["Hybrid_Gen"] != "NA"

            # record plant id
            network_data["hybrid"][id]["Hybrid_Gen_Plant_ID"] = plant_id_str

            # add Hybrid Gen to Plant database 
            network_data["plant"][plant_id_str] = deepcopy(network_data["hybrid"][id]["Component_gen_tech_data"])
            network_data["plant"][plant_id_str]["original_bus_i"] = network_data["hybrid"][id]["original_bus_i"]
            network_data["plant"][plant_id_str]["bus_i"] = network_data["hybrid"][id]["bus_i"]
            network_data["plant"][plant_id_str]["bus_idx"] = network_data["hybrid"][id]["bus_idx"]
            network_data["plant"][plant_id_str]["region_mapping_info"] = network_data["hybrid"][id]["region_mapping_info"]
            network_data["plant"][plant_id_str]["ES_MWh"] = 0.0
            network_data["plant"][plant_id_str]["Tech_Type"] = "Existing"
            network_data["plant"][plant_id_str]["bypass_parameter_check"] = true # does not perform check_and_update_plant_data!
            network_data["plant"][plant_id_str]["hybrid_type"] = "GEN"
            network_data["plant"][plant_id_str]["hybrid_ID"] = id

            # add common info
            for info in common_info_list
                network_data["plant"][plant_id_str][info] = network_data["hybrid"][id][info]
            end

            plant_id += 1
            plant_id_str = string(plant_id)
        end

        if network_data["hybrid"][id]["Hybrid_ES"] != "NA"

            # record plant id
            network_data["hybrid"][id]["Hybrid_ES_Plant_ID"] = plant_id_str

            # add Hybrid Gen to Plant database 
            network_data["plant"][plant_id_str] = deepcopy(network_data["hybrid"][id]["Component_ES_tech_data"])
            network_data["plant"][plant_id_str]["original_bus_i"] = network_data["hybrid"][id]["original_bus_i"]
            network_data["plant"][plant_id_str]["bus_i"] = network_data["hybrid"][id]["bus_i"]
            network_data["plant"][plant_id_str]["bus_idx"] = network_data["hybrid"][id]["bus_idx"]
            network_data["plant"][plant_id_str]["region_mapping_info"] = network_data["hybrid"][id]["region_mapping_info"]
            network_data["plant"][plant_id_str]["ES_MWh"] = network_data["hybrid"][id]["Component_ES_tech_data"]["Charge_CAP"] * network_data["hybrid"][id]["Component_ES_tech_data"]["STOHR_MAX"]
            network_data["plant"][plant_id_str]["Tech_Type"] = "Existing"
            network_data["plant"][plant_id_str]["bypass_parameter_check"] = true # does not perform check_and_update_plant_data!
            network_data["plant"][plant_id_str]["hybrid_type"] = "ES"
            network_data["plant"][plant_id_str]["hybrid_ID"] = id

            # add common info
            for info in common_info_list
                network_data["plant"][plant_id_str][info] = network_data["hybrid"][id][info]
            end

            plant_id += 1
            plant_id_str = string(plant_id)

        end
    end
end


function add_hybrid_plant_to_LFL!(network_data::Dict{String,<:Any})

    function find_hybrid_component_tech_idx(data::Dict{String,<:Any}, hybrid_id::String, LFL_id::String)
        
        hybrid_tech_idx = ""
        for idx in keys(data)
            if (data[idx]["UNITGROUP"] == hybrid_id) 
                hybrid_tech_idx = idx
            end
        end

        if (hybrid_tech_idx == "") 
            terminate_with_error(; msg="Failed to find hybrid component data for LFL id $LFL_id")
        else
            return hybrid_tech_idx
        end
    end
        
    for id in keys(network_data["demand"])
        if network_data["demand"][id]["Hybrid_Gen"] != "NA"
            # get components idx for generator
            hybrid_tech_idx = find_hybrid_component_tech_idx(network_data["gen_technology"], network_data["demand"][id]["Hybrid_Gen"], id)
            network_data["demand"][id]["Component_gen_tech_idx"] = hybrid_tech_idx
            network_data["demand"][id]["Component_gen_tech_data"] = deepcopy(network_data["gen_technology"][hybrid_tech_idx])
            network_data["demand"][id]["Component_gen_tech_data"]["CAP"] =  deepcopy(network_data["demand"][id]["Hybrid_Gen_CAP"])
        end 
        
        if network_data["demand"][id]["Hybrid_ES"] != "NA"
            # get components idx for storage
            hybrid_tech_idx = find_hybrid_component_tech_idx(network_data["gen_technology"], network_data["demand"][id]["Hybrid_ES"], id)
            network_data["demand"][id]["Component_ES_tech_idx"] = hybrid_tech_idx
            network_data["demand"][id]["Component_ES_tech_data"] = deepcopy(network_data["gen_technology"][hybrid_tech_idx])
            network_data["demand"][id]["Component_ES_tech_data"]["CAP"] =  deepcopy(network_data["demand"][id]["Hybrid_ES_CAP"])
            network_data["demand"][id]["Component_ES_tech_data"]["Charge_CAP"] =  deepcopy(network_data["demand"][id]["Hybrid_ES_CAP"])
        end
    end
end


function update_gen_technology_data!(data_location::String, ALEAF_setting::Dict{String,<:Any}, network_data::Dict{String,<:Any}, scenario::Int)
    
    # Read ATB data
    common_data_location = joinpath(pwd(), "data", "common")
    atb_year = ALEAF_setting["Simulation Configuration"][string(scenario)]["ATB_Year"]
    atb_file_name = string("ATB_", atb_year, ".csv")
    ATB_data_raw = ALEAF_setting["ATB_data_raw"] = DataFrame(CSV.File(joinpath(common_data_location, atb_file_name)))
    detailed_log = haskey(ALEAF_setting, "Simulation Setting") &&
                   (lowercase(string(get(ALEAF_setting["Simulation Setting"], "logging_level_value", "simple"))) == "detailed")
    
    # Read ESGC data
    ESGC_ID = ALEAF_setting["Simulation Configuration"][string(scenario)]["ESGC_Setting_ID"]
    delete!(ALEAF_setting["Storage Cost and Performance"], "1") # delete: description
    delete!(ALEAF_setting["Storage Cost and Performance"], "2") # delete: units

    
    function find_ESGC_data(unit_group_id::String, data_label::String)
        storage_data = ALEAF_setting["Storage Cost and Performance"] 
        for idx in keys(storage_data)
            entry = storage_data[idx]  
            if entry["ESGC_Setting_ID"] == ESGC_ID && entry["UNITGROUP"] == unit_group_id
                value = get(entry, data_label, nothing)
                if isnothing(value)
                    @aleaf_warn "Missing ESGC field '$data_label' for UNITGROUP $unit_group_id with ESGC Setting ID $ESGC_ID"
                    terminate_with_error(; msg="Failed to find ESGC field '$data_label' for UNITGROUP $unit_group_id with ESGC Setting ID $ESGC_ID. Please check your ESGC settings.")
                end
                return value
            end
        end

        @aleaf_warn "Failed to find ESGC data for UNITGROUP $unit_group_id with ESGC Setting ID $ESGC_ID"
        terminate_with_error(; msg="Failed to find ESGC data for UNITGROUP $unit_group_id with ESGC Setting ID $ESGC_ID. Please check your ESGC settings.")
    end

    function find_policy_data(unit_group_id::String, data_label::String)
        data = ALEAF_setting[data_label]
        for idx in keys(data)
            if data[idx]["UNITGROUP"] == unit_group_id 
                return data[idx]
            end
        end

        @aleaf_warn "Failed to find policy data '$data_label' for UNITGROUP $unit_group_id in scenario $scenario"
        terminate_with_error(; msg="Failed to find policy data '$data_label' for UNITGROUP $unit_group_id in scenario $scenario. Please check the simulation setting policy input data.")
    end

    # 2) update Gen Technology data using ATB setting and ATB data
    for id in keys(network_data["gen_technology"])

        atb_setting = get_atb_setting(ALEAF_setting, ALEAF_setting["Simulation Configuration"][string(scenario)]["ATB_Setting_ID"], network_data["gen_technology"][id]["UNITGROUP"])
                
        # terminate only when this technology requires ATB data
        tech_data = network_data["gen_technology"][id]
        requires_atb = (tech_data["CRP"] == "ATB") ||
                       (tech_data["CAPEX"] == "ATB") ||
                       (tech_data["FCR"] == "ATB") ||
                       (tech_data["FOM"] == "ATB") ||
                       (tech_data["VOM"] == "ATB") ||
                       (tech_data["FC"] == "ATB") ||
                       (tech_data["STO_CAPEX"] == "ATB")

        if isempty(atb_setting) && requires_atb
            @aleaf_warn "No matching ATB setting found for technology $(tech_data["UNITGROUP"]) in scenario $scenario, but this technology requires ATB data."
            terminate_with_error(; msg="Failed to find ATB setting for technology $(tech_data["UNITGROUP"]) in scenario $scenario. Please check your ATB settings.")
        end

        if requires_atb
            required_atb_keys = ["Case", "CRP", "Tech", "TechDetail", "Scenario", "ATB_Year", "CAPEX_Scale"]
            missing_atb_keys = [k for k in required_atb_keys if !haskey(atb_setting, k)]
            if !isempty(missing_atb_keys)
                missing_keys_str = join(missing_atb_keys, ", ")
                @aleaf_warn "ATB setting is missing required keys for technology $(tech_data["UNITGROUP"]) in scenario $scenario: $missing_keys_str"
                terminate_with_error(; msg="ATB setting is missing required keys for technology $(tech_data["UNITGROUP"]) in scenario $scenario: $missing_keys_str. Please check your ATB settings.")
            end
        end
        
        # Add CRP years
        if network_data["gen_technology"][id]["CRP"] == "ATB"
            network_data["gen_technology"][id]["crpyears"] = atb_setting["CRP"]
        else 
            network_data["gen_technology"][id]["crpyears"] = network_data["gen_technology"][id]["CRP"]
        end

        # investment cost
        network_data["gen_technology"][id]["INVC"] = Dict{String,Any}()
        network_data["gen_technology"][id]["STO_INVC"] = Dict{String,Any}()
        network_data["gen_technology"][id]["Annual_FCR"] = Dict{String,Any}()
        network_data["gen_technology"][id]["Annual_CRF"] = Dict{String,Any}()
        network_data["gen_technology"][id]["Annual_INVC"] = Dict{String,Any}()
        network_data["gen_technology"][id]["Annual_STO_INV"] = Dict{String,Any}()
        network_data["gen_technology"][id]["Annual_FOM"] = Dict{String,Any}()
        network_data["gen_technology"][id]["Annual_VOM"] = Dict{String,Any}()
        network_data["gen_technology"][id]["Annual_FC"] = Dict{String,Any}()

        network_data["gen_technology"][id]["ITC"] = find_policy_data(network_data["gen_technology"][id]["UNITGROUP"], "ITC")
        network_data["gen_technology"][id]["PTC"] = find_policy_data(network_data["gen_technology"][id]["UNITGROUP"], "PTC")
        warned_missing_atb_vom = false
        
        for y in keys(network_data["planning_stages"])

            # update year value (y -> year)
            atb_setting["Year"] = network_data["planning_stages"][y]["year"]
            if atb_setting["Year"] > 2050
                atb_setting["Year"] = 2050
            end

            # update ESGC year  (if year > 2050, ESGC year = 2050)
            ESGC_year = network_data["planning_stages"][y]["year"]
            if ESGC_year > 2050
                ESGC_year = 2050
            end

            # update ITC/PTC year 
            itc_ptc_year = network_data["planning_stages"][y]["year"]
            if itc_ptc_year > 2050 # (if year > 2050, ITC/PTC year = 2050)
                itc_ptc_year = 2050
            elseif itc_ptc_year < 2025 # (if year < 2025, ITC/PTC year = 2025)
                itc_ptc_year = 2025
            end

            # CAPEX
            capex = 0.0
            if network_data["gen_technology"][id]["CAPEX"] == "ATB"
                capex = get_atb_value(ATB_data_raw, atb_setting, "CAPEX") * atb_setting["CAPEX_Scale"]
            elseif network_data["gen_technology"][id]["CAPEX"] == "ESGC"
                label = string("Total Capital Cost_kW_", ESGC_year)
                capex = find_ESGC_data(network_data["gen_technology"][id]["UNITGROUP"], label) * find_ESGC_data(network_data["gen_technology"][id]["UNITGROUP"], "CAPEX_Scale")
            else
                capex = network_data["gen_technology"][id]["CAPEX"] 
            end

            # Apply ITC to CAPEX
            if ALEAF_setting["Simulation Configuration"][string(scenario)]["ITC_Flag"] == true
                if network_data["gen_technology"][id]["ITC Flag"] == true
                    if !haskey(network_data["gen_technology"][id]["ITC"], string(itc_ptc_year))
                        @aleaf_warn "Missing ITC value for UNITGROUP $(network_data["gen_technology"][id]["UNITGROUP"]) in scenario $scenario at year $(itc_ptc_year)"
                        terminate_with_error(; msg="Failed to find ITC value for UNITGROUP $(network_data["gen_technology"][id]["UNITGROUP"]) in scenario $scenario at year $(itc_ptc_year). Please check ITC policy data.")
                    end
                    itc = network_data["gen_technology"][id]["ITC"][string(itc_ptc_year)] # percent
                    capex = capex * (1 - itc)
                end
            end

            # FCR
            if network_data["gen_technology"][id]["FCR"] == "ATB"
                fcr = get_atb_value_FCR(ATB_data_raw, atb_setting, "FCR")
                network_data["gen_technology"][id]["Annual_FCR"][y] = fcr
            else
                fcr = network_data["gen_technology"][id]["FCR"]
                network_data["gen_technology"][id]["Annual_FCR"][y] = fcr
            end

            # CRF
            if network_data["gen_technology"][id]["FCR"] == "ATB"
                life = network_data["gen_technology"][id]["Life"]
                network_data["gen_technology"][id]["Annual_CRF"][y] = get_atb_value_CRF(ATB_data_raw, atb_setting, life)
            else
                network_data["gen_technology"][id]["Annual_CRF"][y] = 0.0
            end

            if fcr == "NA"
                network_data["gen_technology"][id]["INVC"][y] = capex 
                network_data["gen_technology"][id]["Annual_INVC"][y] = get_annual_investment_cost(capex, network_data, ALEAF_setting, id, ATB_data_raw, atb_setting) 
            else
                network_data["gen_technology"][id]["INVC"][y] = capex 
                network_data["gen_technology"][id]["Annual_INVC"][y] = capex * fcr 
            end

            # storage duration investment cost
            if network_data["gen_technology"][id]["UNIT_CATEGORY"] == "STORAGE"
                if network_data["gen_technology"][id]["STO_CAPEX"] == "ESGC"
                    label = string("Total Capital Cost_kWh_", ESGC_year)
                    storage_CAPEX = find_ESGC_data(network_data["gen_technology"][id]["UNITGROUP"], label) * find_ESGC_data(network_data["gen_technology"][id]["UNITGROUP"], "CAPEX_Scale")
                elseif network_data["gen_technology"][id]["STO_CAPEX"] == "ATB"
                    storage_CAPEX = 0.0 # ATB does not report $/kWh for storage duration investment cost
                else
                    storage_CAPEX = network_data["gen_technology"][id]["STO_CAPEX"] * network_data["gen_technology"][id]["CAPEX_Scale"]
                end

                network_data["gen_technology"][id]["STO_INVC"][y] = 0.0 

                # Apply ITC to CAPEX
                if ALEAF_setting["Simulation Configuration"][string(scenario)]["ITC_Flag"] == true
                    if network_data["gen_technology"][id]["ITC Flag"] == true
                        if !haskey(network_data["gen_technology"][id]["ITC"], string(itc_ptc_year))
                            @aleaf_warn "Missing ITC value for UNITGROUP $(network_data["gen_technology"][id]["UNITGROUP"]) in scenario $scenario at year $(itc_ptc_year)"
                            terminate_with_error(; msg="Failed to find ITC value for UNITGROUP $(network_data["gen_technology"][id]["UNITGROUP"]) in scenario $scenario at year $(itc_ptc_year). Please check ITC policy data.")
                        end
                        itc =network_data["gen_technology"][id]["ITC"][string(itc_ptc_year)]
                        storage_CAPEX = storage_CAPEX * (1 - itc)
                    end
                end

                if (fcr == "TECH_LIFE") || (fcr == "NA" )
                    network_data["gen_technology"][id]["STO_INVC"][y] = storage_CAPEX 
                    network_data["gen_technology"][id]["Annual_STO_INV"][y] = get_annual_investment_cost(storage_CAPEX, network_data, ALEAF_setting, id, ATB_data_raw, atb_setting) 
                else
                    network_data["gen_technology"][id]["Annual_STO_INV"][y] = storage_CAPEX * fcr 
                end
            else
                network_data["gen_technology"][id]["Annual_STO_INV"][y] = 0.0
            end

            # FOM & VOM
            if network_data["gen_technology"][id]["FOM"] == "ATB"
                network_data["gen_technology"][id]["Annual_FOM"][y] = get_atb_value(ATB_data_raw, atb_setting, "Fixed O&M") 
            elseif network_data["gen_technology"][id]["FOM"] == "ESGC"
                network_data["gen_technology"][id]["Annual_FOM"][y] = find_ESGC_data(network_data["gen_technology"][id]["UNITGROUP"], "FOM")
            else
                network_data["gen_technology"][id]["Annual_FOM"][y] = network_data["gen_technology"][id]["FOM"]
            end
    
            if network_data["gen_technology"][id]["VOM"] == "ATB"
                # ATB 2025 renamed several technologies; a name absent here silently zeroes VOM instead of reading it
                if atb_setting["Tech"] in ATB_TECHS_WITH_VOM
                    network_data["gen_technology"][id]["Annual_VOM"][y] = get_atb_value(ATB_data_raw, atb_setting, "Variable O&M")
                else
                    if !warned_missing_atb_vom && detailed_log
                        @aleaf_warn "No VOM for technology $(network_data["gen_technology"][id]["UNITGROUP"]) in ATB. VOM will be set to zero"
                        warned_missing_atb_vom = true
                    end
                    network_data["gen_technology"][id]["Annual_VOM"][y] = 0.0
                end
            elseif network_data["gen_technology"][id]["VOM"] == "ESGC"
                network_data["gen_technology"][id]["Annual_VOM"][y] = find_ESGC_data(network_data["gen_technology"][id]["UNITGROUP"], "VOM") 
            else
                network_data["gen_technology"][id]["Annual_VOM"][y] = network_data["gen_technology"][id]["VOM"] 
            end

            # Fuel cost
            if network_data["gen_technology"][id]["FC"] in ["Fuel", "Fuel-Regional"]    # even if the fuel cost is regional, we will use the system-wide fuel cost here. The regional cost will be updated later.
                network_data["gen_technology"][id]["Annual_FC"][y] = 0.0 # fuel price will be updated later in the update_regional_fule_prices_With_pu function.
            elseif network_data["gen_technology"][id]["FC"] == "Water Value"
                network_data["gen_technology"][id]["Annual_FC"][y] = 0.0 # water value will be updated in the objective function
            elseif network_data["gen_technology"][id]["FC"] == "ATB"
                @aleaf_warn "FC='ATB' is not supported for UNITGROUP $(network_data["gen_technology"][id]["UNITGROUP"]) in scenario $scenario"
                terminate_with_error(; msg="Fuel cost source FC='ATB' is not implemented for UNITGROUP $(network_data["gen_technology"][id]["UNITGROUP"]) in scenario $scenario. Please update FC input or implement ATB-based fuel cost lookup.")
            else
                network_data["gen_technology"][id]["Annual_FC"][y] = network_data["gen_technology"][id]["FC"]
            end       

        end
    end

    # 3) calculate unit marginal generation costs ($/MWh)
    
    for id in keys(network_data["gen_technology"])

        scenario = network_data["case_id"]
        network_data["gen_technology"][id]["Annual_MC"] = Dict{String,Any}()
        
        for y in keys(network_data["planning_stages"])
            network_data["gen_technology"][id]["Annual_MC"][y] = 0.0 # MC value will be updated in the update_regional_fule_prices_With_pu function.
        end
                
    end

    # impoundment hydro option update
    hydro_option = ALEAF_setting["Simulation Configuration"][scenario]["Reservoir_Hydro_Operation_Option"]
    for id in keys(network_data["gen_technology"])
        if network_data["gen_technology"][id]["Hydro_Flag"] == "IMPOUNDMENT"
            network_data["gen_technology"][id]["FUEL_LIMIT"] = hydro_option
        end
    end

end


function check_and_update_plant_data!(am, ALEAF_setting::Dict{String,<:Any}, plant_data, gen_tech_info, scenario::Int)
    
    plant_data_backup = deepcopy(plant_data)
    
    pu_power_base = ALEAF_setting["Simulation Setting"]["per_unit_base_value"]
    pu_econ_base = ALEAF_setting["Simulation Setting"]["per_unit_econ_base_value"] / pu_power_base    # per unit economic base
    detailed_log = haskey(ALEAF_setting, "Simulation Setting") &&
                   (lowercase(string(get(ALEAF_setting["Simulation Setting"], "logging_level_value", "simple"))) == "detailed")
    
    # 1) read ATB data
    ATB_data_raw = ALEAF_setting["ATB_data_raw"]

    function find_ESGC_data(unit_group_id::String, data_label::String)
        ESGC_ID = ALEAF_setting["Simulation Configuration"][string(scenario)]["ESGC_Setting_ID"]
        storage_data = ALEAF_setting["Storage Cost and Performance"] 
        for idx in keys(storage_data)
            entry = storage_data[idx]  
            if entry["ESGC_Setting_ID"] == ESGC_ID && entry["UNITGROUP"] == unit_group_id
                return get(entry, data_label, nothing)
            end
        end
    
        return nothing
    end

    function find_fuel_cost(fuel_data, fuel_scenario::String, fuel_type::String, year::Int)
        fuel_cost = 0.0
        found_fuel_flag = false
        for idx in keys(fuel_data)
            if (fuel_data[idx]["Scenario"] == fuel_scenario) && (fuel_data[idx]["FUEL"] == fuel_type) && (fuel_data[idx]["Region"] == "System-wide")
                fuel_cost = fuel_data[idx][string(year)]
                found_fuel_flag = true
            end
        end

        if found_fuel_flag == true
            return fuel_cost
        else
            if detailed_log
                @aleaf_warn "Failed to load fuel price for type $fuel_type in scenario $fuel_scenario. Fuel price will be set to zero"
            end
            return fuel_cost
        end
    end

    # ATB setting 
    atb_setting = get_atb_setting(ALEAF_setting, ALEAF_setting["Simulation Configuration"][string(scenario)]["ATB_Setting_ID"], plant_data["UNITGROUP"])

    # udpate MC flag
    update_plant_parameters_flag = false
    update_MC_flag = false
    
    # Check and update Storage fisrt using ESGC
    if plant_data["FOM"] == "Gen_Tech"
        plant_data["FOM"] = gen_tech_info["FOM"]
    else 
        update_plant_parameters_flag = true
        FOM_value = plant_data["FOM"]
        if FOM_value == "ESGC"
            FOM_value = find_ESGC_data(plant_data["UNITGROUP"], "FOM")
        elseif FOM_value == "ATB"
            FOM_value = get_atb_value(ATB_data_raw, atb_setting, "Fixed O&M")
        end
        for y in keys(plant_data["Annual_FOM"])
            plant_data["Annual_FOM"][y] = FOM_value / pu_econ_base  # apply per unit economic base
        end
    end

    # BATEFF 
    BATEFF_gen_tech_reference = gen_tech_info["BATEFF"]
    if haskey(gen_tech_info, "BATEFF_reference")
        BATEFF_gen_tech_reference = gen_tech_info["BATEFF_reference"]
    end
    
    if plant_data["BATEFF"] == "Gen_Tech"
        plant_data["BATEFF"] = gen_tech_info["BATEFF"]
    else 
        if plant_data["BATEFF"] != BATEFF_gen_tech_reference
            update_plant_parameters_flag = true
            if plant_data["BATEFF"] == "ESGC"
                plant_data["BATEFF"] = find_ESGC_data(plant_data["UNITGROUP"], "RTE")
            end
        else
            plant_data["BATEFF"] = gen_tech_info["BATEFF"]
        end
    end

    # AET 
    AET_gen_tech_reference = gen_tech_info["AET"]
    if haskey(gen_tech_info, "AET_reference")
        AET_gen_tech_reference = gen_tech_info["AET_reference"]
    end

    if plant_data["AET"] == "Gen_Tech"
        plant_data["AET"] = gen_tech_info["AET"]
    else 
        if plant_data["AET"] != AET_gen_tech_reference
            update_plant_parameters_flag = true
            if plant_data["AET"] == "ESGC"
                plant_data["AET"] = find_ESGC_data(plant_data["UNITGROUP"], "AET") / pu_power_base  # apply per unit power base
            end
        else
            plant_data["AET"] = gen_tech_info["AET"] # the gen_tech_info data is already in per unit base
        end   
    end
                     
    # VOM
    if plant_data["VOM"] == "Gen_Tech"
        plant_data["VOM"] = gen_tech_info["VOM"]
    else 
        update_plant_parameters_flag = true
        update_MC_flag = true
        VOM_value = plant_data["VOM"]
        if plant_data["VOM"] == "ESGC"
            VOM_value = find_ESGC_data(plant_data["UNITGROUP"], "VOM")
        elseif plant_data["VOM"] == "ATB"
            # ATB 2025 renamed several technologies; a name absent here silently zeroes VOM instead of reading it
            if atb_setting["Tech"] in ATB_TECHS_WITH_VOM
                VOM_value = get_atb_value(ATB_data_raw, atb_setting, "Variable O&M") 
            else
                @aleaf_warn "No VOM for technology $(plant_data["UNITGROUP"]) in ATB. VOM will be set to zero"
                VOM_value = 0.0
            end
        end
        for y in keys(plant_data["Annual_VOM"])
            plant_data["Annual_VOM"][y] = VOM_value / pu_econ_base  # apply per unit economic base
        end
    end

    # FC
    if plant_data["FC"] == "Gen_Tech"
        plant_data["FC"] = gen_tech_info["FC"]
    else 
        update_plant_parameters_flag = true
        update_MC_flag = true
        
        fuel_cost_value = plant_data["FC"]
        
        if fuel_cost_value in ["Fuel", "Fuel-Regional"]    
            for y in keys(plant_data["Annual_FC"])
                plant_data["Annual_FC"][y] = 0.0 # fuel cost value will be updated later 
            end
            
        elseif fuel_cost_value == "Water Value"
            for y in keys(plant_data["Annual_FC"])
                plant_data["Annual_FC"][y] = 0.0 # water value will be updated in the objective function
            end
        else
            for y in keys(plant_data["Annual_FC"])
                plant_data["Annual_FC"][y] = plant_data["FC"] / pu_econ_base    # apply per unit economic base
            end
        end
    end

    if update_MC_flag == true 
        for y in keys(plant_data["Annual_MC"])
            plant_data["Annual_MC"][y] = 0.0 # MC value will be updated later
        end
    end

    # impoundment hydro option update
    hydro_option = ALEAF_setting["Simulation Configuration"][string(scenario)]["Reservoir_Hydro_Operation_Option"]
    if gen_tech_info["Hydro_Flag"] == "IMPOUNDMENT"
        plant_data["FUEL_LIMIT"] = hydro_option
    end

    # update remaining fields and apply per unit base
    for field in keys(plant_data)

        if ismissing(plant_data[field]) || plant_data[field] == "missing"
            # do nothing (skip missing values)

        elseif plant_data[field] == "Gen_Tech"
            plant_data[field] = gen_tech_info[field]
            # already in per-unit base

        else
            if field in ["NLC", "SUC", "SDC", "DECC"]
                plant_data[field] /= pu_econ_base

            elseif field in ["reg_cost", "spin_cost", "nspin_cost", "flex_cost"] && ALEAF_setting["Planning Design"]["reserve_cost_type_flag"] == "absolute"
                plant_data[field] /= pu_econ_base

            elseif field in ["Emission_CO2", "Emission_1", "Emission_2", "Emission_3"]
                plant_data[field] *= pu_power_base
            end
        end
    end

    # add unit type
    plant_data["Tech_Type"] = "Existing"

end


function get_annual_investment_cost(capex, network_data::Dict{String,<:Any}, ALEAF_setting::Dict{String,<:Any}, tech_id, ATB_data_raw, atb_setting)

    WACC = ALEAF_setting["Planning Design"]["WACC_value"]
    cost_recovery_period = network_data["gen_technology"][tech_id]["crpyears"]
    capital_recovery_factor = (WACC * ((1+WACC)^cost_recovery_period)) / ((1+WACC)^cost_recovery_period - 1)
    
    project_finance_factor = ALEAF_setting["Planning Design"]["project_finance_factor_value"]
    fixed_charge_rate = capital_recovery_factor * project_finance_factor

    annual_cost = capex * fixed_charge_rate

    return annual_cost
end


function get_network_data!(network_data::Dict{String, <:Any}, ALEAF_setting::Dict{String, <:Any})

    # O(1) column indexes over sub_area_mapping, built once (keys-order preserved => keep-first ==
    # original first-match). Keyed on the raw cell value to reproduce == semantics exactly.
    _firstrow_by_col = Dict{String, Dict{Any, Any}}()
    function _get_firstrow_index(sub_area_mapping, col)
        get!(_firstrow_by_col, col) do
            idx = Dict{Any, Any}()
            for key in keys(sub_area_mapping)
                v = sub_area_mapping[key][col]
                v === missing && continue
                get!(idx, v, sub_area_mapping[key])
            end
            idx
        end
    end
    # collect-all member index for a (parent_col, child_col) pair, in keys order (dup-preserving to
    # match the original push! loops; callers keep their existing unique()).
    _members_by_col = Dict{Tuple{String, String}, Dict{Any, Vector{Any}}}()
    function _get_members_index(sub_area_mapping, parent_col, child_col)
        get!(_members_by_col, (parent_col, child_col)) do
            idx = Dict{Any, Vector{Any}}()
            for key in keys(sub_area_mapping)
                pv = sub_area_mapping[key][parent_col]
                pv === missing && continue
                push!(get!(idx, pv, Vector{Any}()), sub_area_mapping[key][child_col])
            end
            idx
        end
    end

    # aggregated-region member -> aggregated bus, built after the bus swap (first wins).
    _member2aggbus = Dict{String, String}()
    function find_new_bus_i_using_agg_resolution(original_id)
        original_id = string(original_id)
        bus_idx = get(_member2aggbus, original_id, "")
        if bus_idx == ""
            return 0, 0, 0
        end
        bus = network_data["bus"][bus_idx]
        lookup_type = bus["region_config"]["region_lookup_type"]
        regional_resolution = bus["region_config"][lookup_type]
        region_mapping_info = bus["region_mapping_info"]
        return region_mapping_info[regional_resolution], bus_idx, region_mapping_info
    end

    function find_new_bus_i(sub_area_mapping, original_id, sub_area_key::String)
        reference_resolution_id = ALEAF_setting["network_resolution_level"]["1"]["Network_Resolution_ID"]
        row = get(_get_firstrow_index(sub_area_mapping, reference_resolution_id), original_id, nothing)
        return row === nothing ? nothing : row[sub_area_key]
    end

    function find_region_id(sub_area_mapping, original_id, reduction_level, region_key)
        row = get(_get_firstrow_index(sub_area_mapping, reduction_level), original_id, nothing)
        return row === nothing ? nothing : row[region_key]
    end

    function find_region_mapping_info(sub_area_mapping, original_id, reduction_level)
        return get(_get_firstrow_index(sub_area_mapping, reduction_level), original_id, nothing)
    end

    function find_sub_region_mapping_info(sub_area_mapping, parent_region_id, parent_region_resolution, sub_resgionl_reduction_level)
        members = get(_get_members_index(sub_area_mapping, parent_region_resolution, sub_resgionl_reduction_level), parent_region_id, nothing)
        return members === nothing ? Any[] : copy(members)
    end

    function collect_inter_tie_capacity(from_bus::String, to_bus::String)
        inter_tie_capacity = 0.0
        total_length = 0.0
        max_length = 0.0
        max_inter_tie_capacity = 0.0
        expansion_flag = false
        merged_line_list = []
        merged_line_uids = []
        sum_inv_x = 0.0   # accumulate parallel conductance Σ(1/x_i) for DC equivalent reactance

        # accumulate parallel reactance 
        function accumulate_reactance!(branch_k)
            x_i = get(branch_k, "br_x_pu", nothing)
            if x_i === nothing
                return
            end
            x_i = Float64(x_i)
            if !isfinite(x_i) || abs(x_i) < 1e-9
                # missing/zero/near-zero reactance: skip to avoid divide-by-zero blow-up
                return
            end
            sum_inv_x += 1.0 / x_i
        end

        for k in keys(network_data["branch"])

            if network_data["branch"][k]["model_flag"] == true

                if (string(network_data["branch"][k]["f_bus"]) == from_bus) & (string(network_data["branch"][k]["t_bus"]) == to_bus)
                    inter_tie_capacity += network_data["branch"][k]["rate_a"]
                    if network_data["branch"][k]["expansion_flag"] == true
                        max_inter_tie_capacity += network_data["branch"][k]["max_rate_a"]
                        expansion_flag = true
                    end

                    push!(merged_line_list, k)
                    push!(merged_line_uids, network_data["branch"][k]["UID"])
                    accumulate_reactance!(network_data["branch"][k])

                    total_length += network_data["branch"][k]["Length"]
                    if network_data["branch"][k]["Length"] >= max_length
                        max_length = network_data["branch"][k]["Length"]
                    end

                elseif (string(network_data["branch"][k]["t_bus"]) == from_bus) & (string(network_data["branch"][k]["f_bus"]) == to_bus)
                    inter_tie_capacity += network_data["branch"][k]["rate_a"]
                    if network_data["branch"][k]["expansion_flag"] == true
                        max_inter_tie_capacity += network_data["branch"][k]["max_rate_a"]
                        expansion_flag = true
                    end

                    push!(merged_line_list, k)
                    push!(merged_line_uids, network_data["branch"][k]["UID"])
                    accumulate_reactance!(network_data["branch"][k])

                    total_length += network_data["branch"][k]["Length"]
                    if network_data["branch"][k]["Length"] >= max_length
                        max_length = network_data["branch"][k]["Length"]
                    end
                end
            end
        end

        # DC equivalent reactance of the (parallel) merged constituent lines: 1 / Σ(1/x_i).
        br_x_pu_eq = sum_inv_x != 0.0 ? 1.0 / sum_inv_x : 0.01

        return inter_tie_capacity, total_length, max_inter_tie_capacity, merged_line_list, merged_line_uids, expansion_flag, br_x_pu_eq
    end


    # Memoized set of modeled bus_i (each bus's region_mapping_info at its own lookup resolution), built
    # once from network_data["bus"]. Replaces the former per-call O(buses) scan so intertie_checker is O(1);
    # only called during branch aggregation (after the bus swap) where network_data["bus"] is stable.
    _intertie_modeled = Ref{Union{Nothing, Set{Any}}}(nothing)
    function intertie_checker(from_bus, to_bus) # check whether the line (from_bus, to_bus) is an inter-tie line or not

        if _intertie_modeled[] === nothing
            s = Set{Any}()
            for key in keys(network_data["bus"])
                lookup_type = network_data["bus"][key]["region_config"]["region_lookup_type"]
                regional_resolution = network_data["bus"][key]["region_config"][lookup_type]
                push!(s, network_data["bus"][key]["region_mapping_info"][regional_resolution])
            end
            _intertie_modeled[] = s
        end
        modeled = _intertie_modeled[]

        # A bus is "external" iff its id is not among the modeled bus_i. f_bus_flag/t_bus_flag in the original
        # were set to from_bus/to_bus on match (so "flags equal" <=> from_bus == to_bus for two modeled buses).
        f_external = !(from_bus in modeled)
        t_external = !(to_bus in modeled)

        external_flag = false
        inter_tie_flag = false
        intra_line_flag = true

        if f_external && t_external
            external_flag = true
            inter_tie_flag = false
            intra_line_flag = false
        elseif f_external || t_external
            external_flag = false
            inter_tie_flag = true
            intra_line_flag = false
        elseif from_bus == to_bus
            external_flag = false
            inter_tie_flag = false
            intra_line_flag = true
        else
            external_flag = false
            inter_tie_flag = false
            intra_line_flag = false
        end

        return external_flag, inter_tie_flag, intra_line_flag
    end


    function aggregate_branches_using_database!(ALEAF_setting, network_data)

        # Capture the FULL nodal AC network (incl. intra-region lines that get deleted below) and the
        # nodal->region map, so the corridor reactance estimator can rebuild the nodal susceptance Laplacian.
        # Ids are canonicalized to String (network_data stores bus ids as strings).
        ptdf_nodal_ac = Tuple{String,String,Float64}[]   # (nodal_f, nodal_t, x) for AC model branches
        ptdf_bus2region = Dict{String,String}()          # nodal bus -> region (aggregated) bus id

        # Complete capture pre-pass over a key SNAPSHOT (collect): the remap loop below deletes while
        # iterating keys(), which can skip branches — an incomplete bus2region drops corridors and
        # disconnects the region graph. This read-only pass leaves the remap loop byte-identical.
        for k in collect(keys(network_data["branch"]))
            b = network_data["branch"][k]
            get(b, "model_flag", true) == true || continue
            of = b["f_bus"]; ot = b["t_bus"]
            nf, _bi1, _rm1 = find_new_bus_i_using_agg_resolution(of)
            nt, _bi2, _rm2 = find_new_bus_i_using_agg_resolution(ot)
            ext, _it, _il = intertie_checker(nf, nt)
            ext == true && continue
            ptdf_bus2region[string(of)] = string(nf)
            ptdf_bus2region[string(ot)] = string(nt)
            (get(b, "dc_line", false) == true) && continue
            xi = get(b, "br_x_pu", nothing); xi === nothing && continue
            xv = Float64(xi); (isfinite(xv) && abs(xv) >= 1e-9) || continue
            push!(ptdf_nodal_ac, (string(of), string(ot), xv))
        end

        # Update branch from and to id & delete internal lines
        for k in keys(network_data["branch"])

            original_f_bus_i = deepcopy(network_data["branch"][k]["f_bus"])
            original_t_bus_i = deepcopy(network_data["branch"][k]["t_bus"])

            new_f_bus_i, bus_idx, region_mapping_info = find_new_bus_i_using_agg_resolution(original_f_bus_i)
            network_data["branch"][k]["original_f_bus"] = original_f_bus_i
            network_data["branch"][k]["f_bus"] = new_f_bus_i

            new_t_bus_i, bus_idx, region_mapping_info = find_new_bus_i_using_agg_resolution(original_t_bus_i)
            network_data["branch"][k]["original_t_bus"] = original_t_bus_i
            network_data["branch"][k]["t_bus"] = new_t_bus_i

            external_flag, inter_tie_flag, intra_line_flag = intertie_checker(new_f_bus_i, new_t_bus_i)

            if external_flag == true
                delete!(network_data["branch"], k)
            elseif intra_line_flag == true
                delete!(network_data["branch"], k)
            end
        end

        # One pass over surviving branches, grouped by normalized bus_i pair: replicates
        # collect_inter_tie_capacity's merge in O(branches), in keys(branch) order.
        # DC lines are keyed separately from AC so a DC tie parallel to an AC line is never
        # merged into (and demoted by) the AC corridor.
        # Buses modeled at the finest (Bus) resolution. A corridor whose BOTH endpoints are such buses
        # keeps every physical line SEPARATE (full nodal: each line its own rate_a/br_x_pu/thermal
        # limit + shared-angle KVL) instead of being parallel-merged. Any corridor touching an
        # aggregated (zonal) bus still merges. This preserves the multi-resolution option.
        nodal_busi = Set{String}()
        for key in keys(network_data["bus"])
            rc = network_data["bus"][key]["region_config"]
            regional_resolution = rc[rc["region_lookup_type"]]        # e.g. "Bus" / "County" / "BA"
            if regional_resolution == "Bus"
                push!(nodal_busi, string(network_data["bus"][key]["region_mapping_info"][regional_resolution]))
            end
        end

        # pair_key 4th slot: "" for merged (zonal-touching) corridors; string(k) for nodal↔nodal lines
        # so each physical line becomes its own single-line corridor (emitted with native rate_a/x).
        inter_tie_accum = Dict{Tuple{String, String, Bool, String}, Dict{String, Any}}()
        for k in keys(network_data["branch"])
            branch_k = network_data["branch"][k]
            branch_k["model_flag"] == true || continue

            f_bus_i = string(branch_k["f_bus"])
            t_bus_i = string(branch_k["t_bus"])
            is_dc = get(branch_k, "dc_line", false) == true
            keep_separate = (f_bus_i in nodal_busi) && (t_bus_i in nodal_busi)
            disambig = keep_separate ? string(k) : ""
            # Normalize the (from_bus_i, to_bus_i) key so both line directions map to the same pair,
            # matching collect_inter_tie_capacity's forward/reverse handling.
            pair_key = f_bus_i <= t_bus_i ? (f_bus_i, t_bus_i, is_dc, disambig) : (t_bus_i, f_bus_i, is_dc, disambig)

            acc = get!(inter_tie_accum, pair_key) do
                Dict{String, Any}(
                    "inter_tie_capacity" => 0.0,
                    "total_length" => 0.0,
                    "weighted_length_capacity_sum" => 0.0,
                    "max_inter_tie_capacity" => 0.0,
                    "expansion_flag" => false,
                    "merged_line_list" => [],
                    "merged_line_uids" => [],
                    "sum_inv_x" => 0.0,
                    # min_i(rate_i * x_i) for the DC-consistent corridor limit rate_a = b0 * min(rate_i*x_i).
                    "min_rate_x" => Inf,
                    # AND-folded across constituents so a merged corridor is DC only if every line is a DC tie.
                    "dc_line" => true,
                    # true when both endpoints are Bus-resolution: this corridor is a single kept line.
                    "nodal_kept" => keep_separate,
                )
            end

            acc["inter_tie_capacity"] += branch_k["rate_a"]
            if branch_k["expansion_flag"] == true
                acc["max_inter_tie_capacity"] += branch_k["max_rate_a"]
                acc["expansion_flag"] = true
            end

            push!(acc["merged_line_list"], k)
            push!(acc["merged_line_uids"], branch_k["UID"])

            # Fail-safe: corridor stays DC only while every constituent is a DC tie.
            acc["dc_line"] = acc["dc_line"] && (get(branch_k, "dc_line", false) == true)

            # accumulate parallel conductance Σ(1/x_i); mirror accumulate_reactance! guards exactly
            x_i = get(branch_k, "br_x_pu", nothing)
            if x_i !== nothing
                x_i = Float64(x_i)
                if isfinite(x_i) && abs(x_i) >= 1e-9
                    acc["sum_inv_x"] += 1.0 / x_i
                    # rate_i*x_i = rate_i / b_i; the min over legs gives the DC-consistent corridor limit.
                    acc["min_rate_x"] = min(acc["min_rate_x"], branch_k["rate_a"] * x_i)
                end
            end

            # Length is stored uniformly in km for every row; convert km -> miles here for the
            # $/MW-mile expansion cost.
            acc["total_length"] += branch_k["Length"] * 0.621371
            # Capacity-weighted numerator for N>1 corridors (see emission site below).
            acc["weighted_length_capacity_sum"] += branch_k["Length"] * 0.621371 * branch_k["rate_a"]
        end

        # create new branch (inter-tie only)
        # Emit by iterating the corridor accumulator (O(corridors)) but in the SAME visitation order as
        # the original nested bus^2 loop so branch IDs stay byte-identical: outer = numerically-smaller
        # bus key (pos in keys order), inner = larger key, AC (is_dc=false) before DC.
        new_branch = Dict{String, Any}()
        new_branch_id = 1
        primal_branch_id = 1

        pos = Dict{String, Int}(k => i for (i, k) in enumerate(keys(network_data["bus"])))
        busi2key = Dict{String, String}(string(network_data["bus"][k]["bus_i"]) => k for k in keys(network_data["bus"]))

        emit_order = Vector{Tuple{Int, Int, Int, String, String, String, Bool, Dict{String, Any}}}()
        for (pair_key, acc) in inter_tie_accum
            a_busi, b_busi, is_dc, disambig = pair_key
            ka = get(busi2key, a_busi, nothing)
            kb = get(busi2key, b_busi, nothing)
            # Skip corridors whose endpoint bus_i is not a real modeled bus (e.g. unmatched bus_i 0);
            # the original nested bus^2 loop never visited these, so they must not be emitted.
            (ka === nothing || kb === nothing) && continue
            na = parse(Int, ka); nb = parse(Int, kb)
            na == nb && continue  # original nested loop requires to_bus > from_bus (distinct keys)
            from_key, to_key = na <= nb ? (na, nb) : (nb, na)
            # disambig (per-line id for nodal↔nodal kept lines) is the final tiebreak so parallel
            # kept lines on the same bus-pair get a deterministic, reproducible branch-id order.
            push!(emit_order, (pos[string(from_key)], pos[string(to_key)], is_dc ? 1 : 0, disambig, string(from_key), string(to_key), is_dc, acc))
        end
        sort!(emit_order; by = t -> (t[1], t[2], t[3], t[4]))

        for (_, _, _, _disambig, from_key, to_key, is_dc, acc) in emit_order
            from_bus = parse(Int, from_key)
            to_bus = parse(Int, to_key)

            inter_tie_capacity = acc["inter_tie_capacity"]
            # Capacity-weighted length for N>1 corridors (e.g. multi-pole DC ties)
            # avoids pricing parallel circuits at N^2 instead of N; N=1 is unaffected.
            total_length = length(acc["merged_line_list"]) > 1 ? acc["weighted_length_capacity_sum"] / inter_tie_capacity : acc["total_length"]
            # DC equivalent reactance of the (parallel) merged lines: 1 / Σ(1/x_i).
            b0 = acc["sum_inv_x"]
            br_x_pu_eq = b0 != 0.0 ? 1.0 / b0 : 0.01

            # AC corridor rating = summed thermal capacity. The DC-consistent rate_a = b0*min_i(rate_i*x_i)
            # form under-rates spatially-aggregated corridors: it assumes the constituent legs are electrically
            # parallel (share endpoints) so the weakest leg caps the corridor, but region corridors bundle
            # lines with DIFFERENT nodal endpoints whose true deliverability is a max-flow ~= Σ rate_i.
            # Enabling it collapsed multi-leg imports (e.g. Harris->Brazoria 9 legs -> 38%) and caused TWh-scale
            # ENS, so it stays off.
            rate_a_corr = inter_tie_capacity
            max_rate_a_corr = acc["max_inter_tie_capacity"]

            if inter_tie_capacity!= 0
                if total_length != 0

                    new_branch[string(new_branch_id)] = Dict{String, Any}(
                        "f_bus" => from_bus,
                        "t_bus" => to_bus,
                        "f_bus_id" => network_data["bus"][from_key]["bus_i"],
                        "t_bus_id" => network_data["bus"][to_key]["bus_i"],
                        "rate_a" => rate_a_corr,
                        "max_rate_a" => max_rate_a_corr,
                        "model_flag" => true,
                        "expansion_flag" => acc["expansion_flag"],
                        "dc_line" => acc["dc_line"],
                        "br_x_pu" => br_x_pu_eq,
                        "length" => total_length,
                        # Fresh arrays per emission (original built new lists each call); avoids
                        # shared references if multiple bus pairs map to the same bus_i pair.
                        "merged_line_list" => copy(acc["merged_line_list"]),
                        "merged_line_uids" => copy(acc["merged_line_uids"]),
                        "nodal_kept" => get(acc, "nodal_kept", false),
                        "primal_branch_id" => primal_branch_id
                    )
                    new_branch_id += 1
                    primal_branch_id += 1

                end
            end
        end

        # ---- corridor reactance for spatially-aggregated AC corridors (snapshot estimator) ----
        # br_x_pu is estimated from synthetic snapshots on the fine network so the aggregated network
        # reproduces its cross-border flows; nodal_kept lines enter as fixed edges. DC ties, nodal_kept
        # lines and pure-nodal / 1:1 runs keep their values; the selection-only build has no branches.
        aggregation_occurred = any((length(br["merged_line_list"]) > 1) && !get(br, "nodal_kept", false)
                                   for br in values(new_branch))
        pf_mode = get(ALEAF_setting["Simulation Setting"], "power_flow_mode_flag", "")
        selection_build = get(network_data, "__selection_agg_override", nothing) !== nothing
        if aggregation_occurred && pf_mode in ("B-theta", "PTDF") && !selection_build
            try
                # corridor key = sorted pair of region ids of the first member line's fine endpoints
                corr_key(br) = begin
                    of = string(network_data["branch"][br["merged_line_list"][1]]["original_f_bus"])
                    ot = string(network_data["branch"][br["merged_line_list"][1]]["original_t_bus"])
                    ca = get(ptdf_bus2region, of, ""); cb = get(ptdf_bus2region, ot, "")
                    ca <= cb ? (ca, cb) : (cb, ca)
                end
                ptdf_corridors = PTDFCorridor[]
                kept_lines = KeptLine[]
                for bk in sort!(collect(keys(new_branch)), by = k -> parse(Int, k))
                    br = new_branch[bk]
                    (get(br, "dc_line", false) == true) && continue
                    isempty(br["merged_line_list"]) && continue
                    if get(br, "nodal_kept", false)
                        lb = network_data["branch"][br["merged_line_list"][1]]
                        xv = Float64(br["br_x_pu"]); (isfinite(xv) && abs(xv) >= 1e-9) || continue
                        push!(kept_lines, KeptLine(bk, string(lb["original_f_bus"]), string(lb["original_t_bus"]), 1.0/abs(xv)))
                        continue
                    end
                    legs = Tuple{String,String,Float64}[]
                    for lk in br["merged_line_list"]
                        haskey(network_data["branch"], lk) || continue
                        lb = network_data["branch"][lk]
                        xi = get(lb, "br_x_pu", nothing); xi === nothing && continue
                        xv = Float64(xi); (isfinite(xv) && abs(xv) >= 1e-9) || continue
                        push!(legs, (string(lb["original_f_bus"]), string(lb["original_t_bus"]), 1.0/abs(xv)))
                    end
                    isempty(legs) && continue
                    k = corr_key(br); (k[1]=="" || k[2]=="" || k[1]==k[2]) && continue
                    push!(ptdf_corridors, PTDFCorridor(k[1], k[2], legs))
                end
                # snapshot inputs: fine-bus peak loads and plant capacities. VRE_Flag / Profile_Type are still
                # "Gen_Tech" placeholders here, so plants are classified by UNIT_CATEGORY; hybrid and demand
                # sheets are not used.
                fine_load = Dict{String,Float64}()
                for bus in values(network_data["bus"])
                    for (b, v) in get(get(bus, "aggregation_info", Dict()), "original_load_(bus_i, MW)", Dict())
                        fine_load[string(b)] = _mw(v)
                    end
                end
                fine_plants = [(string(p["original_bus_i"]), string(get(p, "UNIT_CATEGORY", "")), get(p, "CAP", 0.0))
                               for p in values(network_data["plant"]) if haskey(p, "original_bus_i")]
                est = reduce_network_by_snapshots(ptdf_nodal_ac, ptdf_bus2region, ptdf_corridors, kept_lines,
                                                  fine_injection_data(fine_load, fine_plants))
                if haskey(ALEAF_setting, "Simulation Setting") && lowercase(string(get(ALEAF_setting["Simulation Setting"], "logging_level_value", "simple"))) == "detailed"
                    @aleaf_info "Network reduction (snapshot): $(est.n_islands) island(s), $(est.n_corridors) corridors, $(est.n_kept) kept lines, $(est.n_fallback) non-positive, $(est.n_clamped) at bound, $(est.n_island_fallback) island fallback, held-out flow err $(round(est.err_parallel, digits=3)) -> $(round(est.err_final, digits=3)), leakage $(round(est.leakage, sigdigits=2)), coverage $(round(est.coverage, digits=3)), checksum $(string(est.checksum, base=16)), $(round(est.elapsed, digits=1))s"
                end
                est.coverage < 0.95 && @aleaf_warn "Network reduction: only $(round(100*est.coverage, digits=1))% of load/capacity sits on reduced-region buses (check bus id mapping)"
                est.max_components > 1 && @aleaf_warn "Network reduction: aggregated AC network has $(est.max_components) components in one fine island"
                for br in values(new_branch)
                    (get(br, "dc_line", false) == true) && continue
                    get(br, "nodal_kept", false) && continue
                    isempty(br["merged_line_list"]) && continue
                    key = corr_key(br)
                    if haskey(est.x_by_pair, key) && isfinite(est.x_by_pair[key]) && est.x_by_pair[key] > 0
                        br["br_x_pu"] = est.x_by_pair[key]
                    end
                end
            catch err
                @aleaf_warn "Network reduction failed ($(sprint(showerror, err))) -> keeping parallel-combine reactances"
            end
        end
        # PTDF mode: rows from the final branch list (f_bus -> t_bus, br_x_pu), so PTDF == B-theta
        if pf_mode == "PTDF" && !selection_build
            n_comp = attach_branch_ptdf!(new_branch, length(network_data["bus"]))
            n_comp > 1 && @aleaf_warn "PTDF mode: AC network has $n_comp components; PTDF balances system-wide while B-theta balances each component"
        end

        delete!(network_data, "branch")
        network_data["branch"] = deepcopy(new_branch)
        new_branch = nothing
    end

    network_file_location = ALEAF_setting["network_data_file_location"]
    config_file_location = get(ALEAF_setting, "network_config_file_location", nothing)

    # read data from the network workbook (open once, read every category from the shared handle)
    data_category_list = ["bus", "plant", "branch", "demand", "hybrid"]
    XLSX.openxlsx(network_file_location) do network_xf
        for data_category in data_category_list
            network_data[data_category] = read_xlsx_return_dict_string_any(network_xf, data_category; first_row_value = 2)
        end
    end

    # Clean up the plant data
    for g in collect(keys(network_data["plant"]))
        if network_data["plant"][g]["CAP"] == 0
            delete!(network_data["plant"], g)
        end
    end
    
    network_data["supplement_network_data"] = Dict{String, Any}()
    network_data["additional_network_data"] = network_data["supplement_network_data"]
    if isnothing(config_file_location)
        bus_idx_by_bus_id = Dict{String, String}()
        ALEAF_setting["Network Setting"] = Dict{String, Any}()
        ALEAF_setting["network_resolution_level"] = Dict{String, Any}()
        ALEAF_setting["sub_area_mapping"] = Dict{String, Any}()
        # Without a config workbook, keep the original network topology and synthesize minimal region metadata.
        for (bus_idx, bus_data) in network_data["bus"]
            bus_id = string(bus_data["bus_i"])
            bus_idx_by_bus_id[bus_id] = bus_idx
            bus_name = string(get(bus_data, "bus_name", bus_id))
            bus_data["original_bus_i"] = bus_id
            bus_data["region_mapping_info"] = Dict{String, Any}(
                "system" => "system",
                "System" => "system",
                "bus_i" => bus_data["bus_i"],
                "bus_name" => bus_name,
            )
            bus_data["region_config"] = Dict{String, Any}(
                "region_type" => "original",
                "aggregation_resolution_type" => "bus_i",
                "aggregation_resolution_value" => bus_data["bus_i"],
                "subregional_aggregation_type" => "NA",
                "subregional_aggregation_value" => "NA",
                "region_lookup_type" => "aggregation_resolution_type",
                "region_lookup_value" => "aggregation_resolution_value",
                "parent_region_id" => bus_data["bus_i"],
                "region_name" => bus_name,
                "parent_bus_name" => bus_name,
            )
            bus_data["aggregation_info"] = Dict{String, Any}(
                "aggregated_regions_bus_i" => [bus_id],
                "aggregated_regions_bus_idx" => [parse(Int, bus_idx)],
                "aggregated_regions_bus_(idx, name)" => [(bus_idx, bus_name)],
                "original_load_(bus_i, MW)" => Dict{String, Float64}(bus_id => Float64(get(bus_data, "MW load", 0.0))),
                "local_gen_idx" => String[],
                "new_local_gen_idx" => Int[],
            )
            bus_data["RA_ELCC_Calculation_Flag"] = true
        end

        # Plants, hybrids, and flexible demand still need bus and region references even when no reconfiguration is applied.
        for (idx, plant_data) in network_data["plant"]
            plant_data["original_bus_i"] = plant_data["bus_ID"]
            plant_data["bus_i"] = plant_data["bus_ID"]
            plant_data["bus_idx"] = parse(Int, get(bus_idx_by_bus_id, string(plant_data["bus_ID"]), "0"))
            plant_data["bypass_parameter_check"] = get(plant_data, "bypass_parameter_check", false)
            plant_data["hybrid_type"] = get(plant_data, "hybrid_type", "NA")
            plant_data["region_mapping_info"] = Dict{String, Any}(
                "system" => "system",
                "System" => "system",
                "bus_i" => plant_data["bus_ID"],
                "bus_ID" => plant_data["bus_ID"],
                "bus_name" => plant_data["bus_name"],
            )
            bus_idx = get(bus_idx_by_bus_id, string(plant_data["bus_ID"]), nothing)
            if !isnothing(bus_idx)
                push!(network_data["bus"][bus_idx]["aggregation_info"]["local_gen_idx"], idx)
            end
        end

        for hybrid_data in values(network_data["hybrid"])
            hybrid_data["original_bus_i"] = hybrid_data["bus_ID"]
            hybrid_data["bus_i"] = hybrid_data["bus_ID"]
            hybrid_data["bus_idx"] = parse(Int, get(bus_idx_by_bus_id, string(hybrid_data["bus_ID"]), "0"))
            hybrid_data["region_mapping_info"] = Dict{String, Any}(
                "system" => "system",
                "System" => "system",
                "bus_i" => hybrid_data["bus_ID"],
                "bus_ID" => hybrid_data["bus_ID"],
                "bus_name" => hybrid_data["bus_name"],
            )
            hybrid_data["bypass_parameter_check"] = true
            hybrid_data["hybrid_type"] = get(hybrid_data, "hybrid_type", "NA")
            hybrid_data["hybrid_ID"] = get(hybrid_data, "hybrid_ID", 0)
        end

        for demand_data in values(network_data["demand"])
            demand_data["original_bus_i"] = demand_data["bus_ID"]
            demand_data["bus_i"] = demand_data["bus_ID"]
            demand_data["bus_idx"] = parse(Int, get(bus_idx_by_bus_id, string(demand_data["bus_ID"]), "0"))
            demand_data["region_mapping_info"] = Dict{String, Any}(
                "system" => "system",
                "System" => "system",
                "bus_i" => demand_data["bus_ID"],
                "bus_ID" => demand_data["bus_ID"],
                "bus_name" => demand_data["bus_name"],
            )
        end
        return
    end

    # Config sheets drive all aggregation, zoning, and region-level policy/resource data.
    # Open the config workbook once and read every sheet from the shared handle.
    sub_area_list = nothing
    XLSX.openxlsx(config_file_location) do config_xf
        ALEAF_setting["Network Setting"] = read_setting_xlsx_return_dict_string_any(config_xf, "Network Setting")
        sub_area_list = read_xlsx_return_dict_string_any(config_xf, "sub_area_list"; first_row_value = 2)
        ALEAF_setting["network_resolution_level"] = read_xlsx_return_dict_string_any(config_xf, "network_resolution_level"; first_row_value = 2)
        ALEAF_setting["sub_area_mapping"] = read_xlsx_return_dict_string_any(config_xf, "sub_area_mapping"; first_row_value = 2)

        additional_network_data_list = ["CEGT", "CERT", "RPS"]
        for level in keys(ALEAF_setting["network_resolution_level"])
            push!(additional_network_data_list, string("Network Data Level ", level))
        end

        for data_category in additional_network_data_list
            network_data["supplement_network_data"][data_category] = read_xlsx_return_dict_string_any(config_xf, data_category; first_row_value = 2)
        end
    end

    # Obtain network configuration settings
    # 1) network boundary type and level
    network_boundary_type = ALEAF_setting["Network Setting"]["network_boundary_type"]   # regional_network_level_value
    network_boundary_level = ALEAF_setting["Network Setting"]["network_boundary_level"]   
    network_data["supplement_network_data"]["network_boundary_type"] = network_boundary_type
    network_boundary_data = deepcopy(network_data["supplement_network_data"][string("Network Data Level ", network_boundary_level)])

    # 2) regional network aggregation resolution type and level; keep using "regional_aggregation_resolution_type" and "regional_aggregation_resolution_level" for the regional network level
    aggregation_resolution_type = ALEAF_setting["Network Setting"]["regional_aggregation_resolution_type"]   # network_reduction_level; this is the regional network level
    aggregation_resolution_level = ALEAF_setting["Network Setting"]["regional_aggregation_resolution_level"]

    # Selection-only build: coarsen the regional aggregation to the repday SELECTION resolution so the
    # bus set collapses (fewer buses => cheap build) while every finest region is still enumerated via
    # aggregation_info. The system aggregate is invariant to this grouping, so rep days are unchanged.
    selection_override = get(network_data, "__selection_agg_override", nothing)
    selection_collapse_to_one = false
    if selection_override !== nothing
        aggregation_resolution_type = selection_override[1]
        aggregation_resolution_level = selection_override[2]
        selection_collapse_to_one = length(selection_override) >= 3 && selection_override[3] === true
    end
    network_data["supplement_network_data"]["aggregation_resolution_level"] = aggregation_resolution_level
    network_aggregation_data = deepcopy(network_data["supplement_network_data"][string("Network Data Level ", aggregation_resolution_level)])

    # 3) sub-regional network aggregation resolution type and level
    subregional_aggregation_type = ALEAF_setting["Network Setting"]["subregional_aggregation_resolution_type"]   # sub-regional network level
    subregional_aggregation_level = ALEAF_setting["Network Setting"]["subregional_aggregation_resolution_level"]
    if selection_override !== nothing
        # No sub-regional split for the selection build: keep one flat coarse grouping.
        subregional_aggregation_type = aggregation_resolution_type
        subregional_aggregation_level = aggregation_resolution_level
    end
    network_data["supplement_network_data"]["subregional_aggregation_level"] = subregional_aggregation_level
    subregional_network_aggregation_data = deepcopy(network_data["supplement_network_data"][string("Network Data Level ", subregional_aggregation_level)])

    # Obtain selected regions based on the network boundary type and data
    selected_regions_within_boundary_list = []
    for id in keys(network_boundary_data)
        if network_boundary_data[id]["Model"] == true
            push!(selected_regions_within_boundary_list, network_boundary_data[id]["Region_ID"])
        end
    end

    # Obtain sub-regions within the selected regions
    selected_subregions_within_boundary_list = []
    selected_outerregions_within_boundary_list = []
    # Check if the aggregation resolution type is valid
    if aggregation_resolution_type == subregional_aggregation_type
        
        for id in keys(network_aggregation_data)
            if network_aggregation_data[id]["Model"] == true 
                push!(selected_outerregions_within_boundary_list, network_aggregation_data[id]["Region_ID"])
            end
        end

    else

        @aleaf_info "Sub-regional aggregation will be be applied."

        for id in keys(network_aggregation_data)
            if network_aggregation_data[id]["Model"] == true 
                if network_aggregation_data[id]["Subregional_Aggregation"] == true
                    push!(selected_subregions_within_boundary_list, network_aggregation_data[id]["Region_ID"])
                else
                    push!(selected_outerregions_within_boundary_list, network_aggregation_data[id]["Region_ID"])
                end
            end
        end

    end

    # Update bus id
    for n in keys(network_data["bus"])
        original_bus_i = deepcopy(network_data["bus"][n]["bus_i"])
        network_data["bus"][n]["original_bus_i"] = string(original_bus_i)
        network_data["bus"][n]["bus_i"] = find_new_bus_i(ALEAF_setting["sub_area_mapping"], original_bus_i, aggregation_resolution_type)
    end

    # clean-up the sub_area_list dict (the sub_area_list will be converted to a new bus dictionary)
    for n in keys(sub_area_list)
        if sub_area_list[n][aggregation_resolution_type] === missing
            delete!(sub_area_list, n)
        else
            for id in keys(sub_area_list[n])
                if id != aggregation_resolution_type
                    delete!(sub_area_list[n], id)
                end
            end
        end
    end

    # create a new bus list including all sub-regions within the selected regions
    new_bus_list = Dict{String, Any}()
    bus_idx = 1

    for id in keys(sub_area_list)
       
        sub_area_id = deepcopy(sub_area_list[id][aggregation_resolution_type])
        
        # check if the sub_area_id is within the selected regions
        if sub_area_id in selected_subregions_within_boundary_list

            sub_region_list= find_sub_region_mapping_info(ALEAF_setting["sub_area_mapping"], sub_area_id, aggregation_resolution_type, subregional_aggregation_type)

            for sub_region_idx in sub_region_list
                new_bus_list[string(bus_idx)] = Dict{String, Any}()
                new_bus_list[string(bus_idx)]["region_config"] = Dict{String, Any}(
                    "region_type" => "sub_regional",
                    "aggregation_resolution_type" => aggregation_resolution_type,
                    "aggregation_resolution_value" => sub_area_id,
                    "subregional_aggregation_type" => subregional_aggregation_type,
                    "subregional_aggregation_value" => sub_region_idx,
                    "region_lookup_type" => "subregional_aggregation_type",
                    "region_lookup_value" => "subregional_aggregation_value",
                    "parent_region_id" => sub_area_id
                )
                bus_idx += 1
            end
            
        elseif sub_area_id in selected_outerregions_within_boundary_list
            new_bus_list[string(bus_idx)] = Dict{String, Any}()
            new_bus_list[string(bus_idx)]["region_config"] = Dict{String, Any}(
                "region_type" => "regional",
                "aggregation_resolution_type" => aggregation_resolution_type,
                "aggregation_resolution_value" => sub_area_id,
                "subregional_aggregation_type" => "NA",
                "subregional_aggregation_value" => "NA",
                "region_lookup_type" => "aggregation_resolution_type",
                "region_lookup_value" => "aggregation_resolution_value",
                "parent_region_id" => sub_area_id
            )
            bus_idx += 1
        end
    end

    # remove invalid regions and add region mapping info
    for n in keys(new_bus_list)

        lookup_value = new_bus_list[n]["region_config"]["region_lookup_value"]
        lookup_type = new_bus_list[n]["region_config"]["region_lookup_type"]
        sub_area_id = deepcopy(new_bus_list[n]["region_config"][lookup_value])

        # Find the region mapping info based on the sub_area_id and lookup type
        region_info = find_region_mapping_info(ALEAF_setting["sub_area_mapping"], sub_area_id, new_bus_list[n]["region_config"][lookup_type])

        if region_info === nothing
            delete!(new_bus_list, n)
        elseif !(region_info[network_boundary_type] in selected_regions_within_boundary_list)
            delete!(new_bus_list, n)
        else
            new_bus_list[n]["region_mapping_info"] = region_info
        end
    end

    # Sort sub_area_list and assign to a new dictionary with sequential keys
    sorted_new_bus_listt = Dict{String, Any}(string(i) => deepcopy(value) for (i, (key, value)) in enumerate(new_bus_list))
    new_bus_list = deepcopy(sorted_new_bus_listt)
    sorted_new_bus_listt = nothing

    # Aggregate network data by incorporating the most detailed network information
    highest_network_resolution_level_ID = ALEAF_setting["network_resolution_level"]["1"]["Network_Resolution_ID"]

    # finest bus original_id -> its (finest, pre-swap) bus dict key, for the 2171-2178 inlined scan.
    orig2idx = Dict{String, String}()
    for bus_idx in keys(network_data["bus"])
        get!(orig2idx, string(network_data["bus"][bus_idx]["original_bus_i"]), bus_idx)
    end

    for n in keys(new_bus_list)

        new_bus_list[n]["aggregation_info"] = Dict{String,Any}()
        
        aggregated_regions = []
        aggregated_regions_idx = []
        aggregated_regions_bus_name = []

        lookup_value = new_bus_list[n]["region_config"]["region_lookup_value"]
        lookup_type = new_bus_list[n]["region_config"]["region_lookup_type"]
        region_lookup_level = new_bus_list[n]["region_config"][lookup_type]
        sub_area_id = deepcopy(new_bus_list[n]["region_config"][lookup_value])
                
        # add aggregated regions info
        for member in get(_get_members_index(ALEAF_setting["sub_area_mapping"], region_lookup_level, highest_network_resolution_level_ID), sub_area_id, Any[])
            push!(aggregated_regions, string(member))
        end
        new_bus_list[n]["aggregation_info"]["aggregated_regions_bus_i"] = unique(aggregated_regions)

        # add idx of aggregated regions
        for agg_bus_id in new_bus_list[n]["aggregation_info"]["aggregated_regions_bus_i"]
            bus_idx = get(orig2idx, agg_bus_id, nothing)
            if bus_idx !== nothing
                push!(aggregated_regions_idx, bus_idx)
                push!(aggregated_regions_bus_name, (bus_idx, network_data["bus"][bus_idx]["bus_name"]))
            end
        end
        new_bus_list[n]["aggregation_info"]["aggregated_regions_bus_idx"] = aggregated_regions_idx
        new_bus_list[n]["aggregation_info"]["aggregated_regions_bus_(idx, name)"] = aggregated_regions_bus_name

        # prepare list of generators located in this sub_area
        new_bus_list[n]["aggregation_info"]["local_gen_idx"] = []
        new_bus_list[n]["aggregation_info"]["new_local_gen_idx"] = []  # this will be updated in the "add_ref" function
    
    end

    # bus peak load update and add bus_i
    for n in keys(new_bus_list)
        total_peak_load = 0
        original_load_MW = Dict()

        for nn in new_bus_list[n]["aggregation_info"]["aggregated_regions_bus_idx"]
            total_peak_load += network_data["bus"][nn]["MW load"]
            original_load_MW[network_data["bus"][nn]["original_bus_i"]] = network_data["bus"][nn]["MW load"]
        end
        
        # add peak load MW 
        new_bus_list[n]["MW load"] = total_peak_load
        new_bus_list[n]["aggregation_info"]["original_load_(bus_i, MW)"] = original_load_MW

        # add bus_i
        lookup_value = new_bus_list[n]["region_config"]["region_lookup_value"]
        new_bus_list[n]["bus_i"] = new_bus_list[n]["region_config"][lookup_value]
        # delete!(new_bus_list[n], aggregation_resolution_type)
    end

    # Region_ID -> entry indexes (LAST wins, matching the original non-breaking scans that keep the
    # final match). subregional and regional (network_aggregation_data) kept separate.
    nad_by_region = Dict{Any, Any}()
    for key in keys(network_aggregation_data)
        nad_by_region[network_aggregation_data[key]["Region_ID"]] = network_aggregation_data[key]
    end
    subnad_by_region = Dict{Any, Any}()
    for key in keys(subregional_network_aggregation_data)
        subnad_by_region[subregional_network_aggregation_data[key]["Region_ID"]] = subregional_network_aggregation_data[key]
    end

    # add RA ELCC flag to bus data & Region Name
    for n in keys(new_bus_list)

        if new_bus_list[n]["region_config"]["region_type"] == "sub_regional"
            entry = get(subnad_by_region, new_bus_list[n]["bus_i"], nothing)
            entry === nothing || (new_bus_list[n]["RA_ELCC_Calculation_Flag"] = entry["RA_ELCC_Calculation_Flag"])
        else
            entry = get(nad_by_region, new_bus_list[n]["bus_i"], nothing)
            entry === nothing || (new_bus_list[n]["RA_ELCC_Calculation_Flag"] = entry["RA_ELCC_Calculation_Flag"])
        end
    end

    # add Region Name based on network_boundary_type
    for n in keys(new_bus_list)

        new_bus_list[n]["region_config"]["region_name"] = new_bus_list[n]["region_mapping_info"][network_boundary_type]

        if new_bus_list[n]["region_config"]["region_type"] == "sub_regional"
            entry = get(nad_by_region, new_bus_list[n]["region_config"]["parent_region_id"], nothing)
            entry === nothing || (new_bus_list[n]["region_config"]["parent_bus_name"] = entry["Region_Name"])
        else
            entry = get(nad_by_region, new_bus_list[n]["bus_i"], nothing)
            entry === nothing || (new_bus_list[n]["region_config"]["parent_bus_name"] = entry["Region_Name"])
        end

    end

    # "system" selection: merge the modeled boundary buses into ONE, unioning their finest-region
    # membership. The selection metric iterates finest regions (original_load_(bus_i, MW) keys), so the
    # merged set is identical to the per-boundary buses => byte-identical rep days, fewer buses built.
    # Merging happens AFTER the modeled-region filter, so no unmodeled region can be pulled in (no straddle).
    if selection_collapse_to_one && length(new_bus_list) > 1
        merged = deepcopy(new_bus_list[first(keys(new_bus_list))])
        merged_bus_i = String[]
        merged_bus_idx = String[]
        merged_bus_name = Tuple{String,String}[]
        merged_load = Dict{Any,Any}()
        merged_peak = 0.0
        for n in keys(new_bus_list)
            info = new_bus_list[n]["aggregation_info"]
            append!(merged_bus_i, info["aggregated_regions_bus_i"])
            append!(merged_bus_idx, info["aggregated_regions_bus_idx"])
            append!(merged_bus_name, info["aggregated_regions_bus_(idx, name)"])
            merge!(merged_load, info["original_load_(bus_i, MW)"])
            merged_peak += new_bus_list[n]["MW load"]
        end
        merged["aggregation_info"]["aggregated_regions_bus_i"] = unique(merged_bus_i)
        merged["aggregation_info"]["aggregated_regions_bus_idx"] = unique(merged_bus_idx)
        merged["aggregation_info"]["aggregated_regions_bus_(idx, name)"] = unique(merged_bus_name)
        merged["aggregation_info"]["original_load_(bus_i, MW)"] = merged_load
        merged["aggregation_info"]["local_gen_idx"] = []
        merged["aggregation_info"]["new_local_gen_idx"] = []
        merged["MW load"] = merged_peak
        new_bus_list = Dict{String, Any}("1" => merged)
    end

    # replace bus dict
    delete!(network_data, "bus")
    network_data["bus"] = new_bus_list

    # member finest-region -> aggregated bus key (first wins), powering find_new_bus_i_using_agg_resolution.
    for busidx in keys(network_data["bus"])
        for member in network_data["bus"][busidx]["aggregation_info"]["aggregated_regions_bus_i"]
            get!(_member2aggbus, string(member), busidx)
        end
    end

    ###
    # Update plant data
    ###

    # Remove unlocated plants 
    for g in collect(keys(network_data["plant"]))
        original_bus_i = deepcopy(network_data["plant"][g]["bus_ID"])
        new_bus_info = find_new_bus_i_using_agg_resolution(string(original_bus_i))  # return bus_i and bus_idx
        if new_bus_info[1] == 0
           delete!(network_data["plant"], g)
        end
    end

    # Reassign plant IDs 
    old_plants = values(network_data["plant"])
    new_plant_dict = Dict{String, Any}()  
    i = 1
    for plant in old_plants
        new_plant_dict[string(i)] = plant
        i += 1
    end
    network_data["plant"] = new_plant_dict

    # Update geneator bus id and add region mapping info
    for g in collect(keys(network_data["plant"]))
        original_bus_i = deepcopy(network_data["plant"][g]["bus_ID"])

        new_bus_info = find_new_bus_i_using_agg_resolution(string(original_bus_i))  # return bus_i and bus_idx

        network_data["plant"][g]["original_bus_i"] = original_bus_i
        network_data["plant"][g]["bus_i"] = new_bus_info[1]
        network_data["plant"][g]["bus_idx"] = new_bus_info[2]

        network_data["plant"][g]["region_mapping_info"] = deepcopy(new_bus_info[3])  # deepcopy is used to prevent the original "region_mapping_info" data from being updated
        network_data["plant"][g]["region_mapping_info"]["bus_ID"] = network_data["plant"][g]["bus_ID"]   # update the mapping info to prevent confusion with the original data
        network_data["plant"][g]["region_mapping_info"]["bus_name"] = network_data["plant"][g]["bus_name"]    # update the mapping info to prevent confusion with the original data

        # add hybrid parameters
        network_data["plant"][g]["bypass_parameter_check"] = false # false => perform check_and_update_plant_data!
        network_data["plant"][g]["hybrid_type"] = "NA"
        network_data["plant"][g]["hybrid_ID"] = 0
    end

    ###
    # Update hybrid plant data 
    ###

    # Remove unlocated plants 
    for g in collect(keys(network_data["hybrid"]))
        original_bus_i = deepcopy(network_data["hybrid"][g]["bus_ID"])
        new_bus_info = find_new_bus_i_using_agg_resolution(string(original_bus_i))  # return bus_i and bus_idx
        if new_bus_info[1] == 0
           delete!(network_data["hybrid"], g)
        end
    end
    
    # Reassign plant IDs 
    old_hybrid_plant = values(network_data["hybrid"])
    new_hybrid_plant_dict = Dict{String, Any}()  
    i = 1
    for plant in old_hybrid_plant
        new_hybrid_plant_dict[string(i)] = plant
        i += 1
    end
    network_data["hybrid"] = new_hybrid_plant_dict

    # Update geneator bus id and add region mapping info
    for g in collect(keys(network_data["hybrid"]))
        original_bus_i = deepcopy(network_data["hybrid"][g]["bus_ID"])

        new_bus_info = find_new_bus_i_using_agg_resolution(string(original_bus_i))  # return bus_i and bus_idx

        network_data["hybrid"][g]["original_bus_i"] = original_bus_i
        network_data["hybrid"][g]["bus_i"] = new_bus_info[1]
        network_data["hybrid"][g]["bus_idx"] = new_bus_info[2]

        network_data["hybrid"][g]["region_mapping_info"] = deepcopy(new_bus_info[3])  # deepcopy is used to prevent the original "region_mapping_info" data from being updated
        network_data["hybrid"][g]["region_mapping_info"]["bus_ID"] = network_data["hybrid"][g]["bus_ID"]   # update the mapping info to prevent confusion with the original data
        network_data["hybrid"][g]["region_mapping_info"]["bus_name"] = network_data["hybrid"][g]["bus_name"]    # update the mapping info to prevent confusion with the original data

        # add hybrid parameters
        network_data["hybrid"][g]["bypass_parameter_check"] = true # false => perform check_and_update_plant_data!
        network_data["hybrid"][g]["hybrid_type"] = "NA" # will be updated in the update_hybrid_plant_technology_data! function
        network_data["hybrid"][g]["hybrid_ID"] = 0  # will be updated in the update_hybrid_plant_technology_data! function
    end

    ###
    # Update demand data (large load, demand response)
    ###

    for g in collect(keys(network_data["demand"]))
        original_bus_i = deepcopy(network_data["demand"][g]["bus_ID"])
        new_bus_info = find_new_bus_i_using_agg_resolution(string(original_bus_i))  # return bus_i and bus_idx
        if new_bus_info[1] == 0
           delete!(network_data["demand"], g)
        end
    end

    # Reassign plant IDs 
    old_demand = values(network_data["demand"])
    new_demand_dict = Dict{String, Any}()  
    i = 1
    for plant in old_demand
        new_demand_dict[string(i)] = plant
        i += 1
    end
    network_data["demand"] = new_demand_dict

    # Update geneator bus id and add region mapping info
    for g in collect(keys(network_data["demand"]))
        original_bus_i = deepcopy(network_data["demand"][g]["bus_ID"])

        new_bus_info = find_new_bus_i_using_agg_resolution(string(original_bus_i))  # return bus_i and bus_idx

        network_data["demand"][g]["original_bus_i"] = original_bus_i
        network_data["demand"][g]["bus_i"] = new_bus_info[1]
        network_data["demand"][g]["bus_idx"] = new_bus_info[2]

        network_data["demand"][g]["region_mapping_info"] = deepcopy(new_bus_info[3])  # deepcopy is used to prevent the original "region_mapping_info" data from being updated
        network_data["demand"][g]["region_mapping_info"]["bus_ID"] = network_data["demand"][g]["bus_ID"]   # update the mapping info to prevent confusion with the original data
        network_data["demand"][g]["region_mapping_info"]["bus_name"] = network_data["demand"][g]["bus_name"]    # update the mapping info to prevent confusion with the original data
    end
    
    ###
    # Update branch data
    ###

    if length(network_data["bus"]) > 1
        try
            aggregate_branches_using_database!(ALEAF_setting, network_data)
        catch
            # a partial failure leaves region-remapped but unmerged branches; never run on that network
            @aleaf_warn "Failed generating transmission network!"
            rethrow()
        end
    else
        network_data["branch"] = Dict()
    end

end


function define_zone!(network_data::Dict{String, <:Any}, ALEAF_setting::Dict{String, <:Any}, zone_type::String, zone_boundary_type::String, data_resolution_type::String, data_resolution_level::Int)

    _firstrow_by_col = Dict{String, Dict{Any, Any}}()
    function find_region_mapping_info(sub_area_mapping, original_id, reduction_level)
        idx = get!(_firstrow_by_col, reduction_level) do
            d = Dict{Any, Any}()
            for key in keys(sub_area_mapping)
                v = sub_area_mapping[key][reduction_level]
                v === missing && continue
                get!(d, v, sub_area_mapping[key])
            end
            d
        end
        return get(idx, original_id, nothing)
    end

    function find_data_region_info(sub_area_mapping, zone_id, zone_boundary_type, data_resolution_type)
        sub_area_list = []
        for key in keys(sub_area_mapping)
            if sub_area_mapping[key][zone_boundary_type] == zone_id
                push!(sub_area_list, string(sub_area_mapping[key][data_resolution_type]))
            end
        end
        return unique(sub_area_list)
    end
    
    network_data["zone"][zone_type] = Dict{String, Any}()

    # list of sub-areas based on the reserve zone boundary level & remove duplicates
    zone_list = []
    for id in keys(network_data["bus"])
        push!(zone_list, network_data["bus"][id]["region_mapping_info"][zone_boundary_type])
    end
    zone_list = unique(zone_list)  

    # list of sub-areas based on data resolution
    zone_data_area_list = []
    sub_area_list_dict = read_xlsx_return_dict_string_any(ALEAF_setting["network_config_file_location"], "sub_area_list")
    for n in keys(sub_area_list_dict)
        if sub_area_list_dict[n][data_resolution_type] !== missing
            push!(zone_data_area_list, sub_area_list_dict[n][data_resolution_type])
        end
    end
    zone_data_area_list = unique(zone_data_area_list)  

    # update zone structure
    for n in eachindex(zone_list)
        network_data["zone"][zone_type][string(n)] = Dict{String, Any}()

        # zone identifier 
        network_data["zone"][zone_type][string(n)]["zone_id"] = zone_list[n]

        # zone data resolution level
        network_data["zone"][zone_type][string(n)]["zone_data_resolution_type"] = data_resolution_type
        network_data["zone"][zone_type][string(n)]["zone_data_resolution_level"] = data_resolution_level

        # Zone regional mapping info 
        network_data["zone"][zone_type][string(n)]["zone_regional_mapping_info"] = find_region_mapping_info(ALEAF_setting["sub_area_mapping"], zone_list[n], zone_boundary_type)

        # update reserve data list  & remove duplicates
        network_data["zone"][zone_type][string(n)]["data_area_list"] = []
        for id in keys(network_data["bus"])
            if network_data["bus"][id]["region_mapping_info"][zone_boundary_type] == zone_list[n]
                push!(network_data["zone"][zone_type][string(n)]["data_area_list"], network_data["bus"][id]["region_mapping_info"][data_resolution_type])
            end
        end
        network_data["zone"][zone_type][string(n)]["data_area_list"] = unique(network_data["zone"][zone_type][string(n)]["data_area_list"])

        # add aggregated bus info
        network_data["zone"][zone_type][string(n)]["aggregation_info"] = Dict{String, Any}()
        network_data["zone"][zone_type][string(n)]["aggregation_info"]["aggregated_regions_bus_i"] = []  
        network_data["zone"][zone_type][string(n)]["aggregation_info"]["aggregated_regions_bus_idx"] = []  
        network_data["zone"][zone_type][string(n)]["aggregation_info"]["zone_bus_idx"] = []  
        network_data["zone"][zone_type][string(n)]["aggregation_info"]["local_gen_idx"] = []  
        network_data["zone"][zone_type][string(n)]["aggregation_info"]["new_local_gen_idx"] = []    # this will be updated in the "add_ref" function
        if zone_type == "reserve"
            network_data["zone"][zone_type][string(n)]["aggregation_info"]["technology_list"] = String[]
            network_data["zone"][zone_type][string(n)]["aggregation_info"]["gen_technology_groups"] = Dict{String, Any}()
        end

        # Concatenate the list
        for key in keys(network_data["bus"])
            if network_data["bus"][key]["region_mapping_info"][zone_boundary_type] == zone_list[n]
                
                push!(network_data["zone"][zone_type][string(n)]["aggregation_info"]["zone_bus_idx"], parse(Int64, key))  
                
                append!(network_data["zone"][zone_type][string(n)]["aggregation_info"]["aggregated_regions_bus_i"], 
                        network_data["bus"][key]["aggregation_info"]["aggregated_regions_bus_i"])  

                append!(network_data["zone"][zone_type][string(n)]["aggregation_info"]["aggregated_regions_bus_idx"], 
                        network_data["bus"][key]["aggregation_info"]["aggregated_regions_bus_idx"])  

                append!(network_data["zone"][zone_type][string(n)]["aggregation_info"]["local_gen_idx"], 
                        network_data["bus"][key]["aggregation_info"]["local_gen_idx"])  

            end
        end

        network_data["zone"][zone_type][string(n)]["aggregation_info"]["aggregated_data_regions_id"] = find_data_region_info(ALEAF_setting["sub_area_mapping"], zone_list[n], zone_boundary_type, data_resolution_type)
        
    end

end


function define_zone_wo_data!(network_data::Dict{String, <:Any}, ALEAF_setting::Dict{String, <:Any}, zone_type::String, zone_boundary_type::String)

    _firstrow_by_col = Dict{String, Dict{Any, Any}}()
    function find_region_mapping_info(sub_area_mapping, original_id, reduction_level)
        idx = get!(_firstrow_by_col, reduction_level) do
            d = Dict{Any, Any}()
            for key in keys(sub_area_mapping)
                v = sub_area_mapping[key][reduction_level]
                v === missing && continue
                get!(d, v, sub_area_mapping[key])
            end
            d
        end
        return get(idx, original_id, nothing)
    end

    network_data["zone"][zone_type] = Dict{String, Any}()

    # list of sub-areas based on the reserve zone boundary level & remove duplicates
    zone_list = []
    for id in keys(network_data["bus"])
        push!(zone_list, network_data["bus"][id]["region_mapping_info"][zone_boundary_type])
    end
    zone_list = unique(zone_list)  

    # update zone structure
    for n in eachindex(zone_list)
        network_data["zone"][zone_type][string(n)] = Dict{String, Any}()

        # zone identifier 
        network_data["zone"][zone_type][string(n)]["zone_id"] = zone_list[n]

        # Zone regional mapping info 
        network_data["zone"][zone_type][string(n)]["zone_regional_mapping_info"] = find_region_mapping_info(ALEAF_setting["sub_area_mapping"], zone_list[n], zone_boundary_type)

        # add aggregated bus info
        network_data["zone"][zone_type][string(n)]["aggregation_info"] = Dict{String, Any}()
        network_data["zone"][zone_type][string(n)]["aggregation_info"]["aggregated_regions_bus_i"] = []  
        network_data["zone"][zone_type][string(n)]["aggregation_info"]["aggregated_regions_bus_idx"] = []  
        network_data["zone"][zone_type][string(n)]["aggregation_info"]["zone_bus_idx"] = []  
        network_data["zone"][zone_type][string(n)]["aggregation_info"]["local_gen_idx"] = []  
        network_data["zone"][zone_type][string(n)]["aggregation_info"]["new_local_gen_idx"] = []    # this will be updated in the "add_ref" function

        # Concatenate the list
        for key in keys(network_data["bus"])
            if network_data["bus"][key]["region_mapping_info"][zone_boundary_type] == zone_list[n]
                
                push!(network_data["zone"][zone_type][string(n)]["aggregation_info"]["zone_bus_idx"], parse(Int64, key))  
                
                append!(network_data["zone"][zone_type][string(n)]["aggregation_info"]["aggregated_regions_bus_i"], 
                        network_data["bus"][key]["aggregation_info"]["aggregated_regions_bus_i"])  

                append!(network_data["zone"][zone_type][string(n)]["aggregation_info"]["aggregated_regions_bus_idx"], 
                        network_data["bus"][key]["aggregation_info"]["aggregated_regions_bus_idx"])  

                append!(network_data["zone"][zone_type][string(n)]["aggregation_info"]["local_gen_idx"], 
                        network_data["bus"][key]["aggregation_info"]["local_gen_idx"])  

            end
        end
        
    end

end


function get_policy_targets!(network_data::Dict{String, <:Any}, ALEAF_setting::Dict{String, <:Any}, case_id, target_policy, global_target_value)

    function get_policy_target(policy_label::String, region_id, year_id)
        for key in keys(network_data["supplement_network_data"][policy_label])
            if network_data["supplement_network_data"][policy_label][key]["Region"] == region_id
                return network_data["supplement_network_data"][policy_label][key][string(year_id)]
            end
        end        
    end

    for n in keys(network_data["zone"]["policy"])
        network_data["zone"]["policy"][n][target_policy] = Dict{String, Any}()

        network_data["zone"]["policy"][n][target_policy]["target"] = Dict()
        
        for y in keys(network_data["planning_stages"])
            year_id = network_data["planning_stages"][y]["year"]
            
            if global_target_value > 0
                network_data["zone"]["policy"][n][target_policy]["target"][y] = global_target_value
            else
                target_value = 0.0
                for data_region_id in network_data["zone"]["policy"][n]["data_area_list"]   # we should have only one data region for RPS
                    target_value += get_policy_target(target_policy, data_region_id, year_id)
                end
                network_data["zone"]["policy"][n][target_policy]["target"][y] = target_value
            end
        end

        # add reference emission level for CERT (carbon emission reduction target)
        if target_policy == "CERT"
            network_data["zone"]["policy"][n][target_policy]["reference"] = 0.0
            for data_region_id in network_data["zone"]["policy"][n]["data_area_list"]   
                network_data["zone"]["policy"][n][target_policy]["reference"] += get_policy_target(target_policy, data_region_id, "Ref_Emission_m_ton")
            end
        end

    end
end


# Function to find the value from the data dictionary
function get_resource_supply_curve!(network_data::Dict{String, <:Any}, ALEAF_setting::Dict{String, <:Any}, case_id)

    _firstrow_by_col = Dict{String, Dict{Any, Any}}()
    function find_region_mapping_info(sub_area_mapping, original_id, reduction_level)
        idx = get!(_firstrow_by_col, reduction_level) do
            d = Dict{Any, Any}()
            for key in keys(sub_area_mapping)
                v = sub_area_mapping[key][reduction_level]
                v === missing && continue
                get!(d, v, sub_area_mapping[key])
            end
            d
        end
        return get(idx, original_id, nothing)
    end
    function find_data_region_info(sub_area_mapping, zone_id, zone_boundary_type, data_resolution_type)
        sub_area_list = []
        for key in keys(sub_area_mapping)
            if sub_area_mapping[key][zone_boundary_type] == zone_id
                push!(sub_area_list, string(sub_area_mapping[key][data_resolution_type]))
            end
        end
        return unique(sub_area_list)
    end
    
    # Get the list of supply curve resources
    supply_curve_resource_list = []
    for gen_tech_id in keys(ALEAF_setting["Gen Technology"])
        if ALEAF_setting["Gen Technology"][gen_tech_id]["Resource_Limit_Flag"] == true
            resource_id = ALEAF_setting["Gen Technology"][gen_tech_id]["Resource_Limit_ID"]
            if !(resource_id in supply_curve_resource_list)
                push!(supply_curve_resource_list, resource_id)
            end
        end
    end
    network_data["zone"]["resource_supply_curve"] = Dict{String, Any}()
    regional_resource_supply_curve_resolution_type = ALEAF_setting["Network Setting"]["regional_resource_supply_curve_resolution_type"]   # regional resource supply curve resolution type
    regional_resource_supply_curve_resolution_level = ALEAF_setting["Network Setting"]["regional_resource_supply_curve_resolution_level"]   # regional resource supply curve resolution level
    data_dict = network_data["supplement_network_data"][string("Network Data Level ", regional_resource_supply_curve_resolution_level)]
    network_data["zone"]["resource_supply_curve"]["data_resolution_type"] = data_resolution_type = ALEAF_setting["Network Setting"]["regional_resource_supply_curve_resolution_type"]
    network_data["zone"]["resource_supply_curve"]["supply_curve_resource_list"] = supply_curve_resource_list
    network_data["zone"]["resource_supply_curve"]["supply_curve_region_list"] = []
    
    # list of sub-areas based on the reserve zone boundary level & remove duplicates
    zone_list = []
    for id in keys(network_data["bus"])
        push!(zone_list, network_data["bus"][id]["region_mapping_info"][regional_resource_supply_curve_resolution_type])
    end
    zone_list = unique(zone_list)  
    # list of sub-areas based on data resolution
    sub_area_list_dict = read_xlsx_return_dict_string_any(ALEAF_setting["network_config_file_location"], "sub_area_list")
    zone_data_area_list = []
    for n in keys(sub_area_list_dict)
        if sub_area_list_dict[n][data_resolution_type] !== missing
            push!(zone_data_area_list, sub_area_list_dict[n][data_resolution_type])
        end
    end
    zone_data_area_list = unique(zone_data_area_list)  
    
    # update zone structure
    for n in eachindex(zone_list)
        network_data["zone"]["resource_supply_curve"][string(n)] = Dict{String, Any}()
        # zone identifier 
        zone_id = zone_list[n]
        network_data["zone"]["resource_supply_curve"][string(n)]["zone_id"] = zone_id
        # zone data resolution level
        network_data["zone"]["resource_supply_curve"][string(n)]["zone_data_resolution_type"] = regional_resource_supply_curve_resolution_type
        network_data["zone"]["resource_supply_curve"][string(n)]["zone_data_resolution_level"] = regional_resource_supply_curve_resolution_level
        # Zone regional mapping info 
        network_data["zone"]["resource_supply_curve"][string(n)]["zone_regional_mapping_info"] = find_region_mapping_info(ALEAF_setting["sub_area_mapping"], zone_list[n], regional_resource_supply_curve_resolution_type)
        # Add data 
        network_data["zone"]["resource_supply_curve"][string(n)]["data"] = Dict()
        for tech_id in supply_curve_resource_list
            
            network_data["zone"]["resource_supply_curve"][string(n)]["data"][tech_id] = Dict()
            cap_level = "highcap"
            if ALEAF_setting["Simulation Configuration"][string(case_id)]["Resource_limit_level_value"] == "Low"
                cap_level = "lowcap"
            end
            # tag 
            tag = string("Resource_Limit_", tech_id, "_", cap_level)
            # data idx
            data_keys = [k for (k, v) in data_dict if v["Region_ID"] == zone_id]
            network_data["zone"]["resource_supply_curve"][string(n)]["data"][tech_id]["Resource_Limit_MW"] = data_dict[data_keys[1]][tag]
        end
        
        # update reserve data list  & remove duplicates
        network_data["zone"]["resource_supply_curve"][string(n)]["data_area_list"] = []
        for id in keys(network_data["bus"])
            if network_data["bus"][id]["region_mapping_info"][regional_resource_supply_curve_resolution_type] == zone_list[n]
                push!(network_data["zone"]["resource_supply_curve"][string(n)]["data_area_list"], network_data["bus"][id]["region_mapping_info"][regional_resource_supply_curve_resolution_type])
            end
        end
        network_data["zone"]["resource_supply_curve"][string(n)]["data_area_list"] = unique(network_data["zone"]["resource_supply_curve"][string(n)]["data_area_list"])
        # add aggregated bus info
        network_data["zone"]["resource_supply_curve"][string(n)]["aggregation_info"] = Dict{String, Any}()
        network_data["zone"]["resource_supply_curve"][string(n)]["aggregation_info"]["aggregated_regions_bus_i"] = []  
        network_data["zone"]["resource_supply_curve"][string(n)]["aggregation_info"]["aggregated_regions_bus_idx"] = []  
        network_data["zone"]["resource_supply_curve"][string(n)]["aggregation_info"]["zone_bus_idx"] = []  
        network_data["zone"]["resource_supply_curve"][string(n)]["aggregation_info"]["local_gen_idx"] = []    
        network_data["zone"]["resource_supply_curve"][string(n)]["aggregation_info"]["new_local_gen_idx"] = []    # this will be updated in the "add_ref" function
        # Concatenate the list
        for key in keys(network_data["bus"])
            if network_data["bus"][key]["region_mapping_info"][regional_resource_supply_curve_resolution_type] == zone_list[n]
                
                push!(network_data["zone"]["resource_supply_curve"][string(n)]["aggregation_info"]["zone_bus_idx"], parse(Int64, key))  
                
                append!(network_data["zone"]["resource_supply_curve"][string(n)]["aggregation_info"]["aggregated_regions_bus_i"], 
                        network_data["bus"][key]["aggregation_info"]["aggregated_regions_bus_i"])  
                append!(network_data["zone"]["resource_supply_curve"][string(n)]["aggregation_info"]["aggregated_regions_bus_idx"], 
                        network_data["bus"][key]["aggregation_info"]["aggregated_regions_bus_idx"])  
                append!(network_data["zone"]["resource_supply_curve"][string(n)]["aggregation_info"]["local_gen_idx"],
                        network_data["bus"][key]["aggregation_info"]["local_gen_idx"])
            end
        end
        network_data["zone"]["resource_supply_curve"][string(n)]["aggregation_info"]["aggregated_data_regions_id"] = find_data_region_info(ALEAF_setting["sub_area_mapping"], zone_list[n], regional_resource_supply_curve_resolution_type, data_resolution_type)
        push!(network_data["zone"]["resource_supply_curve"]["supply_curve_region_list"], n)
    end
    
end


function get_resource_cost_variations!(network_data::Dict{String, <:Any}, ALEAF_setting::Dict{String, <:Any}, case_id)
    
    # Get the list of supply curve resources
    resource_list = []
    for gen_tech_id in keys(ALEAF_setting["Gen Technology"])
        if ALEAF_setting["Gen Technology"][gen_tech_id]["Locational_Scaling_Flag"] == true
            if !(gen_tech_id in resource_list)
                push!(resource_list, gen_tech_id)
            end
        end
    end
        
    network_data["zone"]["cost_scaling"] = Dict{String, Any}()
    regional_resource_cost_scaling_resolution_level = ALEAF_setting["Network Setting"]["regional_resource_cost_scaling_resolution_level"]   # regional resource cost variation resolution level
    data_dict = network_data["supplement_network_data"][string("Network Data Level ", regional_resource_cost_scaling_resolution_level)]
    network_data["zone"]["cost_scaling"]["data_identifier"] = data_identifier = ALEAF_setting["Network Setting"]["regional_resource_cost_scaling_resolution_type"]

    # Update resource cost variation data; We need to define resource cost variation for each bus (average value)
    for n in keys(data_dict)

        region_ID = data_dict[n]["Region_ID"]
        region_ID = region_ID isa AbstractString ? region_ID : string(region_ID)

        network_data["zone"]["cost_scaling"][region_ID] = Dict()
        
        # gen technology
        for gen_tech_id in resource_list

            string_id = ALEAF_setting["Gen Technology"][gen_tech_id]["Tech_ID"]
            unit_group = ALEAF_setting["Gen Technology"][gen_tech_id]["UNITGROUP"]
            Resource_Limit_ID = ALEAF_setting["Gen Technology"][gen_tech_id]["Resource_Limit_ID"]

            network_data["zone"]["cost_scaling"][region_ID][string_id] = Dict()

            network_data["zone"]["cost_scaling"][region_ID][string_id]["UNITGROUP"] = unit_group

            network_data["zone"]["cost_scaling"][region_ID][string_id]["CAPAX_scale"] = (1 + data_dict[n][string("CAPAX_scale_", Resource_Limit_ID)])
        end
    end

end


function get_resource_capacity_credits!(network_data::Dict{String, <:Any}, ALEAF_setting::Dict{String, <:Any}, case_id)
    
    # Get the list of supply curve resources
    resource_list = []
    for gen_tech_id in keys(ALEAF_setting["Gen Technology"])
        if ALEAF_setting["Gen Technology"][gen_tech_id]["CAPCRED"] isa String
            if !(ALEAF_setting["Gen Technology"][gen_tech_id]["UNITGROUP"] in resource_list)
                push!(resource_list, ALEAF_setting["Gen Technology"][gen_tech_id]["CAPCRED"])
            end
        end
    end
        
    network_data["zone"]["capacity_credit"] = Dict{String, Any}()
    regional_capacity_credits_resolution_level = ALEAF_setting["Network Setting"]["regional_capacity_credits_resolution_level"]   
    data_dict = network_data["supplement_network_data"][string("Network Data Level ", regional_capacity_credits_resolution_level)]
    network_data["zone"]["capacity_credit"]["data_identifier"] = data_identifier = ALEAF_setting["Network Setting"]["regional_capacity_credits_resolution_type"]

    # Update resource capacity ccredit data; We need to define this for each bus (average value)
    for n in keys(data_dict)

        region_ID = data_dict[n]["Region_ID"]
        region_ID = region_ID isa AbstractString ? region_ID : string(region_ID)

        network_data["zone"]["capacity_credit"][region_ID] = Dict()
        
        # gen technology
        for gen_tech_id in resource_list

            tag = string("CAPCRED_", gen_tech_id)

            network_data["zone"]["capacity_credit"][region_ID][gen_tech_id] = data_dict[n][tag]
        end
        
    end

end


function check_manual_rep_day_file_existence!(network_data::Dict{String, <:Any}, ALEAF_setting::Dict{String,<:Any}, case_id::Int, num_rep_daygroups::Int, num_days_in_group::Int, model_type::String)

    try
        file_name = string("Manual_", model_type, "_repDays_", ALEAF_setting["Simulation Configuration"][string(case_id)]["Repday_File_ID"], ".csv")
        repday_groups_data = read_csv_return_dict_string_any_per_round(joinpath(ALEAF_setting["data_location_timeseries"], file_name), network_data["simulation_round_idx"])
        
        # check num of repday groups from the manual input file
        if (length(repday_groups_data) == num_rep_daygroups) && ((repday_groups_data["1"]["End_Day_Id"] - repday_groups_data["1"]["Start_Day_Id"] + 1) == num_days_in_group)
            network_data["repday_groups"] = repday_groups_data
        else
            @aleaf_info "[ALEAF LC_GTEP Network Generation]:\tERROR! The manual repdays info does not match with the repdays setting. Switch to Scenario Reduction method"
            network_data["repday_groups"] = "Scenario Reduction"
        end
    
    catch
        @aleaf_info "[ALEAF LC_GTEP Network Generation]:\tFailed to collect the specified manual repdays info. Switch to Scenario Reduction method"
        network_data["repday_groups"] = "Scenario Reduction"
    end

end


function get_repday_groups_data!(network_data::Dict{String, <:Any}, ALEAF_setting::Dict{String, <:Any}, case_id::Int, LC_GTEP_model_type::String, year_id; print_output_flag::Bool=true, recorded_investment_decisions::Dict{String,<:Any} = Dict{String,Any}(), GTEP_multi_round_info_data::Dict{Any,<:Any} = Dict{Any,Any}(), target_day_group_id::Int = 0, precomputed_repday_groups::AbstractDict = Dict{String,Any}(), precomputed_repdays::AbstractDict = Dict{String,Any}(), build_hourly_data::Bool = true)

    # Single-day-group worker path -------------------------------------------------------------
    # When the master has already run scenario reduction and passes the resulting assignment
    # (precomputed_repday_groups + per-repday metadata), a distributed OP worker materializes ONLY
    # its target day-group. This skips the ~3-min clustering AND keeps the full-horizon repday/VRE
    # tensors from ever being built on the worker. Index keys are reused verbatim from the master's
    # assignment so worker solution keys line up with the master reporting reference.
    if target_day_group_id > 0 && !isempty(precomputed_repday_groups)
        tgid = string(target_day_group_id)
        haskey(precomputed_repday_groups, tgid) || terminate_with_error(; msg="[ALEAF LC_GTEP Network Generation]:\tprecomputed_repday_groups has no key $tgid; cannot build single day-group")
        network_data["repday_groups"] = Dict{String,Any}(tgid => deepcopy(precomputed_repday_groups[tgid]))
        stage1_repdays = Dict{String,Any}()
        for d in precomputed_repday_groups[tgid]["Day_Idx_List"]
            stage1_repdays[string(d)] = deepcopy(precomputed_repdays[string(d)])
        end
        network_data["planning_stages"]["1"]["repdays"] = stage1_repdays
        # Fall through to the shared tail below (copy to remaining stages, allocate hourly ["data"]
        # for the retained repday(s), hydro budget/water value, and the flat repdays reference).

    elseif target_day_group_id == 0 && !isempty(precomputed_repday_groups)

        # Full serial assignment from a coarse-resolution selection: the day assignment (all day-groups
        # + per-repday metadata) was already picked on a lightweight selection network, so reuse it
        # verbatim and skip selection/probability-normalization/day-list construction (those fields are
        # already baked into the snapshot). Preserves stochastic Scenario_ID metadata. Falls through to
        # the shared tail (hourly ["data"], hydro, flat reference, stochastic timeseries expansion).
        network_data["repday_groups"] = deepcopy(precomputed_repday_groups)
        stage1_repdays = Dict{String,Any}()
        for (d, meta) in precomputed_repdays
            stage1_repdays[string(d)] = deepcopy(meta)
        end
        network_data["planning_stages"]["1"]["repdays"] = stage1_repdays

    else

    # define number of repday groups and days in a group
    num_rep_daygroups = 0
    num_days_in_group = 0
    rep_day_selection_method = "Manual"
    rep_day_model_type = "EXP"
    stochastic_EXP_flag = false
    if LC_GTEP_model_type == "expansion"
        num_rep_daygroups = ALEAF_setting["Simulation Configuration"][string(case_id)]["NDAY_Groups"]
        num_days_in_group = ALEAF_setting["Simulation Configuration"][string(case_id)]["NDAYS_in_Single_Group"]
        rep_day_selection_method = ALEAF_setting["Simulation Configuration"][string(case_id)]["Repday_Selection_Mode_EXP"]
        rep_day_model_type = "EXP"
        if ALEAF_setting["Simulation Configuration"][string(case_id)]["Stochastic_Expansion_Flag"] == true
            stochastic_EXP_flag = true
        end
    elseif LC_GTEP_model_type == "operation"
        num_rep_daygroups = ALEAF_setting["Simulation Configuration"][string(case_id)]["NDAY_Groups_OP"]
        num_days_in_group = ALEAF_setting["Simulation Configuration"][string(case_id)]["NDAYS_in_Single_Group_OP"]
        rep_day_selection_method = ALEAF_setting["Simulation Configuration"][string(case_id)]["Repday_Selection_Mode_OP"]
        rep_day_model_type = "OP"
    elseif LC_GTEP_model_type == "RA"
        num_rep_daygroups = ALEAF_setting["Simulation Configuration"][string(case_id)]["NDAY_Groups_RA"]
        num_days_in_group = ALEAF_setting["Simulation Configuration"][string(case_id)]["NDAYS_in_Single_Group_RA"]
        rep_day_selection_method = ALEAF_setting["Simulation Configuration"][string(case_id)]["Repday_Selection_Mode_RA"]
        rep_day_model_type = "RA"
    end

    # Read manual repday file
    if rep_day_selection_method == "Manual"
        check_manual_rep_day_file_existence!(network_data, ALEAF_setting, case_id, num_rep_daygroups, num_days_in_group, rep_day_model_type)
    end 

    # Perform scenario reduction 
    if rep_day_selection_method == "Scenario Reduction" || network_data["repday_groups"] == "Scenario Reduction"    # we also run scenario reduction method if manual repday file does not exist

        try
            rep_days_file_name = string("repDays_", rep_day_model_type, "_G", num_rep_daygroups, "_D",  num_days_in_group, ".csv")

            # generate repdays; the result df is threaded back in memory (no transient CSV in the data folder)
            rep_day_input_df = generate_rep_days_LCO_GTEP(ALEAF_setting, network_data, num_rep_daygroups, num_days_in_group, case_id, rep_days_file_name; LC_GTEP_model_type, GTEP_multi_round_info_data, recorded_investment_decisions, print_output_flag)

            # CSV round-trip through an IOBuffer so the parsed types match the former file-read path exactly
            io = IOBuffer()
            CSV.write(io, rep_day_input_df)
            seekstart(io)
            network_data["repday_groups"] = convert_dataFrame_to_dict_string_any(CSV.read(io, DataFrame))

        catch
            terminate_with_error(; msg="[ALEAF LC_GTEP Network Generation]:\tFailed to perform Scenario Reduction method; terminate the program")
        end
    end

    # Calculate number of groups
    num_of_day_groups = div(365, ALEAF_setting["Simulation Configuration"][string(case_id)]["NDAY_Groups"])
    total_day_list = []
    
    # Adjust probability to make sure that the sum is 1
    current_prob = sum(network_data["repday_groups"][idx]["Probability"] for idx in keys(network_data["repday_groups"]))
    if current_prob != 1.0
        distribute_difference = (1 - current_prob) / ALEAF_setting["Simulation Configuration"][string(case_id)]["NDAY_Groups"]   
        for idx in keys(network_data["repday_groups"])
            network_data["repday_groups"][idx]["Probability"] += distribute_difference
        end
    end

    # Add Month Info 
    function convert_day_to_month(day::Int)
        month_starts = [1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366] 
        for m in 1:12
            if day >= month_starts[m] && day < month_starts[m+1]
                return m
            end
        end
    end
    
    # update planning stages data
    stage = "1"
    total_day_list_idx = 1
    
    if stochastic_EXP_flag == true   
        
        # update repday groups first

        new_repday_group_id = 1
        new_repday_group = Dict{String, Any}()

        for sce_id in 1:ALEAF_setting["Simulation Configuration"][string(case_id)]["Num_Sto_Scenarios"]
        
            scenario_probability = 1 / ALEAF_setting["Simulation Configuration"][string(case_id)]["Num_Sto_Scenarios"]    
            
            for id in 1:length(keys(network_data["repday_groups"]))


                id = string(id)

                NumDays_of_day_group = 365 * network_data["repday_groups"][id]["Probability"] 
                NumDays_of_each_day_in_a_day_group = NumDays_of_day_group / num_days_in_group

                day_list = []
                day_idx_list = []
                start_day = network_data["repday_groups"][id]["Start_Day_Id"]
                end_day = network_data["repday_groups"][id]["End_Day_Id"]
                
                for day_idx in start_day:end_day
                    append!(day_list, day_idx)
                    append!(total_day_list, day_idx)
                    append!(day_idx_list, total_day_list_idx)

                    # update repdays dict
                    network_data["planning_stages"][stage]["repdays"][string(total_day_list_idx)] = Dict()
                    network_data["planning_stages"][stage]["repdays"][string(total_day_list_idx)]["Day_Group_ID"] = string(new_repday_group_id)
                    network_data["planning_stages"][stage]["repdays"][string(total_day_list_idx)]["Day_Group"] = network_data["repday_groups"][id]["Day_Group"]
                    network_data["planning_stages"][stage]["repdays"][string(total_day_list_idx)]["NumDays_Group"] = NumDays_of_day_group
                    network_data["planning_stages"][stage]["repdays"][string(total_day_list_idx)]["NumDays"] = NumDays_of_each_day_in_a_day_group
                    network_data["planning_stages"][stage]["repdays"][string(total_day_list_idx)]["Day"] = day_idx
                    network_data["planning_stages"][stage]["repdays"][string(total_day_list_idx)]["Scenario_ID"] = sce_id
                    network_data["planning_stages"][stage]["repdays"][string(total_day_list_idx)]["Month"] = convert_day_to_month(day_idx)

                    total_day_list_idx += 1
                end

                # update repday_groups dict
                new_repday_group[string(new_repday_group_id)] = deepcopy(network_data["repday_groups"][id])
                new_repday_group[string(new_repday_group_id)]["NumDays_Group"] = NumDays_of_day_group
                new_repday_group[string(new_repday_group_id)]["Day_List"] = day_list
                new_repday_group[string(new_repday_group_id)]["Day_Idx_List"] = day_idx_list
                new_repday_group[string(new_repday_group_id)]["Scenario_ID"] = sce_id

                new_repday_group_id += 1
            end

        end

        # replace repday_group 
        network_data["repday_groups"] = new_repday_group

    else 

        for id in 1:length(keys(network_data["repday_groups"]))

            id = string(id)
    
            NumDays_of_day_group = 365 * network_data["repday_groups"][id]["Probability"] 
            NumDays_of_each_day_in_a_day_group = NumDays_of_day_group / num_days_in_group
    
            day_list = []
            day_idx_list = []
            start_day = network_data["repday_groups"][id]["Start_Day_Id"]
            end_day = network_data["repday_groups"][id]["End_Day_Id"]
            
            for day_idx in start_day:end_day
                append!(day_list, day_idx)
                append!(total_day_list, day_idx)
                append!(day_idx_list, total_day_list_idx)
    
                # update repdays dict
                network_data["planning_stages"][stage]["repdays"][string(total_day_list_idx)] = Dict()
                network_data["planning_stages"][stage]["repdays"][string(total_day_list_idx)]["Day_Group_ID"] = id
                network_data["planning_stages"][stage]["repdays"][string(total_day_list_idx)]["Day_Group"] = network_data["repday_groups"][id]["Day_Group"]
                network_data["planning_stages"][stage]["repdays"][string(total_day_list_idx)]["NumDays_Group"] = NumDays_of_day_group
                network_data["planning_stages"][stage]["repdays"][string(total_day_list_idx)]["NumDays"] = NumDays_of_each_day_in_a_day_group
                network_data["planning_stages"][stage]["repdays"][string(total_day_list_idx)]["Day"] = day_idx
                network_data["planning_stages"][stage]["repdays"][string(total_day_list_idx)]["Month"] = convert_day_to_month(day_idx)
    
                total_day_list_idx += 1
            end
    
            # update repday_groups dict
            network_data["repday_groups"][id]["NumDays_Group"] = NumDays_of_day_group
            network_data["repday_groups"][id]["Day_List"] = day_list
            network_data["repday_groups"][id]["Day_Idx_List"] = day_idx_list
        end

    end

    end  # end of full-vs-precomputed repday assignment

    # copy repday data to remaining planning stages
    if length(keys(network_data["planning_stages"])) > 1
        for stage in 2:length(keys(network_data["planning_stages"]))
            network_data["planning_stages"][string(stage)]["repdays"] = deepcopy(network_data["planning_stages"]["1"]["repdays"])
        end
    end

    # Collect additional time series data (populates the per-repday hourly ["data"] sub-dicts).
    # build_hourly_data=false (OP-distributed master light reference): skip this — the master only needs
    # repday METADATA + full-year time_series_data DataFrames for reporting; workers rebuild their own
    # single-day ["data"]. This is the O(day-groups × region) memory the master otherwise carries.
    if build_hourly_data
        time_series_data_list = ["load", "wind_ons", "wind_ofs", "pv", "hydro", "rtpv", "csp", "dc", "reg_up", "reg_down", "spin", "nspin", "flex_up", "flex_down", "reg_up_signal", "reg_down_signal"]
        load_timeseries_data!(ALEAF_setting["data_location"], ALEAF_setting, network_data, time_series_data_list, case_id)
    end

    # Read hydropower budget and water value data
    hydro_budget_scenario = ALEAF_setting["Simulation Configuration"][string(case_id)]["Hydro_Budget_File_ID"]
    hydro_budget_file_name = string("HYDRO/timeseries_hydro_daily_budget_" , hydro_budget_scenario, ".csv")
    if hydro_budget_scenario != "Base"
        hydro_budget_file_name = string("0_additional_scenarios/HYDRO/timeseries_hydro_daily_budget_" , hydro_budget_scenario, ".csv")
    end
    network_data["time_series_data"]["hydro_budget_df"] = DataFrame(CSV.File(joinpath(ALEAF_setting["data_location_timeseries"], hydro_budget_file_name)))
    
    water_value_scenario = ALEAF_setting["Simulation Configuration"][string(case_id)]["Hydro_Value_File_ID"]
    water_value_file_name = string("HYDRO/timeseries_hydro_daily_water_value_" , water_value_scenario, ".csv")
    if water_value_scenario != "Base"
        water_value_file_name = string("0_additional_scenarios/HYDRO/timeseries_hydro_daily_water_value_" , water_value_scenario, ".csv")
    end
    network_data["time_series_data"]["hydro_water_value_df"] = DataFrame(CSV.File(joinpath(ALEAF_setting["data_location_timeseries"], water_value_file_name)))

    # Add repday reference for easier access (without data)
    network_data["repdays"] = deepcopy(network_data["planning_stages"]["1"]["repdays"])

    # check and update stochastic timeseries data 
    if (LC_GTEP_model_type == "expansion")
        replace_stochastic_time_series_data!(ALEAF_setting, network_data, case_id, ALEAF_setting["data_location"])
    end

end


function determine_opt_model_type!(network_data::Dict{String, <:Any}, ALEAF_setting::Dict{String, <:Any}, case_id, LC_GTEP_model_type)

    model_type_EXP = "LP"
    EXP_decison_type = "LP"
    model_type_OP = "LP"
    gen_tech_data = ALEAF_setting["Gen Technology"]

    EXP_decison_type = "LP"   # default decision type for expansion model
    EXP_storage_commitment = false   # default storage commitment for expansion model
    EXP_unit_commitment = false
    EXP_model_type = "LP"   # default model type for expansion model

    OP_storage_commitment = false   # default storage commitment for operation model
    OP_unit_commitment = false
    OP_model_type = "LP"   # default model type for operation model

    if LC_GTEP_model_type == "expansion"
        
        # check integrality of investment decisions
        any(gen_tech_data[gen]["Integrality"] for gen in keys(gen_tech_data)) ? EXP_decison_type = "MILP" : EXP_decison_type = "LP" # check integrality of investment decisions

        # check storage commitment
        any(get(gen_tech_data[gen], "Storage Commitment", false) == true for gen in keys(gen_tech_data)) ? EXP_storage_commitment = true : EXP_storage_commitment = false   # check storage commitment; commitment can include "NA"
        
        # check unit commitment
        ALEAF_setting["Simulation Configuration"][string(case_id)]["Dispatch_Mode_in_EXP"] == "Unit Commitment" ? EXP_unit_commitment = true : EXP_unit_commitment = false   # check unit commitment
        
        # update model type
        if EXP_decison_type == "MILP" || EXP_storage_commitment || EXP_unit_commitment
            model_type_EXP = "MILP"
        else
            model_type_EXP = "LP"
        end
    end

    if LC_GTEP_model_type == "operation"

        # check storage commitment
        any(get(gen_tech_data[gen], "Storage Commitment", false) == true for gen in keys(gen_tech_data)) ? OP_storage_commitment = true : OP_storage_commitment = false   # check storage commitment; commitment can include "NA"
        
        # check unit commitment
        ALEAF_setting["Simulation Configuration"][string(case_id)]["Dispatch_Mode_in_OP"] == "Unit Commitment" ? OP_unit_commitment = true : OP_unit_commitment = false   # check unit commitment
        
        # update model type
        if OP_unit_commitment || OP_storage_commitment
            model_type_OP = "MILP"
        else
            model_type_OP = "LP"
        end
    end

    network_data["model_type"] = Dict{String, Any}()
    network_data["model_type"]["model_type_EXP"] = model_type_EXP
    network_data["model_type"]["model_type_OP"] = model_type_OP
    network_data["model_type"]["EXP_decison_type"] = EXP_decison_type

end


"""
Resolve the representative-day SELECTION resolution into `(type_string, level_int, collapse_to_one)`,
or `nothing` when the coarse-decoupled selection should be skipped (fall back to in-line selection).

Driven by the GLOBAL `Scenario Reduction Setting["repday_selection_resolution"]` (case-insensitive):
  * "regional" => build the selection network at the network boundary resolution.
  * "system"   => same boundary resolution, then collapse the modeled boundary buses into a single bus.
Absent/blank/unrecognized => `nothing` (in-line selection; zero behavior change for setting files that
lack the key). The clustering metric is a system-wide sum over the SAME finest-region set for both
options, so both yield byte-identical rep days; only the selection-bus count differs. Both options stay
at the boundary level: aggregating coarser would let one coarse bus straddle the modeled/unmodeled
boundary filter and pull in unmodeled regions, so `system` is realized by a post-filter merge, not a
coarser resolution level.
"""
function resolve_repday_selection_resolution(ALEAF_setting::Dict{String, <:Any}, case_id)

    sr_setting = get(ALEAF_setting, "Scenario Reduction Setting", Dict{String, Any}())
    requested = get(sr_setting, "repday_selection_resolution", nothing)
    requested_str = requested === nothing ? "" : lowercase(strip(string(requested)))

    isempty(requested_str) && return nothing

    if requested_str != "regional" && requested_str != "system"
        @aleaf_info "[ALEAF LC_GTEP]: repday_selection_resolution='$requested' is unrecognized (expected 'regional' or 'system'); using in-line rep-day selection at the run resolution."
        return nothing
    end

    network_setting = get(ALEAF_setting, "Network Setting", Dict{String, Any}())
    boundary_type = get(network_setting, "network_boundary_type", nothing)
    boundary_type === nothing && return nothing

    # Both options select at the network boundary level; "system" additionally merges the modeled
    # boundary buses into one after the modeled-region filter (safe, no straddle).
    resolution_levels = get(ALEAF_setting, "network_resolution_level", Dict{String, Any}())
    boundary_level = nothing
    for entry in values(resolution_levels)
        rid = get(entry, "Network_Resolution_ID", nothing)
        lvl = get(entry, "Network_Resolution_Level", nothing)
        (rid === nothing || lvl === nothing) && continue
        if string(rid) == string(boundary_type)
            boundary_level = Int(lvl)
            break
        end
    end
    boundary_level === nothing && return nothing

    return (string(boundary_type), boundary_level, requested_str == "system")
end


"""
Lever A: build a per-bus map from each finest region id to its data-region at a (possibly coarser)
per-profile resolution, so a nodal net can read one shape per data-region and broadcast it to buses.

For each profile in {load, pv, wind_ons, wind_ofs, rtpv, csp, hydro, load_growth} the resolution type
is read from Network Setting as `<profile>_data_resolution_type`; when absent it defaults to the finest
resolution, which makes the map the identity and preserves byte-identical behavior with existing databases.
"""
function build_profile_data_region_maps!(network_data::Dict{String, <:Any}, ALEAF_setting::Dict{String, <:Any})

    profile_types = ["load", "pv", "wind_ons", "wind_ofs", "rtpv", "csp", "hydro", "load_growth"]

    finest_resolution_id = get(ALEAF_setting, "network_resolution_level", nothing) === nothing ? nothing :
        get(get(ALEAF_setting["network_resolution_level"], "1", Dict{String, Any}()), "Network_Resolution_ID", nothing)
    network_setting = get(ALEAF_setting, "Network Setting", Dict{String, Any}())
    sub_area_mapping = get(ALEAF_setting, "sub_area_mapping", Dict{String, Any}())

    # data-resolution type per profile (finest by default => identity map => unchanged behavior)
    resolution_type_by_profile = Dict{String, Any}()
    for profile in profile_types
        resolution_type_by_profile[profile] = get(network_setting, string(profile, "_data_resolution_type"), finest_resolution_id)
    end

    # finest region id => data-region id, per resolution type; only build for coarser (non-identity) resolutions
    data_region_by_finest = Dict{Any, Dict{String, String}}()
    for resolution_type in unique(values(resolution_type_by_profile))
        (resolution_type === nothing || resolution_type == finest_resolution_id) && continue
        lookup = Dict{String, String}()
        for key in keys(sub_area_mapping)
            finest_id = get(sub_area_mapping[key], finest_resolution_id, nothing)
            data_id = get(sub_area_mapping[key], resolution_type, nothing)
            (finest_id === nothing || data_id === nothing) && continue
            lookup[string(finest_id)] = string(data_id)
        end
        data_region_by_finest[resolution_type] = lookup
    end

    # global map (finest region id independent of bus): profile_type => (finest id => data-region id).
    # An empty per-profile map means identity, i.e. unchanged (finest-resolution) behavior.
    profile_data_region_map = Dict{String, Dict{String, String}}()
    for profile in profile_types
        resolution_type = resolution_type_by_profile[profile]
        if resolution_type !== nothing && resolution_type != finest_resolution_id && haskey(data_region_by_finest, resolution_type)
            profile_data_region_map[profile] = data_region_by_finest[resolution_type]
        else
            profile_data_region_map[profile] = Dict{String, String}()
        end
    end
    network_data["profile_data_region_map"] = profile_data_region_map

end


function get_network_zonal_data!(network_data::Dict{String, <:Any}, ALEAF_setting::Dict{String, <:Any}, case_id)

    network_data["zone"] = Dict{String, Any}()
    build_profile_data_region_maps!(network_data, ALEAF_setting)
    if !get(ALEAF_setting, "network_config_available", false)
        # When no config workbook is present, treat the entire system as one planning/reserve/policy region.
        all_bus_idx = sort(parse.(Int, collect(keys(network_data["bus"]))))
        all_gen_idx = collect(keys(network_data["plant"]))
        aggregated_bus_i = [string(network_data["bus"][string(idx)]["bus_i"]) for idx in all_bus_idx]
        system_data_region_ids = unique(aggregated_bus_i)

        system_aggregation = Dict{String, Any}(
            "aggregated_regions_bus_i" => aggregated_bus_i,
            "aggregated_regions_bus_idx" => all_bus_idx,
            "zone_bus_idx" => all_bus_idx,
            "aggregated_data_regions_id" => system_data_region_ids,
            "local_gen_idx" => all_gen_idx,
            "new_local_gen_idx" => Int[],
            "technology_list" => String[],
            "gen_technology_groups" => Dict{String, Any}(),
        )
        system_mapping = Dict{String, Any}("system" => "system", "System" => "system")

        network_data["zone"]["planning_reserve"] = Dict("1" => Dict(
            "zone_id" => "system",
            "aggregation_info" => deepcopy(system_aggregation),
            "zone_regional_mapping_info" => deepcopy(system_mapping),
        ))
        network_data["zone"]["reserve"] = Dict("1" => Dict(
            "zone_id" => "system",
            "aggregation_info" => deepcopy(system_aggregation),
            "zone_regional_mapping_info" => deepcopy(system_mapping),
            "data_area_list" => deepcopy(system_data_region_ids),
        ))
        network_data["zone"]["policy"] = Dict("1" => Dict(
            "zone_id" => "system",
            "aggregation_info" => deepcopy(system_aggregation),
            "zone_regional_mapping_info" => deepcopy(system_mapping),
            "data_area_list" => deepcopy(system_data_region_ids),
            "vre_aggregated_data" => Dict{String, Any}(),
        ))

        for target_policy in ["RPS", "CEGT", "CERT"]
            network_data["zone"]["policy"]["1"][target_policy] = Dict{String, Any}()
            network_data["zone"]["policy"]["1"][target_policy]["target"] = Dict{String, Any}()
            for y in keys(network_data["planning_stages"])
                if target_policy == "RPS"
                    network_data["zone"]["policy"]["1"][target_policy]["target"][y] = ALEAF_setting["Simulation Configuration"][string(case_id)]["RPS_Global_Target_Value"]
                elseif target_policy == "CEGT"
                    network_data["zone"]["policy"]["1"][target_policy]["target"][y] = ALEAF_setting["Simulation Configuration"][string(case_id)]["Clean_Energy_Generation_Global_Target_Value"]
                else
                    network_data["zone"]["policy"]["1"][target_policy]["target"][y] = ALEAF_setting["Simulation Configuration"][string(case_id)]["Carbon_Emission_Reduction_Global_Target_Value"]
                end
            end
        end
        network_data["zone"]["policy"]["1"]["CERT"]["reference"] = 0.0
        network_data["zone"]["policy"]["1"]["Intertia"] = Dict("target" => 0.0)

        # Resource, cost-scaling, and capacity-credit lookups still expect zone containers, so provide neutral system-wide defaults.
        network_data["zone"]["resource_supply_curve"] = Dict{String, Any}()
        network_data["zone"]["resource_supply_curve"]["data_resolution_type"] = "system"
        network_data["zone"]["resource_supply_curve"]["supply_curve_resource_list"] = String[]
        network_data["zone"]["resource_supply_curve"]["supply_curve_region_list"] = Int[]

        network_data["zone"]["cost_scaling"] = Dict{String, Any}()
        network_data["zone"]["cost_scaling"]["data_identifier"] = "system"
        network_data["zone"]["cost_scaling"]["system"] = Dict{String, Any}()
        for gen_tech_id in keys(ALEAF_setting["Gen Technology"])
            network_data["zone"]["cost_scaling"]["system"][gen_tech_id] = Dict(
                "UNITGROUP" => ALEAF_setting["Gen Technology"][gen_tech_id]["UNITGROUP"],
                "CAPAX_scale" => 1.0,
            )
        end

        network_data["zone"]["capacity_credit"] = Dict{String, Any}()
        network_data["zone"]["capacity_credit"]["data_identifier"] = "system"
        network_data["zone"]["capacity_credit"]["system"] = Dict{String, Any}()
        for gen_tech_id in keys(ALEAF_setting["Gen Technology"])
            unitgroup = ALEAF_setting["Gen Technology"][gen_tech_id]["UNITGROUP"]
            network_data["zone"]["capacity_credit"]["system"][gen_tech_id] = 1.0    # set 1.0 if zonal value is not available
            network_data["zone"]["capacity_credit"]["system"][unitgroup] = 1.0   # set 1.0 if zonal value is not available
        end

        return
    end

    # 1) Planning Reserve Zone 
    zone_boundary_type = ALEAF_setting["Network Setting"]["planning_reserve_zone_boundary_type"]   # reserve zone boundary type
    data_resolution_type = "NA"   # reserve data resolution type
    data_resolution_level = "NA"   # reserve data resolution type
    define_zone_wo_data!(network_data, ALEAF_setting, "planning_reserve", zone_boundary_type)

    # 2) Operating Reserve Zone 
    zone_boundary_type = ALEAF_setting["Network Setting"]["operating_reserve_zone_boundary_type"]   # reserve zone boundary type
    data_resolution_type = ALEAF_setting["Network Setting"]["operating_reserve_data_resolution_type"]   # reserve data resolution type
    data_resolution_level = ALEAF_setting["Network Setting"]["operating_reserve_data_resolution_level"]   # reserve data resolution type
    define_zone!(network_data, ALEAF_setting, "reserve", zone_boundary_type, data_resolution_type, data_resolution_level)

    # 2) Policy Zone
    zone_boundary_type = ALEAF_setting["Network Setting"]["policy_zone_boundary_type"]   # policy zone boundary type
    data_resolution_type = ALEAF_setting["Network Setting"]["policy_data_resolution_type"]   # policy data resolution type
    data_resolution_level = ALEAF_setting["Network Setting"]["policy_data_resolution_level"]   # policy data resolution type
    define_zone!(network_data, ALEAF_setting, "policy", zone_boundary_type, data_resolution_type, data_resolution_level)

    # 2-1) Policy Targets : RPS
    get_policy_targets!(network_data, ALEAF_setting, case_id, "RPS", ALEAF_setting["Simulation Configuration"][string(case_id)]["RPS_Global_Target_Value"])

    # 2-2) Policy Targets : CEGT
    get_policy_targets!(network_data, ALEAF_setting, case_id, "CEGT", ALEAF_setting["Simulation Configuration"][string(case_id)]["Clean_Energy_Generation_Global_Target_Value"])

    # 2-3) Policy Targets : CERT 
    get_policy_targets!(network_data, ALEAF_setting, case_id, "CERT", ALEAF_setting["Simulation Configuration"][string(case_id)]["Carbon_Emission_Reduction_Global_Target_Value"])

    # 3) Resource supply zone 
    get_resource_supply_curve!(network_data, ALEAF_setting, case_id)
    
    # 4) Regional Cost Variations 
    get_resource_cost_variations!(network_data, ALEAF_setting, case_id)

    # 5) Regional Capacity Credits
    get_resource_capacity_credits!(network_data, ALEAF_setting, case_id)

    # 6) Inertia 
    # Hard-coded values for inertia targets 
    for n in keys(network_data["zone"]["policy"])
        network_data["zone"]["policy"][n]["Intertia"] = Dict{String, Any}()
        network_data["zone"]["policy"][n]["Intertia"]["target"] = 0.0

        try 
            inter_connection_id = network_data["zone"]["policy"][n]["zone_regional_mapping_info"]["Interconnection"]
            for id in keys(network_data["additional_network_data"]["Network Data Level 5"])
                if network_data["additional_network_data"]["Network Data Level 5"][id]["Region_ID"] == inter_connection_id
                    network_data["zone"]["policy"][n]["Intertia"]["target"] += network_data["additional_network_data"]["Network Data Level 5"][id]["minimum_system_inertia_MVA*s"]
                end
            end

        catch
            # do nothing if the inertia target is not available
        end
    end
    
end





