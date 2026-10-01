# Copyright © 2026, UChicago Argonne, LLC. All Rights Reserved.
# SPDX-License-Identifier: BSD-3-Clause  (see LICENSE)

# Generic, model-agnostic utilities and helpers shared across ALEAF modules.

using LinearAlgebra
using JSON

function terminate_with_error(; msg::String="debugging")
    error("ERROR: " * msg)
end

function get_index(am::Abstract_ALEAF_Model, nw::Int, key1::Symbol, key2::Symbol)
    return [(key1["index"], key2["index"]) for (i,key1) in am.ref[:nw][nw][key1] for (j, key2) in am.ref[:nw][nw][key2]]
end


function get_index(am::Abstract_ALEAF_Model, nw::Int, key1::Symbol)
    return [(key1["index"]) for (i,key1) in am.ref[:nw][nw][key1]]
end


function get_index(am::Abstract_ALEAF_Model, key1::Symbol; nw::Int=am.cnw)
    return keys(am.ref[:nw][nw][key1])
end


function get_index(am::Abstract_ALEAF_Model, key1::Symbol, nw::Int)
    return keys(am.ref[:nw][nw][key1])
end

function get_index(data::Dict{Int64, Any}, key1::Symbol, nw::Int)
    try
        return keys(data[nw][key1])
    catch
        return keys(data[string(nw)][string(key1)])
    end
end

function get_index(data::Dict{String, Any}, key1::Symbol, nw::Int)
    try
        return keys(data[string(nw)][key1])
    catch
        return keys(data[string(nw)][string(key1)])
    end
end


function get_index(am::Abstract_ALEAF_Model, key1::Symbol, nw::Int, key2::Int)
    return keys(am.ref[:nw][nw][key1][key2])
end


function get_index(am::Abstract_ALEAF_Model, key1::Symbol, nw::Int, key2::String)
    return [parse(Int, k) for k in keys(am.ref[:nw][nw][key1][key2])]
end


function get_index(am::Abstract_ALEAF_Model, nw::Int, key1::Symbol, idx)
    return am.ref[:nw][nw][key1][idx]
end


# One slack bus per synchronous island (connected components of the AC-only graph; DC ties
# excluded so async interconnections each get their own anchor). Reference = min bus id.
function get_ac_reference_buses(am::Abstract_ALEAF_Model; nw::Int=0)
    bus_ids = collect(get_index(am, :bus, nw))
    parent = Dict{Int,Int}(b => b for b in bus_ids)
    function find(x)
        while parent[x] != x
            parent[x] = parent[parent[x]]   # path halving
            x = parent[x]
        end
        return x
    end
    function unite(a, b)
        ra, rb = find(a), find(b)
        ra == rb && return
        ra < rb ? (parent[rb] = ra) : (parent[ra] = rb)   # union by min id
    end
    for k in get_index(am, :branch, nw)
        br = am.ref[:nw][nw][:branch][k]
        get(br, "model_flag", true) == true || continue
        get(br, "dc_line", false) == true && continue      # AC-only graph
        f = br["f_bus"]; t = br["t_bus"]
        (haskey(parent, f) && haskey(parent, t)) || continue
        unite(f, t)
    end
    refs = Set{Int}()
    for b in bus_ids
        push!(refs, find(b))
    end
    return refs
end


function variable(am::Abstract_ALEAF_Model, nw::Int, key1::Symbol, idx)::JuMP.VariableRef
    return am.var[:nw][nw][key1][idx]
end


function variable(am::Abstract_ALEAF_Model, nw::Int, day_id::Int, key1::Symbol, idx)::JuMP.VariableRef
    return am.var[:nw][nw][day_id][key1][idx]
end



function variable(am::Abstract_ALEAF_Model, nw::Int, key1::Symbol)
    return am.var[:nw][nw][key1]
end


function parameter(am::Abstract_ALEAF_Model, nw::Int, key1::Symbol, key2::String, idx)
    return am.ref[:nw][nw][key1][idx][key2]
end


function parameter(data::Dict{Int64, Any}, nw::Int, key1::Symbol)
    try
        return data[nw][key1]
    catch
        return data[string(nw)][string(key1)]
    end
end


function parameter(data::Dict{Int64, Any}, nw::Int, key1::Symbol, key2::String, idx)
    try
        return data[nw][key1][idx][key2]
    catch
        return data[string(nw)][string(key1)][string(idx)][key2]
    end
end

function parameter(data::Dict{String, Any}, nw::Int, key1::Symbol, key2::String, idx)
    try
        return data[string(nw)][key1][idx][key2]
    catch
        return data[string(nw)][string(key1)][string(idx)][key2]
    end        
end

function parameter(am::Abstract_ALEAF_Model, nw::Int, key1::Symbol, idx::Int, key2::String)
    return am.ref[:nw][nw][key1][idx][key2]
end

function parameter(data::Dict{Int64, Any}, nw::Int, key1::Symbol, idx::Int, key2::String)
    try
        return data[nw][key1][idx][key2]
    catch
        return data[string(nw)][string(key1)][string(idx)][key2]
    end
end

function parameter(data::Dict{String, Any}, nw::Int, key1::Symbol, idx::Int, key2::String)
    try
        return data[string(nw)][key1][idx][key2]
    catch
        return data[string(nw)][string(key1)][string(idx)][key2]
    end
end


function parameter(am::Abstract_ALEAF_Model, nw::Int, key1::Symbol, key2::String, idx1, idx2)
    return am.ref[:nw][nw][key1][idx1][key2][idx2]
end


function parameter(am::Abstract_ALEAF_Model, nw::Int, key1::Symbol, idx1, idx2, key2::String)
    return am.ref[:nw][nw][key1][idx1][idx2][key2]
end


function parameter(am::Abstract_ALEAF_Model, nw::Int, key1::Symbol, key2::String, key3::String, d::Int, h::Int, t::Int)
    return am.ref[:nw][nw][key1][d][key2][string(h)][string(t)][key3]
end


function parameter(am::Abstract_ALEAF_Model, nw::Int, key1::Symbol, key2::String, key3::String, key4::String, d::Int, h::Int, t::Int, y::Int)
    return am.ref[:nw][nw][key1][y][key2][string(d)][key3][string(h)][string(t)][key4]
end


function parameter(am::Abstract_ALEAF_Model, nw::Int, key1::Symbol, bus::Int, line::Int)
    return am.ref[:nw][nw][key1][bus][line]
end


function parameter(am::Abstract_ALEAF_Model, nw::Int, key1::Symbol, key2::String)
    return am.ref[:nw][nw][key1][key2]
end


function parameter(am::Abstract_ALEAF_Model, nw::Int, key1::Symbol, idx::Int)
    return am.ref[:nw][nw][key1][idx]
end


function parameter(am::Abstract_ALEAF_Model, nw::Int, key1::Symbol)
    return am.ref[:nw][nw][key1]
end


function parameter(data::Dict{String, Any}, nw::Int, key1::Symbol)
    try
        return data[nw][key1]
    catch
        return data[string(nw)][string(key1)]
    end
end



function parameter(am::Abstract_ALEAF_Model, key1::String)
    return am.setting[key1]
end


function parameter(am::Abstract_ALEAF_Model, key1::String, key2::String)
    return am.setting[key1][key2]
end


function get_load_growth_factor(stage_data::Dict{<:Any,<:Any}, region_id::AbstractString)
    return stage_data["load_growth_by_region"][string(region_id)]
end


function get_load_growth_factor(network_data::Dict{<:Any,<:Any}, stage, region_id::AbstractString)
    return get_load_growth_factor(network_data["planning_stages"][string(stage)], region_id)
end


function get_load_growth_factor(network_data::Dict{Symbol,<:Any}, stage, region_id::AbstractString)
    return network_data[:planning_stages][stage]["load_growth_by_region"][string(region_id)]
end


function get_load_growth_factor(am::Abstract_ALEAF_Model, y::Int, region_id::AbstractString; nw::Int=0)
    return get_load_growth_factor(am.ref[:nw][nw][:planning_stages][y], region_id)
end


# Lever A: resolve the data-region key for a finest region under a profile's data resolution.
# `region_map` is the global profile_type => (finest id => data-region id) map; falls back to the
# region id itself (finest / identity) so behavior is unchanged when unconfigured or at finest resolution.
function profile_data_region(region_map, profile_type::String, region_id)
    region_map === nothing && return region_id
    profile_map = get(region_map, profile_type, nothing)
    (profile_map === nothing || isempty(profile_map)) && return region_id
    return get(profile_map, string(region_id), region_id)
end


function get_bus_demand_with_growth(am::Abstract_ALEAF_Model, n::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=0, system_peak_scale::Float64=1.0)
    demand = 0.0
    load_shape = parameter(am, nw, :planning_stages, "repdays", "data", "load", d, h, t, y)
    bus_data = am.ref[:nw][nw][:bus][n]
    region_map = get(am.ref[:nw][nw], :profile_data_region_map, nothing)
    original_load = bus_data["aggregation_info"]["original_load_(bus_i, MW)"]
    for region_id in bus_data["aggregation_info"]["aggregated_regions_bus_i"]
        demand += load_shape[profile_data_region(region_map, "load", region_id)] * original_load[region_id] * get_load_growth_factor(am, y, region_id; nw) * system_peak_scale
    end
    return demand
end


function apply_tnd_loss(am::Abstract_ALEAF_Model, demand::Float64)
    if am.setting["Planning Design"]["enforce_transmission_loss_flag"] == true
        return demand * (1 + am.setting["Planning Design"]["transmission_loss_percent_value"] * 0.01)
    end
    return demand
end


function get_bus_demand_with_growth(data::Dict, n::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=0, system_peak_scale::Float64=1.0)
    demand = 0.0
    planning_stage_data = data[nw][:planning_stages][y]
    load_shape = planning_stage_data["repdays"][string(d)]["data"][string(h)][string(t)]["load"]
    bus_data = data[nw][:bus][n]
    region_map = get(data[nw], :profile_data_region_map, nothing)
    original_load = bus_data["aggregation_info"]["original_load_(bus_i, MW)"]
    for region_id in bus_data["aggregation_info"]["aggregated_regions_bus_i"]
        demand += load_shape[profile_data_region(region_map, "load", region_id)] * original_load[region_id] * get_load_growth_factor(data[nw], y, region_id) * system_peak_scale
    end
    return demand
end


function get_bus_peak_demand_with_growth(am::Abstract_ALEAF_Model, n::Int, y::Int; nw::Int=0)
    peak_demand = 0.0
    original_load = am.ref[:nw][nw][:bus][n]["aggregation_info"]["original_load_(bus_i, MW)"]
    for region_id in am.ref[:nw][nw][:bus][n]["aggregation_info"]["aggregated_regions_bus_i"]
        peak_demand += original_load[region_id] * get_load_growth_factor(am, y, region_id; nw)
    end
    return peak_demand
end


function get_zone_peak_demand_with_growth(am::Abstract_ALEAF_Model, zone_bus_ids, y::Int; nw::Int=0)
    peak_demand = 0.0
    for bus_id in zone_bus_ids
        peak_demand += get_bus_peak_demand_with_growth(am, bus_id, y; nw)
    end
    return peak_demand
end


function get_zone_peak_demand_with_growth(network_data::Dict{<:Any,<:Any}, zone_bus_ids, y)
    peak_demand = 0.0
    for bus_id in zone_bus_ids
        original_load = network_data["bus"][bus_id]["aggregation_info"]["original_load_(bus_i, MW)"]
        for region_id in network_data["bus"][bus_id]["aggregation_info"]["aggregated_regions_bus_i"]
            peak_demand += original_load[region_id] * get_load_growth_factor(network_data, y, region_id)
        end
    end
    return peak_demand
end


function get_zone_peak_demand_with_growth(network_data::Dict{Symbol,<:Any}, zone_bus_ids, y)
    peak_demand = 0.0
    for bus_id in zone_bus_ids
        original_load = network_data[:bus][bus_id]["aggregation_info"]["original_load_(bus_i, MW)"]
        for region_id in network_data[:bus][bus_id]["aggregation_info"]["aggregated_regions_bus_i"]
            peak_demand += original_load[region_id] * get_load_growth_factor(network_data, y, region_id)
        end
    end
    return peak_demand
end


function get_zone_annual_demand_with_growth(am::Abstract_ALEAF_Model, zone_bus_ids, y::Int; nw::Int=0)
    demand_df = am.ref[:nw][nw][:time_series_data]["load"]
    region_map = get(am.ref[:nw][nw], :profile_data_region_map, nothing)
    annual_demand = 0.0
    for bus_id in zone_bus_ids
        bus_data = am.ref[:nw][nw][:bus][bus_id]
        original_load = bus_data["aggregation_info"]["original_load_(bus_i, MW)"]
        for region_id in bus_data["aggregation_info"]["aggregated_regions_bus_i"]
            annual_demand += sum(demand_df[!, Symbol(profile_data_region(region_map, "load", region_id))]) * original_load[region_id] * get_load_growth_factor(am, y, region_id; nw)
        end
    end
    return annual_demand
end


function get_annual_load_growth_factor(network_data::Dict{<:Any,<:Any}, year, region_id::AbstractString)
    return network_data["annual_load_growth_by_region"][string(year)][string(region_id)]
end


function get_annual_load_growth_factor(network_data::Dict{Symbol,<:Any}, year, region_id::AbstractString)
    annual_growth = network_data[:annual_load_growth_by_region]
    if haskey(annual_growth, year)
        return annual_growth[year][string(region_id)]
    end
    return annual_growth[string(year)][string(region_id)]
end


function get_annual_peak_demand_with_growth(network_data::Dict{<:Any,<:Any}, year)
    peak_demand = 0.0
    for bus_data in values(network_data["bus"])
        original_load = bus_data["aggregation_info"]["original_load_(bus_i, MW)"]
        for region_id in bus_data["aggregation_info"]["aggregated_regions_bus_i"]
            peak_demand += original_load[region_id] * get_annual_load_growth_factor(network_data, year, region_id)
        end
    end
    return peak_demand
end


function get_annual_peak_demand_with_growth(network_data::Dict{Symbol,<:Any}, year)
    peak_demand = 0.0
    for bus_data in values(network_data[:bus])
        original_load = bus_data["aggregation_info"]["original_load_(bus_i, MW)"]
        for region_id in bus_data["aggregation_info"]["aggregated_regions_bus_i"]
            peak_demand += original_load[region_id] * get_annual_load_growth_factor(network_data, year, region_id)
        end
    end
    return peak_demand
end


# Coincident peak: sum the per-region hourly load shapes first, then take the max over the
# full year. Avoids overstating the peak by adding non-simultaneous regional peaks.
function get_zone_coincident_peak_demand_with_growth(am::Abstract_ALEAF_Model, zone_bus_ids, y::Int; nw::Int=0)
    demand_df = am.ref[:nw][nw][:time_series_data]["load"]
    region_map = get(am.ref[:nw][nw], :profile_data_region_map, nothing)
    hourly = nothing
    for bus_id in zone_bus_ids
        bus_data = am.ref[:nw][nw][:bus][bus_id]
        original_load = bus_data["aggregation_info"]["original_load_(bus_i, MW)"]
        for region_id in bus_data["aggregation_info"]["aggregated_regions_bus_i"]
            w = original_load[region_id] * get_load_growth_factor(am, y, region_id; nw)
            col = demand_df[!, Symbol(profile_data_region(region_map, "load", region_id))]
            hourly === nothing && (hourly = zeros(Float64, length(col)))
            hourly .+= w .* col
        end
    end
    return hourly === nothing ? 0.0 : maximum(hourly)
end


# System-wide coincident peak for the PRM report (annual growth basis), mirroring
# get_annual_peak_demand_with_growth but summing region shapes per hour before the max.
function get_annual_coincident_peak_demand_with_growth(network_data::Dict{Symbol,<:Any}, year)
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
    return hourly === nothing ? 0.0 : maximum(hourly)
end



function get_solution_value(result::Dict{String, Any}, key1::String, key2::String, idx::Int)
    return result["solution"][key1][string(idx)][key2]
end


function get_solution_value(result::Dict{String, Any}, key1::String, key2::String, idx::String)
    return result["solution"][key1][idx][key2]
end


function get_solution_value(result::Dict{String, Any}, solution_category::String, key1::String, key2::String, idx::String)
    return result[solution_category]["solution"][key1][idx][key2]
end


function add_sol_component(aim::Abstract_ALEAF_Model, nw::Int, comp_name::Symbol, field_name::Symbol, comp_ids, variables, day_id::Int)
    for i in comp_ids
        d = dict_pointer(aim.sol[:nw][nw][day_id], comp_name, i)   # resolve nested dict once
        @assert !haskey(d, field_name)
        d[field_name] = variables[i]
    end
end


function add_sol_component(aim::Abstract_ALEAF_Model, nw::Int, comp_name::Symbol, field_name::Symbol, comp_ids, variables)
    for i in comp_ids
        d = dict_pointer(aim.sol[:nw][nw], comp_name, i)   # resolve nested dict once
        @assert !haskey(d, field_name)
        d[field_name] = variables[i]
    end
end


function add_con_component(aim::Abstract_ALEAF_Model, nw::Int, day_id::Int, comp_name::Symbol, comp_ids, constraints)
    if !haskey(aim.con[:nw][nw][day_id], comp_name)
        aim.con[:nw][nw][day_id][comp_name] = Dict()
    end
    aim.con[:nw][nw][day_id][comp_name][comp_ids] = constraints
end


function add_con_component(aim::Abstract_ALEAF_Model, nw::Int, comp_name::Symbol, comp_ids, constraints)
    if !haskey(aim.con[:nw][nw], comp_name)
        aim.con[:nw][nw][comp_name] = Dict()
    end
    aim.con[:nw][nw][comp_name][comp_ids] = constraints
end


function dict_pointer(dict::Dict, args...)
    for arg in args
        if haskey(dict, arg)
            dict = dict[arg]
        else
            dict = dict[arg] = Dict()
        end
    end
    return dict
end


function collect_solution_values(var::Dict)
    sol = Dict{String,Any}()
    for (key, val) in var
        sol[string(key)] = collect_solution_values(val)
    end
    return sol
end


function collect_solution_values(var::Array{<:Any,1})
    return [collect_solution_values(val) for val in var]
end


function collect_solution_values(var::Array{<:Any,2})
    return [collect_solution_values(var[i,j]) for i in 1:size(var,1), j in 1:size(var,2)]
end


function collect_solution_values(var::Number)
    return var
end


function collect_solution_values(var)
    
    try
        return JuMP.value(var)
    catch
        @aleaf_warn "Collect_solution_values found unknown type"
        return var
    end
end


function get_atb_setting(ALEAF_setting::Dict{String,<:Any}, atb_setting_id::String, unit_group::String)
    atb_setting = Dict{String, Any}()
    for id in keys(ALEAF_setting["ATB Setting"])
        if (ALEAF_setting["ATB Setting"][id]["ATB_Setting_ID"] == atb_setting_id) & (ALEAF_setting["ATB Setting"][id]["UNITGROUP"] == unit_group)
            atb_setting["Case"] = ALEAF_setting["ATB Setting"][id]["Case"]
            atb_setting["CRP"] = ALEAF_setting["ATB Setting"][id]["CRP"]
            atb_setting["Tech"] = ALEAF_setting["ATB Setting"][id]["Tech"]
            atb_setting["TechDetail"] = ALEAF_setting["ATB Setting"][id]["TechDetail"]
            atb_setting["Scenario"] = ALEAF_setting["ATB Setting"][id]["Scenario"]
            atb_setting["ATB_Year"] = ALEAF_setting["ATB Setting"][id]["ATB Year"]
            atb_setting["CAPEX_Scale"] = ALEAF_setting["ATB Setting"][id]["CAPEX_Scale"]
        end
    end
    return atb_setting
end


function is_atb_2025_schema(ATB_data_raw)
    return ("Parameter" in names(ATB_data_raw)) && ("Technology" in names(ATB_data_raw)) && ("DisplayName" in names(ATB_data_raw))
end


function normalize_atb_token(value)
    if value isa Missing
        return ""
    end
    return lowercase(replace(string(value), r"[^A-Za-z0-9]" => ""))
end


function normalize_atb_numeric_value(value)
    if value isa Number
        return Float64(value)
    end
    if value isa Missing
        return 0.0
    end
    return parse(Float64, strip(string(value)))
end


function get_atb_value_2025(ATB_data_raw, atb_setting::Dict{String,<:Any}, atb_value_type::String; use_display_name::Bool=true)
    parameter_name = atb_value_type
    core_metric_case = atb_setting["Case"]
    crpyears = atb_setting["CRP"]
    technology = atb_setting["Tech"]
    # ATB 2025 DisplayName is the TechDetail verbatim; prefixing Technology never matches.
    display_name = atb_setting["TechDetail"]
    scenario = atb_setting["Scenario"]
    variable_year = atb_setting["Year"]

    if technology == "Nuclear"
        variable_year = max(2030, atb_setting["Year"])
    end

    atb_rows = ATB_data_raw[(ATB_data_raw.Parameter .== parameter_name) .&
                            (ATB_data_raw.Case .== core_metric_case) .&
                            (ATB_data_raw.CRPYears .== crpyears) .&
                            (ATB_data_raw.Technology .== technology) .&
                            (ATB_data_raw.Scenario .== scenario) .&
                            (ATB_data_raw.variable .== variable_year), :]

    if nrow(atb_rows) == 0
        return "NA"
    end

    # A DisplayName that matches nothing means the requested variant does not exist: report it
    # rather than silently falling back to whichever variant happens to sort first.
    if use_display_name
        normalized_display_name = normalize_atb_token(display_name)
        display_mask = [normalize_atb_token(row_display_name) == normalized_display_name for row_display_name in atb_rows.DisplayName]
    else
        display_mask = atb_rows.DisplayName .== "*"
    end
    atb_rows = atb_rows[display_mask, :]

    if nrow(atb_rows) == 0
        return "NA"
    end

    return normalize_atb_numeric_value(atb_rows.value[1])
end


function get_atb_value(ATB_data_raw, atb_setting::Dict{String,<:Any}, atb_value_type::String)
    if is_atb_2025_schema(ATB_data_raw)
        atb_value = get_atb_value_2025(ATB_data_raw, atb_setting, atb_value_type; use_display_name=true)
        if atb_value != "NA"
            return atb_value
        else
            @aleaf_warn "ATB value $atb_value_type not found for $(atb_setting["TechDetail"]), return 0.0"
            return 0.0
        end
    end

    core_metric_parameter = atb_value_type  
    atb_year = atb_setting["ATB_Year"]  
    core_metric_case = atb_setting["Case"] 
    crpyears = atb_setting["CRP"]
    technology = atb_setting["Tech"]
    techdetail = atb_setting["TechDetail"]
    scenario = atb_setting["Scenario"]
    core_metric_variable = atb_setting["Year"]

    # update core_metric_variable for Nuclear (min = 2030)
    if (technology == "Nuclear") && (core_metric_variable < 2030)
        core_metric_variable = 2030
    end

    atb_value = 0.0
    try
        atb_value = normalize_atb_numeric_value(ATB_data_raw[(ATB_data_raw.core_metric_parameter .== core_metric_parameter) .& 
        (ATB_data_raw.atb_year .== atb_year) .& (ATB_data_raw.core_metric_case .== core_metric_case) .& 
        (ATB_data_raw.crpyears .== crpyears) .& (ATB_data_raw.technology_alias .== technology) .& 
        (ATB_data_raw.techdetail .== techdetail) .& (ATB_data_raw.scenario .== scenario) .& 
        (ATB_data_raw.core_metric_variable .== core_metric_variable), :].value[1])
    catch
        atb_value = "NA"
    end

    if atb_value != "NA"
        return atb_value
    else
        @aleaf_warn "ATB value $atb_value_type not found (tech=$technology, techdetail=$techdetail, year=$core_metric_variable, case=$core_metric_case, crp=$crpyears, scenario=$scenario); returning 0.0"
        return 0.0
    end
end


function get_atb_value_FCR(ATB_data_raw, atb_setting::Dict{String,<:Any}, atb_value_type::String)
    if is_atb_2025_schema(ATB_data_raw)
        return get_atb_value_2025(ATB_data_raw, atb_setting, atb_value_type; use_display_name=false)
    end

    core_metric_parameter = atb_value_type
    atb_year = atb_setting["ATB_Year"]
    core_metric_case = atb_setting["Case"]
    crpyears = atb_setting["CRP"]
    technology = atb_setting["Tech"]
    techdetail = "*"
    scenario = atb_setting["Scenario"]
    core_metric_variable = atb_setting["Year"]

    atb_value = 0.0
    try
        atb_value = normalize_atb_numeric_value(ATB_data_raw[(ATB_data_raw.core_metric_parameter .== core_metric_parameter) .& (ATB_data_raw.atb_year .== atb_year) .& (ATB_data_raw.core_metric_case .== core_metric_case) .& (ATB_data_raw.crpyears .== crpyears) .& (ATB_data_raw.technology_alias .== technology) .& (ATB_data_raw.techdetail .== techdetail) .& (ATB_data_raw.scenario .== scenario) .& (ATB_data_raw.core_metric_variable .== core_metric_variable), :].value[1])
    catch
        atb_value = "NA"
    end

    return atb_value
end


function get_atb_value_WACC(ATB_data_raw, atb_setting::Dict{String,<:Any}, atb_value_type::String)
    if is_atb_2025_schema(ATB_data_raw)
        return get_atb_value_2025(ATB_data_raw, atb_setting, atb_value_type; use_display_name=false)
    end

    core_metric_parameter = atb_value_type
    atb_year = atb_setting["ATB_Year"]
    core_metric_case = atb_setting["Case"]
    crpyears = "*"
    technology = atb_setting["Tech"]
    techdetail = "*"
    scenario = atb_setting["Scenario"]
    core_metric_variable = atb_setting["Year"]

    atb_value = 0.0
    try
        atb_value = normalize_atb_numeric_value(ATB_data_raw[(ATB_data_raw.core_metric_parameter .== core_metric_parameter) .& (ATB_data_raw.atb_year .== atb_year) .& (ATB_data_raw.core_metric_case .== core_metric_case) .& (ATB_data_raw.crpyears .== crpyears) .& (ATB_data_raw.technology_alias .== technology) .& (ATB_data_raw.techdetail .== techdetail) .& (ATB_data_raw.scenario .== scenario) .& (ATB_data_raw.core_metric_variable .== core_metric_variable), :].value[1])
    catch
        atb_value = "NA"
    end

    return atb_value
end


function get_atb_value_CRF(ATB_data_raw, atb_setting::Dict{String,<:Any}, life)
    
    WACC = get_atb_value_WACC(ATB_data_raw, atb_setting, "WACC Real")
    
    if WACC != "NA"
        capital_recovery_factor = (WACC) / (1 - (1 + WACC)^(-life))
    else
        capital_recovery_factor = "NA"
    end

    return capital_recovery_factor
end


# Detect a single-round LCO_GTEP expansion-result record
# (produced by report_LCO_GTEP_result_EXP / export_LC_GTEP_result_EXP).
function is_single_round_expansion_record(d::Dict)
    return haskey(d, "expansion model result") && haskey(d, "expansion model system reference")
end


# Convert a single-round LCO_GTEP expansion-result record into the multi-round-info shape
# expected by RA loaders. Mirrors record_decision_year_solutions! from a multi-round run.
function convert_single_round_record_to_multi_round_info(record::Dict, round_id)
    expansion       = record["expansion model result"]["1"]["solution"]["expansion"]
    expansion_line  = record["expansion model result"]["1"]["solution"]["expansion_line"]

    # GTEP serializes u_ESE_iy in MWh, but RA expects PU in gen_bus["ES_MWh"];
    # convert back to PU using the record's own per-unit base.
    pu_power_base = record["setting"]["Simulation Setting"]["per_unit_base_value"]

    parse_iy(key) = begin
        m = match(r"\(\s*(\d+)\s*,\s*(\d+)\s*\)", key)
        (parse(Int, m.captures[1]), parse(Int, m.captures[2]))
    end

    U_G_i        = Dict{String, Float64}()
    U_G_i_year   = Dict{String, Int}()      # latest year stored per unit, so U_G_i tracks the decision-year stock
    U_NEW_G_i    = Dict{String, Float64}()
    U_RET_G_i    = Dict{String, Float64}()
    U_ESE_i      = Dict{String, Float64}()
    U_NEW_ESH_i  = Dict{String, Float64}()
    u_new_G_iy   = Dict{String, Float64}()
    u_ret_G_iy   = Dict{String, Float64}()
    u_new_ESH_iy = Dict{String, Float64}()
    U_NEW_T_k    = Dict{String, Float64}()

    years_seen = Set{Int}()

    for (key, vals) in expansion
        (i, y) = parse_iy(key)
        push!(years_seen, y)
        i_str = string(i)
        u_new = get(vals, "u_new_G_iy", 0.0)
        u_ret = get(vals, "u_ret_G_iy", 0.0)
        u_G   = get(vals, "u_G_iy", 0.0)

        u_new_G_iy[key] = u_new
        u_ret_G_iy[key] = u_ret

        if u_new > 0.001
            U_NEW_G_i[i_str] = get(U_NEW_G_i, i_str, 0.0) + u_new
        else
            U_NEW_G_i[i_str] = get(U_NEW_G_i, i_str, 0.0)
        end
        if u_ret > 0.001
            U_RET_G_i[i_str] = get(U_RET_G_i, i_str, 0.0) + u_ret
        else
            U_RET_G_i[i_str] = get(U_RET_G_i, i_str, 0.0)
        end
        # U_G_i is the decision-year stock (NOT accumulated). `expansion` is an unordered JSON Dict,
        # so keep the value from the latest year seen per unit rather than whichever record lands last.
        # Mirrors record_decision_year_solutions!, which reads u_G_iy at the decision year.
        if !haskey(U_G_i, i_str) || y >= U_G_i_year[i_str]
            U_G_i[i_str] = u_G
            U_G_i_year[i_str] = y
        end

        if haskey(vals, "u_ESE_iy")
            u_ESE = vals["u_ESE_iy"] / pu_power_base   # JSON saved in MWh; convert back to PU
            U_ESE_i[i_str] = u_ESE
        end
        if haskey(vals, "u_new_ESH_iy")
            u_new_ESH = vals["u_new_ESH_iy"]
            u_new_ESH_iy[key] = u_new_ESH
            if u_new_ESH > 0.001
                U_NEW_ESH_i[i_str] = get(U_NEW_ESH_i, i_str, 0.0) + u_new_ESH
            else
                U_NEW_ESH_i[i_str] = get(U_NEW_ESH_i, i_str, 0.0)
            end
        end
    end

    for (key, vals) in expansion_line
        (k, _y) = parse_iy(key)
        k_str = string(k)
        u_new_T = get(vals, "u_new_T_ky", 0.0)
        U_NEW_T_k[k_str] = get(U_NEW_T_k, k_str, 0.0) + u_new_T
    end

    last_year = isempty(years_seen) ? 1 : maximum(years_seen)
    first_year = isempty(years_seen) ? 1 : minimum(years_seen)

    updated_investment_decisions = Dict{String, Any}(
        "U_G_i"        => U_G_i,
        "U_NEW_G_i"    => U_NEW_G_i,
        "U_RET_G_i"    => U_RET_G_i,
        "U_ESE_i"      => U_ESE_i,
        "U_NEW_ESH_i"  => U_NEW_ESH_i,
        "u_new_G_iy"   => u_new_G_iy,
        "u_ret_G_iy"   => u_ret_G_iy,
        "u_new_ESH_iy" => u_new_ESH_iy,
        "U_NEW_T_k"    => U_NEW_T_k,
        # Aggregate VRE/storage capacity trackers (retained; the RA lookup-table capacity-credit path that consumed them was removed)
        "Wind_TotalMW"          => Dict{String, Float64}(),
        "PV_TotalMW"            => Dict{String, Float64}(),
        "RTPV_TotalMW"          => Dict{String, Float64}(),
        "4hr_Storage_TotalMW"   => Dict{String, Float64}(),
        "8hr_Storage_TotalMW"   => Dict{String, Float64}(),
        "10hr_Storage_TotalMW"  => Dict{String, Float64}(),
        "20hr_Storage_TotalMW"  => Dict{String, Float64}(),
    )

    return Dict{String, Any}(
        string(round_id) => Dict{String, Any}(
            "updated_investment_decisions" => updated_investment_decisions,
            "round_ids_y_decision"         => [last_year],
            "round_first_y"                => first_year,
        ),
        "status" => "completed",
    )
end

# Map each operating-reserve product (and its zone requirement / scarcity product) to the
# per-case ON/OFF flag in the `Simulation Configuration` sheet.
const _RESERVE_PRODUCT_FLAG = Dict{Symbol, String}(
    :reg_up  => "regulation_reserve_flag",
    :reg_dn  => "regulation_reserve_flag",
    :spin    => "spinning_reserve_flag",
    :flex_up => "flexibility_reserve_flag",
    :flex_dn => "flexibility_reserve_flag",
    :nonspin => "nonspin_reserve_flag",
)

"""
    reserve_enabled(am, product::Symbol) -> Bool

Return whether an operating-reserve `product` is enabled for the current case.

`product` is one of `:reg_up`, `:reg_dn`, `:spin`, `:flex_up`, `:flex_dn`,
`:nonspin`. Reads the matching `*_reserve_flag` from the per-case
`Simulation Configuration` settings. Missing keys default to `true`, so older
workbooks (without these flags) reproduce today's behaviour. Only meaningful on
the EXP/OP reserve-modeling path; RA never calls reserve provision builders or
zone requirement constraints.
"""
function reserve_enabled(am::Abstract_ALEAF_Model, product::Symbol)::Bool
    flag_key = get(_RESERVE_PRODUCT_FLAG, product, nothing)
    flag_key === nothing && return true   # unknown product -> behave as today
    sim_config = am.setting["Simulation Configuration"]
    haskey(sim_config, flag_key) || return true   # backward-compatible default
    val = sim_config[flag_key]
    if val isa Bool
        return val
    elseif val isa AbstractString
        return uppercase(strip(val)) in ("TRUE", "1", "YES")
    else
        return val == true
    end
end

"""
    up_reserve_eligible(am, i::Int; nw::Int=0) -> Bool

Return whether generation unit `i` is eligible to provide *up* operating-reserve
products (`:reg_up`, `:flex_up`, `:spin`).

Eligibility mirrors the long-standing fix-to-0 partition: only `Dispatchable`
units (per the `Dispatch` parameter on the unit's `gen_bus` technology row) can
carry up-reserves; `Non-dispatchable` units (VRE on a fixed profile) cannot.
This predicate is the single source of truth used for BOTH up-reserve variable
creation (skip-creation for non-eligible units) and every up-reserve variable
reference, so no reference can reach a variable that was never built.

Down-reserve products (`:reg_dn`, `:flex_dn`) and `:nonspin` have no eligibility
restriction (all units), matching current behaviour; do not gate them with this.
"""
function up_reserve_eligible(am::Abstract_ALEAF_Model, i::Int; nw::Int=0)::Bool
    # Memoize this static per-generator predicate (referenced per (i,d,h,t,y)). Caching a pure
    # predicate is transparent: it does not couple to the formulation, so constraint edits need no change here.
    cache = get!(() -> Dict{Int,Bool}(), am.ref[:nw][nw], :_up_elig_cache)::Dict{Int,Bool}
    return get!(cache, i) do
        bus_idx  = parameter(am, nw, :gen_index, "bus_idx", i)
        tech_idx = parameter(am, nw, :gen_index, "genco_tech_id", i)
        parameter(am, bus_idx, :gen_bus, tech_idx, "Dispatch") == "Dispatchable"
    end
end


# ref/var are Dict{Symbol,Any}; these generic readers type the RESULT (a function barrier) so caller
# arithmetic is type-stable, cutting build-time dynamic dispatch. No per-parameter enumeration — plain parameter()/variable() stay available.
param_f64(am::Abstract_ALEAF_Model, nw::Int, key1::Symbol, idx, key2::String)::Float64 = am.ref[:nw][nw][key1][idx][key2]
param_bool(am::Abstract_ALEAF_Model, nw::Int, key1::Symbol, idx, key2::String)::Bool = am.ref[:nw][nw][key1][idx][key2] == true

const _EMPTY_INT_VEC = Int[]

struct GTEPBusIncidence
    to_branches::Dict{Int,Vector{Int}}    # branches k with t_bus == n (inject to node)
    from_branches::Dict{Int,Vector{Int}}  # branches k with f_bus == n (withdraw from node)
    lfl_at_bus::Dict{Int,Vector{Int}}     # large flexible loads whose bus_idx == n
end

# Precompute per-bus branch/LFL incidence once instead of rescanning per LoadBalance call.
# Visit order matches the per-call scan (branch keys, model_flag-filtered) so term order — the model — is unchanged.
function gtep_bus_incidence(am::Abstract_ALEAF_Model)::GTEPBusIncidence
    cached = get(am.ref[:nw][0], :gtep_bus_incidence, nothing)
    cached === nothing || return cached::GTEPBusIncidence

    to_b   = Dict{Int,Vector{Int}}()
    from_b = Dict{Int,Vector{Int}}()
    for k in get_index(am, :branch, 0)
        parameter(am, 0, :branch, "model_flag", k) == true || continue
        push!(get!(to_b,   parameter(am, 0, :branch, k, "t_bus"), Int[]), k)
        push!(get!(from_b, parameter(am, 0, :branch, k, "f_bus"), Int[]), k)
    end
    lfl = Dict{Int,Vector{Int}}()
    for l in get_index(am, :demand, 0)
        push!(get!(lfl, parse(Int, parameter(am, 0, :demand, "bus_idx", l)), Int[]), l)
    end
    inc = GTEPBusIncidence(to_b, from_b, lfl)
    am.ref[:nw][0][:gtep_bus_incidence] = inc
    return inc
end


# --- Enhanced hybrid transmission expansion (B-theta) ------------------------------
# Split-flow McCormick relaxation of the bilinear expansion increment f_exp = b0·u·Δθ.
# Base flow keeps KVL at a fixed rate_a cap (so |Δθ| ≤ rate_a/b0 holds for free); the
# increment f_exp is fenced by the convex hull of the bilinear surface over u ∈ [0, ū].
# Exact at u = 0, at u = ū, and wherever the base line is at its thermal limit. LP-only.

"""
    hybrid_tx_enabled(am) -> Bool

True when the enhanced-hybrid transmission-expansion formulation is active: the
`transmission_expansion_hybrid_flag` (Simulation Setting sheet) is set AND the power-flow
mode is `B-theta` (the only mode with an angle law to relax). Network_Flow needs no such
relaxation and PTDF is a deferred follow-up, so both fall back to the standard formulation.
"""
function hybrid_tx_enabled(am::Abstract_ALEAF_Model)::Bool
    get(am.setting["Simulation Setting"], "transmission_expansion_hybrid_flag", false) == true &&
        am.setting["Simulation Setting"]["power_flow_mode_flag"] == "B-theta"
end

"""
    hybrid_exp_branch(am, k; nw=0) -> Bool

Per-corridor gate for the split-flow formulation: hybrid enabled, branch expandable
(`expansion_flag`), and AC (not a `dc_line` tie — an angle-decoupled tie is already exact
under a plain bounded transfer, so it keeps the standard formulation). Memoized: a pure
static predicate referenced per (k,d,h,t,y), transparent to formulation edits.
"""
function hybrid_exp_branch(am::Abstract_ALEAF_Model, k::Int; nw::Int=0)::Bool
    cache = get!(() -> Dict{Int,Bool}(), am.ref[:nw][nw], :_hybrid_exp_cache)::Dict{Int,Bool}
    return get!(cache, k) do
        hybrid_tx_enabled(am) &&
            parameter(am, nw, :branch, "expansion_flag", k) == true &&
            get(am.ref[:nw][nw][:branch][k], "dc_line", false) != true
    end
end

"""
    tx_expansion_ub(am, k; nw=0) -> Float64

Per-corridor expansion cap ū_k = max_rate_a/rate_a − 1, identical to the `u_T_ky` upper
bound set in `constraint_u_newT_ky_real` (including the `transmission_expansion_limit_value`
floor). This is the width of the McCormick box; sharing the derivation guarantees the
envelope box and the variable bound can never diverge.
"""
function tx_expansion_ub(am::Abstract_ALEAF_Model, k::Int; nw::Int=0)::Float64
    rate_a = parameter(am, nw, :branch, "rate_a", k)
    max_rate_a = parameter(am, nw, :branch, "max_rate_a", k)
    lim = am.setting["Planning Design"]["transmission_expansion_limit_value"]
    if rate_a * (1 + lim) >= max_rate_a
        max_rate_a = rate_a * (1 + lim)
    end
    return rate_a > 0.0 ? max_rate_a / rate_a - 1.0 : 0.0
end
