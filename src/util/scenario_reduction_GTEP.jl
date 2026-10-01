# Copyright © 2026, UChicago Argonne, LLC. All Rights Reserved.
# SPDX-License-Identifier: BSD-3-Clause  (see LICENSE)

# Scenario reduction for GTEP: reduce temporal/scenario sets for tractability.
# Core based on https://gitlab.com/supsi-dacd-isaac/scenred


using CSV
using DataFrames
using Statistics
using LinearAlgebra
using Distances
using Printf


function run_scenario_reduction(; kwargs...)
    setting = Dict{String, Any}(
        "time_resolution" => "Hourly",
        "num_scenarios_list" => [60],
        "fixing_extreme_days_flag" => true,
        "generate_input_data_flag" => true,
        "allow_repday_overlap_flag" => false,
        "type_of_data_set" => ["load_shape", "wind_shape", "solar_shape"],
        "extreme_day_set" => ["peak_net_demand", "peak_pv_generation", "peak_wind_generation"],
        "fixed_extreme_day_list" => Int[],
        "data_location_timeseries" => "",  # required: caller must override with the run's timeseries path
        "data_input_path" => "Data_Input",
        "data_output_path" => "Data_Output",
        "windCapacity" => 1.0,
        "solarCapacity" => 1.0,
        "peakDemand" => 1.0,
        "num_days_in_group" => "num_days_in_group",
        "print_output" => true,
        "file_name" => "scenario_reduction_input.csv",
        "in_memory_timeseries" => nothing,
    )

    for (key, value) in kwargs
        setting[String(key)] = value
    end

    if minimum(setting["num_scenarios_list"]) == 1
        if !setting["allow_repday_overlap_flag"]
            @aleaf_info "[ALEAF Scenario Reduction] Only one scenario requested. Automatically enabling 'allow_repday_overlap_flag' to ensure sufficient day groups."
            setting["allow_repday_overlap_flag"] = true
        end
    end

    if (setting["num_days_in_group"] isa Number) && (setting["num_days_in_group"] > 365)
        println("Warning: 'num_days_in_group' = $(setting["num_days_in_group"]) is too large. Adjusting to 365.")
        setting["num_days_in_group"] = 365
    end

    if !isdir(setting["data_input_path"])
        mkpath(setting["data_input_path"])
    end
    if !isdir(setting["data_output_path"])
        mkpath(setting["data_output_path"])
    end

    if setting["allow_repday_overlap_flag"] == true
        num_of_day_groups = 365 - setting["num_days_in_group"] + 1
    else
        num_of_day_groups = div(365, setting["num_days_in_group"])
    end

    setting["num_data_set"] = length(setting["type_of_data_set"])

    # Build the seven raw-scenario DataFrames in memory (no CSV round-trip), in the order
    # read_input_data_* used to produce.
    local load_shape, wind_shape, solar_shape, load_MWh, wind_MWh, solar_MWh, net_load_MWh
    if setting["generate_input_data_flag"] == true
        if setting["time_resolution"] == "Hourly"
            if setting["allow_repday_overlap_flag"] == true
                load_shape, wind_shape, solar_shape, load_MWh, wind_MWh, solar_MWh, net_load_MWh =
                    generate_input_data_hourly_with_overlap(
                        setting["data_input_path"],
                        setting["data_location_timeseries"],
                        setting["windCapacity"],
                        setting["solarCapacity"],
                        setting["peakDemand"],
                        num_of_day_groups,
                        setting["num_days_in_group"];
                        in_memory_timeseries = setting["in_memory_timeseries"],
                    )
            else
                load_shape, wind_shape, solar_shape, load_MWh, wind_MWh, solar_MWh, net_load_MWh =
                    generate_input_data_hourly(
                        setting["data_input_path"],
                        setting["data_location_timeseries"],
                        setting["windCapacity"],
                        setting["solarCapacity"],
                        setting["peakDemand"],
                        num_of_day_groups,
                        setting["num_days_in_group"];
                        in_memory_timeseries = setting["in_memory_timeseries"],
                    )
            end
        elseif setting["time_resolution"] == "Five-min"
            load_shape, wind_shape, solar_shape, load_MWh, wind_MWh, solar_MWh, net_load_MWh =
                generate_input_data_5min(
                    setting["data_input_path"],
                    setting["data_location_timeseries"],
                    setting["windCapacity"],
                    setting["solarCapacity"],
                    setting["peakDemand"];
                    in_memory_timeseries = setting["in_memory_timeseries"],
                )
        end
    else
        # Fallback when generation disabled: read pre-existing Input_Raw_Scenarios_* CSVs.
        if setting["time_resolution"] == "Hourly"
            load_shape, wind_shape, solar_shape, load_MWh, wind_MWh, solar_MWh, net_load_MWh = read_input_data_Hourly(setting["data_input_path"])
        elseif setting["time_resolution"] == "Five-min"
            load_shape, wind_shape, solar_shape, load_MWh, wind_MWh, solar_MWh, net_load_MWh = read_input_data_5min(setting["data_input_path"])
        end
    end

    extreme_scenarios, extreme_datapoint = identify_extreme_points(
        load_shape,
        wind_shape,
        solar_shape,
        load_MWh,
        wind_MWh,
        solar_MWh,
        net_load_MWh,
        setting["extreme_day_set"],
        setting["fixed_extreme_day_list"],
        setting["num_days_in_group"],
    )

    extreme_scenarios = unique(extreme_scenarios)

    load_shape_mat = Matrix(load_shape[1:end-1, 2:end])
    num_timestep = size(load_shape_mat, 1)
    num_total_scenario = size(load_shape_mat, 2)

    data = zeros(num_timestep, num_total_scenario, setting["num_data_set"])

    flag_load_shape = false
    flag_wind_shape = false
    flag_solar_shape = false
    flag_load_MWh = false
    flag_wind_MWh = false
    flag_solar_MWh = false
    flag_net_load_MWh = false

    for idx in 1:setting["num_data_set"]
        if ("load_shape" in setting["type_of_data_set"]) && !flag_load_shape
            data[:, :, idx] = load_shape_mat
            flag_load_shape = true
        elseif ("wind_shape" in setting["type_of_data_set"]) && !flag_wind_shape
            data[:, :, idx] = Matrix(wind_shape[1:end-1, 2:end])
            flag_wind_shape = true
        elseif ("solar_shape" in setting["type_of_data_set"]) && !flag_solar_shape
            data[:, :, idx] = Matrix(solar_shape[1:end-1, 2:end])
            flag_solar_shape = true
        elseif ("load_MWh" in setting["type_of_data_set"]) && !flag_load_MWh
            data[:, :, idx] = Matrix(load_MWh[1:end-1, 2:end])
            flag_load_MWh = true
        elseif ("wind_MWh" in setting["type_of_data_set"]) && !flag_wind_MWh
            data[:, :, idx] = Matrix(wind_MWh[1:end-1, 2:end])
            flag_wind_MWh = true
        elseif ("solar_MWh" in setting["type_of_data_set"]) && !flag_solar_MWh
            data[:, :, idx] = Matrix(solar_MWh[1:end-1, 2:end])
            flag_solar_MWh = true
        elseif ("net_load_MWh" in setting["type_of_data_set"]) && !flag_net_load_MWh
            data[:, :, idx] = Matrix(net_load_MWh[1:end-1, 2:end])
            flag_net_load_MWh = true
        end
    end

    idx_case = 1
    rep_day_input_df = DataFrame()

    for sce_list in setting["num_scenarios_list"]
        num_scenarios = sce_list
        num_days_in_group = setting["num_days_in_group"]
        if setting["print_output"] == true
            @aleaf_info "[ALEAF Scenario Reduction]: [$(num_scenarios) day groups, $(num_days_in_group) day(s) in each group]"            
        end
        idx_case += 1

        if setting["fixing_extreme_days_flag"] == true
            if num_scenarios < length(extreme_scenarios)
                if setting["print_output"] == true
                    println("The extreme cases will be ignored because the given number of scenarios to select is less than the number of extreme cases")
                end
                setting["fixing_extreme_days_flag"] = false
            end
        end

        extreme_set = Int[]
        for sce in extreme_scenarios
            col_idx = findfirst(==(sce), names(load_shape))
            if col_idx !== nothing
                push!(extreme_set, col_idx - 1)
            end
        end
        sort!(extreme_set)

        selected_scenarios, selected_scenarios_prob, scenarios_tree_structure = scenario_reduction_core(
            setting["fixing_extreme_days_flag"],
            extreme_set,
            copy(data);
            nodes=fill(num_scenarios, num_timestep),
        )

        s_prob = selected_scenarios_prob[1, :]
        sce_reduction_result = scenarios_tree_structure[1, :]

        selected_sce = Int[]
        all_sce_prob = Float64[]
        selected_sce_prob = Float64[]
        rep_day_input_data = Vector{Any}[]
        prob_idx = 0

        for i in 1:size(scenarios_tree_structure, 2)
            if sce_reduction_result[i] == true
                push!(selected_sce, i)
                push!(all_sce_prob, s_prob[prob_idx + 1])
                push!(selected_sce_prob, s_prob[prob_idx + 1])

                if setting["allow_repday_overlap_flag"] == true
                    start_day = i
                    end_day = i + setting["num_days_in_group"] - 1
                else
                    start_day = setting["num_days_in_group"] * i - setting["num_days_in_group"] + 1
                    end_day = setting["num_days_in_group"] * i
                end

                day_list = collect(start_day:end_day)
                push!(rep_day_input_data, Any[prob_idx + 1, i, s_prob[prob_idx + 1], start_day, end_day, day_list])
                prob_idx += 1
            else
                push!(all_sce_prob, 0.0)
            end
        end

        rep_day_input_df = DataFrame(
            index = [row[1] for row in rep_day_input_data],
            Day_Group = [row[2] for row in rep_day_input_data],
            Probability = [row[3] for row in rep_day_input_data],
            Start_Day_Id = [row[4] for row in rep_day_input_data],
            End_Day_Id = [row[5] for row in rep_day_input_data],
            Day_List = [row[6] for row in rep_day_input_data],
        )
        # Returned in memory (see get_repday_groups_data!); no CSV write.
    end

    return rep_day_input_df
end
function get_dist(X, metric::AbstractString)
    if metric == "euclidean"
        return pairwise(Euclidean(), X; dims=1)
    end
    error("Unsupported metric: $(metric)")
end


function scenario_reduction_core(fixing_extreme_days, extreme_set, samples; nodes=nothing, tol=10.0, metric="euclidean")
    T, n_obs, n_data = size(samples)

    if n_obs == 1
        J = trues(T, 1)
        P = ones(T, 1)
        return samples, P, J
    end

    if nodes === nothing
        nodes = ones(T)
    end

    X_list = Matrix{Float64}[]
    for i in 1:n_data
        V = samples[:, :, i]
        row_mean = mean(V, dims=2)
        row_std = std(V, dims=2) .+ 1e-6
        V_norm = (V .- row_mean) ./ row_std
        push!(X_list, V_norm)
    end

    X = reduce(vcat, X_list)
    X = permutedims(X)

    D = get_dist(X, metric)
    D[diagind(D)] .+= 1 + maximum(D)
    infty = 1e12

    if size(D, 1) >= 2
        D[2, :] .= infty
        D[:, 2] .= infty
    end

    default_nodes = ones(T)
    nodes_vec = vec(nodes)
    if all(nodes_vec .== default_nodes)
        exponents = min.(T .- collect(1:T) .+ 1, 300)
        Tol = reverse(tol ./ (1.5 .^ exponents))
        Tol[1] = infty
    else
        Tol = fill(infty, T)
    end

    J = trues(T, n_obs)
    L = zeros(n_obs, n_obs)
    P = ones(T, n_obs) ./ n_obs
    branches = n_obs

    for i in T:-1:1
        delta_rel = 0.0
        delta_p = 0.0
        D_i = D
        delta_rel_2 = 0.0
        delta_p_2 = 0.0

        basic_idx = vcat(trues(i), falses(T - i))
        sel_idx = repeat(basic_idx, n_data)
        X_filt = X[J[i, :], :]
        X_filt = X_filt[:, sel_idx]
        D_j = get_dist(X_filt, metric)
        D_j[diagind(D_j)] .= 0
        delta_max = minimum(vec(sum(D_j, dims=1)))

        while (delta_rel < Tol[i]) && (branches > nodes_vec[i])
            D_i[.!J[i, :], :] .= infty
            D_i[:, .!J[i, :]] .= infty
            d_s = sort(D_i, dims=1)
            z = vec(d_s[1, :]) .* vec(P[i, :])
            z[.!J[i, :]] .= infty
            if fixing_extreme_days == true
                z[extreme_set] .= infty
            end
            idx_rem = argmin(z)
            dp_min = minimum(z)

            idx_aug = argmin(view(D_i, :, idx_rem))
            J[i, idx_rem] = false
            P[i, idx_aug] = P[i, idx_rem] + P[i, idx_aug]
            P[i, idx_rem] = 0.0
            branches = sum(P[i, :] .> 0)
            L[idx_aug, idx_rem] = 1
            L[idx_aug, L[idx_rem, :] .> 0] .= 1
            samples[1:i, idx_rem, :] = samples[1:i, idx_aug, :]
            to_merge_idx = findall(L[idx_rem, :] .> 0)
            for j in to_merge_idx
                samples[1:i, j, :] = samples[1:i, idx_aug, :]
            end

            if Tol[i] != infty
                delta_p += dp_min
            end
            if delta_max != 0
                delta_rel = delta_p / delta_max
            else
                delta_rel = 0.0
            end

            delta_p_2 += dp_min
            if delta_max != 0
                delta_rel_2 = delta_p_2 / delta_max
            else
                delta_rel_2 = 0.0
            end
        end

        if i > 1
            J[i - 1, :] = J[i, :]
            P[i - 1, :] = P[i, :]
            D[.!J[i, :], .!J[i, :]] .= infty
        end
    end

    S = samples[:, J[end, :] .> 0, :]
    P = P[:, J[end, :] .> 0]

    return S, P, J
end
function read_input_data_5min(data_input_path)
    load_shape = CSV.read(joinpath(data_input_path, "Input_Raw_Scenarios_5min_load_Shape.csv"), DataFrame)
    wind_shape = CSV.read(joinpath(data_input_path, "Input_Raw_Scenarios_5min_wind_Shape.csv"), DataFrame)
    solar_shape = CSV.read(joinpath(data_input_path, "Input_Raw_Scenarios_5min_solar_Shape.csv"), DataFrame)

    load_MWh = CSV.read(joinpath(data_input_path, "Input_Raw_Scenarios_5min_load_MWh.csv"), DataFrame)
    wind_MWh = CSV.read(joinpath(data_input_path, "Input_Raw_Scenarios_5min_wind_MWh.csv"), DataFrame)
    solar_MWh = CSV.read(joinpath(data_input_path, "Input_Raw_Scenarios_5min_solar_MWh.csv"), DataFrame)

    net_load_MWh = CSV.read(joinpath(data_input_path, "Input_Raw_Scenarios_5min_netload.csv"), DataFrame)

    return load_shape, wind_shape, solar_shape, load_MWh, wind_MWh, solar_MWh, net_load_MWh
end


function read_input_data_Hourly(data_input_path)
    load_shape = CSV.read(joinpath(data_input_path, "Input_Raw_Scenarios_60min_load_Shape.csv"), DataFrame)
    wind_shape = CSV.read(joinpath(data_input_path, "Input_Raw_Scenarios_60min_wind_Shape.csv"), DataFrame)
    solar_shape = CSV.read(joinpath(data_input_path, "Input_Raw_Scenarios_60min_solar_Shape.csv"), DataFrame)

    load_MWh = CSV.read(joinpath(data_input_path, "Input_Raw_Scenarios_60min_load_MWh.csv"), DataFrame)
    wind_MWh = CSV.read(joinpath(data_input_path, "Input_Raw_Scenarios_60min_wind_MWh.csv"), DataFrame)
    solar_MWh = CSV.read(joinpath(data_input_path, "Input_Raw_Scenarios_60min_solar_MWh.csv"), DataFrame)

    net_load_MWh = CSV.read(joinpath(data_input_path, "Input_Raw_Scenarios_60min_netload.csv"), DataFrame)

    return load_shape, wind_shape, solar_shape, load_MWh, wind_MWh, solar_MWh, net_load_MWh
end


function init_raw_scenarios(nrows, ncols)
    df = DataFrame()
    time_col = Vector{Any}(undef, nrows)
    for i in 1:(nrows - 1)
        time_col[i] = i
    end
    time_col[nrows] = "Prob"
    df[!, "Time"] = time_col
    for i in 1:ncols
        df[!, "Scenario$(i)"] = Vector{Union{Missing, Float64}}(missing, nrows)
    end
    return df
end


function generate_input_data_5min(data_input_path, data_location_timeseries, windCapacity, solarCapacity, peakDemand; in_memory_timeseries::Union{Nothing,Dict{String,DataFrame}}=nothing)
    # 5-min variant scales shapes by capacity at read time, so it always sources the 5-min
    # CSVs (in_memory_timeseries ignored).
    time_series_base_path = isdir(joinpath(data_location_timeseries, "timeseries_data_files")) ? joinpath(data_location_timeseries, "timeseries_data_files") : data_location_timeseries
    timeSeriesParams_load = CSV.read(joinpath(time_series_base_path, "Load", "timeseries_load_5mins.csv"), DataFrame)
    timeSeriesParams_wind = CSV.read(joinpath(time_series_base_path, "WIND", "timeseries_wind_5mins.csv"), DataFrame)
    timeSeriesParams_solar = CSV.read(joinpath(time_series_base_path, "PV", "timeseries_pv_5mins.csv"), DataFrame)

    nrows = 24 * 12 + 1
    ncols = 365

    net_load_MWh = init_raw_scenarios(nrows, ncols)
    load_MWh = init_raw_scenarios(nrows, ncols)
    wind_MWh = init_raw_scenarios(nrows, ncols)
    solar_MWh = init_raw_scenarios(nrows, ncols)
    load_shape = init_raw_scenarios(nrows, ncols)
    wind_shape = init_raw_scenarios(nrows, ncols)
    solar_shape = init_raw_scenarios(nrows, ncols)

    j = 0
    k = 0
    for i in 1:ncols
        k = j + nrows - 1
        net_load_MWh[1:(nrows - 1), "Scenario$(i)"] = timeSeriesParams_load[(j + 1):k, 2] .* peakDemand .-
                                                     timeSeriesParams_wind[(j + 1):k, 2] .* windCapacity .-
                                                     timeSeriesParams_solar[(j + 1):k, 2] .* solarCapacity
        load_MWh[1:(nrows - 1), "Scenario$(i)"] = timeSeriesParams_load[(j + 1):k, 2] .* peakDemand
        wind_MWh[1:(nrows - 1), "Scenario$(i)"] = timeSeriesParams_wind[(j + 1):k, 2] .* windCapacity
        solar_MWh[1:(nrows - 1), "Scenario$(i)"] = timeSeriesParams_solar[(j + 1):k, 2] .* solarCapacity
        load_shape[1:(nrows - 1), "Scenario$(i)"] = timeSeriesParams_load[(j + 1):k, 2]
        wind_shape[1:(nrows - 1), "Scenario$(i)"] = timeSeriesParams_wind[(j + 1):k, 2]
        solar_shape[1:(nrows - 1), "Scenario$(i)"] = timeSeriesParams_solar[(j + 1):k, 2]
        j = k
    end

    return load_shape, wind_shape, solar_shape, load_MWh, wind_MWh, solar_MWh, net_load_MWh
end
function generate_input_data_hourly(data_input_path, data_location_timeseries, windCapacity, solarCapacity, peakDemand, num_of_day_groups, num_days_in_group; in_memory_timeseries::Union{Nothing,Dict{String,DataFrame}}=nothing)
    if in_memory_timeseries === nothing
        time_series_base_path = isdir(joinpath(data_location_timeseries, "timeseries_data_files")) ? joinpath(data_location_timeseries, "timeseries_data_files") : data_location_timeseries
        timeSeriesParams_load = CSV.read(joinpath(time_series_base_path, "Load", "timeseries_load_MW_hourly_scenario_reduction.csv"), DataFrame)
        timeSeriesParams_wind = CSV.read(joinpath(time_series_base_path, "WIND", "timeseries_wind_MW_hourly_scenario_reduction.csv"), DataFrame)
        timeSeriesParams_solar = CSV.read(joinpath(time_series_base_path, "PV", "timeseries_pv_MW_hourly_scenario_reduction.csv"), DataFrame)

        timeSeriesParams_load_shape = CSV.read(joinpath(time_series_base_path, "Load", "timeseries_load_Shape_hourly_scenario_reduction.csv"), DataFrame)
        timeSeriesParams_wind_shape = CSV.read(joinpath(time_series_base_path, "WIND", "timeseries_wind_Shape_hourly_scenario_reduction.csv"), DataFrame)
        timeSeriesParams_solar_shape = CSV.read(joinpath(time_series_base_path, "PV", "timeseries_pv_Shape_hourly_scenario_reduction.csv"), DataFrame)
    else
        timeSeriesParams_load = in_memory_timeseries["load_MW"]
        timeSeriesParams_wind = in_memory_timeseries["wind_MW"]
        timeSeriesParams_solar = in_memory_timeseries["solar_MW"]
        timeSeriesParams_load_shape = in_memory_timeseries["load_shape"]
        timeSeriesParams_wind_shape = in_memory_timeseries["wind_shape"]
        timeSeriesParams_solar_shape = in_memory_timeseries["solar_shape"]
    end

    nrows = 24 * num_days_in_group + 1
    ncols = num_of_day_groups

    net_load_MWh = init_raw_scenarios(nrows, ncols)
    load_MWh = init_raw_scenarios(nrows, ncols)
    wind_MWh = init_raw_scenarios(nrows, ncols)
    solar_MWh = init_raw_scenarios(nrows, ncols)
    load_shape = init_raw_scenarios(nrows, ncols)
    wind_shape = init_raw_scenarios(nrows, ncols)
    solar_shape = init_raw_scenarios(nrows, ncols)

    j = 0
    k = 0
    for i in 1:ncols
        k = j + nrows - 1
        net_load_MWh[1:(nrows - 1), "Scenario$(i)"] = timeSeriesParams_load[(j + 1):k, 1] .-
                                                     timeSeriesParams_wind[(j + 1):k, 1] .-
                                                     timeSeriesParams_solar[(j + 1):k, 1]
        load_MWh[1:(nrows - 1), "Scenario$(i)"] = timeSeriesParams_load[(j + 1):k, 1]
        wind_MWh[1:(nrows - 1), "Scenario$(i)"] = timeSeriesParams_wind[(j + 1):k, 1]
        solar_MWh[1:(nrows - 1), "Scenario$(i)"] = timeSeriesParams_solar[(j + 1):k, 1]
        load_shape[1:(nrows - 1), "Scenario$(i)"] = timeSeriesParams_load_shape[(j + 1):k, 1]
        wind_shape[1:(nrows - 1), "Scenario$(i)"] = timeSeriesParams_wind_shape[(j + 1):k, 1]
        solar_shape[1:(nrows - 1), "Scenario$(i)"] = timeSeriesParams_solar_shape[(j + 1):k, 1]
        j = k
    end

    return load_shape, wind_shape, solar_shape, load_MWh, wind_MWh, solar_MWh, net_load_MWh
end


function generate_input_data_hourly_with_overlap(data_input_path, data_location_timeseries, windCapacity, solarCapacity, peakDemand, num_of_day_groups, num_days_in_group; in_memory_timeseries::Union{Nothing,Dict{String,DataFrame}}=nothing)
    if in_memory_timeseries === nothing
        time_series_base_path = isdir(joinpath(data_location_timeseries, "timeseries_data_files")) ? joinpath(data_location_timeseries, "timeseries_data_files") : data_location_timeseries
        timeSeriesParams_load = CSV.read(joinpath(time_series_base_path, "Load", "timeseries_load_MW_hourly_scenario_reduction.csv"), DataFrame)
        timeSeriesParams_wind = CSV.read(joinpath(time_series_base_path, "WIND", "timeseries_wind_MW_hourly_scenario_reduction.csv"), DataFrame)
        timeSeriesParams_solar = CSV.read(joinpath(time_series_base_path, "PV", "timeseries_pv_MW_hourly_scenario_reduction.csv"), DataFrame)

        timeSeriesParams_load_shape = CSV.read(joinpath(time_series_base_path, "Load", "timeseries_load_Shape_hourly_scenario_reduction.csv"), DataFrame)
        timeSeriesParams_wind_shape = CSV.read(joinpath(time_series_base_path, "WIND", "timeseries_wind_Shape_hourly_scenario_reduction.csv"), DataFrame)
        timeSeriesParams_solar_shape = CSV.read(joinpath(time_series_base_path, "PV", "timeseries_pv_Shape_hourly_scenario_reduction.csv"), DataFrame)
    else
        timeSeriesParams_load = in_memory_timeseries["load_MW"]
        timeSeriesParams_wind = in_memory_timeseries["wind_MW"]
        timeSeriesParams_solar = in_memory_timeseries["solar_MW"]
        timeSeriesParams_load_shape = in_memory_timeseries["load_shape"]
        timeSeriesParams_wind_shape = in_memory_timeseries["wind_shape"]
        timeSeriesParams_solar_shape = in_memory_timeseries["solar_shape"]
    end

    nrows = 24 * num_days_in_group + 1
    ncols = 365 - num_days_in_group + 1

    net_load_MWh = init_raw_scenarios(nrows, ncols)
    load_MWh = init_raw_scenarios(nrows, ncols)
    wind_MWh = init_raw_scenarios(nrows, ncols)
    solar_MWh = init_raw_scenarios(nrows, ncols)
    load_shape = init_raw_scenarios(nrows, ncols)
    wind_shape = init_raw_scenarios(nrows, ncols)
    solar_shape = init_raw_scenarios(nrows, ncols)

    j = 0
    k = 0
    for i in 1:ncols
        k = j + nrows - 1
        net_load_MWh[1:(nrows - 1), "Scenario$(i)"] = timeSeriesParams_load[(j + 1):k, 1] .-
                                                     timeSeriesParams_wind[(j + 1):k, 1] .-
                                                     timeSeriesParams_solar[(j + 1):k, 1]
        load_MWh[1:(nrows - 1), "Scenario$(i)"] = timeSeriesParams_load[(j + 1):k, 1]
        wind_MWh[1:(nrows - 1), "Scenario$(i)"] = timeSeriesParams_wind[(j + 1):k, 1]
        solar_MWh[1:(nrows - 1), "Scenario$(i)"] = timeSeriesParams_solar[(j + 1):k, 1]
        load_shape[1:(nrows - 1), "Scenario$(i)"] = timeSeriesParams_load_shape[(j + 1):k, 1]
        wind_shape[1:(nrows - 1), "Scenario$(i)"] = timeSeriesParams_wind_shape[(j + 1):k, 1]
        solar_shape[1:(nrows - 1), "Scenario$(i)"] = timeSeriesParams_solar_shape[(j + 1):k, 1]
        j += 24
    end

    return load_shape, wind_shape, solar_shape, load_MWh, wind_MWh, solar_MWh, net_load_MWh
end
function identify_extreme_points(load, wind, solar, load_MWh, wind_MWh, solar_MWh, net_load_MWh, extreme_day_set, fixed_extreme_day_list, num_days_in_group; kwargs...)
    extreme_scenarios = String[]
    extreme_datapoint = Dict{Int, Dict{String, Any}}()
    idx = 0

    for day in fixed_extreme_day_list
        scenario_id = div((day - 1), num_days_in_group) + 1
        push!(extreme_scenarios, "Scenario$(scenario_id)")
        extreme_datapoint[idx] = Dict("type" => "fixed extreme day", "value" => 0, "scenario" => day)
        idx += 1
    end

    if "peak_demand" in extreme_day_set
        data = Matrix(load_MWh[1:end-1, 2:end])
        value = maximum(data)
        col_max = vec(maximum(data, dims=1))
        sce_idx = argmax(col_max)
        sce_id = names(load_MWh)[sce_idx + 1]
        push!(extreme_scenarios, sce_id)
        extreme_datapoint[idx] = Dict("type" => "peak demand", "value" => value, "scenario" => sce_id)
        idx += 1
    end

    if "peak_net_demand" in extreme_day_set
        data = Matrix(net_load_MWh[1:end-1, 2:end])
        value = maximum(data)
        col_max = vec(maximum(data, dims=1))
        sce_idx = argmax(col_max)
        sce_id = names(net_load_MWh)[sce_idx + 1]
        push!(extreme_scenarios, sce_id)
        extreme_datapoint[idx] = Dict("type" => "peak net demand", "value" => value, "scenario" => sce_id)
        idx += 1
    end

    if "peak_pv_generation" in extreme_day_set
        data = Matrix(solar[1:end-1, 2:end])
        value = maximum(sum(data, dims=1))
        col_sum = vec(sum(data, dims=1))
        sce_idx = argmax(col_sum)
        sce_id = names(solar)[sce_idx + 1]
        push!(extreme_scenarios, sce_id)
        extreme_datapoint[idx] = Dict("type" => "peak daily solar shape", "value" => value, "scenario" => sce_id)
        idx += 1
    end

    if "peak_wind_generation" in extreme_day_set
        data = Matrix(wind[1:end-1, 2:end])
        value = maximum(sum(data, dims=1))
        col_sum = vec(sum(data, dims=1))
        sce_idx = argmax(col_sum)
        sce_id = names(wind)[sce_idx + 1]
        push!(extreme_scenarios, sce_id)
        extreme_datapoint[idx] = Dict("type" => "peak daily wind shape", "value" => value, "scenario" => sce_id)
        idx += 1
    end

    if "least_pv_generation" in extreme_day_set
        data = Matrix(solar[1:end-1, 2:end])
        value = minimum(sum(data, dims=1))
        col_sum = vec(sum(data, dims=1))
        sce_idx = argmin(col_sum)
        sce_id = names(solar)[sce_idx + 1]
        push!(extreme_scenarios, sce_id)
        extreme_datapoint[idx] = Dict("type" => "lowest daily solar shape", "value" => value, "scenario" => sce_id)
        idx += 1
    end

    if "least_wind_generation" in extreme_day_set
        data = Matrix(wind[1:end-1, 2:end])
        value = minimum(sum(data, dims=1))
        col_sum = vec(sum(data, dims=1))
        sce_idx = argmin(col_sum)
        sce_id = names(wind)[sce_idx + 1]
        push!(extreme_scenarios, sce_id)
        extreme_datapoint[idx] = Dict("type" => "lowest daily wind shape", "value" => value, "scenario" => sce_id)
        idx += 1
    end

    return extreme_scenarios, extreme_datapoint
end


