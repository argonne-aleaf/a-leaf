# Copyright © 2026, UChicago Argonne, LLC. All Rights Reserved.
# SPDX-License-Identifier: BSD-3-Clause  (see LICENSE)

# Purpose: Outage scenario generation support utilities.
# Scope: Creates, samples, and prepares outage scenario data.
# Notes: Shared by RA workflows and outage cache generation.

using CSV
using DataFrames
using Statistics
using Random
using Distributions
using Distributed

function generate_outage_scenarios(am_ref, am_setting, simulation_list, reference_temp, num_risk_scenario, RA_setting; output_option="Status")

    # Separate scenarios to respect the open-file-descriptor limit.
    if RA_setting["distributed_run_flag"] == false
        outage_scenario_report = Dict{Int, Any}()
        risk_scenario_list = [i for i in 1:num_risk_scenario]
        outage_scenario_report, filtered_outage_scenario_list_report = calculate_and_filter_outage_rate(am_ref, am_setting, simulation_list, reference_temp, risk_scenario_list, RA_setting, 0; output_option)

    else    # distributed_run_flag == true

        @aleaf_info "[ALEAF RA Model]: Perfrom Outage Sampling Using $(length(workers())) Workers"
    
        outage_scenario_report = Dict{Int64, Any}()
        temp_outage_scenario_report = Dict{Int64, Any}()
        temp_filtered_outage_scenario_list = Dict{Int64, Any}()
        total_risk_scenario_list = [i for i in 1:num_risk_scenario]

        num_risk_scenario_segment = length(workers())
        # Ceiling division so segments fully cover num_risk_scenario in the parallel
        # @sync block (last segment may be partial; empty segments are skipped).
        num_scenario_for_single_run = max(1, cld(num_risk_scenario, num_risk_scenario_segment))

        for seg_idx in 1:num_risk_scenario_segment
            temp_outage_scenario_report[seg_idx] = Dict{Int64, Any}()
            temp_filtered_outage_scenario_list[seg_idx] = Dict{Int64, Any}()
        end

        np = nprocs()
        i = 1
        nextidx() = (idx=i; i+=1; idx)

        @sync begin
            for p in workers()
                if p != myid() || np == 1
                    @async begin
                        while true
                            seg_idx = nextidx()
                            if seg_idx > num_risk_scenario_segment
                                break
                            end

                            lo = num_scenario_for_single_run*(seg_idx-1) + 1
                            hi = min(num_scenario_for_single_run*seg_idx, num_risk_scenario)
                            if lo > hi
                                continue  # more workers than scenarios
                            end
                            risk_scenario_list = [i for i in lo:hi]

                            temp_outage_scenario_report[seg_idx], temp_filtered_outage_scenario_list[seg_idx] = remotecall_fetch(calculate_and_filter_outage_rate, p, am_ref, am_setting, simulation_list, reference_temp, risk_scenario_list, RA_setting, seg_idx; output_option)
                        end
                    end
                end
            end
        end

        outage_scenario_report = Dict{Int64, Any}()
        for day_group_id in keys(am_ref[:nw][0][:repday_groups])
            outage_scenario_report[day_group_id] = Dict{Int64, Any}()
            for risk_scenario_id in 1:num_risk_scenario
                outage_scenario_report[day_group_id][risk_scenario_id] = []
            end
        end

        # Serial merge (NOT @spawn): under `using Distributed`, bare @spawn would
        # serialize outage_scenario_report to a worker and lose mutations on master.
        for (seg_idx, seg_data) in temp_outage_scenario_report
            for (day_group_id, day_data) in seg_data
                if !haskey(outage_scenario_report, day_group_id)
                    continue
                end
                target = outage_scenario_report[day_group_id]
                for (risk_scenario_id, val) in day_data
                    target[risk_scenario_id] = val
                end
            end
        end

        # Sanity check: nonzero count = the merge dropped writes (regression of the @spawn bug above).
        unpopulated_pairs = 0
        for (day_group_id, day_map) in outage_scenario_report
            for (risk_scenario_id, val) in day_map
                if val isa AbstractVector && isempty(val)
                    unpopulated_pairs += 1
                end
            end
        end
        if unpopulated_pairs > 0
            @aleaf_warn "[ALEAF RA Model]: outage merge left $(unpopulated_pairs) (day_group, risk_scenario) pairs unpopulated; check upstream worker chunking."
        end

        filtered_outage_scenario_list_report = Dict{Int64, Any}()
        for day_group_id in keys(am_ref[:nw][0][:repday_groups])
            filtered_outage_scenario_list_report[day_group_id] = []
            for seg_idx in keys(temp_outage_scenario_report)
                filtered_outage_scenario_list_report[day_group_id] = unique(vcat(filtered_outage_scenario_list_report[day_group_id], temp_filtered_outage_scenario_list[seg_idx][day_group_id]), dims=1)
            end
        end
    end

    return outage_scenario_report, filtered_outage_scenario_list_report
end

function generate_outage_samples(am_ref, am_setting, simulation_list, reference_temp, num_risk_scenario, RA_setting; output_option="Status")

    if RA_setting["distributed_run_flag"] == false
        outage_scenario_report = Dict{Int, Any}()
        risk_scenario_list = [i for i in 1:num_risk_scenario]
        outage_scenario_report, _ = calculate_and_filter_outage_rate(am_ref, am_setting, simulation_list, reference_temp, risk_scenario_list, RA_setting, 0; output_option, apply_filter=false)
    else
        @aleaf_info "[ALEAF RA Model]: Perfrom Outage Sampling Using $(length(workers())) Workers"

        outage_scenario_report = Dict{Int64, Any}()
        temp_outage_scenario_report = Dict{Int64, Any}()
        total_risk_scenario_list = [i for i in 1:num_risk_scenario]

        num_risk_scenario_segment = length(workers())
        # See generate_outage_scenarios: ceiling division keeps all chunks in the
        # @sync block; empty chunks (more workers than scenarios) are skipped.
        num_scenario_for_single_run = max(1, cld(num_risk_scenario, num_risk_scenario_segment))

        for seg_idx in 1:num_risk_scenario_segment
            temp_outage_scenario_report[seg_idx] = Dict{Int64, Any}()
        end

        np = nprocs()
        i = 1
        nextidx() = (idx=i; i+=1; idx)

        @sync begin
            for p in workers()
                if p != myid() || np == 1
                    @async begin
                        while true
                            seg_idx = nextidx()
                            if seg_idx > num_risk_scenario_segment
                                break
                            end

                            lo = num_scenario_for_single_run*(seg_idx-1) + 1
                            hi = min(num_scenario_for_single_run*seg_idx, num_risk_scenario)
                            if lo > hi
                                continue
                            end
                            risk_scenario_list = [i for i in lo:hi]

                            temp_outage_scenario_report[seg_idx], _ = remotecall_fetch(calculate_and_filter_outage_rate, p, am_ref, am_setting, simulation_list, reference_temp, risk_scenario_list, RA_setting, seg_idx; output_option, apply_filter=false)
                        end
                    end
                end
            end
        end

        outage_scenario_report = Dict{Int64, Any}()
        for day_group_id in keys(am_ref[:nw][0][:repday_groups])
            outage_scenario_report[day_group_id] = Dict{Int64, Any}()
            for risk_scenario_id in 1:num_risk_scenario
                outage_scenario_report[day_group_id][risk_scenario_id] = []
            end
        end

        # Serial merge — see generate_outage_scenarios for rationale (the bare @spawn
        # pattern is both expensive and semantically unsafe under `using Distributed`).
        for (seg_idx, seg_data) in temp_outage_scenario_report
            for (day_group_id, day_data) in seg_data
                if !haskey(outage_scenario_report, day_group_id)
                    continue
                end
                target = outage_scenario_report[day_group_id]
                for (risk_scenario_id, val) in day_data
                    target[risk_scenario_id] = val
                end
            end
        end

        unpopulated_pairs = 0
        for (day_group_id, day_map) in outage_scenario_report
            for (risk_scenario_id, val) in day_map
                if val isa AbstractVector && isempty(val)
                    unpopulated_pairs += 1
                end
            end
        end
        if unpopulated_pairs > 0
            @aleaf_warn "[ALEAF RA Model]: outage sample merge left $(unpopulated_pairs) (day_group, risk_scenario) pairs unpopulated; check upstream worker chunking."
        end
    end

    return outage_scenario_report
end


function calculate_and_filter_outage_rate(am_ref, am_setting, simulation_list, reference_temp, risk_scenario_list, RA_setting, seed_value=1; output_option="Status", apply_filter::Bool=true)
    
    num_risk_scenario = length(risk_scenario_list)
    # Geometric repair model: P(T <= k) = 1 - (1 - p_repair)^k, with p_repair = 1/MTTR (hours)
    repair_quantile_from_mttr(mttr_hours, q) = begin
        mttr_hours = max(mttr_hours, 1.0)
        p_repair = min(1.0, 1 / mttr_hours)
        if p_repair >= 1.0
            1
        else
            max(1, Int(ceil(log(1 - q) / log(1 - p_repair))))
        end
    end
    
    outage_scenario = Dict{Int64, Any}()
    filtered_outage_scenario_list = Dict{Int64, Any}()
    for day_group_id in keys(am_ref[:nw][0][:repday_groups])
        outage_scenario[day_group_id] = Dict{Int64, Any}()
        filtered_outage_scenario_list[day_group_id] = []
    end

    regressionAD = am_ref[:nw][0][:regressionAD]
    regressionDD = am_ref[:nw][0][:regressionDD]

    num_gens = length(am_ref[:nw][0][:RA_gen_index])

    # Warm-start over the trailing half of the horizon for better initial gen status.
    num_warm_start_days = Int(floor(length(simulation_list) / 2))
    new_simulation_list = simulation_list[end-num_warm_start_days:end]

    for day_group_idx in eachindex(simulation_list)

        day_group_id = simulation_list[day_group_idx]
        start_day = am_ref[:nw][0][:repday_groups][day_group_id]["Start_Day_Id"]
        end_day = am_ref[:nw][0][:repday_groups][day_group_id]["End_Day_Id"]
        num_hours = (end_day - start_day + 1) * 24

        outage_scenario[day_group_id] = Dict{Int64, Any}()
        for sce in risk_scenario_list
            outage_scenario[day_group_id][sce] = trues(num_hours, num_gens)
        end
    end

    temperature_cache = Dict{Int, Matrix{Float32}}()
    get_local_average_temperature(day_group_id) = get!(temperature_cache, day_group_id) do
        start_day = am_ref[:nw][0][:repday_groups][day_group_id]["Start_Day_Id"]
        end_day = am_ref[:nw][0][:repday_groups][day_group_id]["End_Day_Id"]
        day_list = [i for i in start_day:end_day]
        num_hours = (end_day - start_day + 1) * 24

        temperatureData = am_ref[:nw][0][:temperatureData]
        local_average_temperature = zeros(Float32, num_hours, length(am_ref[:nw][0][:bus]))
        for bus_idx in keys(am_ref[:nw][0][:bus])
            local_temperature = filter(:Day_ID => in(day_list), temperatureData)[:, am_ref[:nw][0][:bus][bus_idx]["aggregation_info"]["aggregated_regions_bus_i"]]
            local_average_temperature[:, bus_idx] = [mean(c) for c in eachrow(local_temperature)]
        end
        local_average_temperature
    end

    for day_group_idx in eachindex(new_simulation_list)

        day_group_id = new_simulation_list[day_group_idx]
        start_day = am_ref[:nw][0][:repday_groups][day_group_id]["Start_Day_Id"]
        end_day = am_ref[:nw][0][:repday_groups][day_group_id]["End_Day_Id"]
        day_list = [i for i in start_day:end_day]
        num_hours = (end_day - start_day + 1) * 24

        local_average_temperature = get_local_average_temperature(day_group_id)

        gen_indices = collect(keys(am_ref[:nw][0][:RA_gen_index]))
        Threads.@threads for i in eachindex(gen_indices)
            gen_idx = gen_indices[i]
            rng = Random.Xoshiro(seed_value + gen_idx + day_group_id * 100000)

            UnitprobabilityType = am_ref[:nw][0][:RA_gen_index][gen_idx]["RA_FOR"]
            bus_idx = am_ref[:nw][0][:RA_gen_index][gen_idx]["bus_idx"]
            ICAP = am_ref[:nw][0][:RA_gen_index][gen_idx]["ICAP"]

            rateAD = zeros(1, num_hours)
            rateDD = zeros(1, num_hours)

                if typeof(UnitprobabilityType) == String

                    UnitprobabilityTypeAD = string(UnitprobabilityType, "AD")
                    UnitprobabilityTypeDD = string(UnitprobabilityType, "DD")

                    for t in [t for t in 1:num_hours]

                        if local_average_temperature[t, bus_idx] > reference_temp
                            # regression cols AD=[Hi,Lo,tmpH,tmpH2,tmpL,tmpL2]; DD appends [rt_avg,rt_95,rt_90]
                            tmp = local_average_temperature[t, bus_idx] - reference_temp
                            rateAD[t] = 1 / (1 + exp(-1 * (regressionAD[:, UnitprobabilityTypeAD][1] + regressionAD[:, UnitprobabilityTypeAD][3] * tmp + regressionAD[:, UnitprobabilityTypeAD][4] * tmp^2)))
                            rateDD[t] = 1 / (1 + exp(-1 * (regressionDD[:, UnitprobabilityTypeDD][1] + regressionDD[:, UnitprobabilityTypeDD][3] * tmp + regressionDD[:, UnitprobabilityTypeDD][4] * tmp^2)))
                        else
                            tmp = reference_temp - local_average_temperature[t, bus_idx]
                            rateAD[t] = 1 / (1 + exp(-1 * (regressionAD[:, UnitprobabilityTypeAD][2] + regressionAD[:, UnitprobabilityTypeAD][5] * tmp + regressionAD[:, UnitprobabilityTypeAD][6] * tmp^2)))
                            rateDD[t] = 1 / (1 + exp(-1 * (regressionDD[:, UnitprobabilityTypeDD][2] + regressionDD[:, UnitprobabilityTypeDD][5] * tmp + regressionDD[:, UnitprobabilityTypeDD][6] * tmp^2)))
                        end

                    end
                else
                    # Numeric UnitprobabilityType is FOR. Derive AD/DD from FOR + MTTR.
                    # MTTR must be in hours.
                    mttr_hours = get(RA_setting, "repair_time_average_hours", 32)
                    mttr_hours = max(mttr_hours, 1.0) # avoid divide-by-zero or invalid probabilities
                    
                    for t in 1:num_hours
                        # rateDD = P(stay down) per hour, geometric repair
                        rateDD[t] = 1 - (1 / mttr_hours)
                        # rateAD = P(fail) per hour derived from FOR
                        for_val = UnitprobabilityType
                        rateAD[t] = for_val / (mttr_hours * max(1e-6, (1 - for_val)))
                        rateAD[t] = clamp(rateAD[t], 1e-6, 1 - 1e-6)
                        rateDD[t] = clamp(rateDD[t], 1e-6, 1 - 1e-6)
                    end
                end

                max_repair_hours = 10000   # default; no bound
                if typeof(UnitprobabilityType) == String

                    UnitprobabilityTypeDD = string(UnitprobabilityType, "DD")
                    
                    if RA_setting["repair_time_bound"] == "average"
                        max_repair_hours = regressionDD[:, UnitprobabilityTypeDD][7]    #repair_time_average
                    elseif RA_setting["repair_time_bound"] == "95 percentile"
                        max_repair_hours = regressionDD[:, UnitprobabilityTypeDD][8]    #repair_time_95
                    elseif RA_setting["repair_time_bound"] == "90 percentile"
                        max_repair_hours = regressionDD[:, UnitprobabilityTypeDD][9]    #repair_time_90
                    end
                else
                    mttr_hours = get(RA_setting, "repair_time_average_hours", 32)
                    mttr_hours = max(mttr_hours, 1.0)
                    if RA_setting["repair_time_bound"] == "average"
                        max_repair_hours = Int(ceil(mttr_hours))
                    elseif RA_setting["repair_time_bound"] == "95 percentile"
                        max_repair_hours = repair_quantile_from_mttr(mttr_hours, 0.95)
                    elseif RA_setting["repair_time_bound"] == "90 percentile"
                        max_repair_hours = repair_quantile_from_mttr(mttr_hours, 0.90)
                    end
                end

            @inbounds for sce_id in risk_scenario_list

                unit_status = outage_scenario[day_group_id][sce_id][:, gen_idx]

                    previous_status = unit_status[1]
                    current_status = true
                    @inbounds for t in 1:num_hours

                        # carry outage state across day-group boundary (t==1 uses prior day's last hour)
                        if (day_group_id > 1) && (t == 1)
                            if outage_scenario[day_group_id - 1][sce_id][num_hours, gen_idx] == true  # if it was on, check AD probability
                                rand(rng, Bernoulli(rateAD[1])) == 1 ? current_status = false : current_status = true
                            end
                        else
                            if previous_status == true  # if currently on, check AD probability
                                rand(rng, Bernoulli(rateAD[t])) == 1 ? current_status = false : current_status = true
                            end
                        end

                        if (previous_status == true) && (current_status == false)  # if unit is off, check when it will be on again
                            repair_time = min(rand(rng, Geometric(1 - rateDD[t])), max_repair_hours)
                            repair_time = max(1, Int(ceil(repair_time)))

                            if (num_hours - t + 1) >= repair_time   
                                for t2 in t:(t+repair_time-1)
                                    unit_status[t2] = false
                                end
                            else

                                # fill the current day first
                                for t2 in t:num_hours
                                    unit_status[t2] = false
                                end

                                # check remaining repair hours
                                remaining_repair_hours = (repair_time - (num_hours - t + 1))
                                num_remaining_repair_days = div(remaining_repair_hours, num_hours)

                                if num_remaining_repair_days != 0
                                    
                                    for next_day_count in 1:num_remaining_repair_days
                                        if (day_group_id + next_day_count) <= last(new_simulation_list)
                                            outage_scenario[day_group_id + next_day_count][sce_id][:, gen_idx] = falses(num_hours)
                                        end
                                    end
                                end
                                
                                # fill the last day
                                final_remaining_repair_hours = remaining_repair_hours - num_remaining_repair_days * num_hours
                                if (day_group_id + num_remaining_repair_days + 1) <= last(new_simulation_list)
                                    for t3 in 1:final_remaining_repair_hours
                                        outage_scenario[day_group_id + num_remaining_repair_days + 1][sce_id][t3, gen_idx] = false
                                    end
                                end
                            end
                        end
                    
                        previous_status = unit_status[t]
                    end

                outage_scenario[day_group_id][sce_id][:, gen_idx] = unit_status
                
            end
        end
    end

    # Seed the recalculation pass with the last warm-start day's outage state.
    last_day_outage_data = deepcopy(outage_scenario[maximum(new_simulation_list)])

    for day_group_idx in eachindex(simulation_list)

        day_group_id = simulation_list[day_group_idx]
        start_day = am_ref[:nw][0][:repday_groups][day_group_id]["Start_Day_Id"]
        end_day = am_ref[:nw][0][:repday_groups][day_group_id]["End_Day_Id"]
        num_hours = (end_day - start_day + 1) * 24

        outage_scenario[day_group_id] = Dict{Int64, Any}()
        for sce in risk_scenario_list
            outage_scenario[day_group_id][sce] = trues(num_hours, num_gens)
        end
    end

    # Recalculate carrying prior-day state forward from the warm-start draws.
    for day_group_idx in eachindex(simulation_list)

        day_group_id = simulation_list[day_group_idx]
        start_day = am_ref[:nw][0][:repday_groups][day_group_id]["Start_Day_Id"]
        end_day = am_ref[:nw][0][:repday_groups][day_group_id]["End_Day_Id"]
        day_list = [i for i in start_day:end_day]
        num_hours = (end_day - start_day + 1) * 24

        local_average_temperature = get_local_average_temperature(day_group_id)

        gen_indices = collect(keys(am_ref[:nw][0][:RA_gen_index]))
        Threads.@threads for i in eachindex(gen_indices)
            gen_idx = gen_indices[i]
            rng = Random.Xoshiro(seed_value + gen_idx + day_group_id * 100000)

            UnitprobabilityType = am_ref[:nw][0][:RA_gen_index][gen_idx]["RA_FOR"]
            bus_idx = am_ref[:nw][0][:RA_gen_index][gen_idx]["bus_idx"]
            ICAP = am_ref[:nw][0][:RA_gen_index][gen_idx]["ICAP"]

            rateAD = zeros(1, num_hours)
            rateDD = zeros(1, num_hours)

                if typeof(UnitprobabilityType) == String

                    UnitprobabilityTypeAD = string(UnitprobabilityType, "AD")
                    UnitprobabilityTypeDD = string(UnitprobabilityType, "DD")

                    for t in [t for t in 1:num_hours]

                        if local_average_temperature[t, bus_idx] > reference_temp
                            # regression cols AD=[Hi,Lo,tmpH,tmpH2,tmpL,tmpL2]; DD appends [rt_avg,rt_95,rt_90]
                            tmp = local_average_temperature[t, bus_idx] - reference_temp
                            rateAD[t] = 1 / (1 + exp(-1 * (regressionAD[:, UnitprobabilityTypeAD][1] + regressionAD[:, UnitprobabilityTypeAD][3] * tmp + regressionAD[:, UnitprobabilityTypeAD][4] * tmp^2)))
                            rateDD[t] = 1 / (1 + exp(-1 * (regressionDD[:, UnitprobabilityTypeDD][1] + regressionDD[:, UnitprobabilityTypeDD][3] * tmp + regressionDD[:, UnitprobabilityTypeDD][4] * tmp^2)))
                        else
                            tmp = reference_temp - local_average_temperature[t, bus_idx]
                            rateAD[t] = 1 / (1 + exp(-1 * (regressionAD[:, UnitprobabilityTypeAD][2] + regressionAD[:, UnitprobabilityTypeAD][5] * tmp + regressionAD[:, UnitprobabilityTypeAD][6] * tmp^2)))
                            rateDD[t] = 1 / (1 + exp(-1 * (regressionDD[:, UnitprobabilityTypeDD][2] + regressionDD[:, UnitprobabilityTypeDD][5] * tmp + regressionDD[:, UnitprobabilityTypeDD][6] * tmp^2)))
                        end

                    end
                else
                    # Numeric UnitprobabilityType is FOR. Derive AD/DD from FOR + MTTR.
                    # MTTR must be in hours.
                    mttr_hours = get(RA_setting, "repair_time_average_hours", 32)
                    mttr_hours = max(mttr_hours, 1.0) # avoid divide-by-zero or invalid probabilities
                    
                    for t in 1:num_hours
                        # rateDD = P(stay down) per hour, geometric repair
                        rateDD[t] = 1 - (1 / mttr_hours)
                        # rateAD = P(fail) per hour derived from FOR
                        for_val = UnitprobabilityType
                        rateAD[t] = for_val / (mttr_hours * max(1e-6, (1 - for_val)))
                        rateAD[t] = clamp(rateAD[t], 1e-6, 1 - 1e-6)
                        rateDD[t] = clamp(rateDD[t], 1e-6, 1 - 1e-6)
                    end
                end

                max_repair_hours = 10000   # default; no bound
                if typeof(UnitprobabilityType) == String

                    UnitprobabilityTypeDD = string(UnitprobabilityType, "DD")
                    
                    if RA_setting["repair_time_bound"] == "average"
                        max_repair_hours = regressionDD[:, UnitprobabilityTypeDD][7]    #repair_time_average
                    elseif RA_setting["repair_time_bound"] == "95 percentile"
                        max_repair_hours = regressionDD[:, UnitprobabilityTypeDD][8]    #repair_time_95
                    elseif RA_setting["repair_time_bound"] == "90 percentile"
                        max_repair_hours = regressionDD[:, UnitprobabilityTypeDD][9]    #repair_time_90
                    elseif RA_setting["repair_time_bound"] == "NA"
                        max_repair_hours = 10000   # no bounds
                    end
                else
                    mttr_hours = get(RA_setting, "repair_time_average_hours", 32)
                    mttr_hours = max(mttr_hours, 1.0)
                    if RA_setting["repair_time_bound"] == "average"
                        max_repair_hours = Int(ceil(mttr_hours))
                    elseif RA_setting["repair_time_bound"] == "95 percentile"
                        max_repair_hours = repair_quantile_from_mttr(mttr_hours, 0.95)
                    elseif RA_setting["repair_time_bound"] == "90 percentile"
                        max_repair_hours = repair_quantile_from_mttr(mttr_hours, 0.90)
                    elseif RA_setting["repair_time_bound"] == "NA"
                        max_repair_hours = 10000   # no bounds
                    end
                end

            @inbounds for sce_id in risk_scenario_list

                unit_status = outage_scenario[day_group_id][sce_id][:, gen_idx]

                    previous_status = true
                    current_status = true

                    if day_group_id == 1
                        previous_status = last_day_outage_data[sce_id][num_hours, gen_idx]
                    end

                    @inbounds for t in 1:num_hours

                        # carry outage state across day-group boundary (t==1 uses prior day's last hour)
                        if (day_group_id > 1) && (t == 1)
                            if outage_scenario[day_group_id - 1][sce_id][num_hours, gen_idx] == true  # if it was on, check AD probability
                                rand(rng, Bernoulli(rateAD[1])) == 1 ? current_status = false : current_status = true
                            end
                        else
                            if previous_status == true  # if currently on, check AD probability
                                rand(rng, Bernoulli(rateAD[t])) == 1 ? current_status = false : current_status = true
                            end
                        end

                        if (previous_status == true) && (current_status == false)  # if unit is off, check when it will be on again
                            repair_time = min(rand(rng, Geometric(1 - rateDD[t])), max_repair_hours)
                            repair_time = max(1, Int(ceil(repair_time)))

                            if (num_hours - t + 1) >= repair_time
                                for t2 in t:(t+repair_time-1)
                                    unit_status[t2] = false
                                end
                            else

                                # fill the current day first
                                for t2 in t:num_hours
                                    unit_status[t2] = false
                                end

                                # check remaining repair hours
                                remaining_repair_hours = (repair_time - (num_hours - t + 1))
                                num_remaining_repair_days = div(remaining_repair_hours, num_hours)

                                if num_remaining_repair_days != 0
                                    
                                    for next_day_count in 1:num_remaining_repair_days
                                        if (day_group_id + next_day_count) <= last(simulation_list)
                                            outage_scenario[day_group_id + next_day_count][sce_id][:, gen_idx] = falses(num_hours)
                                        end
                                    end
                                end
                                
                                # fill the last day
                                final_remaining_repair_hours = remaining_repair_hours - num_remaining_repair_days * num_hours
                                if (day_group_id + num_remaining_repair_days + 1) <= last(simulation_list)
                                    for t3 in 1:final_remaining_repair_hours
                                        outage_scenario[day_group_id + num_remaining_repair_days + 1][sce_id][t3, gen_idx] = false
                                    end
                                end
                            end
                        end
                    
                        previous_status = unit_status[t]
                    end

                outage_scenario[day_group_id][sce_id][:, gen_idx] = unit_status
                
            end
        end
        
    end

    if apply_filter == false
        return outage_scenario, filtered_outage_scenario_list
    end

    # filter risk scenarios 
    for day_group_idx in eachindex(simulation_list)
        
        if RA_setting["risk_filtering_flag"] == false

            # no filtering 
            if length(RA_setting["preselected_days_list"]) > 0
            
                # selected days without filtering
                day_group_id = simulation_list[day_group_idx]
                if day_group_id in RA_setting["preselected_days_list"]
                    risk_key_list = collect(keys(outage_scenario[day_group_id]))
                    for scenario_id in 1:length(outage_scenario[day_group_id])
                        push!(filtered_outage_scenario_list[day_group_id], risk_key_list[scenario_id])
                    end
                else
                    filtered_outage_scenario_list[day_group_id] = []
                end

            else
                # full day without filtering 
                day_group_id = simulation_list[day_group_idx]
                risk_key_list = collect(keys(outage_scenario[day_group_id]))
                for scenario_id in 1:length(outage_scenario[day_group_id])
                    push!(filtered_outage_scenario_list[day_group_id], risk_key_list[scenario_id])
                end
            end

        else
        
            if length(RA_setting["preselected_days_list"]) > 0
            
                # selected days filtering
                day_group_id = simulation_list[day_group_idx]
                if day_group_id in RA_setting["preselected_days_list"]
                    filtered_outage_scenario_list[day_group_id] = filter_risk_scenarios(am_ref, am_setting, day_group_id, RA_setting, outage_scenario[day_group_id])
                else
                    filtered_outage_scenario_list[day_group_id] = []
                end

            else
                # full day filtering
                day_group_id = simulation_list[day_group_idx]
                filtered_outage_scenario_list[day_group_id] = filter_risk_scenarios(am_ref, am_setting, day_group_id, RA_setting, outage_scenario[day_group_id])
            end

        end

        
    
        for risk_id in keys(outage_scenario[day_group_id])
            if !(risk_id in filtered_outage_scenario_list[day_group_id])
                outage_scenario[day_group_id][risk_id] = []
            end
        end

    end

    

    return outage_scenario, filtered_outage_scenario_list
end



function generate_outage_scenarios_for_single_unit(am_ref, simulation_list, reference_temp, num_risk_scenario, gen_info, RA_setting; output_option="Status")

    outage_scenario_report = Dict{Int, Any}()
    risk_scenario_list = [i for i in 1:num_risk_scenario]
    outage_scenario_report = calculate_outage_rate_of_single_unit(am_ref, simulation_list, gen_info, reference_temp, risk_scenario_list, RA_setting; output_option)

    return outage_scenario_report
end



function calculate_outage_rate_of_single_unit(am_ref, simulation_list, gen_info, reference_temp, risk_scenario_list, RA_setting, seed_value=1; output_option="Status")

    # Geometric repair model: P(T <= k) = 1 - (1 - p_repair)^k, with p_repair = 1/MTTR (hours)
    repair_quantile_from_mttr(mttr_hours, q) = begin
        mttr_hours = max(mttr_hours, 1.0)
        p_repair = min(1.0, 1 / mttr_hours)
        if p_repair >= 1.0
            1
        else
            max(1, Int(ceil(log(1 - q) / log(1 - p_repair))))
        end
    end
    
    outage_scenario = Dict{Int64, Any}()
    for day_group_id in keys(am_ref[:repday_groups])
        start_day = am_ref[:repday_groups][day_group_id]["Start_Day_Id"]
        end_day = am_ref[:repday_groups][day_group_id]["End_Day_Id"]
        num_hours = (end_day - start_day + 1) * 24

        outage_scenario[day_group_id] = Dict{Int64, Any}()
        for sce in risk_scenario_list
            outage_scenario[day_group_id][sce] = trues(num_hours, 1)
        end
    end

    regressionAD = am_ref[:regressionAD]
    regressionDD = am_ref[:regressionDD]

    num_gens = 1

    temperature_cache = Dict{Int, Matrix{Float32}}()
    get_local_average_temperature(day_group_id) = get!(temperature_cache, day_group_id) do
        start_day = am_ref[:repday_groups][day_group_id]["Start_Day_Id"]
        end_day = am_ref[:repday_groups][day_group_id]["End_Day_Id"]
        day_list = [i for i in start_day:end_day]
        num_hours = (end_day - start_day + 1) * 24

        temperatureData = am_ref[:temperatureData]
        local_average_temperature = zeros(Float32, num_hours, length(am_ref[:bus]))
        for bus_idx in keys(am_ref[:bus])
            local_temperature = filter(:Day_ID => in(day_list), temperatureData)[:, am_ref[:bus][bus_idx]["aggregation_info"]["aggregated_regions_bus_i"]]
            local_average_temperature[:, bus_idx] = [mean(c) for c in eachrow(local_temperature)]
        end
        local_average_temperature
    end

    for day_group_id in simulation_list

        start_day = am_ref[:repday_groups][day_group_id]["Start_Day_Id"]
        end_day = am_ref[:repday_groups][day_group_id]["End_Day_Id"]
        day_list = [i for i in start_day:end_day]
        num_hours = (end_day - start_day + 1) * 24

        local_average_temperature = get_local_average_temperature(day_group_id)

        UnitprobabilityType = gen_info["RA_FOR"]
        bus_idx = gen_info["bus_idx"]
        ICAP = gen_info["ICAP"]
        gen_idx = 1

        if ICAP == 0
            # Dummy generators (zero capacity) stay online
            continue
        end
        
        rateAD = zeros(num_hours)
        rateDD = zeros(num_hours)

        if typeof(UnitprobabilityType) == String

            UnitprobabilityTypeAD = string(UnitprobabilityType, "AD")
            UnitprobabilityTypeDD = string(UnitprobabilityType, "DD")

            for t in 1:num_hours

                if local_average_temperature[t, bus_idx] > reference_temp
                    # regression cols = [Hi, Lo, tmpH, tmpH2, tmpL, tmpL2]
                    tmp = local_average_temperature[t, bus_idx] - reference_temp
                    rateAD[t] = 1 / (1 + exp(-1 * (regressionAD[:, UnitprobabilityTypeAD][1] + regressionAD[:, UnitprobabilityTypeAD][3] * tmp + regressionAD[:, UnitprobabilityTypeAD][4] * tmp^2)))
                    rateDD[t] = 1 / (1 + exp(-1 * (regressionDD[:, UnitprobabilityTypeDD][1] + regressionDD[:, UnitprobabilityTypeDD][3] * tmp + regressionDD[:, UnitprobabilityTypeDD][4] * tmp^2)))
                else
                    tmp = reference_temp - local_average_temperature[t, bus_idx]
                    rateAD[t] = 1 / (1 + exp(-1 * (regressionAD[:, UnitprobabilityTypeAD][2] + regressionAD[:, UnitprobabilityTypeAD][5] * tmp + regressionAD[:, UnitprobabilityTypeAD][6] * tmp^2)))
                    rateDD[t] = 1 / (1 + exp(-1 * (regressionDD[:, UnitprobabilityTypeDD][2] + regressionDD[:, UnitprobabilityTypeDD][5] * tmp + regressionDD[:, UnitprobabilityTypeDD][6] * tmp^2)))
                end

            end
        else
            # Numeric UnitprobabilityType is FOR. Derive AD/DD from FOR + MTTR.
            # MTTR must be in hours.
            mttr_hours = get(RA_setting, "repair_time_average_hours", get(RA_setting, "repair_time_average_hours_value", 32))
            mttr_hours = max(mttr_hours, 1.0) # avoid divide-by-zero or invalid probabilities
            
            for t in 1:num_hours
                # rateDD = P(stay down) per hour, geometric repair
                rateDD[t] = 1 - (1 / mttr_hours)
                # rateAD = P(fail) per hour derived from FOR
                for_val = UnitprobabilityType
                rateAD[t] = for_val / (mttr_hours * max(1e-6, (1 - for_val)))
                rateAD[t] = clamp(rateAD[t], 1e-6, 1 - 1e-6)
                rateDD[t] = clamp(rateDD[t], 1e-6, 1 - 1e-6)
            end
        end

        max_repair_hours = 10000   # default; no bound
        if typeof(UnitprobabilityType) == String

            UnitprobabilityTypeDD = string(UnitprobabilityType, "DD")
            
            if RA_setting["repair_time_bound"] == "average"
                max_repair_hours = regressionDD[:, UnitprobabilityTypeDD][7]    #repair_time_average
            elseif RA_setting["repair_time_bound"] == "95 percentile"
                max_repair_hours = regressionDD[:, UnitprobabilityTypeDD][8]    #repair_time_95
            elseif RA_setting["repair_time_bound"] == "90 percentile"
                max_repair_hours = regressionDD[:, UnitprobabilityTypeDD][9]    #repair_time_90
            end
        else
            mttr_hours = get(RA_setting, "repair_time_average_hours", 32)
            mttr_hours = max(mttr_hours, 1.0)
            if RA_setting["repair_time_bound"] == "average"
                max_repair_hours = Int(ceil(mttr_hours))
            elseif RA_setting["repair_time_bound"] == "95 percentile"
                max_repair_hours = repair_quantile_from_mttr(mttr_hours, 0.95)
            elseif RA_setting["repair_time_bound"] == "90 percentile"
                max_repair_hours = repair_quantile_from_mttr(mttr_hours, 0.90)
            end
        end

        Threads.@threads for i in eachindex(risk_scenario_list)
            sce_id = risk_scenario_list[i]
            rng = Random.Xoshiro(seed_value + gen_idx + day_group_id * 100000 + sce_id)

            unit_status = outage_scenario[day_group_id][sce_id]

            previous_status = unit_status[1]
            current_status = true
            @inbounds for t in 1:num_hours

                # carry outage state across day-group boundary (t==1 uses prior day's last hour)
                if (day_group_id > 1) && (t == 1)
                    if outage_scenario[day_group_id - 1][sce_id][num_hours, gen_idx] == true  # if it was on, check AD probability
                        rand(rng, Bernoulli(rateAD[1])) == 1 ? current_status = false : current_status = true
                    end
                else
                    if previous_status == true  # if currently on, check AD probability
                        rand(rng, Bernoulli(rateAD[t])) == 1 ? current_status = false : current_status = true
                    end
                end

                if (previous_status == true) && (current_status == false)  # if unit is off, check when it will be on again
                    repair_time = min(rand(rng, Geometric(1 - rateDD[t])), max_repair_hours)
                    repair_time = max(1, Int(ceil(repair_time)))

                    if (num_hours - t + 1) >= repair_time   
                        for t2 in t:(t+repair_time-1)
                            unit_status[t2] = false
                        end
                    else

                        # fill the current day first
                        for t2 in t:num_hours
                            unit_status[t2] = false
                        end

                        # check remaining repair hours
                        remaining_repair_hours = (repair_time - (num_hours - t + 1))
                        num_remaining_repair_days = div(remaining_repair_hours, num_hours)

                        if num_remaining_repair_days != 0
                            
                            for next_day_count in 1:num_remaining_repair_days
                                if (day_group_id + next_day_count) <= last(simulation_list)
                                    outage_scenario[day_group_id + next_day_count][sce_id][:, gen_idx] = falses(num_hours)
                                end
                            end
                        end
                        
                        # fill the last day
                        final_remaining_repair_hours = remaining_repair_hours - num_remaining_repair_days * num_hours
                        if (day_group_id + num_remaining_repair_days + 1) <= last(simulation_list)
                            for t3 in 1:final_remaining_repair_hours
                                outage_scenario[day_group_id + num_remaining_repair_days + 1][sce_id][t3, gen_idx] = false
                            end
                        end
                    end
                end
            
                previous_status = unit_status[t]
            end

            outage_scenario[day_group_id][sce_id] = unit_status
        end
    end

    return outage_scenario
end

