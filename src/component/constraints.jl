# Copyright © 2026, UChicago Argonne, LLC. All Rights Reserved.
# SPDX-License-Identifier: BSD-3-Clause  (see LICENSE)

# Purpose: Shared JuMP constraint generator functions for ALEAF models.
# Scope: Centralized constraint builders reused by multiple models.
# Notes: Keep signatures stable for caller compatibility.


function constraint_carbon_emission_reduction_target_ny(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, n::Int, y::Int, ids_d; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=true)

    ids_local_i = am.ref[:nw][0][:zone]["policy"][string(n)]["aggregation_info"]["new_local_gen_idx"]
    ids_dht = [(d,h,t) for d in ids_d for h in am.setting["run_H"] for t in am.setting["run_T"]]
    
    reference_carbon_emission_level_value = am.ref[:nw][0][:zone]["policy"][string(n)]["CERT"]["reference"] * 1000 # million ton -> M ton
    Carbon_Emission_Reduction_Target = (1 - am.ref[:nw][0][:zone]["policy"][string(n)]["CERT"]["target"][string(y)])

    sum_emission = JuMP.AffExpr(0.0)
    for i in ids_local_i
        bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
        tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)     
        
        for (d,h,t) in ids_dht
            Emission = parameter(am, bus_idx, :gen_bus, tech_idx, "Emission_CO2") # metric tonnes CO2 per MWh -> M ton
            if Emission > 0.0001
                JuMP.add_to_expression!(sum_emission, parameter(am, 0, :repdays, "NumDays", d) * (1/1000) * Emission, variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y)))
            end
        end
    end

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($n,$y"))) end    
    expr = JuMP.@expression(JuMP_model,  
        sum_emission - reference_carbon_emission_level_value * Carbon_Emission_Reduction_Target   
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($n,$y)")) end    
end

function constraint_max_ENS_MWh_cap_y(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=true)


    Max_ENS_MWh = am.setting["Simulation Configuration"]["Max_ENS_MWh"]

    max_ENS_MWh_y = variable(am, nw, decomp_group, :max_ENS_MWh_y, (y))

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        max_ENS_MWh_y - Max_ENS_MWh
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($y)")) end    
    
end

function constraint_max_ENS_MWh_y(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, y::Int, ids_n, ids_d; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=true)

    ids_dht = [(d,h,t) for d in ids_d for h in am.setting["run_H"] for t in am.setting["run_T"]]


    max_ENS_MWh_y = variable(am, nw, decomp_group, :max_ENS_MWh_y, (y))

    for n in ids_n
        for (d,h,t) in ids_dht
            
            if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($n,$d,$h,$t,$y)"))) end    
            expr = JuMP.@expression(JuMP_model,  
                max_ENS_MWh_y - variable(am, nw, decomp_group, :ens_ndhty, (n,d,h,t,y))
                )   
            constraint = JuMP.@constraint(JuMP_model, 
                expr >= 0
                )    
            if const_name_flag JuMP.set_name(constraint, string(const_name, "_($n,$d,$h,$t,$y)")) end    
            
        end
    end
        
end

function constraint_total_ENS_MWh_cap_y(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, y::Int, ids_n, ids_d; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=true)

    ids_dht = [(d,h,t) for d in ids_d for h in am.setting["run_H"] for t in am.setting["run_T"]]

    Total_ENS_MWh = am.setting["Simulation Configuration"]["Total_ENS_MWh"]

    sum_ens_idht = JuMP.AffExpr(0.0)
    for n in ids_n
        for (d,h,t) in ids_dht
            JuMP.add_to_expression!(sum_ens_idht, parameter(am, 0, :repdays, "NumDays", d), variable(am, nw, decomp_group, :ens_ndhty, (n,d,h,t,y)))
        end
    end
    
    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        sum_ens_idht - Total_ENS_MWh
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($y)")) end    
end

function constraint_ENS_hours_approx_cap_y_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, y::Int, ids_n, ids_d; nw::Int=am.cnw, update::Bool=false, ens_flag::Bool=true)

    ids_dht = [(d,h,t) for d in ids_d for h in am.setting["run_H"] for t in am.setting["run_T"]]

    ENS_Hours = am.setting["Simulation Configuration"]["ENS_Hours"] * 60 

    sum_ens_hours = JuMP.AffExpr(0.0)
    for n in ids_n

        for (d,h,t) in ids_dht

            Demand = 0.0
            bus_data = am.ref[:nw][0][:bus][n]
            region_map = get(am.ref[:nw][0], :profile_data_region_map, nothing)
            for ba_id in bus_data["aggregation_info"]["aggregated_regions_bus_i"]
                Demand += parameter(am, 0, :planning_stages, "repdays", "data", "load", d, h, t, y)[profile_data_region(region_map, "load", ba_id)] * bus_data["original load"][ba_id]
            end

            if Demand > 0
                JuMP.add_to_expression!(sum_ens_hours, 60 * (1 / Demand) * parameter(am, 0, :repdays, "NumDays", d), variable(am, nw, decomp_group, :ens_ndhty, (n,d,h,t,y)))
            end
        end
    end
        
    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        sum_ens_hours - ENS_Hours
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    JuMP.set_name(constraint, string(const_name, "_($y)"))

end

function constraint_clean_energy_generation_ny(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, n::Int, y::Int, ids_d; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=true)

    ids_local_i = am.ref[:nw][0][:zone]["policy"][string(n)]["aggregation_info"]["new_local_gen_idx"]
    ids_dht = [(d,h,t) for d in ids_d for h in am.setting["run_H"] for t in am.setting["run_T"]]

    Clean_Energy_Generation_Target = am.ref[:nw][0][:zone]["policy"][string(n)]["CEGT"]["target"][string(y)] 

    sum_g_idht = JuMP.AffExpr(0.0)
    sum_clean_energy_g_idht = JuMP.AffExpr(0.0)
    for i in ids_local_i
        bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
        tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)     
        Clean_Energy_Flag = parameter(am, bus_idx, :gen_bus, tech_idx, "Clean_Energy_Flag")

        for (d,h,t) in ids_dht

            if parameter(am, 0, :gen_index, "UNIT_CATEGORY", i) == "STORAGE" 
                if Clean_Energy_Flag == true
                    JuMP.add_to_expression!(sum_g_idht, 1, variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y)))
                end
            else
                JuMP.add_to_expression!(sum_g_idht, 1, variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y)))
            end                       
            
            if Clean_Energy_Flag == true
                JuMP.add_to_expression!(sum_clean_energy_g_idht, 1, variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y)))
            end
        end
    end


    slack_CEG_ny = variable(am, nw, decomp_group, :slack_CEG_ny, (n,y))
    
    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($n,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        sum_clean_energy_g_idht + slack_CEG_ny - Clean_Energy_Generation_Target * sum_g_idht
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr >= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($n,$y)")) end    
end

function constraint_clean_energy_generation_ndy(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, n::Int, d::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=true)

    ids_local_i = am.ref[:nw][0][:zone]["policy"][string(n)]["aggregation_info"]["new_local_gen_idx"]
    ids_ht = [(h,t) for h in am.setting["run_H"] for t in am.setting["run_T"]]

    Clean_Energy_Generation_Target = am.ref[:nw][0][:zone]["policy"][string(n)]["CEGT"]["target"][string(y)] 

    sum_g_idht = JuMP.AffExpr(0.0)
    sum_clean_energy_g_idht = JuMP.AffExpr(0.0)
    for i in ids_local_i
        bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
        tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)     
        Clean_Energy_Flag = parameter(am, bus_idx, :gen_bus, tech_idx, "Clean_Energy_Flag")
        
        for (h,t) in ids_ht
            
            if parameter(am, 0, :gen_index, "UNIT_CATEGORY", i) == "STORAGE" 
                if Clean_Energy_Flag == true
                    JuMP.add_to_expression!(sum_g_idht, 1, variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y)))
                end
            else
                JuMP.add_to_expression!(sum_g_idht, 1, variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y)))
            end
            
            if Clean_Energy_Flag == true
                JuMP.add_to_expression!(sum_clean_energy_g_idht, 1, variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y)))
            end
        end
    end

    slack_CEG_ndy = variable(am, nw, decomp_group, :slack_CEG_ndy, (n,d,y))
    
    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($n,$d,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        sum_clean_energy_g_idht + slack_CEG_ndy - Clean_Energy_Generation_Target * sum_g_idht
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr >= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($n,$d,$y)")) end    
end

function constraint_RPS_regional_ny(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, n::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=true)

    ids_local_i = am.ref[:nw][0][:zone]["policy"][string(n)]["aggregation_info"]["new_local_gen_idx"]

    Wind_Ons_genA = am.ref[:nw][0][:zone]["policy"][string(n)]["vre_aggregated_data"][string(y)]["wind_ons_shape"]
    Wind_Ofs_genA = am.ref[:nw][0][:zone]["policy"][string(n)]["vre_aggregated_data"][string(y)]["wind_ofs_shape"]
    Solar_genA = am.ref[:nw][0][:zone]["policy"][string(n)]["vre_aggregated_data"][string(y)]["pv_shape"]
    Hydro_genA = am.ref[:nw][0][:zone]["policy"][string(n)]["vre_aggregated_data"][string(y)]["hydro_shape"]
    CSP_genA = am.ref[:nw][0][:zone]["policy"][string(n)]["vre_aggregated_data"][string(y)]["csp_shape"]
    RTPV_genA = am.ref[:nw][0][:zone]["policy"][string(n)]["vre_aggregated_data"][string(y)]["rtpv_shape"]

    bus_list = am.ref[:nw][0][:zone]["policy"][string(n)]["aggregation_info"]["zone_bus_idx"]
    regional_annualDemand = get_zone_annual_demand_with_growth(am, bus_list, y; nw=0)
    regional_RPS = am.ref[:nw][0][:zone]["policy"][string(n)]["RPS"]["target"][string(y)]

    sum_VRE = JuMP.AffExpr(0.0)
    for i in ids_local_i
        bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
        tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)     
        CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")
        profile_type = parameter(am, bus_idx, :gen_bus, tech_idx, "Profile_Type")

        if profile_type == "wind_ons"
            JuMP.add_to_expression!(sum_VRE, Wind_Ons_genA * CAP, variable(am, nw, decomp_group, :u_G_iy, (i,y)))
        elseif profile_type == "wind_ofs"
            JuMP.add_to_expression!(sum_VRE, Wind_Ofs_genA * CAP, variable(am, nw, decomp_group, :u_G_iy, (i,y)))
        elseif profile_type == "pv"
            JuMP.add_to_expression!(sum_VRE, Solar_genA * CAP, variable(am, nw, decomp_group, :u_G_iy, (i,y)))
        elseif profile_type == "hydro"
            JuMP.add_to_expression!(sum_VRE, Hydro_genA * CAP, variable(am, nw, decomp_group, :u_G_iy, (i,y)))
        elseif profile_type == "csp"
            JuMP.add_to_expression!(sum_VRE, CSP_genA * CAP, variable(am, nw, decomp_group, :u_G_iy, (i,y)))
        elseif profile_type == "rtpv"
            JuMP.add_to_expression!(sum_VRE, RTPV_genA * CAP, variable(am, nw, decomp_group, :u_G_iy, (i,y)))
        end
    end

    slack_RPS_ny = variable(am, nw, decomp_group, :slack_RPS_ny, (n,y))

    
    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($n,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        sum_VRE + slack_RPS_ny - regional_RPS * regional_annualDemand
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr >= 0
        )    

    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($n,$y)")) end    

    # force slack variable to zero if the flag is false
    if am.setting["Simulation Configuration"]["Allow_Alternative_RPS_Compliance_Flag"] == false
        JuMP.fix(slack_RPS_ny, 0, force=true)
    end

end

function constraint_ES_AET_y(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, y::Int, ids_d; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    ids_dhty = [(d,h,t,y) for d in ids_d for h in am.setting["run_H"] for t in am.setting["run_T"]]

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
    
    AET = parameter(am, bus_idx, :gen_bus, tech_idx, "AET")
    BATEFF = sqrt(parameter(am, bus_idx, :gen_bus, tech_idx, "BATEFF"))

    total_num_days = 0.0
    for d in ids_d
        total_num_days += parameter(am, 0, :repdays, "NumDays", d)
    end
    AET_scale = total_num_days / 365    # this will be 1 if PH is not used.


    if AET > 0
    
        sum_AET = JuMP.AffExpr(0.0)
        for (d, h, t, y) in ids_dhty

            reserve_enabled(am, :reg_up) && up_reserve_eligible(am, i) && JuMP.add_to_expression!(sum_AET, 0.15 * parameter(am, 0, :repdays, "NumDays", d) / BATEFF, variable(am, nw, decomp_group, :reg_up_idhty, (i,d,h,t,y)))
            reserve_enabled(am, :reg_dn) && JuMP.add_to_expression!(sum_AET, 0.15 * parameter(am, 0, :repdays, "NumDays", d) * BATEFF, variable(am, nw, decomp_group, :reg_dn_idhty, (i,d,h,t,y)))
        
            JuMP.add_to_expression!(sum_AET, parameter(am, 0, :repdays, "NumDays", d) * BATEFF, variable(am, nw, decomp_group, :chg_idhty, (i,d,h,t,y)))
            JuMP.add_to_expression!(sum_AET, parameter(am, 0, :repdays, "NumDays", d) / BATEFF, variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y)))
            
        end

        u_i = variable(am, nw, decomp_group, :u_G_iy, (i,y))

        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$y)"))) end    
        expr = JuMP.@expression(JuMP_model,  
            (1.0 / sqrt(2 * AET * AET_scale)) * sum_AET - sqrt(2 * AET * AET_scale) * u_i
            )   
        constraint = JuMP.@constraint(JuMP_model, 
            expr <= 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$y)")) end   

    end
end

function constraint_ES_AET_OP_y(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, y::Int, ids_d; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false, reserve_flag::Bool=true)

    ids_dhty = [(d,h,t,y) for d in ids_d for h in am.setting["run_H"] for t in am.setting["run_T"]]

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
    
    AET = parameter(am, bus_idx, :gen_bus, tech_idx, "AET")
    BATEFF = sqrt(parameter(am, bus_idx, :gen_bus, tech_idx, "BATEFF"))

    total_num_days = 0.0
    for d in ids_d
        total_num_days += parameter(am, 0, :repdays, "NumDays", d)
    end
    AET_scale = total_num_days / 365    # this will be 1 if PH is not used.

    if AET > 0
    
        sum_AET = JuMP.AffExpr(0.0)
        for (d, h, t, y) in ids_dhty
            

            reserve_enabled(am, :reg_up) && up_reserve_eligible(am, i) && JuMP.add_to_expression!(sum_AET, 0.15 * parameter(am, 0, :repdays, "NumDays", d) / BATEFF, variable(am, nw, decomp_group, :reg_up_idhty, (i,d,h,t,y)))
            reserve_enabled(am, :reg_dn) && JuMP.add_to_expression!(sum_AET, 0.15 * parameter(am, 0, :repdays, "NumDays", d) * BATEFF, variable(am, nw, decomp_group, :reg_dn_idhty, (i,d,h,t,y)))
        
            JuMP.add_to_expression!(sum_AET, parameter(am, 0, :repdays, "NumDays", d) * BATEFF, variable(am, nw, decomp_group, :chg_idhty, (i,d,h,t,y)))
            JuMP.add_to_expression!(sum_AET, parameter(am, 0, :repdays, "NumDays", d) / BATEFF, variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y)))
            
        end

        u_i = variable(am, nw, decomp_group, :u_G_iy, (i,y))

        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$y)"))) end    
        expr = JuMP.@expression(JuMP_model,  
            (1.0 / sqrt(2 * AET * AET_scale)) * sum_AET - sqrt(2 * AET * AET_scale) * u_i
            )   
        constraint = JuMP.@constraint(JuMP_model, 
            expr <= 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$y)")) end    

    end
end

function constraint_ES_SOC_Balance_Inter_SubHour_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int, day_group::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    start_day = am.ref[:nw][0][:repday_groups][day_group]["Start_Day_Id"]
    end_day = am.ref[:nw][0][:repday_groups][day_group]["End_Day_Id"]
    start_day_idx = am.ref[:nw][0][:repday_groups][day_group]["Day_Idx_List"][1]
    
    if (h==1) && (t==1)
        if am.ref[:nw][0][:repdays][d]["Day"] == start_day
            bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
            tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
            BATEFF = sqrt(parameter(am, bus_idx, :gen_bus, tech_idx, "BATEFF"))
            card_T = length(am.setting["run_T"])
            STOMIN = parameter(am, bus_idx, :gen_bus, tech_idx, "STOMIN")
            
            soc_idhty = variable(am, nw, decomp_group, :soc_idhty, (i,d,h,t,y))
            chg_idhty = variable(am, nw, decomp_group, :chg_idhty, (i,d,h,t,y))
            g_idhty = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
            
            if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
            expr = JuMP.@expression(JuMP_model,  
                soc_idhty - STOMIN - (BATEFF * chg_idhty - (1/BATEFF) * g_idhty) / card_T
                )   
            constraint = JuMP.@constraint(JuMP_model, 
                expr == 0
                )    
            if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end   

        else
            bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
            tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
            
            BATEFF = sqrt(parameter(am, bus_idx, :gen_bus, tech_idx, "BATEFF"))
            HFREQ = am.setting["Simulation Configuration"]["HFREQ"]
            card_T = length(am.setting["run_T"])
            last_hour = last(am.setting["run_H"])
            last_subhour = last(am.setting["run_T"])
            
            soc_idhty = variable(am, nw, decomp_group, :soc_idhty, (i,d,h,t,y))
            chg_idhty = variable(am, nw, decomp_group, :chg_idhty, (i,d,h,t,y))
            prior_soc_idhty = variable(am, nw, decomp_group, :soc_idhty, (i,d-1,last_hour,last_subhour,y)) 
            g_idhty = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
            
            if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
            expr = JuMP.@expression(JuMP_model,  
                soc_idhty - prior_soc_idhty - (BATEFF * chg_idhty - (1/BATEFF) * g_idhty)/card_T
                )   
            constraint = JuMP.@constraint(JuMP_model, 
                expr == 0
                )    
            if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end   
        end

    elseif (h != 1) && (t == 1)
        bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
        tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
        
        BATEFF = sqrt(parameter(am, bus_idx, :gen_bus, tech_idx, "BATEFF"))
        HFREQ = am.setting["Simulation Configuration"]["HFREQ"]
        card_T = length(am.setting["run_T"])
        last_subhour = last(am.setting["run_T"])

        soc_idhty = variable(am, nw, decomp_group, :soc_idhty, (i,d,h,t,y))
        chg_idhty = variable(am, nw, decomp_group, :chg_idhty, (i,d,h,t,y))
        prior_soc_idhty = variable(am, nw, decomp_group, :soc_idhty, (i,d,h-HFREQ,last_subhour,y)) 
        g_idhty = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
        
        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
        expr = JuMP.@expression(JuMP_model,  
            soc_idhty - prior_soc_idhty - (BATEFF * chg_idhty - (1/BATEFF) * g_idhty)/card_T
            )   
        constraint = JuMP.@constraint(JuMP_model, 
            expr == 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end    
    
    elseif (h != 1) && (t != 1)
        bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
        tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
        
        BATEFF = sqrt(parameter(am, bus_idx, :gen_bus, tech_idx, "BATEFF"))
        HFREQ = am.setting["Simulation Configuration"]["HFREQ"]
        card_T = length(am.setting["run_T"])

        soc_idhty = variable(am, nw, decomp_group, :soc_idhty, (i,d,h,t,y))
        chg_idhty = variable(am, nw, decomp_group, :chg_idhty, (i,d,h,t,y))
        prior_soc_idhty = variable(am, nw, decomp_group, :soc_idhty, (i,d,h,t-1,y)) 
        g_idhty = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
        
        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
        expr = JuMP.@expression(JuMP_model,  
            soc_idhty - prior_soc_idhty - (BATEFF * chg_idhty - (1/BATEFF) * g_idhty)/card_T
            )   
        constraint = JuMP.@constraint(JuMP_model, 
            expr == 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end    
    
    end

end

function constraint_hybrid_ES_SOC_Neutral_idy(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, y::Int, start_day_idx::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  

    last_subhour = last(am.setting["run_T"])
    last_hour = last(am.setting["run_H"])
    
    BATEFF = sqrt(parameter(am, bus_idx, :gen_bus, tech_idx, "BATEFF"))
    card_T = length(am.setting["run_T"])

    hybrid_gen_idx = am.ref[:nw][0][:gen_index][i]["hybrid_main_gen_idx"]

    soc_end = variable(am, nw, decomp_group, :soc_idhty, (i,d,last_hour,last_subhour,y))
    soc_beginning = variable(am, nw, decomp_group, :soc_idhty, (i,start_day_idx,1,1,y))
    chg_beginning = variable(am, nw, decomp_group, :chg_idhty, (i,start_day_idx,1,1,y))
    g_G_ES_beginning = variable(am, nw, decomp_group, :g_G_ES_idhty, (hybrid_gen_idx,start_day_idx,1,1,y))
    g_beginning = variable(am, nw, decomp_group, :g_idhty, (i,start_day_idx,1,1,y))
    
    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        soc_end - soc_beginning + (BATEFF * (chg_beginning + g_G_ES_beginning) - (1/BATEFF) * g_beginning)/card_T
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$y)")) end    

end

function constraint_hybrid_ES_SOC_Balance_Inter_Hour_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int, day_group::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    start_day = am.ref[:nw][0][:repday_groups][day_group]["Start_Day_Id"]
    end_day = am.ref[:nw][0][:repday_groups][day_group]["End_Day_Id"]
    start_day_idx = am.ref[:nw][0][:repday_groups][day_group]["Day_Idx_List"][1]
    hybrid_gen_idx = am.ref[:nw][0][:gen_index][i]["hybrid_main_gen_idx"]
    
    if (h==1) && (t==1)
        if am.ref[:nw][0][:repdays][d]["Day"] == start_day
            bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
            tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
            BATEFF = sqrt(parameter(am, bus_idx, :gen_bus, tech_idx, "BATEFF"))
            card_T = length(am.setting["run_T"])
            STOMIN = parameter(am, bus_idx, :gen_bus, tech_idx, "STOMIN")

            soc_idhty = variable(am, nw, decomp_group, :soc_idhty, (i,d,h,t,y))
            chg_idhty = variable(am, nw, decomp_group, :chg_idhty, (i,d,h,t,y))
            g_G_ES_idhty = variable(am, nw, decomp_group, :g_G_ES_idhty, (hybrid_gen_idx,d,h,t,y))
            g_idhty = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
            
            ESE_MWh_iy = haskey(am.var[:nw][nw][decomp_group], :u_ESE_iy) ?
                variable(am, nw, decomp_group, :u_ESE_iy, (i, y)) :
                parameter(am, bus_idx, :gen_bus, tech_idx, "ES_MWh")

            if am.setting["Simulation Configuration"]["storage initialization option"] == "Minimum"
            
                if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
                expr = JuMP.@expression(JuMP_model,  
                    soc_idhty - STOMIN - (BATEFF * (chg_idhty + g_G_ES_idhty) - (1/BATEFF) * g_idhty) / card_T
                    )   
                constraint = JuMP.@constraint(JuMP_model, 
                    expr == 0
                    )    
                if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end    

            elseif am.setting["Simulation Configuration"]["storage initialization option"] == "Middle"

                if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
                expr = JuMP.@expression(JuMP_model,  
                    soc_idhty - ESE_MWh_iy / 2 - (BATEFF * (chg_idhty + g_G_ES_idhty) - (1/BATEFF) * g_idhty) / card_T
                    )   
                constraint = JuMP.@constraint(JuMP_model, 
                    expr == 0
                    )    
                if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end    

            elseif am.setting["Simulation Configuration"]["storage initialization option"] == "Maximum"

                if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
                expr = JuMP.@expression(JuMP_model,  
                    soc_idhty - ESE_MWh_iy - (BATEFF * (chg_idhty + g_G_ES_idhty) - (1/BATEFF) * g_idhty) / card_T
                    )   
                constraint = JuMP.@constraint(JuMP_model, 
                    expr == 0
                    )    
                if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end    


            end
        else
            bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
            tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
            
            BATEFF = sqrt(parameter(am, bus_idx, :gen_bus, tech_idx, "BATEFF"))
            HFREQ = am.setting["Simulation Configuration"]["HFREQ"]
            card_T = length(am.setting["run_T"])
            last_hour = last(am.setting["run_H"])

            soc_idhty = variable(am, nw, decomp_group, :soc_idhty, (i,d,h,t,y))
            chg_idhty = variable(am, nw, decomp_group, :chg_idhty, (i,d,h,t,y))
            g_G_ES_idhty = variable(am, nw, decomp_group, :g_G_ES_idhty, (hybrid_gen_idx,d,h,t,y))
            prior_soc_idhty = variable(am, nw, decomp_group, :soc_idhty, (i,d-1,last_hour,1,y)) # Last time step of last hour is T=1
            g_idhty = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
            
            if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
            expr = JuMP.@expression(JuMP_model,  
                soc_idhty - prior_soc_idhty - (BATEFF * (chg_idhty + g_G_ES_idhty) - (1/BATEFF) * g_idhty)/card_T
                )   
            constraint = JuMP.@constraint(JuMP_model, 
                expr == 0
                )    
            if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end    
        end

    elseif (h != 1) && (t == 1)
        bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
        tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
        
        BATEFF = sqrt(parameter(am, bus_idx, :gen_bus, tech_idx, "BATEFF"))
        HFREQ = am.setting["Simulation Configuration"]["HFREQ"]
        card_T = length(am.setting["run_T"])

        soc_idhty = variable(am, nw, decomp_group, :soc_idhty, (i,d,h,t,y))
        chg_idhty = variable(am, nw, decomp_group, :chg_idhty, (i,d,h,t,y))
        g_G_ES_idhty = variable(am, nw, decomp_group, :g_G_ES_idhty, (hybrid_gen_idx,d,h,t,y))
        prior_soc_idhty = variable(am, nw, decomp_group, :soc_idhty, (i,d,h-HFREQ,1,y)) # Last time step of last hour is T=1
        g_idhty = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
        
        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
        expr = JuMP.@expression(JuMP_model,  
            soc_idhty - prior_soc_idhty - (BATEFF * (chg_idhty + g_G_ES_idhty) - (1/BATEFF) * g_idhty)/card_T
            )   
        constraint = JuMP.@constraint(JuMP_model, 
            expr == 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end    
    end

end

function constraint_ES_SOC_Balance_Inter_Hour_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int, day_group::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    start_day = am.ref[:nw][0][:repday_groups][day_group]["Start_Day_Id"]
    end_day = am.ref[:nw][0][:repday_groups][day_group]["End_Day_Id"]
    start_day_idx = am.ref[:nw][0][:repday_groups][day_group]["Day_Idx_List"][1]
    
    if (h==1) && (t==1)
        if am.ref[:nw][0][:repdays][d]["Day"] == start_day
            bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
            tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
            BATEFF = sqrt(parameter(am, bus_idx, :gen_bus, tech_idx, "BATEFF"))
            card_T = length(am.setting["run_T"])
            STOMIN = parameter(am, bus_idx, :gen_bus, tech_idx, "STOMIN")

            soc_idhty = variable(am, nw, decomp_group, :soc_idhty, (i,d,h,t,y))
            chg_idhty = variable(am, nw, decomp_group, :chg_idhty, (i,d,h,t,y))
            g_idhty = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
            u_ESE_iy = variable(am, nw, decomp_group, :u_ESE_iy, (i, y))

            if am.setting["Simulation Configuration"]["storage initialization option"] == "Minimum"
            
                if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
                expr = JuMP.@expression(JuMP_model,  
                    soc_idhty - STOMIN - (BATEFF * chg_idhty - (1/BATEFF) * g_idhty) / card_T
                    )   
                constraint = JuMP.@constraint(JuMP_model, 
                    expr == 0
                    )    
                if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end    

            elseif am.setting["Simulation Configuration"]["storage initialization option"] == "Middle"

                if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
                expr = JuMP.@expression(JuMP_model,  
                    soc_idhty - u_ESE_iy / 2 - (BATEFF * chg_idhty - (1/BATEFF) * g_idhty) / card_T
                    )   
                constraint = JuMP.@constraint(JuMP_model, 
                    expr == 0
                    )    
                if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end    

            elseif am.setting["Simulation Configuration"]["storage initialization option"] == "Maximum"

                if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
                expr = JuMP.@expression(JuMP_model,  
                    soc_idhty - u_ESE_iy - (BATEFF * chg_idhty - (1/BATEFF) * g_idhty) / card_T
                    )   
                constraint = JuMP.@constraint(JuMP_model, 
                    expr == 0
                    )    
                if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end    


            end
        else
            bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
            tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
            
            BATEFF = sqrt(parameter(am, bus_idx, :gen_bus, tech_idx, "BATEFF"))
            HFREQ = am.setting["Simulation Configuration"]["HFREQ"]
            card_T = length(am.setting["run_T"])
            last_hour = last(am.setting["run_H"])

            soc_idhty = variable(am, nw, decomp_group, :soc_idhty, (i,d,h,t,y))
            chg_idhty = variable(am, nw, decomp_group, :chg_idhty, (i,d,h,t,y))
            prior_soc_idhty = variable(am, nw, decomp_group, :soc_idhty, (i,d-1,last_hour,1,y)) # Last time step of last hour is T=1
            g_idhty = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
            
            if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
            expr = JuMP.@expression(JuMP_model,  
                soc_idhty - prior_soc_idhty - (BATEFF * chg_idhty - (1/BATEFF) * g_idhty)/card_T
                )   
            constraint = JuMP.@constraint(JuMP_model, 
                expr == 0
                )    
            if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end    
        end

    elseif (h != 1) && (t == 1)
        bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
        tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
        
        BATEFF = sqrt(parameter(am, bus_idx, :gen_bus, tech_idx, "BATEFF"))
        HFREQ = am.setting["Simulation Configuration"]["HFREQ"]
        card_T = length(am.setting["run_T"])

        soc_idhty = variable(am, nw, decomp_group, :soc_idhty, (i,d,h,t,y))
        chg_idhty = variable(am, nw, decomp_group, :chg_idhty, (i,d,h,t,y))
        prior_soc_idhty = variable(am, nw, decomp_group, :soc_idhty, (i,d,h-HFREQ,1,y)) # Last time step of last hour is T=1
        g_idhty = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
        
        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
        expr = JuMP.@expression(JuMP_model,  
            soc_idhty - prior_soc_idhty - (BATEFF * chg_idhty - (1/BATEFF) * g_idhty)/card_T
            )   
        constraint = JuMP.@constraint(JuMP_model, 
            expr == 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end    
    end

end

function constraint_ES_DisCharge_Max_Sto_UC_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  

    if parameter(am, bus_idx, :gen_bus, tech_idx, "Storage Commitment") == true

        CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")
        PMAX = parameter(am, bus_idx, :gen_bus, tech_idx, "PMAX")
        BigM = max(1000, parameter(am, bus_idx, :gen_bus, tech_idx, "MAXINVEST"))
        
        sum_variable = JuMP.AffExpr(0.0)
        JuMP.add_to_expression!(sum_variable, 1, variable(am, nw, decomp_group, :g_idhty, (i, d, h, t, y)))
        if up_reserve_eligible(am, i) && am.setting["Simulation Configuration"]["operating_reserve_modeling_option"] == "individual"# [aggregated, individual]
            reserve_enabled(am, :reg_up)  && JuMP.add_to_expression!(sum_variable, 1, variable(am, nw, decomp_group, :reg_up_idhty, (i, d, h, t, y)))
            reserve_enabled(am, :flex_up) && JuMP.add_to_expression!(sum_variable, 1, variable(am, nw, decomp_group, :flex_up_idhty, (i, d, h, t, y)))
            reserve_enabled(am, :spin)    && JuMP.add_to_expression!(sum_variable, 1, variable(am, nw, decomp_group, :spin_idhty, (i, d, h, t, y)))
        end

        sto_c_idhty = variable(am, nw, decomp_group, :sto_c_idhty, (i,d,h,t,y))
        
        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
        expr = JuMP.@expression(JuMP_model,  
            sum_variable - sto_c_idhty * PMAX * CAP * BigM
            )   
        constraint = JuMP.@constraint(JuMP_model, 
            expr <= 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end    
    end

end

function constraint_ES_DisCharge_Max_Sto_UC_NonReserve_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  

    if parameter(am, bus_idx, :gen_bus, tech_idx, "Storage Commitment") == true

        CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")
        MAXINVEST = parameter(am, bus_idx, :gen_bus, tech_idx, "MAXINVEST")
        Big_M = max(MAXINVEST, 100)
        
        g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
        sto_c_idhty = variable(am, nw, decomp_group, :sto_c_idhty, (i,d,h,t,y))
        
        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
        expr = JuMP.@expression(JuMP_model,  
            g_idht - sto_c_idhty * CAP * Big_M
            )   
        constraint = JuMP.@constraint(JuMP_model, 
            expr <= 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end   
    end

end


function constraint_hybrid_ES_Charge_Max_Sto_UC_NoReserve_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  

    if parameter(am, bus_idx, :gen_bus, tech_idx, "Storage Commitment") == true
        
        CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")
        PMAX = parameter(am, bus_idx, :gen_bus, tech_idx, "PMAX")
        BigM = max(1000, parameter(am, bus_idx, :gen_bus, tech_idx, "MAXINVEST"))

        hybrid_gen_idx = am.ref[:nw][0][:gen_index][i]["hybrid_main_gen_idx"]
        
        chg_idht = variable(am, nw, decomp_group, :chg_idhty, (i,d,h,t,y))
        g_G_ES_idhty = variable(am, nw, decomp_group, :g_G_ES_idhty, (hybrid_gen_idx,d,h,t,y))
            
        sto_c_idhty = variable(am, nw, decomp_group, :sto_c_idhty, (i,d,h,t,y))
    

        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
        expr = JuMP.@expression(JuMP_model,  
            chg_idht + g_G_ES_idhty - (1 - sto_c_idhty) * CAP * PMAX * BigM
            )   
        constraint = JuMP.@constraint(JuMP_model, 
            expr <= 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end   
    end

end

function constraint_hybrid_ES_Charge_Max_Sto_UC_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  

    if parameter(am, bus_idx, :gen_bus, tech_idx, "Storage Commitment") == true
        
        CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")
        PMAX = parameter(am, bus_idx, :gen_bus, tech_idx, "PMAX")
        BigM = max(1000, parameter(am, bus_idx, :gen_bus, tech_idx, "MAXINVEST"))

        hybrid_gen_idx = am.ref[:nw][0][:gen_index][i]["hybrid_main_gen_idx"]
        
        sum_variable = JuMP.AffExpr(0.0)
        JuMP.add_to_expression!(sum_variable, 1, variable(am, nw, decomp_group, :chg_idhty, (i, d, h, t, y)))
        JuMP.add_to_expression!(sum_variable, 1, variable(am, nw, decomp_group, :g_G_ES_idhty, (hybrid_gen_idx, d, h, t, y)))
        if am.setting["Simulation Configuration"]["operating_reserve_modeling_option"] == "individual"# [aggregated, individual]
            reserve_enabled(am, :reg_dn)  && JuMP.add_to_expression!(sum_variable, 1, variable(am, nw, decomp_group, :reg_dn_idhty, (i, d, h, t, y)))
            reserve_enabled(am, :flex_dn) && JuMP.add_to_expression!(sum_variable, 1, variable(am, nw, decomp_group, :flex_dn_idhty, (i, d, h, t, y)))
        end
            
        sto_c_idhty = variable(am, nw, decomp_group, :sto_c_idhty, (i,d,h,t,y))
    
        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
        expr = JuMP.@expression(JuMP_model,  
            sum_variable - (1 - sto_c_idhty) * CAP * PMAX * BigM
            )   
        constraint = JuMP.@constraint(JuMP_model, 
            expr <= 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end   
    end

end

function constraint_ES_Charge_Max_Sto_UC_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  

    if parameter(am, bus_idx, :gen_bus, tech_idx, "Storage Commitment") == true
        
        CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")
        PMAX = parameter(am, bus_idx, :gen_bus, tech_idx, "PMAX")
        BigM = max(1000, parameter(am, bus_idx, :gen_bus, tech_idx, "MAXINVEST"))
        
        sum_variable = JuMP.AffExpr(0.0)
        JuMP.add_to_expression!(sum_variable, 1, variable(am, nw, decomp_group, :chg_idhty, (i, d, h, t, y)))
        if am.setting["Simulation Configuration"]["operating_reserve_modeling_option"] == "individual"# [aggregated, individual]
            reserve_enabled(am, :reg_dn)  && JuMP.add_to_expression!(sum_variable, 1, variable(am, nw, decomp_group, :reg_dn_idhty, (i, d, h, t, y)))
            reserve_enabled(am, :flex_dn) && JuMP.add_to_expression!(sum_variable, 1, variable(am, nw, decomp_group, :flex_dn_idhty, (i, d, h, t, y)))
        end

        sto_c_idhty = variable(am, nw, decomp_group, :sto_c_idhty, (i,d,h,t,y))
    

        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
        expr = JuMP.@expression(JuMP_model,  
            sum_variable - (1 - sto_c_idhty) * CAP * PMAX * BigM
            )   
        constraint = JuMP.@constraint(JuMP_model, 
            expr <= 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end   
    end

end

function constraint_hybrid_ES_Charge_from_Grid_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, ids_d, ids_h, ids_t, ids_y; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
    hybrid_ID = parse(Int, parameter(am, bus_idx, :gen_bus, tech_idx, "hybrid_ID"))
    CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")

    Grid_Charge_flag = am.ref[:nw][0][:hybrid][hybrid_ID]["Grid_Charge"]

    if Grid_Charge_flag == false
        for d in ids_d, h in ids_h, t in ids_t, y in ids_y
            JuMP.fix(variable(am, nw, decomp_group, :chg_idhty, (i,d,h,t,y)), 0.0; force = true)
        end
    end
end

function constraint_hybrid_ES_Charge_Max_NonReserve_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
    CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")
    PMAX = parameter(am, bus_idx, :gen_bus, tech_idx, "PMAX")

    hybrid_gen_idx = am.ref[:nw][0][:gen_index][i]["hybrid_main_gen_idx"]
    
    chg_idht = variable(am, nw, decomp_group, :chg_idhty, (i,d,h,t,y))
    g_G_ES_idhty = variable(am, nw, decomp_group, :g_G_ES_idhty, (hybrid_gen_idx,d,h,t,y))
    u_i = parameter(am, bus_idx, :gen_bus, tech_idx, "EXUNITS")
        
    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        chg_idht + g_G_ES_idhty - u_i * CAP * PMAX
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end    

end

function constraint_hybrid_ES_Charge_Max_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
    CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")
    PMAX = parameter(am, bus_idx, :gen_bus, tech_idx, "PMAX")

    hybrid_gen_idx = am.ref[:nw][0][:gen_index][i]["hybrid_main_gen_idx"]
    
    sum_variable = JuMP.AffExpr(0.0)
    JuMP.add_to_expression!(sum_variable, 1, variable(am, nw, decomp_group, :chg_idhty, (i, d, h, t, y)))
    JuMP.add_to_expression!(sum_variable, 1, variable(am, nw, decomp_group, :g_G_ES_idhty, (hybrid_gen_idx, d, h, t, y)))
    if am.setting["Simulation Configuration"]["operating_reserve_modeling_option"] == "individual"# [aggregated, individual]
        reserve_enabled(am, :reg_dn)  && JuMP.add_to_expression!(sum_variable, 1, variable(am, nw, decomp_group, :reg_dn_idhty, (i, d, h, t, y)))
        reserve_enabled(am, :flex_dn) && JuMP.add_to_expression!(sum_variable, 1, variable(am, nw, decomp_group, :flex_dn_idhty, (i, d, h, t, y)))
    end
    u_i = variable(am, nw, decomp_group, :u_G_iy, (i,y))

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        sum_variable - u_i * CAP * PMAX
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end    

end

function constraint_agg_ES_Charge_Max_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, r::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    group_info = am.ref[:nw][0][:reserve_group_lookup][r]                                         
    unit_group = group_info["unit_group"]                                                         
    z = group_info["zone_idx"]                                                                    
    gen_idx = am.ref[:nw][0][:zone]["reserve"][string(z)]["aggregation_info"]["gen_technology_groups"][unit_group]["gen_idx"]
    
    total_capacity = am.ref[:nw][0][:zone]["reserve"][string(z)]["aggregation_info"]["gen_technology_groups"][unit_group]["total_capacity"]       
    
    sum_chg = JuMP.AffExpr(0.0)
    sum_u = JuMP.AffExpr(0.0)
    for i in gen_idx
        JuMP.add_to_expression!(sum_chg, 1, variable(am, nw, decomp_group, :chg_idhty, (i, d, h, t, y)))
        JuMP.add_to_expression!(sum_u, 1, variable(am, nw, decomp_group, :u_G_iy, (i, y)))

        hybrid_type = am.ref[:nw][0][:gen_index][i]["hybrid_type"] 
        if hybrid_type == "ES"  
            hybrid_gen_idx = am.ref[:nw][0][:gen_index][i]["hybrid_main_gen_idx"]
            JuMP.add_to_expression!(sum_chg, 1, variable(am, nw, decomp_group, :g_G_ES_idhty, (hybrid_gen_idx, d, h, t, y)))
        end
        
    end

    sum_res = JuMP.AffExpr(0.0)
    reserve_enabled(am, :reg_dn)  && JuMP.add_to_expression!(sum_res, 1, variable(am, nw, decomp_group, :reg_dn_idhty, (r, d, h, t, y)))
    reserve_enabled(am, :flex_dn) && JuMP.add_to_expression!(sum_res, 1, variable(am, nw, decomp_group, :flex_dn_idhty, (r, d, h, t, y)))

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($r,$d,$h,$t,$y)"))) end
    expr = JuMP.@expression(JuMP_model,
        sum_chg + sum_res - total_capacity * sum_u
        )
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($r,$d,$h,$t,$y)")) end    

end

function constraint_ES_Charge_Max_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
    CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")
    PMAX = parameter(am, bus_idx, :gen_bus, tech_idx, "PMAX")
    
    sum_variable = JuMP.AffExpr(0.0)
    JuMP.add_to_expression!(sum_variable, 1, variable(am, nw, decomp_group, :chg_idhty, (i, d, h, t, y)))
    if am.setting["Simulation Configuration"]["operating_reserve_modeling_option"] == "individual"# [aggregated, individual]
        reserve_enabled(am, :reg_dn)  && JuMP.add_to_expression!(sum_variable, 1, variable(am, nw, decomp_group, :reg_dn_idhty, (i, d, h, t, y)))
        reserve_enabled(am, :flex_dn) && JuMP.add_to_expression!(sum_variable, 1, variable(am, nw, decomp_group, :flex_dn_idhty, (i, d, h, t, y)))
    end
    u_i = variable(am, nw, decomp_group, :u_G_iy, (i,y))
        
    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        sum_variable - u_i * CAP * PMAX
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end    

end

function constraint_agg_ES_SOC_Min_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, r::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    group_info = am.ref[:nw][0][:reserve_group_lookup][r]                                         
    unit_group = group_info["unit_group"]                                                         
    z = group_info["zone_idx"]                                                                    
    gen_idx = am.ref[:nw][0][:zone]["reserve"][string(z)]["aggregation_info"]["gen_technology_groups"][unit_group]["gen_idx"]
    
    STOMIN = 0.0
    BATEFF = 0.0
    contingency_reserve_min_hr = am.setting["Simulation Configuration"]["contingency_reserve_min_duration_value"]

    sum_soc = JuMP.AffExpr(0.0)
    sum_u_EUE = JuMP.AffExpr(0.0)
    for i in gen_idx
        JuMP.add_to_expression!(sum_soc, 1, variable(am, nw, decomp_group, :soc_idhty, (i, d, h, t, y)))
        JuMP.add_to_expression!(sum_u_EUE, 1, variable(am, nw, decomp_group, :u_ESE_iy, (i, y)))

        bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
        tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
        BATEFF = sqrt(sqrt(parameter(am, bus_idx, :gen_bus, tech_idx, "BATEFF")))
        STOMIN += parameter(am, bus_idx, :gen_bus, tech_idx, "STOMIN")
    end

    reg_up_idht  = reserve_enabled(am, :reg_up)  ? variable(am, nw, decomp_group, :reg_up_idhty, (r,d,h,t,y))  : 0.0
    flex_up_idht = reserve_enabled(am, :flex_up) ? variable(am, nw, decomp_group, :flex_up_idhty, (r,d,h,t,y)) : 0.0
    spin_idht    = reserve_enabled(am, :spin)    ? variable(am, nw, decomp_group, :spin_idhty, (r,d,h,t,y))    : 0.0

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($r,$d,$h,$t,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        sum_soc - (1/BATEFF) * (reg_up_idht + flex_up_idht + contingency_reserve_min_hr * spin_idht) - STOMIN
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr >= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($r,$d,$h,$t,$y)")) end    
end

function constraint_ES_SOC_Min_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    if am.setting["Simulation Configuration"]["operating_reserve_modeling_option"] == "individual"# [aggregated, individual]
        bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
        tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
        STOMIN = parameter(am, bus_idx, :gen_bus, tech_idx, "STOMIN")
        BATEFF = sqrt(parameter(am, bus_idx, :gen_bus, tech_idx, "BATEFF"))

        soc_idhty = variable(am, nw, decomp_group, :soc_idhty, (i,d,h,t,y))

        contingency_reserve_min_hr = am.setting["Simulation Configuration"]["contingency_reserve_min_duration_value"]

        sum_res_up = JuMP.AffExpr(0.0)
        if up_reserve_eligible(am, i)
            reserve_enabled(am, :reg_up)  && JuMP.add_to_expression!(sum_res_up, 1, variable(am, nw, decomp_group, :reg_up_idhty, (i,d,h,t,y)))
            reserve_enabled(am, :flex_up) && JuMP.add_to_expression!(sum_res_up, 1, variable(am, nw, decomp_group, :flex_up_idhty, (i,d,h,t,y)))
            reserve_enabled(am, :spin)    && JuMP.add_to_expression!(sum_res_up, contingency_reserve_min_hr, variable(am, nw, decomp_group, :spin_idhty, (i,d,h,t,y)))
        end

        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end
        expr = JuMP.@expression(JuMP_model,
            soc_idhty - (1/BATEFF) * sum_res_up - STOMIN
            )
        constraint = JuMP.@constraint(JuMP_model, 
            expr >= 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end    
    else
        bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
        tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
        STOMIN = parameter(am, bus_idx, :gen_bus, tech_idx, "STOMIN")
        BATEFF = sqrt(parameter(am, bus_idx, :gen_bus, tech_idx, "BATEFF"))

        if STOMIN > 0

            soc_idhty = variable(am, nw, decomp_group, :soc_idhty, (i,d,h,t,y))
            
            if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
            expr = JuMP.@expression(JuMP_model,  
                soc_idhty - STOMIN
                )   
            constraint = JuMP.@constraint(JuMP_model, 
                expr >= 0
                )    
            if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end    
        end
    end

end

function constraint_agg_ES_SOC_Max_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, r::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    group_info = am.ref[:nw][0][:reserve_group_lookup][r]                                         
    unit_group = group_info["unit_group"]                                                         
    z = group_info["zone_idx"]                                                                    
    gen_idx = am.ref[:nw][0][:zone]["reserve"][string(z)]["aggregation_info"]["gen_technology_groups"][unit_group]["gen_idx"]
    
    total_capacity = am.ref[:nw][0][:zone]["reserve"][string(z)]["aggregation_info"]["gen_technology_groups"][unit_group]["total_capacity"]       
    BATEFF = 0.0

    sum_soc = JuMP.AffExpr(0.0)
    sum_u_EUE = JuMP.AffExpr(0.0)
    for i in gen_idx
        JuMP.add_to_expression!(sum_soc, 1, variable(am, nw, decomp_group, :soc_idhty, (i, d, h, t, y)))
        JuMP.add_to_expression!(sum_u_EUE, 1, variable(am, nw, decomp_group, :u_ESE_iy, (i, y)))

        bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
        tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
        BATEFF = sqrt(sqrt(parameter(am, bus_idx, :gen_bus, tech_idx, "BATEFF")))
    end

    reg_dn_idht  = reserve_enabled(am, :reg_dn)  ? variable(am, nw, decomp_group, :reg_dn_idhty, (r,d,h,t,y))  : 0.0
    flex_dn_idht = reserve_enabled(am, :flex_dn) ? variable(am, nw, decomp_group, :flex_dn_idhty, (r,d,h,t,y)) : 0.0

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($r,$d,$h,$t,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        sum_soc + BATEFF * (reg_dn_idht + flex_dn_idht) - sum_u_EUE
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($r,$d,$h,$t,$y)")) end   
end

function constraint_ES_SOC_Max_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    if am.setting["Simulation Configuration"]["operating_reserve_modeling_option"] == "individual"# [aggregated, individual]
        bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
        tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)
        BATEFF = sqrt(parameter(am, bus_idx, :gen_bus, tech_idx, "BATEFF"))

        soc_idhty = variable(am, nw, decomp_group, :soc_idhty, (i,d,h,t,y))
        u_ESE_iy = variable(am, nw, decomp_group, :u_ESE_iy, (i,y))

        sum_res_dn = JuMP.AffExpr(0.0)
        reserve_enabled(am, :reg_dn)  && JuMP.add_to_expression!(sum_res_dn, 1, variable(am, nw, decomp_group, :reg_dn_idhty, (i,d,h,t,y)))
        reserve_enabled(am, :flex_dn) && JuMP.add_to_expression!(sum_res_dn, 1, variable(am, nw, decomp_group, :flex_dn_idhty, (i,d,h,t,y)))

        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end
        expr = JuMP.@expression(JuMP_model,
            soc_idhty + BATEFF * sum_res_dn - u_ESE_iy
            )
        constraint = JuMP.@constraint(JuMP_model,
            expr <= 0
            )
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end

    else
        bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
        tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)
        BATEFF = sqrt(parameter(am, bus_idx, :gen_bus, tech_idx, "BATEFF"))

        soc_idhty = variable(am, nw, decomp_group, :soc_idhty, (i,d,h,t,y))
        u_ESE_iy = variable(am, nw, decomp_group, :u_ESE_iy, (i,y))

        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
        expr = JuMP.@expression(JuMP_model,  
            soc_idhty - u_ESE_iy
            )   
        constraint = JuMP.@constraint(JuMP_model, 
            expr <= 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end    
    end

end

function constraint_RampDn_Inter_SubHour_idhy(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int, day_group::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  

    if t != 1
        HFREQ = am.setting["Simulation Configuration"]["HFREQ"]
        last_subhour = last(am.setting["run_T"])
        
        RDL = parameter(am, bus_idx, :gen_bus, tech_idx, "RDL") # ramp rate %/min
        CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")
        card_T = length(am.setting["run_T"])
        
        g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
        prior_g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t-1,y))
        u_i = variable(am, nw, decomp_group, :u_G_iy, (i,y))
        
        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
        expr = JuMP.@expression(JuMP_model,  
            prior_g_idht - g_idht - (u_i * RDL * CAP * (HFREQ - 1 + 1/card_T) * 60) 
            )   
        constraint = JuMP.@constraint(JuMP_model, 
            expr <= 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end    
    
    elseif (h != 1) && (t == 1)
        HFREQ = am.setting["Simulation Configuration"]["HFREQ"]
        last_subhour = last(am.setting["run_T"])
        
        RDL = parameter(am, bus_idx, :gen_bus, tech_idx, "RDL") # ramp rate %/min
        CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")
        card_T = length(am.setting["run_T"])
        
        g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
        prior_g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h-HFREQ,last_subhour,y))
        u_i = variable(am, nw, decomp_group, :u_G_iy, (i,y))
        
        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
        expr = JuMP.@expression(JuMP_model,  
            prior_g_idht - g_idht - (u_i * RDL * CAP * (HFREQ - 1 + 1/card_T) * 60) 
            )   
        constraint = JuMP.@constraint(JuMP_model, 
            expr <= 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end   
    
    elseif (h == 1) && (t == 1)
        if am.ref[:nw][0][:repdays][d]["Day"] != am.ref[:nw][0][:repday_groups][day_group]["Start_Day_Id"]
            HFREQ = am.setting["Simulation Configuration"]["HFREQ"]
            last_hour = last(am.setting["run_H"])
            last_subhour = last(am.setting["run_T"])
            
            RDL = parameter(am, bus_idx, :gen_bus, tech_idx, "RDL") # ramp rate %/min
            CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")
            card_T = length(am.setting["run_T"])
            
            g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
            prior_g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d-1,last_hour,last_subhour,y))
            u_i = variable(am, nw, decomp_group, :u_G_iy, (i,y))
            
            if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
            expr = JuMP.@expression(JuMP_model,  
                prior_g_idht - g_idht - (u_i * RDL * CAP * (HFREQ - 1 + 1/card_T) * 60) 
                )   
            constraint = JuMP.@constraint(JuMP_model, 
                expr <= 0
                )    
            if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end    
        end
    end
end

function constraint_RampUp_Inter_SubHour_idhy(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int, day_group::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  

    if t != 1
        HFREQ = am.setting["Simulation Configuration"]["HFREQ"]
        last_subhour = last(am.setting["run_T"])
        
        RUL = parameter(am, bus_idx, :gen_bus, tech_idx, "RUL") # ramp rate %/min
        CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")
        card_T = length(am.setting["run_T"])
        
        g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
        prior_g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t-1,y))
        u_i = variable(am, nw, decomp_group, :u_G_iy, (i,y))
        
        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
        expr = JuMP.@expression(JuMP_model,  
            g_idht - prior_g_idht - (u_i * RUL * CAP * ( HFREQ - 1 + 1/card_T) * 60)
            )   
        constraint = JuMP.@constraint(JuMP_model, 
            expr <= 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)"))  end   
    
    elseif (h != 1) && (t == 1)
        HFREQ = am.setting["Simulation Configuration"]["HFREQ"]
        last_subhour = last(am.setting["run_T"])
        
        RUL = parameter(am, bus_idx, :gen_bus, tech_idx, "RUL") # ramp rate %/min
        CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")
        card_T = length(am.setting["run_T"])
        
        g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
        prior_g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h-HFREQ,last_subhour,y))
        u_i = variable(am, nw, decomp_group, :u_G_iy, (i,y))
        
        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
        expr = JuMP.@expression(JuMP_model,  
            g_idht - prior_g_idht - (u_i * RUL * CAP * ( HFREQ - 1 + 1/card_T) * 60)
            )   
        constraint = JuMP.@constraint(JuMP_model, 
            expr <= 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)"))  end   
    
    elseif (h == 1) && (t == 1)
        if am.ref[:nw][0][:repdays][d]["Day"] != am.ref[:nw][0][:repday_groups][day_group]["Start_Day_Id"]
            HFREQ = am.setting["Simulation Configuration"]["HFREQ"]
            last_hour = last(am.setting["run_H"])
            last_subhour = last(am.setting["run_T"])
            
            RUL = parameter(am, bus_idx, :gen_bus, tech_idx, "RUL") # ramp rate %/min
            CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")
            card_T = length(am.setting["run_T"])
            
            g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
            prior_g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d-1,last_hour,last_subhour,y))
            u_i = variable(am, nw, decomp_group, :u_G_iy, (i,y))
            
            if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
            expr = JuMP.@expression(JuMP_model,  
                g_idht - prior_g_idht - (u_i * RUL * CAP * ( HFREQ - 1 + 1/card_T) * 60)
                )   
            constraint = JuMP.@constraint(JuMP_model, 
                expr <= 0
                )    
            if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end    
        end
    end

end

function constraint_VREBalance_Fixed_Profile_wifh_Flexibilty_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int, type::String; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
    CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")
    Timeseries_Tag = parameter(am, bus_idx, :gen_bus, tech_idx, "Timeseries_Tag")
    PMAX = parameter(am, bus_idx, :gen_bus, tech_idx, "PMAX")    
    bus_id = parameter(am, 0, :bus, "bus_i", bus_idx)

    operational_option = parameter(am, bus_idx, :gen_bus, tech_idx, "Dispatch")
    flexibility_percent = am.setting["Simulation Configuration"]["Hydro_Flexibility_Percent"]

    # VRE shape (available shapes: wind_ons, wind_ofs, pv, rtpv, hydro, csp)
    shape = 0.0
    profile_type = parameter(am, bus_idx, :gen_bus, tech_idx, "Profile_Type")
    if Timeseries_Tag == "LOCAL"
        type_key = profile_type * "_shape"
        shape = get_vre_zdt_shape(am.ref[:nw][0], y, bus_idx, d, h, type_key)
    else
        if profile_type in ("wind_ons", "wind_ofs", "csp")
            shape = parameter(am, 0, :planning_stages, "repdays", "data", profile_type, d, h, t, y)[Timeseries_Tag]
        end
    end

    g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
    u_i = variable(am, nw, decomp_group, :u_G_iy, (i,y))

    # Up-reserves only exist for dispatchable units (skip-creation elsewhere).
    up_res_eligible = up_reserve_eligible(am, i)
    reg_up_idht  = (up_res_eligible && reserve_enabled(am, :reg_up))  ? variable(am, nw, decomp_group, :reg_up_idhty, (i,d,h,t,y))  : 0.0
    flex_up_idht = (up_res_eligible && reserve_enabled(am, :flex_up)) ? variable(am, nw, decomp_group, :flex_up_idhty, (i,d,h,t,y)) : 0.0
    spin_idht    = (up_res_eligible && reserve_enabled(am, :spin))    ? variable(am, nw, decomp_group, :spin_idhty, (i,d,h,t,y))    : 0.0


    if operational_option == "Dispatchable"

        # up direction: fixed + flexibility band
        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "up_($i,$d,$h,$t,$y)"))) end
        expr = JuMP.@expression(JuMP_model,
            g_idht + reg_up_idht + flex_up_idht + spin_idht - shape * u_i * CAP * PMAX * (1 + flexibility_percent)
            )
        constraint = JuMP.@constraint(JuMP_model,
            expr <= 0
            )
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end

        # down direction: fixed - flexibility band

        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "down_($i,$d,$h,$t,$y)"))) end
        expr = JuMP.@expression(JuMP_model,
            g_idht + reg_up_idht + flex_up_idht + spin_idht - shape * u_i * CAP * PMAX * (1 - flexibility_percent)
            )
        constraint = JuMP.@constraint(JuMP_model,
            expr >= 0
            )
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end

    elseif operational_option == "Non-dispatchable"

        # No AS procurement from non-dispatchable resources: up-reserve variables
        # are not created for these units, so nothing to fix to 0 here.

        # up direction: fixed + flexibility band
        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "up_($i,$d,$h,$t,$y)"))) end
        expr = JuMP.@expression(JuMP_model,
            g_idht - shape * u_i * CAP * PMAX * (1 + flexibility_percent)
            )
        constraint = JuMP.@constraint(JuMP_model, 
            expr <= 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)"))  end   

        # down direction: fixed - flexibility band
        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "down_($i,$d,$h,$t,$y)"))) end    
        expr = JuMP.@expression(JuMP_model,  
            g_idht - shape * u_i * CAP * PMAX * (1 - flexibility_percent)
            )   
        constraint = JuMP.@constraint(JuMP_model, 
            expr >= 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)"))  end   

    end

end

function constraint_VREBalance_Fixed_Profile_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int, type::String; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
    CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")
    Timeseries_Tag = parameter(am, bus_idx, :gen_bus, tech_idx, "Timeseries_Tag")
    PMAX = parameter(am, bus_idx, :gen_bus, tech_idx, "PMAX")    
    bus_id = parameter(am, 0, :bus, "bus_i", bus_idx)

    operational_option = parameter(am, bus_idx, :gen_bus, tech_idx, "Dispatch")

    # VRE shape (available shapes: wind_ons, wind_ofs, pv, rtpv, hydro, csp)
    shape = 0.0
    profile_type = parameter(am, bus_idx, :gen_bus, tech_idx, "Profile_Type")
    if Timeseries_Tag == "LOCAL"
        type_key = profile_type * "_shape"
        shape = get_vre_zdt_shape(am.ref[:nw][0], y, bus_idx, d, h, type_key)
    else
        if profile_type in ("wind_ons", "wind_ofs", "csp")
            shape = parameter(am, 0, :planning_stages, "repdays", "data", profile_type, d, h, t, y)[Timeseries_Tag]
        end
    end

    # Up-reserves only exist for dispatchable units (skip-creation elsewhere).
    up_res_eligible = up_reserve_eligible(am, i)

    sum_variable = JuMP.AffExpr(0.0)
    JuMP.add_to_expression!(sum_variable, 1, variable(am, nw, decomp_group, :g_idhty, (i, d, h, t, y)))
    if up_res_eligible && am.setting["Simulation Configuration"]["operating_reserve_modeling_option"] == "individual"# [aggregated, individual]
        reserve_enabled(am, :reg_up)  && JuMP.add_to_expression!(sum_variable, 1, variable(am, nw, decomp_group, :reg_up_idhty, (i, d, h, t, y)))
        reserve_enabled(am, :flex_up) && JuMP.add_to_expression!(sum_variable, 1, variable(am, nw, decomp_group, :flex_up_idhty, (i, d, h, t, y)))
        reserve_enabled(am, :spin)    && JuMP.add_to_expression!(sum_variable, 1, variable(am, nw, decomp_group, :spin_idhty, (i, d, h, t, y)))
    end

    g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
    u_i = variable(am, nw, decomp_group, :u_G_iy, (i,y))

    # Numerical conditioning: near-zero-availability hours give a tiny matrix coefficient
    # (~1e-6), so fix output (and up-reserves) to 0 and skip the ill-scaled row.
    if shape * CAP * PMAX < 1e-4
        JuMP.fix(g_idht, 0.0, force=true)
        # Only fix up-reserves that were actually created (dispatchable units).
        if up_res_eligible && am.setting["Simulation Configuration"]["operating_reserve_modeling_option"] == "individual"
            reserve_enabled(am, :reg_up)  && JuMP.fix(variable(am, nw, decomp_group, :reg_up_idhty, (i, d, h, t, y)), 0.0, force=true)
            reserve_enabled(am, :flex_up) && JuMP.fix(variable(am, nw, decomp_group, :flex_up_idhty, (i, d, h, t, y)), 0.0, force=true)
            reserve_enabled(am, :spin)    && JuMP.fix(variable(am, nw, decomp_group, :spin_idhty, (i, d, h, t, y)), 0.0, force=true)
        end
        return
    end


    if operational_option == "Dispatchable"

        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
        expr = JuMP.@expression(JuMP_model,  
            sum_variable - shape * u_i * CAP * PMAX
            )   
        constraint = JuMP.@constraint(JuMP_model, 
            expr <= 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end    

    elseif operational_option == "Non-dispatchable"

        # No AS procurement from non-dispatchable resources: up-reserve variables
        # are not created for these units, so nothing to fix to 0 here.

        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end
        expr = JuMP.@expression(JuMP_model,
            g_idht - shape * u_i * CAP * PMAX
            )
        constraint = JuMP.@constraint(JuMP_model, 
            expr <= 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)"))  end   

    end

end

function constraint_VREBalance_Budget_iy(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, day_group_id::Int, y::Int, type::String; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    ids_ht = [(h,t) for h in am.setting["run_H"] for t in am.setting["run_T"]]
    
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
    CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")
    PMAX = parameter(am, bus_idx, :gen_bus, tech_idx, "PMAX") 

    # water budget of day_group_id in bus_idx (MWh/MW)
    budget = am.ref[:nw][0][:hydro_budget][bus_idx][day_group_id]

    operational_option = parameter(am, bus_idx, :gen_bus, tech_idx, "Dispatch")

    sum_generation = JuMP.AffExpr(0.0)
    u_i = variable(am, nw, decomp_group, :u_G_iy, (i,y))

    if operational_option == "Dispatchable"
        
        for d in am.ref[:nw][0][:repday_groups][day_group_id]["Day_Idx_List"]
            for (h, t) in ids_ht
                JuMP.add_to_expression!(sum_generation, 1.0, variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y)))
                reserve_enabled(am, :reg_up)  && JuMP.add_to_expression!(sum_generation, 1.0, variable(am, nw, decomp_group, :reg_up_idhty, (i,d,h,t,y)))
                reserve_enabled(am, :flex_up) && JuMP.add_to_expression!(sum_generation, 1.0, variable(am, nw, decomp_group, :flex_up_idhty, (i,d,h,t,y)))
                reserve_enabled(am, :spin)    && JuMP.add_to_expression!(sum_generation, 1.0, variable(am, nw, decomp_group, :spin_idhty, (i,d,h,t,y)))
            end
        end


    elseif operational_option == "Non-dispatchable"

        for d in am.ref[:nw][0][:repday_groups][day_group_id]["Day_Idx_List"]
            for (h, t) in ids_ht
                JuMP.add_to_expression!(sum_generation, 1.0, variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y)))
            end
        end

    end

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        sum_generation - budget * u_i * CAP * PMAX
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$y)"))  end   

end

function constraint_VREBalance_Budget_Annual_iy(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, y::Int, ids_d, type::String; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    ids_dhty = [(d,h,t,y) for d in ids_d for h in am.setting["run_H"] for t in am.setting["run_T"]]
    
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
    CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")
    PMAX = parameter(am, bus_idx, :gen_bus, tech_idx, "PMAX")    

    operational_option = parameter(am, bus_idx, :gen_bus, tech_idx, "Dispatch")
    
    budget = 0.0
    for day_group_id in keys(am.ref[:nw][0][:hydro_budget][bus_idx])
        budget += am.ref[:nw][0][:hydro_budget][bus_idx][day_group_id]
    end
    
    total_num_days = 0.0
    for d in ids_d
        total_num_days += parameter(am, 0, :repdays, "NumDays", d)
    end
    scale = total_num_days / 365    # this will be 1 if PH is not used.

    sum_generation = JuMP.AffExpr(0.0)
    for (d, h, t, y) in ids_dhty
        JuMP.add_to_expression!(sum_generation, parameter(am, 0, :repdays, "NumDays", d), variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y)))
    end

    u_i = variable(am, nw, decomp_group, :u_G_iy, (i,y))

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$y)"))) end
    expr = JuMP.@expression(JuMP_model,
        sum_generation - budget * u_i * CAP * PMAX * scale
        )
    constraint = JuMP.@constraint(JuMP_model,
        expr <= 0
        )
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$y)"))  end

end

function constraint_water_use_to_generation_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    ids_j = am.ref[:nw][0][:water_management]["segment_index"]
    
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
    CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")
    PMAX = parameter(am, bus_idx, :gen_bus, tech_idx, "PMAX")    

    sum_water_use = JuMP.AffExpr(0.0)
    for (j) in ids_j
        JuMP.add_to_expression!(sum_water_use, 1.0, variable(am, nw, decomp_group, :water_use_ijdhty, (i,j,d,h,t,y)))        
    end
    g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        sum_water_use - g_idht
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)"))  end   

end

function constraint_water_use_limit_per_daygroup_idy(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, day_group_id::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    ids_j = am.ref[:nw][0][:water_management]["segment_index"]
    ids_d = am.ref[:nw][0][:repday_groups][day_group_id]["Day_Idx_List"]
    ids_jdhty = [(j,d,h,t,y) for j in ids_j for d in ids_d for h in am.setting["run_H"] for t in am.setting["run_T"]]

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
    CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")

    budget = 0.0    # unit = MWh
    start_day = am.ref[:nw][0][:repday_groups][day_group_id]["Start_Day_Id"] + 1 # add 1 because the first row is annual budget
    end_day = am.ref[:nw][0][:repday_groups][day_group_id]["End_Day_Id"] + 1
    region_map = get(am.ref[:nw][0], :profile_data_region_map, nothing)
    # coarse hydro resolution can fold many finest regions onto one column; count each distinct hydro column once
    hydro_cols = unique(profile_data_region(region_map, "hydro", ba_ids)
        for ba_ids in am.ref[:nw][0][:bus][bus_idx]["aggregation_info"]["aggregated_regions_bus_i"])
    for col in hydro_cols
        budget += sum(am.ref[:nw][0][:time_series_data]["hydro_budget_df"][start_day:end_day, col])  # the first row = annual budget
    end

    shape = budget / parameter(am, bus_idx, :gen_bus, tech_idx, "EXCAPS")

    sum_water_use = JuMP.AffExpr(0.0)
    for (j,d,h,t,y) in ids_jdhty
        JuMP.add_to_expression!(sum_water_use, 1.0, variable(am, nw, decomp_group, :water_use_ijdhty, (i,j,d,h,t,y)))        
    end

    u_i = variable(am, nw, decomp_group, :u_G_iy, (i,y))

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$day_group_id,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        sum_water_use - shape * u_i * CAP
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$day_group_id,$y)")) end    

end

function constraint_water_use_limit_per_segment_idy(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, day_group_id::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)
    ids_j = am.ref[:nw][0][:water_management]["segment_index"]
    ids_d = am.ref[:nw][0][:repday_groups][day_group_id]["Day_Idx_List"]
    ids_dhty = [(d,h,t,y) for d in ids_d for h in am.setting["run_H"] for t in am.setting["run_T"]]
    
    water_value_table = am.ref[:nw][0][:hydro_water_value_df]
    reservoir_level_list = water_value_table[(water_value_table.Day_ID .== 1), "Reservoir_Level"] # all regions and days have the same segment size

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
    CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")

    budget = 0.0    # unit = MWh-year
    start_day = am.ref[:nw][0][:repday_groups][day_group_id]["Start_Day_Id"] + 1 # add 1 because the first row is annual budget
    end_day = am.ref[:nw][0][:repday_groups][day_group_id]["End_Day_Id"] + 1
    region_map = get(am.ref[:nw][0], :profile_data_region_map, nothing)
    # coarse hydro resolution can fold many finest regions onto one column; count each distinct hydro column once
    hydro_cols = unique(profile_data_region(region_map, "hydro", ba_ids)
        for ba_ids in am.ref[:nw][0][:bus][bus_idx]["aggregation_info"]["aggregated_regions_bus_i"])
    for col in hydro_cols
        budget += sum(am.ref[:nw][0][:time_series_data]["hydro_budget_df"][start_day:end_day, col])  # the first row = annual budget
    end

    shape = budget / parameter(am, bus_idx, :gen_bus, tech_idx, "EXCAPS")

    u_i = variable(am, nw, decomp_group, :u_G_iy, (i,y))
    
    for j in ids_j
        
        segment_size = 0.0
        if j == length(ids_j)
            segment_size = 0.01 * (reservoir_level_list[j])
        else
            segment_size = 0.01 * (reservoir_level_list[j] - reservoir_level_list[j+1])
        end
        
        sum_water_use = JuMP.AffExpr(0.0)
        for (d,h,t,y) in ids_dhty
            JuMP.add_to_expression!(sum_water_use, 1.0, variable(am, nw, decomp_group, :water_use_ijdhty, (i,j,d,h,t,y)))           
        end

        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($j,$i,$day_group_id,$y)"))) end    
        expr = JuMP.@expression(JuMP_model,  
            sum_water_use - segment_size * shape * u_i * CAP
            )   
        constraint = JuMP.@constraint(JuMP_model, 
            expr <= 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($j,$i,$day_group_id,$y)"))  end   
    end

end

function constraint_water_use_segment_order_idy(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, day_group_id::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)
    ids_j = am.ref[:nw][0][:water_management]["segment_index"]
    ids_d = sort(am.ref[:nw][0][:repday_groups][day_group_id]["Day_Idx_List"])
    ids_dhty = [(d,h,t,y) for d in ids_d for h in am.setting["run_H"] for t in am.setting["run_T"]]
    ids_hty = [(h,t,y) for h in am.setting["run_H"] for t in am.setting["run_T"]]
    
    water_value_table = am.ref[:nw][0][:hydro_water_value_df]
    reservoir_level_list = water_value_table[(water_value_table.Day_ID .== 1), "Reservoir_Level"] # all regions and days have the same segment size
    
    
    for day_idx in eachindex(ids_d)
        
        if day_idx != 1
            d = ids_d[day_idx]

            for j in ids_j
                
                prior_sum_water_use = JuMP.AffExpr(0.0)
                current_sum_water_use = JuMP.AffExpr(0.0)
                for (h,t,y) in ids_hty
                    JuMP.add_to_expression!(prior_sum_water_use, 1.0, variable(am, nw, decomp_group, :water_use_ijdhty, (i,j,d-1,h,t,y)))           
                    JuMP.add_to_expression!(current_sum_water_use, 1.0, variable(am, nw, decomp_group, :water_use_ijdhty, (i,j,d,h,t,y)))           
                end

                if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($j,$i,$d,$y)"))) end    
                expr = JuMP.@expression(JuMP_model,  
                    current_sum_water_use - prior_sum_water_use
                    )   
                constraint = JuMP.@constraint(JuMP_model, 
                    expr <= 0
                    )    
                if const_name_flag JuMP.set_name(constraint, string(const_name, "_($j,$i,$d,$y)"))  end   
            end
        end
    end

end

function constraint_RampDnMax_UC_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)   
    PMAX = parameter(am, bus_idx, :gen_bus, tech_idx, "PMAX")
    RDL = parameter(am, bus_idx, :gen_bus, tech_idx, "RDL")   # 10-mins ramp rate in %
    CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")
    MAXC = parameter(am, bus_idx, :gen_bus, tech_idx, "MAXC")  # Max contingency reserve percent
    Res_cap = min(RDL, MAXC)

    reg_dn_idht  = reserve_enabled(am, :reg_dn)  ? variable(am, nw, decomp_group, :reg_dn_idhty, (i,d,h,t,y))  : 0.0
    flex_dn_idht = reserve_enabled(am, :flex_dn) ? variable(am, nw, decomp_group, :flex_dn_idhty, (i,d,h,t,y)) : 0.0

    c_idht = variable(am, nw, decomp_group, :c_idhty, (i,d,h,t,y)) 

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        reg_dn_idht + flex_dn_idht - c_idht * PMAX * CAP * Res_cap
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end    

end

function constraint_RegUp_Max_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)
    PMAX = param_f64(am, bus_idx, :gen_bus, tech_idx, "PMAX")
    MAXR = param_f64(am, bus_idx, :gen_bus, tech_idx, "MAXR")   # 5 min ramp limit
    CAP = param_f64(am, bus_idx, :gen_bus, tech_idx, "CAP")

    reg_up_idht = variable(am, nw, decomp_group, :reg_up_idhty, (i,d,h,t,y))
    u_i = variable(am, nw, decomp_group, :u_G_iy, (i,y))

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        reg_up_idht - u_i * PMAX * CAP * MAXR 
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end    
end

function constraint_Agg_RegUp_Max_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, r::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    group_info = am.ref[:nw][0][:reserve_group_lookup][r]                                         
    unit_group = group_info["unit_group"]                                                         
    z = group_info["zone_idx"]                                                                    
    gen_idx = am.ref[:nw][0][:zone]["reserve"][string(z)]["aggregation_info"]["gen_technology_groups"][unit_group]["gen_idx"]

    regulation_fraction = am.ref[:nw][0][:zone]["reserve"][string(z)]["aggregation_info"]["gen_technology_groups"][unit_group]["regulation_fraction"]       
    total_capacity = am.ref[:nw][0][:zone]["reserve"][string(z)]["aggregation_info"]["gen_technology_groups"][unit_group]["total_capacity"]       

    reg_up_idht = variable(am, nw, decomp_group, :reg_up_idhty, (r,d,h,t,y))
    sum_u = JuMP.AffExpr(0.0)
    for i in gen_idx
        JuMP.add_to_expression!(sum_u, 1, variable(am, nw, decomp_group, :u_G_iy, (i, y)))
    end

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($r,$d,$h,$t,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        reg_up_idht - regulation_fraction * total_capacity * sum_u 
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($r,$d,$h,$t,$y)")) end    
end

function constraint_RegDn_Max_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)
    PMAX = param_f64(am, bus_idx, :gen_bus, tech_idx, "PMAX")
    MAXR = param_f64(am, bus_idx, :gen_bus, tech_idx, "MAXR")   # 5 min ramp limit
    CAP = param_f64(am, bus_idx, :gen_bus, tech_idx, "CAP")

    reg_dn_idht = variable(am, nw, decomp_group, :reg_dn_idhty, (i,d,h,t,y))
    u_i = variable(am, nw, decomp_group, :u_G_iy, (i,y))

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        reg_dn_idht - u_i * PMAX * CAP * MAXR 
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)"))  end   

end

function constraint_Agg_RegDn_Max_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, r::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    group_info = am.ref[:nw][0][:reserve_group_lookup][r]                                         
    unit_group = group_info["unit_group"]                                                         
    z = group_info["zone_idx"]                                                                    
    gen_idx = am.ref[:nw][0][:zone]["reserve"][string(z)]["aggregation_info"]["gen_technology_groups"][unit_group]["gen_idx"]

    regulation_fraction = am.ref[:nw][0][:zone]["reserve"][string(z)]["aggregation_info"]["gen_technology_groups"][unit_group]["regulation_fraction"]       
    total_capacity = am.ref[:nw][0][:zone]["reserve"][string(z)]["aggregation_info"]["gen_technology_groups"][unit_group]["total_capacity"]       

    reg_dn_idht = variable(am, nw, decomp_group, :reg_dn_idhty, (r,d,h,t,y))
    
    sum_u = JuMP.AffExpr(0.0)
    for i in gen_idx
        JuMP.add_to_expression!(sum_u, 1, variable(am, nw, decomp_group, :u_G_iy, (i, y)))
    end

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($r,$d,$h,$t,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        reg_dn_idht - regulation_fraction * total_capacity * sum_u
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($r,$d,$h,$t,$y)"))  end   

end

function constraint_RampUpMax_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)
    PMAX = param_f64(am, bus_idx, :gen_bus, tech_idx, "PMAX")
    RUL = param_f64(am, bus_idx, :gen_bus, tech_idx, "RUL")     # 10 min ramp limit
    CAP = param_f64(am, bus_idx, :gen_bus, tech_idx, "CAP")
    MAXC = param_f64(am, bus_idx, :gen_bus, tech_idx, "MAXC")  # Max contingency reserve percent
    Res_cap = min(RUL, MAXC)

    # Up-reserves only exist for dispatchable units.
    up_res_eligible = up_reserve_eligible(am, i)
    reg_up_idht  = (up_res_eligible && reserve_enabled(am, :reg_up))  ? variable(am, nw, decomp_group, :reg_up_idhty, (i,d,h,t,y))  : 0.0
    flex_up_idht = (up_res_eligible && reserve_enabled(am, :flex_up)) ? variable(am, nw, decomp_group, :flex_up_idhty, (i,d,h,t,y)) : 0.0
    spin_idht    = (up_res_eligible && reserve_enabled(am, :spin))    ? variable(am, nw, decomp_group, :spin_idhty, (i,d,h,t,y))    : 0.0

    u_i = variable(am, nw, decomp_group, :u_G_iy, (i,y))

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end
    expr = JuMP.@expression(JuMP_model,
        reg_up_idht + flex_up_idht + spin_idht - u_i * PMAX * CAP * Res_cap
        )
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end    

end

function constraint_Agg_RampUpMax_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, r::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    group_info = am.ref[:nw][0][:reserve_group_lookup][r]                                         
    unit_group = group_info["unit_group"]                                                         
    z = group_info["zone_idx"]                                                                    
    gen_idx = am.ref[:nw][0][:zone]["reserve"][string(z)]["aggregation_info"]["gen_technology_groups"][unit_group]["gen_idx"]

    reserve_up_fraction = am.ref[:nw][0][:zone]["reserve"][string(z)]["aggregation_info"]["gen_technology_groups"][unit_group]["reserve_up_fraction"]       
    total_capacity = am.ref[:nw][0][:zone]["reserve"][string(z)]["aggregation_info"]["gen_technology_groups"][unit_group]["total_capacity"]       

    reg_up_idht  = reserve_enabled(am, :reg_up)  ? variable(am, nw, decomp_group, :reg_up_idhty, (r,d,h,t,y))  : 0.0
    flex_up_idht = reserve_enabled(am, :flex_up) ? variable(am, nw, decomp_group, :flex_up_idhty, (r,d,h,t,y)) : 0.0
    spin_idht    = reserve_enabled(am, :spin)    ? variable(am, nw, decomp_group, :spin_idhty, (r,d,h,t,y))    : 0.0

    sum_u = JuMP.AffExpr(0.0)
    for i in gen_idx
        JuMP.add_to_expression!(sum_u, 1, variable(am, nw, decomp_group, :u_G_iy, (i, y)))
    end

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($r,$d,$h,$t,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        reg_up_idht + flex_up_idht + spin_idht - reserve_up_fraction * total_capacity * sum_u
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($r,$d,$h,$t,$y)")) end    

end

function constraint_RampDnMax_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)
    PMAX = param_f64(am, bus_idx, :gen_bus, tech_idx, "PMAX")
    RDL = param_f64(am, bus_idx, :gen_bus, tech_idx, "RDL")  # 10 min ramp rate %
    CAP = param_f64(am, bus_idx, :gen_bus, tech_idx, "CAP")
    MAXC = param_f64(am, bus_idx, :gen_bus, tech_idx, "MAXC")  # Max contingency reserve percent
    Res_cap = min(RDL, MAXC)

    reg_dn_idht  = reserve_enabled(am, :reg_dn)  ? variable(am, nw, decomp_group, :reg_dn_idhty, (i,d,h,t,y))  : 0.0
    flex_dn_idht = reserve_enabled(am, :flex_dn) ? variable(am, nw, decomp_group, :flex_dn_idhty, (i,d,h,t,y)) : 0.0

    u_i = variable(am, nw, decomp_group, :u_G_iy, (i,y))

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        reg_dn_idht + flex_dn_idht - u_i * PMAX * CAP * Res_cap
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end    

end


function constraint_Agg_RampDnMax_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, r::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    group_info = am.ref[:nw][0][:reserve_group_lookup][r]                                         
    unit_group = group_info["unit_group"]                                                         
    z = group_info["zone_idx"]                                                                    
    gen_idx = am.ref[:nw][0][:zone]["reserve"][string(z)]["aggregation_info"]["gen_technology_groups"][unit_group]["gen_idx"]

    reserve_down_fraction = am.ref[:nw][0][:zone]["reserve"][string(z)]["aggregation_info"]["gen_technology_groups"][unit_group]["reserve_down_fraction"]       
    total_capacity = am.ref[:nw][0][:zone]["reserve"][string(z)]["aggregation_info"]["gen_technology_groups"][unit_group]["total_capacity"]       
    
    reg_dn_idht  = reserve_enabled(am, :reg_dn)  ? variable(am, nw, decomp_group, :reg_dn_idhty, (r,d,h,t,y))  : 0.0
    flex_dn_idht = reserve_enabled(am, :flex_dn) ? variable(am, nw, decomp_group, :flex_dn_idhty, (r,d,h,t,y)) : 0.0

    sum_u = JuMP.AffExpr(0.0)
    for i in gen_idx
        JuMP.add_to_expression!(sum_u, 1, variable(am, nw, decomp_group, :u_G_iy, (i, y)))
    end

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($r,$d,$h,$t,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        reg_dn_idht + flex_dn_idht - total_capacity * sum_u * reserve_down_fraction
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($r,$d,$h,$t,$y)")) end    

end

function constraint_RampUpMax_UC_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)   
    PMAX = parameter(am, bus_idx, :gen_bus, tech_idx, "PMAX")
    RUL = parameter(am, bus_idx, :gen_bus, tech_idx, "RUL")  # 10-mins ramp rate in %
    CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")
    MAXC = parameter(am, bus_idx, :gen_bus, tech_idx, "MAXC")  # Max contingency reserve percent
    Res_cap = min(RUL, MAXC)

    # Up-reserves only exist for dispatchable units.
    up_res_eligible = up_reserve_eligible(am, i)
    reg_up_idht  = (up_res_eligible && reserve_enabled(am, :reg_up))  ? variable(am, nw, decomp_group, :reg_up_idhty, (i,d,h,t,y))  : 0.0
    flex_up_idht = (up_res_eligible && reserve_enabled(am, :flex_up)) ? variable(am, nw, decomp_group, :flex_up_idhty, (i,d,h,t,y)) : 0.0
    spin_idht    = (up_res_eligible && reserve_enabled(am, :spin))    ? variable(am, nw, decomp_group, :spin_idhty, (i,d,h,t,y))    : 0.0

    c_idht = variable(am, nw, decomp_group, :c_idhty, (i,d,h,t,y))

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        reg_up_idht + flex_up_idht + spin_idht - c_idht * PMAX * CAP * Res_cap
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end    

end

function constraint_MustRun_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)
    PMAX = param_f64(am, bus_idx, :gen_bus, tech_idx, "PMAX")
    CAP = param_f64(am, bus_idx, :gen_bus, tech_idx, "CAP")
    MustRun_Flag = param_bool(am, bus_idx, :gen_bus, tech_idx, "Must_Run_Flag")
    MustRun_Level = param_f64(am, bus_idx, :gen_bus, tech_idx, "Must_Run_Level")

    if MustRun_Level > PMAX
        MustRun_Level = PMAX
    end

    if MustRun_Flag == true
    
        g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
        u_i = variable(am, nw, decomp_group, :u_G_iy, (i,y))

        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
        expr = JuMP.@expression(JuMP_model,  
            g_idht - u_i * MustRun_Level * CAP
            )   
        constraint = JuMP.@constraint(JuMP_model, 
            expr >= 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)"))  end   

    end

end

function constraint_TherMax_Agg_Res_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, r::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    group_info = am.ref[:nw][0][:reserve_group_lookup][r]                                         
    unit_group = group_info["unit_group"]                                                         
    z = group_info["zone_idx"]                                                                    
    gen_idx = am.ref[:nw][0][:zone]["reserve"][string(z)]["aggregation_info"]["gen_technology_groups"][unit_group]["gen_idx"]
        
    total_capacity = am.ref[:nw][0][:zone]["reserve"][string(z)]["aggregation_info"]["gen_technology_groups"][unit_group]["total_capacity"]  
    
    sum_g = JuMP.AffExpr(0.0)
    sum_u = JuMP.AffExpr(0.0)
    for i in gen_idx
        JuMP.add_to_expression!(sum_g, 1, variable(am, nw, decomp_group, :g_idhty, (i, d, h, t, y)))
        JuMP.add_to_expression!(sum_u, 1, variable(am, nw, decomp_group, :u_G_iy, (i, y)))
    end

    sum_res = JuMP.AffExpr(0.0)
    reserve_enabled(am, :reg_up)  && JuMP.add_to_expression!(sum_res, 1, variable(am, nw, decomp_group, :reg_up_idhty, (r, d, h, t, y)))
    reserve_enabled(am, :flex_up) && JuMP.add_to_expression!(sum_res, 1, variable(am, nw, decomp_group, :flex_up_idhty, (r, d, h, t, y)))
    reserve_enabled(am, :spin)    && JuMP.add_to_expression!(sum_res, 1, variable(am, nw, decomp_group, :spin_idhty, (r, d, h, t, y)))

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($r,$d,$h,$t,$y)"))) end
    expr = JuMP.@expression(JuMP_model,
        sum_g + sum_res - sum_u * total_capacity
        )
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($r,$d,$h,$t,$y)"))  end   
end


function constraint_TherMaxDispatch_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)
    PMAX = param_f64(am, bus_idx, :gen_bus, tech_idx, "PMAX")
    CAP = param_f64(am, bus_idx, :gen_bus, tech_idx, "CAP")

    sum_variable = JuMP.AffExpr(0.0)
    JuMP.add_to_expression!(sum_variable, 1, variable(am, nw, decomp_group, :g_idhty, (i, d, h, t, y)))
    if up_reserve_eligible(am, i) && am.setting["Simulation Configuration"]["operating_reserve_modeling_option"] == "individual"# [aggregated, individual]
        reserve_enabled(am, :reg_up)  && JuMP.add_to_expression!(sum_variable, 1, variable(am, nw, decomp_group, :reg_up_idhty, (i, d, h, t, y)))
        reserve_enabled(am, :flex_up) && JuMP.add_to_expression!(sum_variable, 1, variable(am, nw, decomp_group, :flex_up_idhty, (i, d, h, t, y)))
        reserve_enabled(am, :spin)    && JuMP.add_to_expression!(sum_variable, 1, variable(am, nw, decomp_group, :spin_idhty, (i, d, h, t, y)))
    end

    u_i = variable(am, nw, decomp_group, :u_G_iy, (i,y))

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end
    expr = JuMP.@expression(JuMP_model,
        sum_variable - u_i * PMAX * CAP
        )
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)"))  end   
end

function constraint_TherMax_Dispatch_UC_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)    
    PMAX = parameter(am, bus_idx, :gen_bus, tech_idx, "PMAX")   # %
    CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP") # MW
    
    # Up-reserves only exist for dispatchable units.
    up_res_eligible = up_reserve_eligible(am, i)
    g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
    reg_up_idht  = (up_res_eligible && reserve_enabled(am, :reg_up))  ? variable(am, nw, decomp_group, :reg_up_idhty, (i,d,h,t,y))  : 0.0
    flex_up_idht = (up_res_eligible && reserve_enabled(am, :flex_up)) ? variable(am, nw, decomp_group, :flex_up_idhty, (i,d,h,t,y)) : 0.0
    spin_idht    = (up_res_eligible && reserve_enabled(am, :spin))    ? variable(am, nw, decomp_group, :spin_idhty, (i,d,h,t,y))    : 0.0
    c_idht = variable(am, nw, decomp_group, :c_idhty, (i,d,h,t,y))

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end
    expr = JuMP.@expression(JuMP_model,
        g_idht + reg_up_idht + flex_up_idht + spin_idht - c_idht * PMAX * CAP
        )
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)"))  end   

end

function constraint_TherMax_Dispatch_UC_NoReserve_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)    
    CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP") # MW
    
    g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
    c_idht = variable(am, nw, decomp_group, :c_idhty, (i,d,h,t,y)) 

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        g_idht - c_idht * CAP
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)"))  end   

end

function constraint_TherMinDispatch_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    # PMIN fixed at 0
    PMIN = 0.0

    g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
    reg_dn_idht  = reserve_enabled(am, :reg_dn)  ? variable(am, nw, decomp_group, :reg_dn_idhty, (i,d,h,t,y))  : 0.0
    flex_dn_idht = reserve_enabled(am, :flex_dn) ? variable(am, nw, decomp_group, :flex_dn_idhty, (i,d,h,t,y)) : 0.0
    
    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        g_idht - reg_dn_idht - flex_dn_idht - PMIN
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr >= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end    
    
end

function constraint_TherMin_Dispatch_UC_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)   
    PMIN = parameter(am, bus_idx, :gen_bus, tech_idx, "PMIN")   # %
    CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP") # MW

    g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
    reg_dn_idht  = reserve_enabled(am, :reg_dn)  ? variable(am, nw, decomp_group, :reg_dn_idhty, (i,d,h,t,y))  : 0.0
    flex_dn_idht = reserve_enabled(am, :flex_dn) ? variable(am, nw, decomp_group, :flex_dn_idhty, (i,d,h,t,y)) : 0.0
    c_idht = variable(am, nw, decomp_group, :c_idhty, (i,d,h,t,y)) 
    
    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        g_idht - reg_dn_idht - flex_dn_idht - PMIN * CAP * c_idht
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr >= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end    
    
end

function constraint_TherMin_Dispatch_UC_NoReserve_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)   
    PMIN = parameter(am, bus_idx, :gen_bus, tech_idx, "PMIN")   # %
    CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP") # MW

    g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
    c_idht = variable(am, nw, decomp_group, :c_idhty, (i,d,h,t,y)) 
    
    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        g_idht - PMIN * CAP * c_idht
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr >= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end    
    
end

function constraint_Commit_Limit_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)


    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)
    u_G_iy = parameter(am, bus_idx, :gen_bus, tech_idx, "EXUNITS")

    c_idht = variable(am, nw, decomp_group, :c_idhty, (i,d,h,t,y))
    su_idht = variable(am, nw, decomp_group, :su_idhty, (i,d,h,t,y))

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "c_($i,$d,$h,$t,$y)"))) end
    expr = JuMP.@expression(JuMP_model,
        c_idht - u_G_iy
        )
    constraint = JuMP.@constraint(JuMP_model,
        expr <= 0
        )
    if const_name_flag JuMP.set_name(constraint, string(const_name, "c_($i,$d,$h,$t,$y)")) end

end

function constraint_flex_Dn_zdhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, z::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=true, rns_flag::Bool=true)

    ids_r = am.ref[:nw][0][:zone]["reserve"][string(z)]["aggregation_info"]["new_local_gen_idx"]
    if am.setting["Simulation Configuration"]["operating_reserve_modeling_option"] == "aggregated"
        ids_r = sort(collect(keys(am.ref[:nw][0][:reserve_group_lookup])))
    end

    ReqR = 0.0
    for ba_id in am.ref[:nw][0][:zone]["reserve"][string(z)]["aggregation_info"]["aggregated_data_regions_id"]
        ReqR += parameter(am, 0, :planning_stages, "repdays", "data", "flex_down", d, h, t, y)[ba_id]
    end

    sum_flex_dn_idht = JuMP.AffExpr(0.0)
    for i in ids_r
        JuMP.add_to_expression!(sum_flex_dn_idht, 1, variable(am, nw, decomp_group, :flex_dn_idhty, (i, d, h, t, y)))
    end
    rns_flex_down_zdht = variable(am, nw, decomp_group, :rns_flex_dn_zdhty, (z, d, h, t, y))
    
    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($z,$d,$h,$t,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        sum_flex_dn_idht + rns_flex_down_zdht - ReqR 
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr >= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($z,$d,$h,$t,$y)")) end    

    if rns_flag == true
        JuMP.set_upper_bound(rns_flex_down_zdht, ReqR)
    end
end

function constraint_flex_Up_zdhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, z::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=true, rns_flag::Bool=true)

    ids_r = am.ref[:nw][0][:zone]["reserve"][string(z)]["aggregation_info"]["new_local_gen_idx"]
    if am.setting["Simulation Configuration"]["operating_reserve_modeling_option"] == "aggregated"
        ids_r = sort(collect(keys(am.ref[:nw][0][:reserve_group_lookup])))
    end

    ReqR = 0.0
    for ba_id in am.ref[:nw][0][:zone]["reserve"][string(z)]["aggregation_info"]["aggregated_data_regions_id"]
        ReqR += parameter(am, 0, :planning_stages, "repdays", "data", "flex_up", d, h, t, y)[ba_id]
    end

    # Individual mode: up-reserves only exist for dispatchable units.
    individual_mode = am.setting["Simulation Configuration"]["operating_reserve_modeling_option"] == "individual"
    sum_flex_up_idht = JuMP.AffExpr(0.0)
    for i in ids_r
        (individual_mode && !up_reserve_eligible(am, i)) && continue
        JuMP.add_to_expression!(sum_flex_up_idht, 1, variable(am, nw, decomp_group, :flex_up_idhty, (i, d, h, t, y)))
    end
    rns_flex_up_zdht = variable(am, nw, decomp_group, :rns_flex_up_zdhty, (z, d, h, t, y))
    
    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($z,$d,$h,$t,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        sum_flex_up_idht + rns_flex_up_zdht - ReqR
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr >= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($z,$d,$h,$t,$y)")) end    

    if rns_flag == true
        JuMP.set_upper_bound(rns_flex_up_zdht, ReqR)
    end
end

function constraint_R_Cont_zdhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, z::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=true, rns_flag::Bool=true)

    ids_r = am.ref[:nw][0][:zone]["reserve"][string(z)]["aggregation_info"]["new_local_gen_idx"]
    if am.setting["Simulation Configuration"]["operating_reserve_modeling_option"] == "aggregated"
        ids_r = sort(collect(keys(am.ref[:nw][0][:reserve_group_lookup])))
    end

    ReqR_spin = 0.0
    ReqR_nspin = 0.0
    for ba_id in am.ref[:nw][0][:zone]["reserve"][string(z)]["aggregation_info"]["aggregated_data_regions_id"]
        ReqR_spin += parameter(am, 0, :planning_stages, "repdays", "data", "spin", d, h, t, y)[ba_id]
        ReqR_nspin += parameter(am, 0, :planning_stages, "repdays", "data", "nspin", d, h, t, y)[ba_id]
    end

    # Individual mode: up-reserves only exist for dispatchable units.
    individual_mode = am.setting["Simulation Configuration"]["operating_reserve_modeling_option"] == "individual"
    sum_cont_idht = JuMP.AffExpr(0.0)
    for i in ids_r
        (individual_mode && !up_reserve_eligible(am, i)) && continue
        JuMP.add_to_expression!(sum_cont_idht, 1, variable(am, nw, decomp_group, :spin_idhty, (i, d, h, t, y)))
    end
    rns_cont_zdht = variable(am, nw, decomp_group, :rns_cont_zdhty, (z, d, h, t, y))

    # Demand-side (Load Resource) contingency reserve provision; off (ub=0) when fraction=0.
    demand_reserve_fraction = get(am.setting["Simulation Configuration"], "demand_reserve_provision_fraction", 0.0)
    demand_reserve_zdht = variable(am, nw, decomp_group, :demand_reserve_zdhty, (z, d, h, t, y))

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($z,$d,$h,$t,$y)"))) end
    expr = JuMP.@expression(JuMP_model,
        sum_cont_idht + demand_reserve_zdht + rns_cont_zdht - (ReqR_spin + ReqR_nspin)
        )
    constraint = JuMP.@constraint(JuMP_model,
        expr >= 0
        )
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($z,$d,$h,$t,$y)")) end

    JuMP.set_upper_bound(demand_reserve_zdht, demand_reserve_fraction * (ReqR_spin + ReqR_nspin))

    if rns_flag == true
        JuMP.set_upper_bound(rns_cont_zdht, (ReqR_spin + ReqR_nspin))
    end
end

function constraint_R_Spin_zdhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, z::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=true, rns_flag::Bool=true)

    ids_r = am.ref[:nw][0][:zone]["reserve"][string(z)]["aggregation_info"]["new_local_gen_idx"]
    if am.setting["Simulation Configuration"]["operating_reserve_modeling_option"] == "aggregated"
        ids_r = sort(collect(keys(am.ref[:nw][0][:reserve_group_lookup])))
    end

    ReqR_spin = 0.0
    for ba_id in am.ref[:nw][0][:zone]["reserve"][string(z)]["aggregation_info"]["aggregated_data_regions_id"]
        ReqR_spin += parameter(am, 0, :planning_stages, "repdays", "data", "spin", d, h, t, y)[ba_id]
    end

    # Individual mode: up-reserves only exist for dispatchable units.
    individual_mode = am.setting["Simulation Configuration"]["operating_reserve_modeling_option"] == "individual"
    sum_cont_idht = JuMP.AffExpr(0.0)
    for i in ids_r
        (individual_mode && !up_reserve_eligible(am, i)) && continue
        JuMP.add_to_expression!(sum_cont_idht, 1, variable(am, nw, decomp_group, :spin_idhty, (i, d, h, t, y)))
    end
    rns_cont_zdht = variable(am, nw, decomp_group, :rns_cont_zdhty, (z, d, h, t, y))

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($z,$d,$h,$t,$y)"))) end
    expr = JuMP.@expression(JuMP_model,  
        sum_cont_idht + rns_cont_zdht - (ReqR_spin)
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr >= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($z,$d,$h,$t,$y)"))  end   

    if rns_flag == true
        JuMP.set_upper_bound(rns_cont_zdht, (ReqR_spin))
    end
end

function constraint_R_NonSpin_zdhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, z::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=true, rns_flag::Bool=true)

    ids_r = am.ref[:nw][0][:zone]["reserve"][string(z)]["aggregation_info"]["new_local_gen_idx"]
    if am.setting["Simulation Configuration"]["operating_reserve_modeling_option"] == "aggregated"
        ids_r = sort(collect(keys(am.ref[:nw][0][:reserve_group_lookup])))
    end

    ReqR_nspin = 0.0
    for ba_id in am.ref[:nw][0][:zone]["reserve"][string(z)]["aggregation_info"]["aggregated_data_regions_id"]
        ReqR_nspin += parameter(am, 0, :planning_stages, "repdays", "data", "nspin", d, h, t, y)[ba_id]
    end

    sum_cont_idht = JuMP.AffExpr(0.0)
    for i in ids_r
        JuMP.add_to_expression!(sum_cont_idht, 1, variable(am, nw, decomp_group, :nonspin_idhty, (i, d, h, t, y)))
    end
    rns_nonspin_zdhty = variable(am, nw, decomp_group, :rns_nonspin_zdhty, (z, d, h, t, y))
    
    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($z,$d,$h,$t,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        sum_cont_idht + rns_nonspin_zdhty - (ReqR_nspin)
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr >= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($z,$d,$h,$t,$y)")) end    

    if rns_flag == true
        JuMP.set_upper_bound(rns_nonspin_zdhty, (ReqR_nspin))
    end
end

function constraint_reg_Dn_zdhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, z::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=true, rns_flag::Bool=true)

    ids_r = am.ref[:nw][0][:zone]["reserve"][string(z)]["aggregation_info"]["new_local_gen_idx"]
    if am.setting["Simulation Configuration"]["operating_reserve_modeling_option"] == "aggregated"
        ids_r = sort(collect(keys(am.ref[:nw][0][:reserve_group_lookup])))
    end

    ReqR = 0.0
    for ba_id in am.ref[:nw][0][:zone]["reserve"][string(z)]["aggregation_info"]["aggregated_data_regions_id"]
        ReqR += parameter(am, 0, :planning_stages, "repdays", "data", "reg_down", d, h, t, y)[ba_id]
    end

    rns_reg_dn_zdhty = variable(am, nw, decomp_group, :rns_reg_dn_zdhty, (z, d, h, t, y))
    sum_reg_dn_idht = JuMP.AffExpr(0.0)
    for i in ids_r
        JuMP.add_to_expression!(sum_reg_dn_idht, 1, variable(am, nw, decomp_group, :reg_dn_idhty, (i, d, h, t, y)))
    end
    
    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($z,$d,$h,$t,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        sum_reg_dn_idht + rns_reg_dn_zdhty - ReqR
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr >= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($z,$d,$h,$t,$y)")) end    

    if rns_flag == true
        JuMP.set_upper_bound(rns_reg_dn_zdhty, (ReqR))
    end
end

function constraint_reg_Up_zdhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, z::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=true, rns_flag::Bool=true)

    ids_r = am.ref[:nw][0][:zone]["reserve"][string(z)]["aggregation_info"]["new_local_gen_idx"]
    if am.setting["Simulation Configuration"]["operating_reserve_modeling_option"] == "aggregated"
        ids_r = sort(collect(keys(am.ref[:nw][0][:reserve_group_lookup])))
    end

    ReqR = 0.0
    for ba_id in am.ref[:nw][0][:zone]["reserve"][string(z)]["aggregation_info"]["aggregated_data_regions_id"]
        ReqR += parameter(am, 0, :planning_stages, "repdays", "data", "reg_up", d, h, t, y)[ba_id] 
    end

    # Individual mode: up-reserves only exist for dispatchable units.
    individual_mode = am.setting["Simulation Configuration"]["operating_reserve_modeling_option"] == "individual"
    rns_reg_up_zdhty = variable(am, nw, decomp_group, :rns_reg_up_zdhty, (z, d, h, t, y))
    sum_reg_up_idht = JuMP.AffExpr(0.0)
    for i in ids_r
        (individual_mode && !up_reserve_eligible(am, i)) && continue
        JuMP.add_to_expression!(sum_reg_up_idht, 1, variable(am, nw, decomp_group, :reg_up_idhty, (i, d, h, t, y)))
    end
    
    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($z,$d,$h,$t,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,
        sum_reg_up_idht + rns_reg_up_zdhty - ReqR
        )
    constraint = JuMP.@constraint(JuMP_model,
        expr >= 0
        )
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($z,$d,$h,$t,$y)")) end    

    if rns_flag == true
        JuMP.set_upper_bound(rns_reg_up_zdhty, (ReqR))
    end
    
end

function constraint_power_flow_min_kdhty_OP(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, k::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false, expansion_record_flag::Bool = false, result::Dict{String,<:Any} = Dict{String,Any}())

    rate_a = parameter(am, 0, :branch, "rate_a", k) 

    total_line_expansion = 0.0
    if expansion_record_flag == true    # only when we have expansion decisions
        total_line_expansion = result["solution"]["expansion"][string("(", k, ", ", y, ")")]["u_T_ky"]
    end

    f_kdhty = variable(am, nw, decomp_group, :f_kdhty, (k, d, h, t, y))

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($k,$d,$h,$t,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        f_kdhty + rate_a * (1 + total_line_expansion)
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr >= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($k,$d,$h,$t,$y)")) end    
end

function constraint_power_flow_max_kdhty_OP(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, k::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false, expansion_record_flag::Bool = false, result::Dict{String,<:Any} = Dict{String,Any}())

    rate_a = parameter(am, 0, :branch, "rate_a", k) 

    total_line_expansion = 0.0
    if expansion_record_flag == true    # only when we have expansion decisions
        total_line_expansion = result["solution"]["expansion"][string("(", k, ", ", y, ")")]["u_T_ky"]
    end

    f_kdhty = variable(am, nw, decomp_group, :f_kdhty, (k, d, h, t, y))

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($k,$d,$h,$t,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        f_kdhty - rate_a * (1 + total_line_expansion)
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($k,$d,$h,$t,$y)")) end    
end

function constraint_power_flow_min_kdhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, k::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false, GTEP_multi_round_info::Dict{Any,<:Any} = Dict{Any,Any}())


    rate_a = parameter(am, 0, :branch, "rate_a", k) 

    f_kdhty = variable(am, nw, decomp_group, :f_kdhty, (k, d, h, t, y))
    u_T_ky = variable(am, nw, decomp_group, :u_T_ky, (k, y))

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($k,$d,$h,$t,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        f_kdhty + rate_a * (1 + u_T_ky)
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr >= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($k,$d,$h,$t,$y)")) end    
end

function constraint_power_flow_max_kdhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, k::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false, GTEP_multi_round_info::Dict{Any,<:Any} = Dict{Any,Any}())


    rate_a = parameter(am, 0, :branch, "rate_a", k) 

    f_kdhty = variable(am, nw, decomp_group, :f_kdhty, (k, d, h, t, y))
    u_T_ky = variable(am, nw, decomp_group, :u_T_ky, (k, y))

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($k,$d,$h,$t,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        f_kdhty - rate_a * (1 + u_T_ky)
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($k,$d,$h,$t,$y)")) end    
end

function constraint_PTDF_Power_injection_ndhty_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, n::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false)

    ids_local_i = am.ref[:nw][0][:bus][n]["aggregation_info"]["new_local_gen_idx"]
        
    sum_g_idht = JuMP.AffExpr(0.0)
    for i in ids_local_i
        JuMP.add_to_expression!(sum_g_idht, 1, variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y)))
        if parameter(am, 0, :gen_index, "UNIT_CATEGORY", i) == "STORAGE"
            JuMP.add_to_expression!(sum_g_idht, -1, variable(am, nw, decomp_group, :chg_idhty, (i,d,h,t,y)))
        end
    end
    
    p_inj_ndhty = variable(am, nw, decomp_group, :p_inj_ndhty, (n,d,h,t,y))
    demand_ndhty = variable(am, nw, decomp_group, :demand_ndhty, (n,d,h,t,y))
        
    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($n,$d,$h,$t,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        p_inj_ndhty - sum_g_idht + demand_ndhty
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
    JuMP.set_name(constraint, string(const_name, "_($n,$d,$h,$t,$y)"))

end

function constraint_PTDF_LoadBalance_ndhty_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, n::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, ens_flag::Bool=true)

    inc = gtep_bus_incidence(am)   # per-bus LFL incidence, precomputed once

    Demand = apply_tnd_loss(am, get_bus_demand_with_growth(am, n, d, h, t, y; nw=0))

    ens_ndht = variable(am, nw, decomp_group, :ens_ndhty, (n,d,h,t,y))
    demand_ndhty = variable(am, nw, decomp_group, :demand_ndhty, (n,d,h,t,y))

    sum_lfl_ldht = JuMP.AffExpr(0.0)
    for lfl in get(inc.lfl_at_bus, n, _EMPTY_INT_VEC)
        JuMP.add_to_expression!(sum_lfl_ldht, 1, variable(am, nw, decomp_group, :lfl_ldhty, (lfl,d,h,t,y)))
    end
        
    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($n,$d,$h,$t,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        demand_ndhty + ens_ndht - sum_lfl_ldht - Demand
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
    JuMP.set_name(constraint, string(const_name, "_($n,$d,$h,$t,$y)"))

    if ens_flag == true
        JuMP.set_upper_bound(ens_ndht, Demand)
    end
end

function constraint_LoadBalance_ndhty_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, n::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, ens_flag::Bool=true)

    ids_local_i = am.ref[:nw][0][:bus][n]["aggregation_info"]["new_local_gen_idx"]
    inc = gtep_bus_incidence(am)   # per-bus branch/LFL incidence, precomputed once

    Demand = apply_tnd_loss(am, get_bus_demand_with_growth(am, n, d, h, t, y; nw=0))

    sum_g_idht = JuMP.AffExpr(0.0)
    for i in ids_local_i
        JuMP.add_to_expression!(sum_g_idht, 1, variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y)))
        if parameter(am, 0, :gen_index, "UNIT_CATEGORY", i) == "STORAGE"
            JuMP.add_to_expression!(sum_g_idht, -1, variable(am, nw, decomp_group, :chg_idhty, (i,d,h,t,y)))
        end
    end

    sum_lfl_ldht = JuMP.AffExpr(0.0)
    for lfl in get(inc.lfl_at_bus, n, _EMPTY_INT_VEC)
        JuMP.add_to_expression!(sum_lfl_ldht, 1, variable(am, nw, decomp_group, :lfl_ldhty, (lfl,d,h,t,y)))
    end

    ens_ndht = variable(am, nw, decomp_group, :ens_ndhty, (n,d,h,t,y))

    # Enhanced-hybrid split flow: eligible AC corridors inject f_kdhty + f_exp. The
    # variable-existence guard keeps this a no-op wherever f_exp was never built (OP model,
    # flag off), so those paths stay byte-identical.
    has_fexp = haskey(am.var[:nw][nw][decomp_group], :f_exp_kdhty)

    sum_f_ktd_to_node = JuMP.AffExpr(0.0)    # injection to the node n
    sum_f_ktd_from_node = JuMP.AffExpr(0.0)  # withdrawals from the node n
    for k in get(inc.to_branches, n, _EMPTY_INT_VEC)
        JuMP.add_to_expression!(sum_f_ktd_to_node, 1, variable(am, nw, decomp_group, :f_kdhty, (k,d,h,t,y)))
        if has_fexp && hybrid_exp_branch(am, k)
            JuMP.add_to_expression!(sum_f_ktd_to_node, 1, variable(am, nw, decomp_group, :f_exp_kdhty, (k,d,h,t,y)))
        end
    end
    for k in get(inc.from_branches, n, _EMPTY_INT_VEC)
        JuMP.add_to_expression!(sum_f_ktd_from_node, 1, variable(am, nw, decomp_group, :f_kdhty, (k,d,h,t,y)))
        if has_fexp && hybrid_exp_branch(am, k)
            JuMP.add_to_expression!(sum_f_ktd_from_node, 1, variable(am, nw, decomp_group, :f_exp_kdhty, (k,d,h,t,y)))
        end
    end

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($n,$d,$h,$t,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,
        - (sum_f_ktd_from_node - sum_f_ktd_to_node) + sum_g_idht + ens_ndht - sum_lfl_ldht - Demand
        )
    constraint = JuMP.@constraint(JuMP_model,
        expr == 0
        )
    # LMP is collected by constraint name, so name unconditionally (like the PTDF twin
    # and reserve constraints) — otherwise const_name_flag=false yields all-zero LMPs.
    JuMP.set_name(constraint, string(const_name, "_($n,$d,$h,$t,$y)"))

    if ens_flag == true
        JuMP.set_upper_bound(ens_ndht, Demand)
    end
end

function constraint_dc_power_flow_max_kdhty_OP(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, k::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false, expansion_record_flag::Bool = false, result::Dict{String,<:Any} = Dict{String,Any}())

    rate_a = parameter(am, 0, :branch, "rate_a", k) 

    total_line_expansion = 0.0
    if expansion_record_flag == true    # only when we have expansion decisions
        total_line_expansion = result["solution"]["expansion"][string("(", k, ", ", y, ")")]["u_T_ky"]
    end

    f_kdhty = variable(am, nw, decomp_group, :f_kdhty, (k, d, h, t, y))

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($k,$d,$h,$t,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        f_kdhty - rate_a * (1 + total_line_expansion)
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($k,$d,$h,$t,$y)"))  end   
end

function constraint_dc_power_flow_min_kdhty_OP(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, k::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false, expansion_record_flag::Bool = false, result::Dict{String,<:Any} = Dict{String,Any}())

    rate_a = parameter(am, 0, :branch, "rate_a", k) 

    total_line_expansion = 0.0
    if expansion_record_flag == true    # only when we have expansion decisions
        total_line_expansion = result["solution"]["expansion"][string("(", k, ", ", y, ")")]["u_T_ky"]
    end

    f_kdhty = variable(am, nw, decomp_group, :f_kdhty, (k, d, h, t, y))

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($k,$d,$h,$t,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        f_kdhty + rate_a * (1 + total_line_expansion)
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr >= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($k,$d,$h,$t,$y)")) end    
end

function constraint_dc_power_flow_max_kdhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, k::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false, GTEP_multi_round_info::Dict{Any,<:Any} = Dict{Any,Any}())


    rate_a = parameter(am, 0, :branch, "rate_a", k) 

    f_kdhty = variable(am, nw, decomp_group, :f_kdhty, (k, d, h, t, y))
    u_T_ky = variable(am, nw, decomp_group, :u_T_ky, (k, y))

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($k,$d,$h,$t,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        f_kdhty - rate_a * (1 + u_T_ky)
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($k,$d,$h,$t,$y)"))  end   

end

function constraint_dc_power_flow_min_kdhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, k::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false, GTEP_multi_round_info::Dict{Any,<:Any} = Dict{Any,Any}())


    rate_a = parameter(am, 0, :branch, "rate_a", k) 

    f_kdhty = variable(am, nw, decomp_group, :f_kdhty, (k, d, h, t, y))
    u_T_ky = variable(am, nw, decomp_group, :u_T_ky, (k, y))

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($k,$d,$h,$t,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        f_kdhty + rate_a * (1 + u_T_ky)
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr >= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($k,$d,$h,$t,$y)")) end    

end

# --- Enhanced hybrid transmission expansion (B-theta) ------------------------------
# For eligible AC corridors these three builders REPLACE constraint_dc_power_flow_max/min.
# f_kdhty is the KVL base flow (still driven by constraint_b_theta_power_flow_kdhty);
# f_exp_kdhty is the expansion increment. Both enter the nodal balance. See
# hybrid_exp_branch / tx_expansion_ub (base_functions.jl) for the gate and the box width ū.

# Base-flow thermal cap held FIXED at rate_a (the (1+u) relaxation moves to f_exp). This
# also pins the angle spread to |Δθ| ≤ rate_a/b0 = D, the premise every envelope row rests on.
function constraint_hybrid_base_cap_kdhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, k::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    rate_a = parameter(am, 0, :branch, "rate_a", k)
    f_kdhty = variable(am, nw, decomp_group, :f_kdhty, (k, d, h, t, y))

    if update
        JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_max_($k,$d,$h,$t,$y)")))
        JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_min_($k,$d,$h,$t,$y)")))
    end
    c_max = JuMP.@constraint(JuMP_model, f_kdhty - rate_a <= 0)
    c_min = JuMP.@constraint(JuMP_model, f_kdhty + rate_a >= 0)
    if const_name_flag
        JuMP.set_name(c_max, string(const_name, "_max_($k,$d,$h,$t,$y)"))
        JuMP.set_name(c_min, string(const_name, "_min_($k,$d,$h,$t,$y)"))
    end
end

# Increment capacity rows |f_exp| ≤ rate_a·u_T_ky (the pure-hybrid pair: no angle term).
function constraint_hybrid_exp_capacity_kdhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, k::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    rate_a = parameter(am, 0, :branch, "rate_a", k)
    f_exp = variable(am, nw, decomp_group, :f_exp_kdhty, (k, d, h, t, y))
    u_T_ky = variable(am, nw, decomp_group, :u_T_ky, (k, y))

    if update
        JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_max_($k,$d,$h,$t,$y)")))
        JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_min_($k,$d,$h,$t,$y)")))
    end
    c_max = JuMP.@constraint(JuMP_model, f_exp - rate_a * u_T_ky <= 0)
    c_min = JuMP.@constraint(JuMP_model, f_exp + rate_a * u_T_ky >= 0)
    if const_name_flag
        JuMP.set_name(c_max, string(const_name, "_max_($k,$d,$h,$t,$y)"))
        JuMP.set_name(c_min, string(const_name, "_min_($k,$d,$h,$t,$y)"))
    end
end

# Angle-coupling rows: ū·f_base − rate_a·(ū−u) ≤ f_exp ≤ ū·f_base + rate_a·(ū−u).
# Written via f_kdhty (≡ b0·Δθ from the base KVL equation) so no raw b0·Δθ product appears —
# coefficients stay at flow/dimensionless scale, no row scaling needed. At u = ū these pin
# f_exp = ū·f_base (exact expanded-line KVL); at u = 0 the capacity rows pin f_exp = 0.
function constraint_hybrid_exp_coupling_kdhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, k::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    rate_a = parameter(am, 0, :branch, "rate_a", k)
    ub = tx_expansion_ub(am, k)
    f_exp = variable(am, nw, decomp_group, :f_exp_kdhty, (k, d, h, t, y))
    f_kdhty = variable(am, nw, decomp_group, :f_kdhty, (k, d, h, t, y))
    u_T_ky = variable(am, nw, decomp_group, :u_T_ky, (k, y))

    if update
        JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_max_($k,$d,$h,$t,$y)")))
        JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_min_($k,$d,$h,$t,$y)")))
    end
    c_max = JuMP.@constraint(JuMP_model, f_exp - ub * f_kdhty - rate_a * ub + rate_a * u_T_ky <= 0)
    c_min = JuMP.@constraint(JuMP_model, f_exp - ub * f_kdhty + rate_a * ub - rate_a * u_T_ky >= 0)
    if const_name_flag
        JuMP.set_name(c_max, string(const_name, "_max_($k,$d,$h,$t,$y)"))
        JuMP.set_name(c_min, string(const_name, "_min_($k,$d,$h,$t,$y)"))
    end
end

function constraint_dc_power_flow_kdhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, k::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    ptdf_k = parameter(am, 0, :branch, k, "ptdf")

    f_ptdf_kdhty = JuMP.AffExpr(0.0)
    for n in get_index(am, :bus, 0)

        if abs(ptdf_k[n]) >= am.setting["Simulation Setting"]["PTDF_threshold_value"]
        
            p_inj_ndhty = variable(am, nw, decomp_group, :p_inj_ndhty, (n,d,h,t,y))
                
            expr = JuMP.@expression(JuMP_model,  
                ptdf_k[n]*(p_inj_ndhty)
            )  

            JuMP.add_to_expression!(f_ptdf_kdhty, expr)

        end


    end

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($k,$d,$h,$t,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        f_ptdf_kdhty - variable(am, nw, decomp_group, :f_kdhty, (k,d,h,t,y))
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($k,$d,$h,$t,$y)")) end   
end

function constraint_PRM_y_p_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, y::Int, n::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=true)

    ids_local_i = am.ref[:nw][0][:zone]["planning_reserve"][string(n)]["aggregation_info"]["new_local_gen_idx"]

    zone_bus_ids = am.ref[:nw][0][:zone]["planning_reserve"][string(n)]["aggregation_info"]["zone_bus_idx"]
    peak_demand = get_zone_coincident_peak_demand_with_growth(am, zone_bus_ids, y; nw=0)

    PRM = am.setting["Simulation Configuration"]["planning_reserve_margin_value"]

    totalucap = JuMP.AffExpr(0.0)
    for i in ids_local_i
        bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
        tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)    
        CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")

        unit_group = parameter(am, 0, :gen_index, "UNIT_GROUP", i)
        CAPCRED = parameter(am, bus_idx, :gen_bus, tech_idx, "CAPCRED")
        if CAPCRED isa String
            data_identifier = am.ref[:nw][0][:zone]["capacity_credit"]["data_identifier"]
            data_region = am.ref[:nw][0][:bus][bus_idx]["region_mapping_info"][data_identifier]
            data_region = data_region isa AbstractString ? data_region : string(data_region)
            CAPCRED = am.ref[:nw][0][:zone]["capacity_credit"][data_region][CAPCRED]
        end

        JuMP.add_to_expression!(totalucap, CAP * CAPCRED, variable(am, nw, decomp_group, :u_G_iy, (i, y)))
    end

    if am.setting["Simulation Configuration"]["planning_reserve_margin_type"] == "maximum"
    
        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($n,$y)"))) end    
        expr = JuMP.@expression(JuMP_model,  
            totalucap - peak_demand * (1+PRM) 
            )
        constraint = JuMP.@constraint(JuMP_model, 
            expr <= 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($n,$y)")) end   

    elseif am.setting["Simulation Configuration"]["planning_reserve_margin_type"] == "minimum"

        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($n,$y)"))) end    
        expr = JuMP.@expression(JuMP_model,  
            totalucap - peak_demand * (1+PRM) 
            )
        constraint = JuMP.@constraint(JuMP_model, 
            expr >= 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($n,$y)")) end   

    end
    
end

function constraint_ESH_investment_min_iy_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)    
    lower_limit = parameter(am, bus_idx, :gen_bus, tech_idx, "STOHR_MIN")
    
    u_new_ESH_iy = variable(am, nw, decomp_group, :u_new_ESH_iy, (i, y))
    u_new_G_iy = variable(am, nw, decomp_group, :u_new_G_iy, (i, y))
    u_ret_G_iy = variable(am, nw, decomp_group, :u_ret_G_iy, (i, y))

    if parameter(am, bus_idx, :gen_bus, tech_idx, "ES_STO_INVEST_FLAG") == true
        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$y)"))) end    
        expr = JuMP.@expression(JuMP_model,  
        u_new_ESH_iy - lower_limit * u_new_G_iy
            )
        constraint = JuMP.@constraint(JuMP_model, 
            expr >= 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$y)")) end    
    else
        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$y)"))) end    
        expr = JuMP.@expression(JuMP_model,  
        u_new_ESH_iy - lower_limit * u_new_G_iy
            )
        constraint = JuMP.@constraint(JuMP_model, 
            expr == 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$y)")) end    
    end
end

function constraint_u_ESE_balance_iy_OP_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)    
    Existing_ES_MWh = parameter(am, bus_idx, :gen_bus, tech_idx, "ES_MWh")
    
    u_ESE_iy = variable(am, nw, decomp_group, :u_ESE_iy, (i, y))

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        u_ESE_iy - Existing_ES_MWh 
        )
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$y)")) end    
end

function constraint_u_ESE_balance_iy_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false, GTEP_multi_round_info::Dict{Any,<:Any} = Dict{Any,Any}())

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)    
    Existing_ES_MWh = parameter(am, bus_idx, :gen_bus, tech_idx, "ES_MWh")
    Charge_CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "Charge_CAP")
    life_time = parameter(am, bus_idx, :gen_bus, tech_idx, "Life")

    recorded_decision = 0.0
    first_active_year = 1
    round_year_list = [i for i in 1:y]
    if am.setting["Planning Design"]["multi_round_solution_process_flag"] == true
        recorded_decision = GTEP_multi_round_info["recorded_investment_decisions"]["U_NEW_ESH_i"][string(i)]
        first_active_year = GTEP_multi_round_info["round_first_y"]
        round_year_list = [i for i in GTEP_multi_round_info["round_first_y"]:y]
    end

    # check max service life of new asset
    stage_lengh = am.ref[:nw][0][:planning_stages][y]["stage_length"]
    max_service_stages_of_new_asset = Int(ceil(life_time / stage_lengh)) + 1
    new_asset_for_retirement_year = y - max_service_stages_of_new_asset + 1 # planning stages
    
    u_ESE_iy = variable(am, nw, decomp_group, :u_ESE_iy, (i, y))
    u_new_ESH_iy = variable(am, nw, decomp_group, :u_new_ESH_iy, (i, y))

    retire_u_new_ESH_iy = 0.0
    if y >= max_service_stages_of_new_asset   # now we need to consider retirement of new assets after lifetime

        # new investment subject to retirement after lifetime
        if new_asset_for_retirement_year in round_year_list
            retire_u_new_ESH_iy = variable(am, nw, decomp_group, :u_new_ESH_iy, (i, new_asset_for_retirement_year))
        else    # in this case, we will use investment decision from previous simulation rounds
            retire_u_new_ESH_iy = GTEP_multi_round_info["recorded_investment_decisions"]["u_new_ESH_iy"][string("(", i, ", ", new_asset_for_retirement_year, ")")]
        end
    end

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$y)"))) end
    expr = JuMP.AffExpr(0.0)
    if y == first_active_year
        expr = JuMP.@expression(JuMP_model,
            u_ESE_iy - Existing_ES_MWh - Charge_CAP * recorded_decision - Charge_CAP * u_new_ESH_iy + Charge_CAP * retire_u_new_ESH_iy
        )
    else
        u_ESE_i_prev_y = variable(am, nw, decomp_group, :u_ESE_iy, (i, y - 1))
        expr = JuMP.@expression(JuMP_model,
            u_ESE_iy - u_ESE_i_prev_y - Charge_CAP * u_new_ESH_iy + Charge_CAP * retire_u_new_ESH_iy
        )
    end
    constraint = JuMP.@constraint(JuMP_model,
        expr == 0
    )
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$y)")) end
end

function constraint_annual_raw_mmaterial_limits_my_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, m::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    pu_base = am.setting["Simulation Setting"]["per_unit_base_value"]

    ids_i = [(i) for (i) in get_index(am, :gen_index, 0)]

    material_type = am.ref[:nw][0][:raw_materials][m]["Material Type"]
    material_limit = am.ref[:nw][0][:raw_materials][m]["Annual Limit"]

    total_consumption = JuMP.AffExpr(0.0)
    for i in ids_i
        bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
        tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)    
        
        if parameter(am, bus_idx, :gen_bus, tech_idx, "Material_Flag") == true
            CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP") * pu_base
            material_consumption = parameter(am, bus_idx, :gen_bus, tech_idx, "GenTech_Raw_Materials")[material_type]
            
            JuMP.add_to_expression!(total_consumption, CAP * material_consumption, variable(am, nw, decomp_group, :u_new_G_iy, (i, y)))
        end
    end

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($m,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        total_consumption - material_limit
        )
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($m,$y)")) end    
end

function constraint_u_G_balance_iy_OP_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)    
    EXUNITS = parameter(am, bus_idx, :gen_bus, tech_idx, "EXUNITS")
    
    u_G_iy = variable(am, nw, decomp_group, :u_G_iy, (i, y))

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        u_G_iy - EXUNITS 
        )
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$y)")) end    
end

function constraint_resource_limit_per_tech(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, y::Int, resource_zone_id::Int, resource_idx::Int, resource_id::String; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false, GTEP_multi_round_info::Dict{Any,<:Any} = Dict{Any,Any}())

    if resource_id == "psh" # psh resource supply limits = new investment only

        EX_INVEST = 0.0
        resource_cap = am.ref[:nw][0][:zone]["resource_supply_curve"][string(resource_zone_id)]["data"][resource_id]["Resource_Limit_MW"]
        resource_id_list = [resource_id, string(resource_id, "_new")]

        total_new_capacity = JuMP.AffExpr(0.0)
        for local_i in am.ref[:nw][0][:zone]["resource_supply_curve"][string(resource_zone_id)]["aggregation_info"]["new_local_gen_idx"]
            bus_idx = parameter(am, 0, :gen_index, "bus_idx", local_i)
            tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", local_i)    
            
            if parameter(am, bus_idx, :gen_bus, tech_idx, "Resource_Limit_Flag") == true
                if parameter(am, bus_idx, :gen_bus, tech_idx, "Resource_Limit_ID") in resource_id_list

                    CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")
                    
                    year_list = [i for i in 1:y]
                    if am.setting["Planning Design"]["multi_round_solution_process_flag"] == true
                        EX_INVEST += CAP * GTEP_multi_round_info["recorded_investment_decisions"]["U_NEW_G_i"][string(local_i)] 
                        year_list = [i for i in GTEP_multi_round_info["round_first_y"]:y]
                    end

                    JuMP.add_to_expression!(total_new_capacity, CAP, variable(am, nw, decomp_group, :u_new_G_iy, (local_i, y)))
                end
            end
        end

        if EX_INVEST > resource_cap
            resource_cap = EX_INVEST
        end

        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "$resource_id_($resource_zone_id,$resource_idx,$y)"))) end    
        expr = JuMP.@expression(JuMP_model,  
            (total_new_capacity + EX_INVEST) - resource_cap
            )
        constraint = JuMP.@constraint(JuMP_model, 
            expr <= 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($resource_zone_id,$resource_idx,$y)")) end    

    else

        EXCAP = 0.0
        resource_cap = am.ref[:nw][0][:zone]["resource_supply_curve"][string(resource_zone_id)]["data"][resource_id]["Resource_Limit_MW"]
        resource_id_list = [resource_id, string(resource_id, "_new")]

        total_capacity = JuMP.AffExpr(0.0)
        for local_i in am.ref[:nw][0][:zone]["resource_supply_curve"][string(resource_zone_id)]["aggregation_info"]["new_local_gen_idx"]
            bus_idx = parameter(am, 0, :gen_index, "bus_idx", local_i)
            tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", local_i)    
            
            if parameter(am, bus_idx, :gen_bus, tech_idx, "Resource_Limit_Flag") == true
                if parameter(am, bus_idx, :gen_bus, tech_idx, "Resource_Limit_ID") in resource_id_list
                    
                    EXUNITS = parameter(am, bus_idx, :gen_bus, tech_idx, "EXUNITS")
                    CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")
                    EXCAP += EXUNITS * CAP
                    
                    year_list = [i for i in 1:y]
                    if am.setting["Planning Design"]["multi_round_solution_process_flag"] == true
                        EXCAP += CAP * GTEP_multi_round_info["recorded_investment_decisions"]["U_NEW_G_i"][string(local_i)] 
                        EXCAP -= CAP * GTEP_multi_round_info["recorded_investment_decisions"]["U_RET_G_i"][string(local_i)]
                        year_list = [i for i in GTEP_multi_round_info["round_first_y"]:y]
                    end

                    JuMP.add_to_expression!(total_capacity, CAP, variable(am, nw, decomp_group, :u_G_iy, (local_i, y)))
                end
            end
        end

        if EXCAP > resource_cap
            resource_cap = EXCAP
        end

        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "$resource_id_($resource_zone_id,$resource_idx,$y)"))) end    
        expr = JuMP.@expression(JuMP_model,  
            total_capacity - resource_cap
            )
        constraint = JuMP.@constraint(JuMP_model, 
            expr <= 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($resource_zone_id,$resource_idx,$y)")) end    

    end
      
end

function constraint_u_newT_ky_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, k::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false, fix_flag::Bool=false, GTEP_multi_round_info::Dict{Any,<:Any} = Dict{Any,Any}())

    rate_a = parameter(am, 0, :branch, "rate_a", k) 
    max_rate_a = parameter(am, 0, :branch, "max_rate_a", k) 
    expansion_flag = parameter(am, 0, :branch, "expansion_flag", k) 

    if rate_a * (1 + am.setting["Planning Design"]["transmission_expansion_limit_value"]) >= max_rate_a
        max_rate_a = rate_a * (1 + am.setting["Planning Design"]["transmission_expansion_limit_value"])
    end

    recorded_decision = 0.0
    first_active_year = 1
    if am.setting["Planning Design"]["multi_round_solution_process_flag"] == true
        recorded_decision = GTEP_multi_round_info["recorded_investment_decisions"]["U_NEW_T_k"][string(k)]
        first_active_year = GTEP_multi_round_info["round_first_y"]
    end

    if expansion_flag == true
             
        u_T_ky = variable(am, nw, decomp_group, :u_T_ky, (k, y))
        u_new_T_ky = variable(am, nw, decomp_group, :u_new_T_ky, (k, y))
        if rate_a > 0.0
            JuMP.set_upper_bound(u_T_ky, max_rate_a / rate_a - 1.0)
        end

        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($k,$y)"))) end    
        expr = JuMP.AffExpr(0.0)
        if y == first_active_year
            expr = JuMP.@expression(JuMP_model,
                u_T_ky - recorded_decision - u_new_T_ky
            )
        else
            u_T_k_prev_y = variable(am, nw, decomp_group, :u_T_ky, (k, y - 1))
            expr = JuMP.@expression(JuMP_model,
                u_T_ky - u_T_k_prev_y - u_new_T_ky
            )
        end
        constraint = JuMP.@constraint(JuMP_model, 
            expr == 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($k,$y)")) end    

    else

        JuMP.fix(variable(am, nw, decomp_group, :u_new_T_ky, (k, y)), 0.0; force = true)
        JuMP.fix(variable(am, nw, decomp_group, :u_T_ky, (k, y)), 0.0; force = true)

    end
    

end

function constraint_investment_external_iy_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, y::Int, investment_MW; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)
    
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)    
    CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP") 

    u_new_G_iy = variable(am, nw, decomp_group, :u_new_G_iy, (i, y))

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$y)"))) end    
    expr = JuMP.@expression(JuMP_model, 
        CAP * u_new_G_iy - (investment_MW)
        )
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$y)")) end   
    
end

function constraint_ret_G_external_iy_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, y::Int, econ_ret_MW; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)
    
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)    
    CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP") 

    year = am.ref[:nw][0][:planning_stages][y]["year"]
    tag = string("Ret_", year)
    planned_retirement = am.ref[:nw][bus_idx][:gen_bus][tech_idx]["Planned_Retirement"][tag] 

    u_ret_G_iy = variable(am, nw, decomp_group, :u_ret_G_iy, (i, y))

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$y)"))) end    
    expr = JuMP.@expression(JuMP_model, 
        CAP * u_ret_G_iy - (planned_retirement + econ_ret_MW)
        )
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$y)")) end    
    
end

function constraint_u_G_balance_iy_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false, GTEP_multi_round_info::Dict{Any,<:Any} = Dict{Any,Any}())

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)    
    EXUNITS = parameter(am, bus_idx, :gen_bus, tech_idx, "EXUNITS")

    recorded_decision = 0.0
    first_active_year = 1
    if am.setting["Planning Design"]["multi_round_solution_process_flag"] == true
        recorded_decision = GTEP_multi_round_info["recorded_investment_decisions"]["U_NEW_G_i"][string(i)] - GTEP_multi_round_info["recorded_investment_decisions"]["U_RET_G_i"][string(i)]
        first_active_year = GTEP_multi_round_info["round_first_y"]
    end
    
    u_G_iy = variable(am, nw, decomp_group, :u_G_iy, (i, y))
    u_new_G_iy = variable(am, nw, decomp_group, :u_new_G_iy, (i, y))
    u_ret_G_iy = variable(am, nw, decomp_group, :u_ret_G_iy, (i, y))

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$y)"))) end    
    expr = JuMP.AffExpr(0.0)
    if y == first_active_year
        expr = JuMP.@expression(JuMP_model,
            u_G_iy - EXUNITS - recorded_decision - u_new_G_iy + u_ret_G_iy
        )
    else
        u_G_i_prev_y = variable(am, nw, decomp_group, :u_G_iy, (i, y - 1))
        expr = JuMP.@expression(JuMP_model,
            u_G_iy - u_G_i_prev_y - u_new_G_iy + u_ret_G_iy
        )
    end
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$y)")) end   
end

function constraint_system_total_investment_tech_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, tech::Int, ids_i, ids_y; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false, GTEP_multi_round_info::Dict{Any,<:Any} = Dict{Any,Any}())


    tech_UNITGROUP = parameter(am, 0, :gen_technology, "UNITGROUP", tech)
    max_investment = parameter(am, 0, :gen_technology, "SYSTEM_MAXINVEST", tech)
    min_investment = parameter(am, 0, :gen_technology, "SYSTEM_MININVEST", tech)

    recorded_decision = 0.0
    total_u_new_G_sum = JuMP.AffExpr(0.0)   

    for i in ids_i
        
        if parameter(am, 0, :gen_index, "UNIT_GROUP", i) == tech_UNITGROUP
            
            if am.setting["Planning Design"]["multi_round_solution_process_flag"] == true
                recorded_decision = GTEP_multi_round_info["recorded_investment_decisions"]["U_NEW_G_i"][string(i)]
            end

            for y in ids_y
                JuMP.add_to_expression!(total_u_new_G_sum, 1, variable(am, nw, decomp_group, :u_new_G_iy, (i, y)))  
            end
        end
    end

    if min_investment > 0
            
        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_Min_($tech)"))) end    
        expr = JuMP.@expression(JuMP_model, 
            total_u_new_G_sum + recorded_decision - min_investment
            )
        constraint = JuMP.@constraint(JuMP_model, 
            expr >= 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($tech)")) end    

    end

    if max_investment > 0

        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_Max_($tech)"))) end    
        expr = JuMP.@expression(JuMP_model, 
            total_u_new_G_sum + recorded_decision - max_investment
            )
        constraint = JuMP.@constraint(JuMP_model, 
            expr <= 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($tech)")) end    

    end

end

function constraint_planned_ret_G_balance_iy_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false, GTEP_multi_round_info::Dict{Any,<:Any} = Dict{Any,Any}())
    
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)    
    EXUNITS = parameter(am, bus_idx, :gen_bus, tech_idx, "EXUNITS")
    CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP") 
    ret_flag = parameter(am, bus_idx, :gen_bus, tech_idx, "RET_FLAG")
    life_time = parameter(am, bus_idx, :gen_bus, tech_idx, "Life")

    recorded_decision = 0.0
    year_list = [i for i in 1:y]
    round_year_list = [i for i in 1:y]
    if am.setting["Planning Design"]["multi_round_solution_process_flag"] == true
        recorded_decision = CAP * GTEP_multi_round_info["recorded_investment_decisions"]["U_RET_G_i"][string(i)]    # accumulated retirement decisions
        round_year_list = [i for i in GTEP_multi_round_info["round_first_y"]:y]
    end

    # check max service life of new asset
    stage_lengh = am.ref[:nw][0][:planning_stages][y]["stage_length"]
    max_service_stages_of_new_asset = Int(ceil(life_time / stage_lengh)) + 1
    new_asset_for_retirement_year = y - max_service_stages_of_new_asset + 1 # planning stages

    if ret_flag == true # economic retirement is considered

        if y >= max_service_stages_of_new_asset   # now we need to consider retirement of new assets after lifetime

            planned_retirement = 0.0
            
            for q in year_list
                year = am.ref[:nw][0][:planning_stages][q]["year"]
                tag = string("Ret_", year)
                planned_retirement += am.ref[:nw][bus_idx][:gen_bus][tech_idx]["Planned_Retirement"][tag] 
            end
    
            u_ret_sum = JuMP.AffExpr(0.0)
            for q in round_year_list
                JuMP.add_to_expression!(u_ret_sum, CAP, variable(am, nw, decomp_group, :u_ret_G_iy, (i, q)))      
            end

            # new investment subject to retirement after lifetime 
            u_new_ret_sum = JuMP.AffExpr(0.0)
            prior_u_new_ret_sum = 0.0
            for q in [i for i in max_service_stages_of_new_asset:y]
                new_asset_retirement_year = q - max_service_stages_of_new_asset + 1

                if new_asset_retirement_year in round_year_list
                    JuMP.add_to_expression!(u_new_ret_sum, CAP, variable(am, nw, decomp_group, :u_new_G_iy, (i, new_asset_retirement_year)))                       
                else    # in this case, we will use investment decision from previous simulation rounds
                    prior_u_new_ret_sum += GTEP_multi_round_info["recorded_investment_decisions"]["u_new_G_iy"][string("(", i, ", ", new_asset_retirement_year, ")")]
                end
            end
    
            if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$y)"))) end    
            expr = JuMP.@expression(JuMP_model, 
                u_ret_sum - (planned_retirement + (u_new_ret_sum + prior_u_new_ret_sum) - recorded_decision)
                )
            constraint = JuMP.@constraint(JuMP_model, 
                expr >= 0
                )    
            if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$y)")) end    

        else # no consideration of new asset retirements

            planned_retirement = 0.0
            
            for q in year_list
                year = am.ref[:nw][0][:planning_stages][q]["year"]
                tag = string("Ret_", year)
                planned_retirement += am.ref[:nw][bus_idx][:gen_bus][tech_idx]["Planned_Retirement"][tag] 
            end
    
            u_ret_sum = JuMP.AffExpr(0.0)
            for q in round_year_list
                JuMP.add_to_expression!(u_ret_sum, CAP, variable(am, nw, decomp_group, :u_ret_G_iy, (i, q)))      
            end
    
            if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$y)"))) end    
            expr = JuMP.@expression(JuMP_model, 
                u_ret_sum - (planned_retirement - recorded_decision)
                )
            constraint = JuMP.@constraint(JuMP_model, 
                expr >= 0
                )    
            if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$y)")) end    

        end


    else    # economic retirement is NOT considered

        year = am.ref[:nw][0][:planning_stages][y]["year"]
        planned_retirement = am.ref[:nw][bus_idx][:gen_bus][tech_idx]["Planned_Retirement"][string("Ret_", year)] 

        u_ret_G_iy = variable(am, nw, decomp_group, :u_ret_G_iy, (i, y))

        if y >= max_service_stages_of_new_asset   # now we need to consider retirement of new assets after lifetime

            # new investment subject to retirement after lifetime
            retire_new_G_iy = 0.0
            if new_asset_for_retirement_year in round_year_list
                retire_new_G_iy = variable(am, nw, decomp_group, :u_new_G_iy, (i, new_asset_for_retirement_year))
            else    # in this case, we will use investment decision from previous simulation rounds
                retire_new_G_iy = GTEP_multi_round_info["recorded_investment_decisions"]["u_new_G_iy"][string("(", i, ", ", new_asset_for_retirement_year, ")")]
            end
    
            if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$y)"))) end    
            expr = JuMP.@expression(JuMP_model, 
                CAP * u_ret_G_iy - (planned_retirement + CAP * retire_new_G_iy)
                )
            constraint = JuMP.@constraint(JuMP_model, 
                expr == 0
                )    
            if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$y)")) end    


        else # no consideration of new asset retirements

            if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$y)"))) end    
            expr = JuMP.@expression(JuMP_model, 
                CAP * u_ret_G_iy - planned_retirement
                )
            constraint = JuMP.@constraint(JuMP_model, 
                expr == 0
                )    
            if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$y)")) end    

        end

    end
    
end

function constraint_inertia_ndhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, n::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    ids_local_i = am.ref[:nw][0][:zone]["policy"][string(n)]["aggregation_info"]["new_local_gen_idx"]

    inertia_requirement = am.ref[:nw][0][:zone]["policy"][n]["Inertia"]["target"]

    if inertia_requirement > 0

        sum_inertia = JuMP.AffExpr(0.0)
        for i in ids_local_i
            bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
            tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)
            CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")
            power_factor = parameter(am, bus_idx, :gen_bus, tech_idx, "Power_Factor")
            inertia_constant = parameter(am, bus_idx, :gen_bus, tech_idx, "Inertia_Constant")

            JuMP.add_to_expression!(sum_inertia, CAP / power_factor * inertia_constant, variable(am, nw, decomp_group, :c_idhty, (i,d,h,t,y)))
        end

        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($n,$d,$h,$t,$y)"))) end # specify indices
        expr = JuMP.@expression(JuMP_model,
            sum_inertia - inertia_requirement
            )
        constraint = JuMP.@constraint(JuMP_model,
            expr >= 0
            )
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($n,$d,$h,$t,$y)")) end # specify indices

    end
end

function constraint_sum_p_injdhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, constant_load=0.0, system_peak_scale=1.0)

    total_p_inj_dhty = JuMP.AffExpr(0.0)

    for n in get_index(am, :bus, 0)
        JuMP.add_to_expression!(total_p_inj_dhty, 1,  variable(am, nw, decomp_group, :p_inj_ndhty, (n,d,h,t,y)))
    end

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($d,$h,$t,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        total_p_inj_dhty 
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
    JuMP.set_name(constraint, string(const_name, "_($d,$h,$t,$y)"))    
end

function constraint_Start_Up_Limit_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)


    c_idht = variable(am, nw, decomp_group, :c_idhty, (i,d,h,t,y))
    su_idht = variable(am, nw, decomp_group, :su_idhty, (i,d,h,t,y))

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "su_($i,$d,$h,$t,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        su_idht - c_idht 
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "su_($i,$d,$h,$t,$y)")) end    

end

function constraint_Start_Up_Status_Dn_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    day_group = parse(Int, am.ref[:nw][0][:repdays][d]["Day_Group_ID"])
    HFREQ = am.setting["Simulation Configuration"]["HFREQ"]

    c_idht = variable(am, nw, decomp_group, :c_idhty, (i,d,h,t,y))
    su_idht = variable(am, nw, decomp_group, :su_idhty, (i,d,h,t,y))

    if (h == 1) 
        if am.ref[:nw][0][:repdays][d]["Day"] != am.ref[:nw][0][:repday_groups][day_group]["Start_Day_Id"]
            last_hour = last(am.setting["run_H"])
            c_idht_prior = variable(am, nw, decomp_group, :c_idhty, (i,d-1,last_hour,t,y))
            
            if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
            expr = JuMP.@expression(JuMP_model,  
                su_idht - c_idht + c_idht_prior
                )   
            constraint = JuMP.@constraint(JuMP_model, 
                expr >= 0
                )    
            if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end    

        elseif am.ref[:nw][0][:repdays][d]["Day"] == am.ref[:nw][0][:repday_groups][day_group]["Start_Day_Id"]

            num_days_in_group = length(am.ref[:nw][0][:repday_groups][day_group]["Day_Idx_List"])
            end_day_id = am.ref[:nw][0][:repday_groups][day_group]["Day_Idx_List"][num_days_in_group]
            last_hour = last(am.setting["run_H"])
            
            c_idht_prior = variable(am, nw, decomp_group, :c_idhty, (i,end_day_id,last_hour,t,y))
            
            if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
            expr = JuMP.@expression(JuMP_model,  
                su_idht - c_idht + c_idht_prior
                )   
            constraint = JuMP.@constraint(JuMP_model, 
                expr >= 0
                )    
            if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end   

        end

    else

        c_idht_prior = variable(am, nw, decomp_group, :c_idhty, (i,d,h-HFREQ,t,y))
    
        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
        expr = JuMP.@expression(JuMP_model,  
            su_idht - c_idht + c_idht_prior
            )   
        constraint = JuMP.@constraint(JuMP_model, 
            expr >= 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end    
    
    end

end

function constraint_Start_Up_Status_Up_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    day_group = parse(Int, am.ref[:nw][0][:repdays][d]["Day_Group_ID"])
    HFREQ = am.setting["Simulation Configuration"]["HFREQ"]

    u_G_iy = variable(am, nw, decomp_group, :u_G_iy, (i,y))
    su_idht = variable(am, nw, decomp_group, :su_idhty, (i,d,h,t,y))

    if (h == 1) 
        if am.ref[:nw][0][:repdays][d]["Day"] != am.ref[:nw][0][:repday_groups][day_group]["Start_Day_Id"]
            last_hour = last(am.setting["run_H"])
            c_idht_prior = variable(am, nw, decomp_group, :c_idhty, (i,d-1,last_hour,t,y))
            
            if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
            expr = JuMP.@expression(JuMP_model,  
                su_idht - u_G_iy + c_idht_prior
                )   
            constraint = JuMP.@constraint(JuMP_model, 
                expr <= 0
                )    
            if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end    

        elseif am.ref[:nw][0][:repdays][d]["Day"] == am.ref[:nw][0][:repday_groups][day_group]["Start_Day_Id"]

            num_days_in_group = length(am.ref[:nw][0][:repday_groups][day_group]["Day_Idx_List"])
            end_day_id = am.ref[:nw][0][:repday_groups][day_group]["Day_Idx_List"][num_days_in_group]
            last_hour = last(am.setting["run_H"])
            
            c_idht_prior = variable(am, nw, decomp_group, :c_idhty, (i,end_day_id,last_hour,t,y))
            
            if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
            expr = JuMP.@expression(JuMP_model,  
                su_idht - u_G_iy + c_idht_prior
                )   
            constraint = JuMP.@constraint(JuMP_model, 
                expr <= 0
                )    
            if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end    

        end

    else

        c_idht_prior = variable(am, nw, decomp_group, :c_idhty, (i,d,h-HFREQ,t,y))
    
        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
        expr = JuMP.@expression(JuMP_model,  
            su_idht - u_G_iy + c_idht_prior
            )   
        constraint = JuMP.@constraint(JuMP_model, 
            expr <= 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end    
    
    end

end

function constraint_RegDnMax_UC_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)   
    PMAX = parameter(am, bus_idx, :gen_bus, tech_idx, "PMAX")
    MAXR = parameter(am, bus_idx, :gen_bus, tech_idx, "MAXR")   # 5 min ramp limit in %
    CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")

    reg_dn_idht = variable(am, nw, decomp_group, :reg_dn_idhty, (i,d,h,t,y))
    c_idht = variable(am, nw, decomp_group, :c_idhty, (i,d,h,t,y)) 

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        reg_dn_idht- c_idht * PMAX * CAP * MAXR
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end    

end

function constraint_RegUpMax_UC_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)   
    PMAX = parameter(am, bus_idx, :gen_bus, tech_idx, "PMAX")
    MAXR = parameter(am, bus_idx, :gen_bus, tech_idx, "MAXR")   # 5 min ramp limit in %
    CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")

    reg_up_idht = variable(am, nw, decomp_group, :reg_up_idhty, (i,d,h,t,y))
    c_idht = variable(am, nw, decomp_group, :c_idhty, (i,d,h,t,y)) 

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        reg_up_idht - c_idht * PMAX * CAP * MAXR
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end    

end

function constraint_NonSpinMax_UC_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)   
    PMAX = parameter(am, bus_idx, :gen_bus, tech_idx, "PMAX")
    RUL = parameter(am, bus_idx, :gen_bus, tech_idx, "RUL")   # 10-mins ramp rate in %
    CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")

    nonspin_idhty = variable(am, nw, decomp_group, :nonspin_idhty, (i,d,h,t,y))
    c_idht = variable(am, nw, decomp_group, :c_idhty, (i,d,h,t,y)) 
    u_G_iy = variable(am, nw, decomp_group, :u_G_iy, (i,y))

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        nonspin_idhty - PMAX * CAP * RUL * (u_G_iy - c_idht)
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end    

end

function constraint_RampUp_InterTemporal_Hour_UC_idhy(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false, apply_PMAX::Bool=true)

    day_group = parse(Int, am.ref[:nw][0][:repdays][d]["Day_Group_ID"])
    HFREQ = am.setting["Simulation Configuration"]["HFREQ"]

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
    CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")
    PMAX = parameter(am, bus_idx, :gen_bus, tech_idx, "PMAX")   # %
    if apply_PMAX == false
        PMAX = 1.0
    end
    Ramp_Rate = min(1.0, 60 * parameter(am, bus_idx, :gen_bus, tech_idx, "Ramp"))   # hourly Ramp Rate (%)
    SU_Ramp_Rate = min(1.0, 60 * parameter(am, bus_idx, :gen_bus, tech_idx, "Ramp"))   # start-up hourly Ramp Rate (%)

    su_idht = variable(am, nw, decomp_group, :su_idhty, (i,d,h,t,y))

    if (h == 1)
        if am.ref[:nw][0][:repdays][d]["Day"] != am.ref[:nw][0][:repday_groups][day_group]["Start_Day_Id"]
            last_hour = last(am.setting["run_H"])
            
            g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
            prior_g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d-1,last_hour,t,y))
            prior_c_idht = variable(am, nw, decomp_group, :c_idhty, (i,d-1,last_hour,t,y))
            
            if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$y)"))) end    
            expr = JuMP.@expression(JuMP_model,  
                g_idht - prior_g_idht - (prior_c_idht * Ramp_Rate * CAP * PMAX) - (su_idht * SU_Ramp_Rate * CAP * PMAX)
                )   
            constraint = JuMP.@constraint(JuMP_model, 
                expr <= 0
                )    
            if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$y)")) end    
        
        elseif am.ref[:nw][0][:repdays][d]["Day"] == am.ref[:nw][0][:repday_groups][day_group]["Start_Day_Id"]
            num_days_in_group = length(am.ref[:nw][0][:repday_groups][day_group]["Day_Idx_List"])
            end_day_id = am.ref[:nw][0][:repday_groups][day_group]["Day_Idx_List"][num_days_in_group]
            last_hour = last(am.setting["run_H"])
            
            g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
            prior_g_idht = variable(am, nw, decomp_group, :g_idhty, (i,end_day_id,last_hour,t,y))
            prior_c_idht = variable(am, nw, decomp_group, :c_idhty, (i,end_day_id,last_hour,t,y))
            
            if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$y)"))) end    
            expr = JuMP.@expression(JuMP_model,  
                g_idht - prior_g_idht - (prior_c_idht * Ramp_Rate * CAP * PMAX) - (su_idht * SU_Ramp_Rate * CAP * PMAX)
                )   
            constraint = JuMP.@constraint(JuMP_model, 
                expr <= 0
                )    
            if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$y)")) end   
        
        end
    elseif (h != 1) 
        HFREQ = am.setting["Simulation Configuration"]["HFREQ"]
        
        g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
        prior_g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h-HFREQ,t,y))
        prior_c_idht = variable(am, nw, decomp_group, :c_idhty, (i,d,h-HFREQ,t,y))
        
        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$y)"))) end    
        expr = JuMP.@expression(JuMP_model,  
            g_idht - prior_g_idht - (prior_c_idht * Ramp_Rate * CAP * PMAX) - (su_idht * SU_Ramp_Rate * CAP * PMAX)
            )   
        constraint = JuMP.@constraint(JuMP_model, 
            expr <= 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$y)")) end    

    end
end

function constraint_RampDn_InterTemporal_Hour_UC_idhy(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false, apply_PMAX::Bool=true)

    day_group = parse(Int, am.ref[:nw][0][:repdays][d]["Day_Group_ID"])
    HFREQ = am.setting["Simulation Configuration"]["HFREQ"]

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
    CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")
    PMAX = parameter(am, bus_idx, :gen_bus, tech_idx, "PMAX")   # %
    if apply_PMAX == false
        PMAX = 1.0
    end
    Ramp_Rate = min(1.0, 60 * parameter(am, bus_idx, :gen_bus, tech_idx, "Ramp"))   # hourly Ramp Rate (%)
    SD_Ramp_Rate = min(1.0, 60 * parameter(am, bus_idx, :gen_bus, tech_idx, "Ramp"))   # shut-down hourly Ramp Rate (%)

    su_idht = variable(am, nw, decomp_group, :su_idhty, (i,d,h,t,y))
    c_idhty = variable(am, nw, decomp_group, :c_idhty, (i,d,h,t,y))

    if (h == 1)
        if am.ref[:nw][0][:repdays][d]["Day"] != am.ref[:nw][0][:repday_groups][day_group]["Start_Day_Id"]
            last_hour = last(am.setting["run_H"])
            
            g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
            prior_g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d-1,last_hour,t,y))
            prior_c_idht = variable(am, nw, decomp_group, :c_idhty, (i,d-1,last_hour,t,y))
            
            if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$y)"))) end    
            expr = JuMP.@expression(JuMP_model,  
                - g_idht + prior_g_idht - (c_idhty * Ramp_Rate * CAP * PMAX) - SD_Ramp_Rate * CAP * PMAX * (su_idht - c_idhty + prior_c_idht)
                )   
            constraint = JuMP.@constraint(JuMP_model, 
                expr <= 0
                )    
            if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$y)")) end    
        
        elseif am.ref[:nw][0][:repdays][d]["Day"] == am.ref[:nw][0][:repday_groups][day_group]["Start_Day_Id"]
            num_days_in_group = length(am.ref[:nw][0][:repday_groups][day_group]["Day_Idx_List"])
            end_day_id = am.ref[:nw][0][:repday_groups][day_group]["Day_Idx_List"][num_days_in_group]
            last_hour = last(am.setting["run_H"])
            
            g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
            prior_g_idht = variable(am, nw, decomp_group, :g_idhty, (i,end_day_id,last_hour,t,y))
            prior_c_idht = variable(am, nw, decomp_group, :c_idhty, (i,end_day_id,last_hour,t,y))
            
            if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$y)"))) end    
            expr = JuMP.@expression(JuMP_model,  
                - g_idht + prior_g_idht - (c_idhty * Ramp_Rate * CAP * PMAX) - SD_Ramp_Rate * CAP * PMAX * (su_idht - c_idhty + prior_c_idht)
                )   
            constraint = JuMP.@constraint(JuMP_model, 
                expr <= 0
                )    
            if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$y)")) end    
        
        end
    elseif (h != 1) 
        HFREQ = am.setting["Simulation Configuration"]["HFREQ"]
        
        g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
        prior_g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h-HFREQ,t,y))
        prior_c_idht = variable(am, nw, decomp_group, :c_idhty, (i,d,h-HFREQ,t,y))
        
        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$y)"))) end    
        expr = JuMP.@expression(JuMP_model,  
            - g_idht + prior_g_idht - (c_idhty * Ramp_Rate * CAP * PMAX) - SD_Ramp_Rate * CAP * PMAX * (su_idht - c_idhty + prior_c_idht)
            )   
        constraint = JuMP.@constraint(JuMP_model, 
            expr <= 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$y)")) end    

    end
end

function constraint_Agg_Ramp_InterTemporal_Hour_idhy(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, r::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false, apply_PMAX::Bool=true)

    day_group = parse(Int, am.ref[:nw][0][:repdays][d]["Day_Group_ID"])
    HFREQ = am.setting["Simulation Configuration"]["HFREQ"]

    group_info = am.ref[:nw][0][:reserve_group_lookup][r]                                         
    unit_group = group_info["unit_group"]                                                         
    z = group_info["zone_idx"]                                                                    
    gen_idx = am.ref[:nw][0][:zone]["reserve"][string(z)]["aggregation_info"]["gen_technology_groups"][unit_group]["gen_idx"]

    if am.ref[:nw][0][:reserve_group_lookup][r]["UNIT_CATEGORY"] in ["THERMAL", "NUCLEAR", "OTHER"]

        sum_g = JuMP.AffExpr(0.0)
        sum_prior_g = JuMP.AffExpr(0.0)
        total_ramp = JuMP.AffExpr(0.0)

        for i in gen_idx

            bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
            tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
            CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")
            PMAX = parameter(am, bus_idx, :gen_bus, tech_idx, "PMAX")   # %
            if apply_PMAX == false
                PMAX = 1.0
            end
            Ramp_Rate = min(1.0, 60 * parameter(am, bus_idx, :gen_bus, tech_idx, "Ramp"))   # hourly Ramp Rate (%)

            JuMP.add_to_expression!(sum_g, 1, variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y)))
            JuMP.add_to_expression!(total_ramp, Ramp_Rate * CAP * PMAX, variable(am, nw, decomp_group, :u_G_iy, (i,y)))
        end

        if (h == 1)
            if am.ref[:nw][0][:repdays][d]["Day"] != am.ref[:nw][0][:repday_groups][day_group]["Start_Day_Id"]
                last_hour = last(am.setting["run_H"])
                
                for i in gen_idx
                    JuMP.add_to_expression!(sum_prior_g, 1, variable(am, nw, decomp_group, :g_idhty, (i,d-1,last_hour,t,y)))
                end
            
            elseif am.ref[:nw][0][:repdays][d]["Day"] == am.ref[:nw][0][:repday_groups][day_group]["Start_Day_Id"]
                
                num_days_in_group = length(am.ref[:nw][0][:repday_groups][day_group]["Day_Idx_List"])
                end_day_id = am.ref[:nw][0][:repday_groups][day_group]["Day_Idx_List"][num_days_in_group]
                last_hour = last(am.setting["run_H"])

                for i in gen_idx
                    JuMP.add_to_expression!(sum_prior_g, 1, variable(am, nw, decomp_group, :g_idhty, (i,end_day_id,last_hour,t,y)))
                end
               
            end
        elseif (h != 1) 
            
            for i in gen_idx
                JuMP.add_to_expression!(sum_prior_g, 1, variable(am, nw, decomp_group, :g_idhty, (i,d,h-HFREQ,t,y)))
            end

        end


        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($r,$d,$h,$y)"))) end    
        expr = JuMP.@expression(JuMP_model,  
            sum_g - sum_prior_g - total_ramp
            )   
        constraint = JuMP.@constraint(JuMP_model, 
            expr <= 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_up_($r,$d,$h,$y)")) end    

        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($r,$d,$h,$y)"))) end    
        expr = JuMP.@expression(JuMP_model,  
            - sum_g + sum_prior_g - total_ramp
            )   
        constraint = JuMP.@constraint(JuMP_model, 
            expr <= 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_dn_($r,$d,$h,$y)")) end    
    end
    

end

function constraint_RampUp_InterTemporal_Hour_idhy(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false, apply_PMAX::Bool=true)

    day_group = parse(Int, am.ref[:nw][0][:repdays][d]["Day_Group_ID"])
    HFREQ = am.setting["Simulation Configuration"]["HFREQ"]

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
    CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")
    PMAX = parameter(am, bus_idx, :gen_bus, tech_idx, "PMAX")   # %
    if apply_PMAX == false
        PMAX = 1.0
    end
    Ramp_Rate = min(1.0, 60 * parameter(am, bus_idx, :gen_bus, tech_idx, "Ramp"))   # hourly Ramp Rate (%)

    u_G_iy = variable(am, nw, decomp_group, :u_G_iy, (i,y))

    if (h == 1)
        if am.ref[:nw][0][:repdays][d]["Day"] != am.ref[:nw][0][:repday_groups][day_group]["Start_Day_Id"]
            last_hour = last(am.setting["run_H"])
            
            g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
            prior_g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d-1,last_hour,t,y))
            
            if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$y)"))) end    
            expr = JuMP.@expression(JuMP_model,  
                g_idht - prior_g_idht - (u_G_iy * Ramp_Rate * CAP * PMAX)
                )   
            constraint = JuMP.@constraint(JuMP_model, 
                expr <= 0
                )    
            if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$y)")) end    
        
        elseif am.ref[:nw][0][:repdays][d]["Day"] == am.ref[:nw][0][:repday_groups][day_group]["Start_Day_Id"]
            
            num_days_in_group = length(am.ref[:nw][0][:repday_groups][day_group]["Day_Idx_List"])
            end_day_id = am.ref[:nw][0][:repday_groups][day_group]["Day_Idx_List"][num_days_in_group]
            last_hour = last(am.setting["run_H"])

            g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
            prior_g_idht = variable(am, nw, decomp_group, :g_idhty, (i,end_day_id,last_hour,t,y))
            
            if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$y)"))) end    
            expr = JuMP.@expression(JuMP_model,  
                g_idht - prior_g_idht - (u_G_iy * Ramp_Rate * CAP * PMAX)
                )   
            constraint = JuMP.@constraint(JuMP_model, 
                expr <= 0
                )    
            if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$y)")) end    
        
        end
    elseif (h != 1) 
        HFREQ = am.setting["Simulation Configuration"]["HFREQ"]
        
        g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
        prior_g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h-HFREQ,t,y))
        
        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$y)"))) end    
        expr = JuMP.@expression(JuMP_model,  
            g_idht - prior_g_idht - (u_G_iy * Ramp_Rate * CAP * PMAX)
            )   
        constraint = JuMP.@constraint(JuMP_model, 
            expr <= 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$y)")) end    

    end
end

function constraint_RampDn_InterTemporal_Hour_idhy(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false, apply_PMAX::Bool=true)

    day_group = parse(Int, am.ref[:nw][0][:repdays][d]["Day_Group_ID"])
    HFREQ = am.setting["Simulation Configuration"]["HFREQ"]

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
    CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")
    PMAX = parameter(am, bus_idx, :gen_bus, tech_idx, "PMAX")   # %
    if apply_PMAX == false
        PMAX = 1.0
    end
    Ramp_Rate = min(1.0, 60 * parameter(am, bus_idx, :gen_bus, tech_idx, "Ramp"))   # hourly Ramp Rate (%)

    u_G_iy = variable(am, nw, decomp_group, :u_G_iy, (i,y))

    if (h == 1)
        if am.ref[:nw][0][:repdays][d]["Day"] != am.ref[:nw][0][:repday_groups][day_group]["Start_Day_Id"]
            last_hour = last(am.setting["run_H"])
            
            g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
            prior_g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d-1,last_hour,t,y))
            
            if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$y)"))) end    
            expr = JuMP.@expression(JuMP_model,  
                - g_idht + prior_g_idht - (u_G_iy * Ramp_Rate * CAP * PMAX)
                )   
            constraint = JuMP.@constraint(JuMP_model, 
                expr <= 0
                )    
            if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$y)")) end    
        
        elseif am.ref[:nw][0][:repdays][d]["Day"] == am.ref[:nw][0][:repday_groups][day_group]["Start_Day_Id"]
            num_days_in_group = length(am.ref[:nw][0][:repday_groups][day_group]["Day_Idx_List"])
            end_day_id = am.ref[:nw][0][:repday_groups][day_group]["Day_Idx_List"][num_days_in_group]
            last_hour = last(am.setting["run_H"])
            
            g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
            prior_g_idht = variable(am, nw, decomp_group, :g_idhty, (i,end_day_id,last_hour,t,y))
            
            if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$y)"))) end    
            expr = JuMP.@expression(JuMP_model,  
                - g_idht + prior_g_idht - (u_G_iy * Ramp_Rate * CAP * PMAX)
                )   
            constraint = JuMP.@constraint(JuMP_model, 
                expr <= 0
                )    
            if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$y)")) end    
        
        end
    elseif (h != 1) 
        HFREQ = am.setting["Simulation Configuration"]["HFREQ"]
        
        g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
        prior_g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h-HFREQ,t,y))
        
        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$y)"))) end    
        expr = JuMP.@expression(JuMP_model,  
            - g_idht + prior_g_idht - (u_G_iy * Ramp_Rate * CAP * PMAX)
            )   
        constraint = JuMP.@constraint(JuMP_model, 
            expr <= 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$y)")) end    

    end
end

function constraint_RampUp_InterTemporal_Hour_idhy_OP(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int, day_group_id; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false, apply_PMAX::Bool=true)

    day_group = day_group_id
    HFREQ = am.setting["Simulation Configuration"]["HFREQ"]

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
    CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")
    PMAX = parameter(am, bus_idx, :gen_bus, tech_idx, "PMAX")   # %
    if apply_PMAX == false
        PMAX = 1.0
    end
    Ramp_Rate = min(1.0, 60 * parameter(am, bus_idx, :gen_bus, tech_idx, "Ramp"))   # hourly Ramp Rate (%)
    
    u_G_iy = variable(am, nw, decomp_group, :u_G_iy, (i,y))

    if (h == 1)
        if am.ref[:nw][0][:repdays][d]["Day"] != am.ref[:nw][0][:repday_groups][day_group]["Start_Day_Id"]
            last_hour = last(am.setting["run_H"])
                
            g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
            prior_g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d-1,last_hour,t,y))
            
            if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$y)"))) end    
            expr = JuMP.@expression(JuMP_model,  
                g_idht - prior_g_idht - (u_G_iy * Ramp_Rate * CAP * PMAX)
                )   
            constraint = JuMP.@constraint(JuMP_model, 
                expr <= 0
                )    
            if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$y)")) end    
            
        
        elseif am.ref[:nw][0][:repdays][d]["Day"] == am.ref[:nw][0][:repday_groups][day_group]["Start_Day_Id"]
            num_days_in_group = length(am.ref[:nw][0][:repday_groups][day_group]["Day_Idx_List"])
            end_day_id = am.ref[:nw][0][:repday_groups][day_group]["Day_Idx_List"][num_days_in_group]
            last_hour = last(am.setting["run_H"])
            
            g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
            prior_g_idht = variable(am, nw, decomp_group, :g_idhty, (i,end_day_id,last_hour,t,y))
            
            if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$y)"))) end    
            expr = JuMP.@expression(JuMP_model,  
                g_idht - prior_g_idht - (u_G_iy * Ramp_Rate * CAP * PMAX)
                )   
            constraint = JuMP.@constraint(JuMP_model, 
                expr <= 0
                )    
            if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$y)")) end    
        
        end
    elseif (h != 1) 
        HFREQ = am.setting["Simulation Configuration"]["HFREQ"]
        
        g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
        prior_g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h-HFREQ,t,y))
        
        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$y)"))) end    
        expr = JuMP.@expression(JuMP_model,  
            g_idht - prior_g_idht - (u_G_iy * Ramp_Rate * CAP * PMAX)
            )   
        constraint = JuMP.@constraint(JuMP_model, 
            expr <= 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$y)")) end    

    end
end

function constraint_RampDn_InterTemporal_Hour_idhy_OP(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int, day_group_id; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false, apply_PMAX::Bool=true)

    day_group = day_group_id
    HFREQ = am.setting["Simulation Configuration"]["HFREQ"]

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
    CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")
    PMAX = parameter(am, bus_idx, :gen_bus, tech_idx, "PMAX")   # %
    if apply_PMAX == false
        PMAX = 1.0
    end
    Ramp_Rate = min(1.0, 60 * parameter(am, bus_idx, :gen_bus, tech_idx, "Ramp"))   # hourly Ramp Rate (%)

    u_G_iy = variable(am, nw, decomp_group, :u_G_iy, (i,y))

    if (h == 1)
        if am.ref[:nw][0][:repdays][d]["Day"] != am.ref[:nw][0][:repday_groups][day_group]["Start_Day_Id"]
            last_hour = last(am.setting["run_H"])
            
            g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
            prior_g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d-1,last_hour,t,y))
            
            if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$y)"))) end    
            expr = JuMP.@expression(JuMP_model,  
                - g_idht + prior_g_idht - (u_G_iy * Ramp_Rate * CAP * PMAX)
                )   
            constraint = JuMP.@constraint(JuMP_model, 
                expr <= 0
                )    
            if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$y)")) end    
        
        elseif am.ref[:nw][0][:repdays][d]["Day"] == am.ref[:nw][0][:repday_groups][day_group]["Start_Day_Id"]
            num_days_in_group = length(am.ref[:nw][0][:repday_groups][day_group]["Day_Idx_List"])
            end_day_id = am.ref[:nw][0][:repday_groups][day_group]["Day_Idx_List"][num_days_in_group]
            last_hour = last(am.setting["run_H"])
            
            g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
            prior_g_idht = variable(am, nw, decomp_group, :g_idhty, (i,end_day_id,last_hour,t,y))
            
            if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$y)"))) end    
            expr = JuMP.@expression(JuMP_model,  
                - g_idht + prior_g_idht - (u_G_iy * Ramp_Rate * CAP * PMAX)
                )   
            constraint = JuMP.@constraint(JuMP_model, 
                expr <= 0
                )    
            if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$y)")) end    
        
        end
    elseif (h != 1) 
        HFREQ = am.setting["Simulation Configuration"]["HFREQ"]
        
        g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
        prior_g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h-HFREQ,t,y))
        
        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$y)"))) end    
        expr = JuMP.@expression(JuMP_model,  
            - g_idht + prior_g_idht - (u_G_iy * Ramp_Rate * CAP * PMAX)
            )   
        constraint = JuMP.@constraint(JuMP_model, 
            expr <= 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$y)")) end    

    end
end

function constraint_NonSpinMax_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    JuMP.fix(variable(am, nw, decomp_group, :nonspin_idhty, (i,d,h,t,y)), 0.0; force = true)

end

function constraint_ES_Sto_UC_Limit_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
    
    if parameter(am, bus_idx, :gen_bus, tech_idx, "Storage Commitment") == true

        u_i = variable(am, nw, decomp_group, :u_G_iy, (i,y))
        sto_c_idhty = variable(am, nw, decomp_group, :sto_c_idhty, (i,d,h,t,y))
        
        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
        expr = JuMP.@expression(JuMP_model,  
            sto_c_idhty - u_i 
            )   
        constraint = JuMP.@constraint(JuMP_model, 
            expr <= 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end    

    end

end

function constraint_b_theta_power_flow_kdhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, k::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    t_bus = parameter(am, 0, :branch, k, "t_bus")
    f_bus = parameter(am, 0, :branch, k, "f_bus")

    b = 1.0 / parameter(am, 0, :branch, k, "br_x_pu")
    # geometric-mean row scaling keeps the ~1e5 susceptance off the matrix ceiling; abs() handles
    # negative-reactance branches (sign stays on the angle term via b/sb). Exact row-scaling, solution-neutral.
    sb = sqrt(abs(b))

    f_kdhty = variable(am, nw, decomp_group, :f_kdhty, (k,d,h,t,y))
    t_bus_angle = variable(am, nw, decomp_group, :bus_angle_ndhty, (t_bus,d,h,t,y))
    f_bus_angle = variable(am, nw, decomp_group, :bus_angle_ndhty, (f_bus,d,h,t,y))

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($k,$d,$h,$t,$y)"))) end
    expr = JuMP.@expression(JuMP_model,
        (1.0 / sb) * f_kdhty - (b / sb) * (f_bus_angle - t_bus_angle)
        )
    constraint = JuMP.@constraint(JuMP_model,
        expr == 0
        )
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($k,$d,$h,$t,$y)")) end
end

# --- Large Flexible Load (LFL) effective per-hour capacity -------------------------------------
# Folds two features into the single point where CAP sizes an LFL's hourly demand:
#   1. Online-year gate — an LFL (e.g. a data center) is inactive in planning years before its
#      Online_Year, mirroring the plant online/retire pattern. Online_Year == "NA" => always on.
#   2. Per-LFL hourly profile — if the LFL names a "Profile_Tag" and a "dc" timeseries is loaded, its
#      demand follows that hourly shape (0..1); otherwise the shape defaults to 1.0 (flat CAP).
# Both default to no-op, so DBs without Online_Year / Profile_Tag / a dc timeseries behave as before.
function lfl_online(am::Abstract_ALEAF_Model, lfl::Int, y::Int)
    oy = am.ref[:nw][0][:demand][lfl]["Online_Year"]
    (oy === nothing || oy == "NA") && return true
    return oy <= am.ref[:nw][0][:planning_stages][y]["year"]
end

function lfl_profile_shape(am::Abstract_ALEAF_Model, lfl::Int, d::Int, h::Int, t::Int, y::Int)
    tag = get(am.ref[:nw][0][:demand][lfl], "Profile_Tag", "NA")
    (tag === nothing || tag == "NA") && return 1.0
    data_ht = am.ref[:nw][0][:planning_stages][y]["repdays"][string(d)]["data"][string(h)][string(t)]
    haskey(data_ht, "dc") || return 1.0
    return get(data_ht["dc"], tag, 1.0)
end

function lfl_effective_cap(am::Abstract_ALEAF_Model, lfl::Int, d::Int, h::Int, t::Int, y::Int)
    lfl_online(am, lfl, y) || return 0.0
    return am.ref[:nw][0][:demand][lfl]["CAP"] * lfl_profile_shape(am, lfl, d, h, t, y)
end

function constraint_LFL_power_balance_ldhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, lfl::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    segment_vector = [s for s in 1:am.ref[:nw][0][:demand][lfl]["Num_DR_Segments"]]

    
    lfl_ldhty = variable(am, nw, decomp_group, :lfl_ldhty, (lfl,d,h,t,y))

    sum_g_idht = JuMP.AffExpr(0.0)
    for s in segment_vector
        JuMP.add_to_expression!(sum_g_idht, 1, variable(am, nw, decomp_group, :lfl_seg_lsdhty, (lfl,s,d,h,t,y)))
    end

    if am.ref[:nw][0][:demand][lfl]["Hybrid_Gen"] != "NA"
        JuMP.add_to_expression!(sum_g_idht, -1, variable(am, nw, decomp_group, :lfl_g_G_LFL_ldhty, (lfl,d,h,t,y)))
    end

    if am.ref[:nw][0][:demand][lfl]["Hybrid_ES"] != "NA"
        JuMP.add_to_expression!(sum_g_idht, -1, variable(am, nw, decomp_group, :lfl_g_ES_LFL_ldhty, (lfl,d,h,t,y)))
    end

    expr = JuMP.@expression(JuMP_model,  
        lfl_ldhty - sum_g_idht
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($lfl,$d,$h,$t,$y)")) end   
end

function constraint_LFL_Limit_ldhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, lfl::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    CAP = lfl_effective_cap(am, lfl, d, h, t, y)   # online-year gate x hourly profile (see helpers above)

    lfl_ldhty = variable(am, nw, decomp_group, :lfl_ldhty, (lfl,d,h,t,y))

    expr = JuMP.@expression(JuMP_model,  
        lfl_ldhty - CAP
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($lfl,$d,$h,$t,$y)")) end   
end

function constraint_LFL_inter_connection_limit_injection_ldhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, lfl::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    add_constraint_flag = false

    INTERCON_LIM = am.ref[:nw][0][:demand][lfl]["INTERCON_LIM"]
    
    sum_g_idht = JuMP.AffExpr(0.0)
    if am.ref[:nw][0][:demand][lfl]["Hybrid_Gen"] != "NA"
        JuMP.add_to_expression!(sum_g_idht, 1, variable(am, nw, decomp_group, :lfl_g_G_Grid_ldhty, (lfl,d,h,t,y)))
        add_constraint_flag = true
    end

    if am.ref[:nw][0][:demand][lfl]["Hybrid_ES"] != "NA"
        JuMP.add_to_expression!(sum_g_idht, 1, variable(am, nw, decomp_group, :lfl_g_ES_Grid_ldhty, (lfl,d,h,t,y)))
        add_constraint_flag = true
    end

    if add_constraint_flag == true
        expr = JuMP.@expression(JuMP_model,  
            sum_g_idht - INTERCON_LIM
            )   
        constraint = JuMP.@constraint(JuMP_model, 
            expr <= 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($lfl,$d,$h,$t,$y)")) end   
    end
end

function constraint_LFL_inter_connection_limit_withdraw_ldhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, lfl::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)


    INTERCON_LIM = am.ref[:nw][0][:demand][lfl]["INTERCON_LIM"]
    
    sum_g_idht = JuMP.AffExpr(0.0)
    JuMP.add_to_expression!(sum_g_idht, 1, variable(am, nw, decomp_group, :lfl_ldhty, (lfl,d,h,t,y)))
       
    if am.ref[:nw][0][:demand][lfl]["Hybrid_ES"] != "NA"
        JuMP.add_to_expression!(sum_g_idht, 1, variable(am, nw, decomp_group, :lfl_chg_Grid_ES_ldhty, (lfl,d,h,t,y)))
    end

    expr = JuMP.@expression(JuMP_model,  
        sum_g_idht - INTERCON_LIM
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($lfl,$d,$h,$t,$y)")) end   
end

function constraint_LFL_segment_bound_lsdhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, lfl::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    segment_vector = [s for s in 1:am.ref[:nw][0][:demand][lfl]["Num_DR_Segments"]]

    for s in segment_vector 

        Pct_MW = am.ref[:nw][0][:demand][lfl][string("Pct_MW_",s)]
        CAP = lfl_effective_cap(am, lfl, d, h, t, y)   # online-year gate x hourly profile
        seg_max = Pct_MW * CAP

        lfl_seg_lsdhty = variable(am, nw, decomp_group, :lfl_seg_lsdhty, (lfl,s,d,h,t,y))
        
        if s == 1

            expr = JuMP.@expression(JuMP_model,  
                lfl_seg_lsdhty - seg_max
                )   
            constraint = JuMP.@constraint(JuMP_model, 
                expr == 0
                )    
            if const_name_flag JuMP.set_name(constraint, string(const_name, "_($lfl,$s,$d,$h,$t,$y)")) end   

        else
        
            lfl_ind_lsdhty = variable(am, nw, decomp_group, :lfl_ind_lsdhty, (lfl,s,d,h,t,y))

            expr = JuMP.@expression(JuMP_model,  
                lfl_seg_lsdhty - seg_max * lfl_ind_lsdhty
                )   
            constraint = JuMP.@constraint(JuMP_model, 
                expr <= 0
                )    
            if const_name_flag JuMP.set_name(constraint, string(const_name, "_($lfl,$s,$d,$h,$t,$y)")) end   

        end
    end
    
end

function constraint_LFL_segment_relation_lsdhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, lfl::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    segment_vector = [s for s in 1:am.ref[:nw][0][:demand][lfl]["Num_DR_Segments"]]

    for s in segment_vector 

        if s == 1
            continue # skip the first segment
        elseif (s+1 <= length(segment_vector))
            next_s = s + 1
            lfl_ind_lsdhty = variable(am, nw, decomp_group, :lfl_ind_lsdhty, (lfl,s,d,h,t,y))
            next_lfl_ind_lsdhty = variable(am, nw, decomp_group, :lfl_ind_lsdhty, (lfl,next_s,d,h,t,y))

            expr = JuMP.@expression(JuMP_model,  
                lfl_ind_lsdhty - next_lfl_ind_lsdhty
                )   
            constraint = JuMP.@constraint(JuMP_model, 
                expr >= 0
                )    
            if const_name_flag JuMP.set_name(constraint, string(const_name, "_($lfl,$s,$d,$h,$t,$y)")) end   

        end
        

    end
    
end

function constraint_lfl_DR_balance_ldhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, lfl::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    segment_vector = [s for s in 1:am.ref[:nw][0][:demand][lfl]["Num_DR_Segments"]]

    total_seg_max = 0.0
    for s in segment_vector
        Pct_MW = am.ref[:nw][0][:demand][lfl][string("Pct_MW_",s)]
        CAP = lfl_effective_cap(am, lfl, d, h, t, y)   # online-year gate x hourly profile (matches segment_bound)
        total_seg_max += Pct_MW * CAP
    end
    
    lfl_ldhty = variable(am, nw, decomp_group, :lfl_ldhty, (lfl,d,h,t,y))
    lfl_DR_ldhty = variable(am, nw, decomp_group, :lfl_DR_ldhty, (lfl,d,h,t,y))

    expr = JuMP.@expression(JuMP_model,  
        lfl_ldhty + lfl_DR_ldhty - total_seg_max
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($lfl,$d,$h,$t,$y)")) end   
end

function constraint_lfl_DR_daily_limit_ldhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, lfl::Int, d::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    ids_ht = [(h, t) for (h) in am.setting["run_H"] for (t) in am.setting["run_T"]]

    Daily_DR_Limit_MWh = am.ref[:nw][0][:demand][lfl]["Daily_DR_Limit_MWh"]
    
    sum_DR_idht = JuMP.AffExpr(0.0)
    for (h, t) in ids_ht
        JuMP.add_to_expression!(sum_DR_idht, 1, variable(am, nw, decomp_group, :lfl_DR_ldhty, (lfl,d,h,t,y)))
    end

    expr = JuMP.@expression(JuMP_model,  
        sum_DR_idht - Daily_DR_Limit_MWh
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($lfl,$d,$y)")) end   
end

function constraint_lfl_onsite_gen_thermal_cap_ldhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, lfl::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    
    bus_idx = parse(Int, parameter(am, 0, :demand, "bus_idx", lfl))
    CAP = am.ref[:nw][0][:demand][lfl]["Component_gen_tech_data"]["CAP"]
    profile_type = am.ref[:nw][0][:demand][lfl]["Component_gen_tech_data"]["Profile_Type"]

    # VRE shape (available shapes: wind_ons, wind_ofs, pv, rtpv, hydro, csp)
    shape = 1.0
    type_key = profile_type == "NA" ? "NA" : profile_type * "_shape"

    if type_key != "NA"
        shape = get_vre_zdt_shape(am.ref[:nw][0], y, bus_idx, d, h, type_key)
    end
    
    lfl_g_G_LFL_ldhty = variable(am, nw, decomp_group, :lfl_g_G_LFL_ldhty, (lfl,d,h,t,y))
    lfl_g_G_Grid_ldhty = variable(am, nw, decomp_group, :lfl_g_G_Grid_ldhty, (lfl,d,h,t,y))
    lfl_g_G_ES_ldhty = variable(am, nw, decomp_group, :lfl_g_G_ES_ldhty, (lfl,d,h,t,y))

    expr = JuMP.@expression(JuMP_model,  
        lfl_g_G_LFL_ldhty + lfl_g_G_Grid_ldhty + lfl_g_G_ES_ldhty - (CAP * shape)
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($lfl,$d,$h,$t,$y)")) end   
end

function constraint_lfl_onsite_ES_charge_cap_ldhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, lfl::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    
    CAP =  am.ref[:nw][0][:demand][lfl]["Component_ES_tech_data"]["Charge_CAP"]

    sum_charge_idht = JuMP.AffExpr(0.0)
    JuMP.add_to_expression!(sum_charge_idht, 1, variable(am, nw, decomp_group, :lfl_chg_Grid_ES_ldhty, (lfl,d,h,t,y)))
    if am.ref[:nw][0][:demand][lfl]["Hybrid_Gen"] != "NA"
        JuMP.add_to_expression!(sum_charge_idht, 1, variable(am, nw, decomp_group, :lfl_g_G_ES_ldhty, (lfl,d,h,t,y)))
    end

    expr = JuMP.@expression(JuMP_model,  
        sum_charge_idht - CAP
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($lfl,$d,$h,$t,$y)")) end   
end

function constraint_lfl_onsite_ES_discharge_cap_ldhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, lfl::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    
    CAP =  am.ref[:nw][0][:demand][lfl]["Component_ES_tech_data"]["CAP"]

    lfl_g_ES_Grid_ldhty = variable(am, nw, decomp_group, :lfl_g_ES_Grid_ldhty, (lfl,d,h,t,y))
    lfl_g_ES_LFL_ldhty = variable(am, nw, decomp_group, :lfl_g_ES_LFL_ldhty, (lfl,d,h,t,y))

    expr = JuMP.@expression(JuMP_model,  
        lfl_g_ES_Grid_ldhty + lfl_g_ES_LFL_ldhty - CAP
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($lfl,$d,$h,$t,$y)")) end   
end

function constraint_lfl_onsite_ES_SOC_cap_ldhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, lfl::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    
    STOMIN =  am.ref[:nw][0][:demand][lfl]["Component_ES_tech_data"]["STOMIN"]
    # Energy capacity (per-unit MWh) = storage power rating x maximum duration, as for hybrid plants.
    STOMAX =  am.ref[:nw][0][:demand][lfl]["Component_ES_tech_data"]["Charge_CAP"] * am.ref[:nw][0][:demand][lfl]["Component_ES_tech_data"]["STOHR_MAX"]

    lfl_soc_ldhty = variable(am, nw, decomp_group, :lfl_soc_ldhty, (lfl,d,h,t,y))
    
    if STOMIN > 0
        expr = JuMP.@expression(JuMP_model,  
            lfl_soc_ldhty - STOMIN
            )   
        constraint = JuMP.@constraint(JuMP_model, 
            expr >= 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_min_($lfl,$d,$h,$t,$y)")) end   
    end
    
    expr = JuMP.@expression(JuMP_model,  
        lfl_soc_ldhty - STOMAX
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_max_($lfl,$d,$h,$t,$y)")) end   
end

function constraint_lfl_onsite_ES_SOC_balance_ldhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, lfl::Int, d::Int, h::Int, t::Int, y::Int, day_group::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    start_day = am.ref[:nw][0][:repday_groups][day_group]["Start_Day_Id"]
    end_day = am.ref[:nw][0][:repday_groups][day_group]["End_Day_Id"]
    start_day_idx = am.ref[:nw][0][:repday_groups][day_group]["Day_Idx_List"][1]
    card_T = length(am.setting["run_T"])
    last_hour = last(am.setting["run_H"])
    HFREQ = am.setting["Simulation Configuration"]["HFREQ"]

    BATEFF = sqrt(am.ref[:nw][0][:demand][lfl]["Component_ES_tech_data"]["BATEFF"])
    STOMIN = am.ref[:nw][0][:demand][lfl]["Component_ES_tech_data"]["STOMIN"]
    # Energy capacity (per-unit MWh) = storage power rating x maximum duration, as for hybrid plants.
    STOMAX =  am.ref[:nw][0][:demand][lfl]["Component_ES_tech_data"]["Charge_CAP"] * am.ref[:nw][0][:demand][lfl]["Component_ES_tech_data"]["STOHR_MAX"]

    INI_SOC = "Optimal"
    if am.setting["Simulation Configuration"]["storage initialization option"] == "Minimum"
        INI_SOC = STOMIN
    elseif am.setting["Simulation Configuration"]["storage initialization option"] == "Middle"
        INI_SOC = STOMAX / 2
    elseif am.setting["Simulation Configuration"]["storage initialization option"] == "Maximum"
        INI_SOC = STOMAX
    end

    if (h==1) && (t==1)
        if am.ref[:nw][0][:repdays][d]["Day"] == start_day

            if INI_SOC != "Optimal"
                lfl_soc_ldhty = variable(am, nw, decomp_group, :lfl_soc_ldhty, (lfl,d,h,t,y))

                sum_charge_idht = JuMP.AffExpr(0.0)
                JuMP.add_to_expression!(sum_charge_idht, 1, variable(am, nw, decomp_group, :lfl_chg_Grid_ES_ldhty, (lfl,d,h,t,y)))
                if am.ref[:nw][0][:demand][lfl]["Hybrid_Gen"] != "NA"
                    JuMP.add_to_expression!(sum_charge_idht, 1, variable(am, nw, decomp_group, :lfl_g_G_ES_ldhty, (lfl,d,h,t,y)))
                end

                lfl_g_ES_Grid_ldhty = variable(am, nw, decomp_group, :lfl_g_ES_Grid_ldhty, (lfl,d,h,t,y))
                lfl_g_ES_LFL_ldhty = variable(am, nw, decomp_group, :lfl_g_ES_LFL_ldhty, (lfl,d,h,t,y))

                if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($lfl,$d,$h,$t,$y)"))) end    
                expr = JuMP.@expression(JuMP_model,  
                    lfl_soc_ldhty - INI_SOC - (BATEFF * sum_charge_idht - (1/BATEFF) * (lfl_g_ES_Grid_ldhty + lfl_g_ES_LFL_ldhty)) / card_T
                    )   
                constraint = JuMP.@constraint(JuMP_model, 
                    expr == 0
                    )    
                if const_name_flag JuMP.set_name(constraint, string(const_name, "_min_($lfl,$d,$h,$t,$y)")) end    

            end
            
        else
            lfl_soc_ldhty = variable(am, nw, decomp_group, :lfl_soc_ldhty, (lfl,d,h,t,y))
            prior_lfl_soc_ldhty = variable(am, nw, decomp_group, :lfl_soc_ldhty, (lfl,d-1,last_hour,1,y)) # Last time step of last hour is T=1

            sum_charge_idht = JuMP.AffExpr(0.0)
            JuMP.add_to_expression!(sum_charge_idht, 1, variable(am, nw, decomp_group, :lfl_chg_Grid_ES_ldhty, (lfl,d,h,t,y)))
            if am.ref[:nw][0][:demand][lfl]["Hybrid_Gen"] != "NA"
                JuMP.add_to_expression!(sum_charge_idht, 1, variable(am, nw, decomp_group, :lfl_g_G_ES_ldhty, (lfl,d,h,t,y)))
            end

            lfl_g_ES_Grid_ldhty = variable(am, nw, decomp_group, :lfl_g_ES_Grid_ldhty, (lfl,d,h,t,y))
            lfl_g_ES_LFL_ldhty = variable(am, nw, decomp_group, :lfl_g_ES_LFL_ldhty, (lfl,d,h,t,y))

            if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($lfl,$d,$h,$t,$y)"))) end    
            expr = JuMP.@expression(JuMP_model,  
                lfl_soc_ldhty - prior_lfl_soc_ldhty - (BATEFF * sum_charge_idht - (1/BATEFF) * (lfl_g_ES_Grid_ldhty + lfl_g_ES_LFL_ldhty)) / card_T
                )   
            constraint = JuMP.@constraint(JuMP_model, 
                expr == 0
                )    
            if const_name_flag JuMP.set_name(constraint, string(const_name, "_min_($lfl,$d,$h,$t,$y)")) end    
        end

    elseif (h != 1) && (t == 1)

        lfl_soc_ldhty = variable(am, nw, decomp_group, :lfl_soc_ldhty, (lfl,d,h,t,y))
        prior_lfl_soc_ldhty = variable(am, nw, decomp_group, :lfl_soc_ldhty, (lfl,d,h-HFREQ,1,y)) # Last time step of last hour is T=1

        sum_charge_idht = JuMP.AffExpr(0.0)
        JuMP.add_to_expression!(sum_charge_idht, 1, variable(am, nw, decomp_group, :lfl_chg_Grid_ES_ldhty, (lfl,d,h,t,y)))
        if am.ref[:nw][0][:demand][lfl]["Hybrid_Gen"] != "NA"
            JuMP.add_to_expression!(sum_charge_idht, 1, variable(am, nw, decomp_group, :lfl_g_G_ES_ldhty, (lfl,d,h,t,y)))
        end
        
        lfl_g_ES_Grid_ldhty = variable(am, nw, decomp_group, :lfl_g_ES_Grid_ldhty, (lfl,d,h,t,y))
        lfl_g_ES_LFL_ldhty = variable(am, nw, decomp_group, :lfl_g_ES_LFL_ldhty, (lfl,d,h,t,y))

        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($lfl,$d,$h,$t,$y)"))) end    
        expr = JuMP.@expression(JuMP_model,  
            lfl_soc_ldhty - prior_lfl_soc_ldhty - (BATEFF * sum_charge_idht - (1/BATEFF) * (lfl_g_ES_Grid_ldhty + lfl_g_ES_LFL_ldhty)) / card_T
            )   
        constraint = JuMP.@constraint(JuMP_model, 
            expr == 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_min_($lfl,$d,$h,$t,$y)")) end    
        
    end

end

function constraint_lfl_onsite_ES_SOC_neutral_ldhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, lfl::Int, d::Int, y::Int, start_day_idx::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    card_T = length(am.setting["run_T"])
    last_subhour = last(am.setting["run_T"])
    last_hour = last(am.setting["run_H"])
    BATEFF = sqrt(am.ref[:nw][0][:demand][lfl]["Component_ES_tech_data"]["BATEFF"])

    end_lfl_soc_ldhty = variable(am, nw, decomp_group, :lfl_soc_ldhty, (lfl,d,last_hour,last_subhour,y))
    begin_lfl_soc_ldhty = variable(am, nw, decomp_group, :lfl_soc_ldhty, (lfl,start_day_idx,1,1,y))

    sum_charge_idht = JuMP.AffExpr(0.0)
    JuMP.add_to_expression!(sum_charge_idht, 1, variable(am, nw, decomp_group, :lfl_chg_Grid_ES_ldhty, (lfl,start_day_idx,1,1,y)))
    if am.ref[:nw][0][:demand][lfl]["Hybrid_Gen"] != "NA"
        JuMP.add_to_expression!(sum_charge_idht, 1, variable(am, nw, decomp_group, :lfl_g_G_ES_ldhty, (lfl,start_day_idx,1,1,y)))
    end
    
    begin_lfl_g_ES_Grid_ldhty = variable(am, nw, decomp_group, :lfl_g_ES_Grid_ldhty, (lfl,start_day_idx,1,1,y))
    begin_lfl_g_ES_LFL_ldhty = variable(am, nw, decomp_group, :lfl_g_ES_LFL_ldhty, (lfl,start_day_idx,1,1,y))
    
    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($lfl,$d,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        end_lfl_soc_ldhty - begin_lfl_soc_ldhty + (BATEFF * sum_charge_idht - (1/BATEFF) * (begin_lfl_g_ES_Grid_ldhty + begin_lfl_g_ES_LFL_ldhty))/card_T
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($lfl,$d,$y)")) end    

end

function constraint_hybrid_onsite_gen_thermal_cap_ldhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false, reserve_flag::Bool=true)

    
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)
    VRE_Flag = parameter(am, bus_idx, :gen_bus, tech_idx, "VRE_Flag")
    CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")
    PMAX = parameter(am, bus_idx, :gen_bus, tech_idx, "PMAX")
    FUEL_LIMIT = parameter(am, bus_idx, :gen_bus, tech_idx, "FUEL_LIMIT")

    # VRE shape (available shapes: wind_ons, wind_ofs, pv, rtpv, hydro, csp)
    shape = 1.0
    profile_type = parameter(am, bus_idx, :gen_bus, tech_idx, "Profile_Type")
    type_key = profile_type == "NA" ? "NA" : profile_type * "_shape"

    if type_key != "NA"
        shape = get_vre_zdt_shape(am.ref[:nw][0], y, bus_idx, d, h, type_key)
    end
    
    sum_gen_idht = JuMP.AffExpr(0.0)

    JuMP.add_to_expression!(sum_gen_idht, 1, variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y)))
    JuMP.add_to_expression!(sum_gen_idht, 1, variable(am, nw, decomp_group, :g_G_ES_idhty, (i,d,h,t,y)))
    if am.setting["Simulation Configuration"]["operating_reserve_modeling_option"] == "individual"# [aggregated, individual]
        if reserve_flag && up_reserve_eligible(am, i)
            reserve_enabled(am, :reg_up)  && JuMP.add_to_expression!(sum_gen_idht, 1, variable(am, nw, decomp_group, :reg_up_idhty, (i,d,h,t,y)))
            reserve_enabled(am, :flex_up) && JuMP.add_to_expression!(sum_gen_idht, 1, variable(am, nw, decomp_group, :flex_up_idhty, (i,d,h,t,y)))
            reserve_enabled(am, :spin)    && JuMP.add_to_expression!(sum_gen_idht, 1, variable(am, nw, decomp_group, :spin_idhty, (i,d,h,t,y)))
        end
    end


    if VRE_Flag == true
        if FUEL_LIMIT == "Fixed Profile"

            expr = JuMP.@expression(JuMP_model,
                sum_gen_idht - (CAP * shape)
                )
            constraint = JuMP.@constraint(JuMP_model,
                expr <= 0
                )
            if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end
        else
            # Budget hydro: per-hour cap is full nameplate; energy budget is enforced by a separate aggregate constraint
            expr = JuMP.@expression(JuMP_model,
                sum_gen_idht - (CAP * PMAX)
                )
            constraint = JuMP.@constraint(JuMP_model,
                expr <= 0
                )
            if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end
        end
    end
end

function constraint_hybrid_inter_connection_limit_injection_ldhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false, reserve_flag::Bool=true)

    hybrid_gen_idx = i
    hybrid_ES_idx = am.ref[:nw][0][:gen_index][i]["hybrid_ES_gen_idx"]

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
    unit_group = parameter(am, 0, :gen_index, "UNIT_GROUP", i)  
    hybrid_ID = parse(Int, parameter(am, bus_idx, :gen_bus, tech_idx, "hybrid_ID"))
    
    INTERCON_LIM = am.ref[:nw][0][:hybrid][hybrid_ID]["INTERCON_LIM"]
    
    sum_injection_idht = JuMP.AffExpr(0.0)

    JuMP.add_to_expression!(sum_injection_idht, 1, variable(am, nw, decomp_group, :g_idhty, (hybrid_gen_idx,d,h,t,y)))
    if am.setting["Simulation Configuration"]["operating_reserve_modeling_option"] == "individual"# [aggregated, individual]
        if reserve_flag && up_reserve_eligible(am, hybrid_gen_idx)
            reserve_enabled(am, :reg_up)  && JuMP.add_to_expression!(sum_injection_idht, 1, variable(am, nw, decomp_group, :reg_up_idhty, (hybrid_gen_idx,d,h,t,y)))
            reserve_enabled(am, :flex_up) && JuMP.add_to_expression!(sum_injection_idht, 1, variable(am, nw, decomp_group, :flex_up_idhty, (hybrid_gen_idx,d,h,t,y)))
            reserve_enabled(am, :spin)    && JuMP.add_to_expression!(sum_injection_idht, 1, variable(am, nw, decomp_group, :spin_idhty, (hybrid_gen_idx,d,h,t,y)))
        end
    end

    JuMP.add_to_expression!(sum_injection_idht, 1, variable(am, nw, decomp_group, :g_idhty, (hybrid_ES_idx,d,h,t,y)))
    if am.setting["Simulation Configuration"]["operating_reserve_modeling_option"] == "individual"# [aggregated, individual]
        if reserve_flag && up_reserve_eligible(am, hybrid_ES_idx)
            reserve_enabled(am, :reg_up)  && JuMP.add_to_expression!(sum_injection_idht, 1, variable(am, nw, decomp_group, :reg_up_idhty, (hybrid_ES_idx,d,h,t,y)))
            reserve_enabled(am, :flex_up) && JuMP.add_to_expression!(sum_injection_idht, 1, variable(am, nw, decomp_group, :flex_up_idhty, (hybrid_ES_idx,d,h,t,y)))
            reserve_enabled(am, :spin)    && JuMP.add_to_expression!(sum_injection_idht, 1, variable(am, nw, decomp_group, :spin_idhty, (hybrid_ES_idx,d,h,t,y)))
        end
    end

    expr = JuMP.@expression(JuMP_model,  
        sum_injection_idht - INTERCON_LIM
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end   
end

function constraint_hybrid_inter_connection_limit_withdraw_ldhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false, reserve_flag::Bool=true)

    hybrid_gen_idx = i
    hybrid_ES_idx = am.ref[:nw][0][:gen_index][i]["hybrid_ES_gen_idx"]
    

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
    unit_group = parameter(am, 0, :gen_index, "UNIT_GROUP", i)  
    hybrid_ID = parse(Int, parameter(am, bus_idx, :gen_bus, tech_idx, "hybrid_ID"))
    
    INTERCON_LIM = am.ref[:nw][0][:hybrid][hybrid_ID]["INTERCON_LIM"]
    
    sum_withdraw_idht = JuMP.AffExpr(0.0)


    JuMP.add_to_expression!(sum_withdraw_idht, 1, variable(am, nw, decomp_group, :chg_idhty, (hybrid_ES_idx,d,h,t,y)))
    if am.setting["Simulation Configuration"]["operating_reserve_modeling_option"] == "individual"# [aggregated, individual]
        if reserve_flag
            reserve_enabled(am, :reg_dn)  && JuMP.add_to_expression!(sum_withdraw_idht, 1, variable(am, nw, decomp_group, :reg_dn_idhty, (hybrid_ES_idx,d,h,t,y)))
            reserve_enabled(am, :flex_dn) && JuMP.add_to_expression!(sum_withdraw_idht, 1, variable(am, nw, decomp_group, :flex_dn_idhty, (hybrid_ES_idx,d,h,t,y)))
        end
    end

    expr = JuMP.@expression(JuMP_model,  
        sum_withdraw_idht - INTERCON_LIM
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end   
end


function constraint_ES_SOC_Balance_Inter_Hour_NoReserve_with_SATA_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int, day_group::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=true)

    start_day = am.ref[:nw][0][:repday_groups][day_group]["Start_Day_Id"]
    end_day = am.ref[:nw][0][:repday_groups][day_group]["End_Day_Id"]
    start_day_idx = am.ref[:nw][0][:repday_groups][day_group]["Day_Idx_List"][1]
    
    if (h==1) && (t==1)
        if am.ref[:nw][0][:repdays][d]["Day"] == start_day
            bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
            tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
            BATEFF = sqrt(parameter(am, bus_idx, :gen_bus, tech_idx, "BATEFF"))
            card_T = length(am.setting["run_T"])
            STOMIN = parameter(am, bus_idx, :gen_bus, tech_idx, "STOMIN")

            dual_use_percent = 0.0  # available SOC to be used in markets (pre-contingency)
            dual_use_flag = false
            if haskey(am.ref[:nw][bus_idx][:gen_bus][tech_idx], "asset_type")   # this is a new battery added as SATA
                if parameter(am, bus_idx, :gen_bus, tech_idx, "asset_type") == "Transmission"
                    dual_use_percent = parameter(am, bus_idx, :gen_bus, tech_idx, "dual_use_percent")
                    dual_use_flag = true
                end
            end

            soc_idhty = variable(am, nw, decomp_group, :soc_idhty, (i,d,h,t,y))
            chg_idhty = variable(am, nw, decomp_group, :chg_idhty, (i,d,h,t,y))
            g_idhty = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
            u_ESE_iy = parameter(am, bus_idx, :gen_bus, tech_idx, "ES_MWh")

            if dual_use_flag == true # this asset is SATA
                

                if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
                expr = JuMP.@expression(JuMP_model,  
                    soc_idhty - u_ESE_iy * (1 - dual_use_percent) - (BATEFF * chg_idhty - (1/BATEFF) * g_idhty) / card_T
                    )   
                constraint = JuMP.@constraint(JuMP_model, 
                    expr == 0
                    )    
                if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end    
            
            else

                if am.setting["Simulation Configuration"]["storage initialization option"] == "Minimum"
            
                    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
                    expr = JuMP.@expression(JuMP_model,  
                        soc_idhty - STOMIN - (BATEFF * chg_idhty - (1/BATEFF) * g_idhty) / card_T
                        )   
                    constraint = JuMP.@constraint(JuMP_model, 
                        expr == 0
                        )    
                    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end    

                elseif am.setting["Simulation Configuration"]["storage initialization option"] == "Middle"

                    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
                    expr = JuMP.@expression(JuMP_model,  
                        soc_idhty - u_ESE_iy / 2 - (BATEFF * chg_idhty - (1/BATEFF) * g_idhty) / card_T
                        )   
                    constraint = JuMP.@constraint(JuMP_model, 
                        expr == 0
                        )    
                    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end    

                elseif am.setting["Simulation Configuration"]["storage initialization option"] == "Maximum"

                    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
                    expr = JuMP.@expression(JuMP_model,  
                        soc_idhty - u_ESE_iy - (BATEFF * chg_idhty - (1/BATEFF) * g_idhty) / card_T
                        )   
                    constraint = JuMP.@constraint(JuMP_model, 
                        expr == 0
                        )    
                    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end    


                end
            end

        else
            bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
            tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
            
            BATEFF = sqrt(parameter(am, bus_idx, :gen_bus, tech_idx, "BATEFF"))
            HFREQ = am.setting["Simulation Configuration"]["HFREQ"]
            card_T = length(am.setting["run_T"])
            last_hour = last(am.setting["run_H"])
            
            soc_idhty = variable(am, nw, decomp_group, :soc_idhty, (i,d,h,t,y))
            chg_idhty = variable(am, nw, decomp_group, :chg_idhty, (i,d,h,t,y))
            prior_soc_idhty = variable(am, nw, decomp_group, :soc_idhty, (i,d-1,last_hour,1,y)) # Last time step of last hour is T=1
            g_idhty = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
                
            if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
            expr = JuMP.@expression(JuMP_model,  
                soc_idhty - prior_soc_idhty - (BATEFF * chg_idhty - (1/BATEFF) * g_idhty)/card_T
                )   
            constraint = JuMP.@constraint(JuMP_model, 
                expr == 0
                )    
            if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end    
        end

    elseif (h != 1) && (t == 1)
        bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
        tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
        
        BATEFF = sqrt(parameter(am, bus_idx, :gen_bus, tech_idx, "BATEFF"))
        HFREQ = am.setting["Simulation Configuration"]["HFREQ"]
        card_T = length(am.setting["run_T"])

        soc_idhty = variable(am, nw, decomp_group, :soc_idhty, (i,d,h,t,y))
        chg_idhty = variable(am, nw, decomp_group, :chg_idhty, (i,d,h,t,y))
        prior_soc_idhty = variable(am, nw, decomp_group, :soc_idhty, (i,d,h-HFREQ,1,y)) # Last time step of last hour is T=1
        g_idhty = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
        
        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
        expr = JuMP.@expression(JuMP_model,  
            soc_idhty - prior_soc_idhty - (BATEFF * chg_idhty - (1/BATEFF) * g_idhty)/card_T
            )   
        constraint = JuMP.@constraint(JuMP_model, 
            expr == 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end    
    end

end

function constraint_ES_SOC_Min_NonReserve_SATA_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=true)

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
    BATEFF = sqrt(sqrt(parameter(am, bus_idx, :gen_bus, tech_idx, "BATEFF")))

    dual_use_percent = 1.0  # available SOC to be used in markets (pre-contingency)
    if haskey(am.ref[:nw][bus_idx][:gen_bus][tech_idx], "asset_type")   # this is a new battery added as SATA
        if parameter(am, bus_idx, :gen_bus, tech_idx, "asset_type") == "Transmission"
            dual_use_percent = parameter(am, bus_idx, :gen_bus, tech_idx, "dual_use_percent")
        end
    end

    if dual_use_percent < 1.0

        soc_idhty = variable(am, nw, decomp_group, :soc_idhty, (i,d,h,t,y))
        u_ESE_iy = parameter(am, bus_idx, :gen_bus, tech_idx, "ES_MWh")

        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
        expr = JuMP.@expression(JuMP_model,  
            soc_idhty - u_ESE_iy * (1 - dual_use_percent)
            )   
        constraint = JuMP.@constraint(JuMP_model, 
            expr >= 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end    

    end
     

end

function _base_ra_bus_demand(am::Abstract_ALEAF_Model, n::Int, d::Int, h::Int, t::Int, y::Int; system_peak_scale=1.0)
    return apply_tnd_loss(am, get_bus_demand_with_growth(am, n, d, h, t, y; nw=0, system_peak_scale))
end

function _constant_load_adjustment_RA(am::Abstract_ALEAF_Model, n::Int, d::Int, h::Int, t::Int, y::Int; constant_load=0.0, constant_load_bus_idx=nothing, constant_load_distribution="single_bus", constant_load_distribution_profile=Dict{Tuple{Int, Int, Int}, Dict{Int, Float64}}(), system_peak_scale=1.0)
    if abs(constant_load) <= 1e-12
        return 0.0
    end

    if constant_load_distribution == "load_weighted_systemwide"
        bus_demand = _base_ra_bus_demand(am, n, d, h, t, y; system_peak_scale)
        system_demand = sum(_base_ra_bus_demand(am, bus_idx, d, h, t, y; system_peak_scale) for bus_idx in get_index(am, :bus, 0))
        return system_demand > 1e-9 ? constant_load * bus_demand / system_demand : 0.0
    elseif constant_load_distribution == "ens_weighted_systemwide"
        if !isempty(constant_load_distribution_profile)
            return constant_load * get(constant_load_distribution_profile, n, 0.0)
        end
        bus_demand = _base_ra_bus_demand(am, n, d, h, t, y; system_peak_scale)
        system_demand = sum(_base_ra_bus_demand(am, bus_idx, d, h, t, y; system_peak_scale) for bus_idx in get_index(am, :bus, 0))
        return system_demand > 1e-9 ? constant_load * bus_demand / system_demand : 0.0
    elseif constant_load_distribution == "single_bus"
        load_bus_idx = constant_load_bus_idx === nothing ? am.ref[:nw][0][:gen_index][maximum(get_index(am, :gen_index, 0))]["bus_idx"] : constant_load_bus_idx
        return n == load_bus_idx ? constant_load : 0.0
    else
        error("Unsupported RA constant_load_distribution: $(constant_load_distribution)")
    end
end

function constraint_PTDF_LoadBalance_ndhty_real_RA(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, n::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, constant_load=0.0, constant_load_bus_idx=nothing, constant_load_distribution="single_bus", constant_load_distribution_profile=Dict{Tuple{Int, Int, Int}, Dict{Int, Float64}}(), system_peak_scale=1.0)

    inc = gtep_bus_incidence(am)   # per-bus LFL incidence, precomputed once

    Demand = _base_ra_bus_demand(am, n, d, h, t, y; system_peak_scale)

    # add constant load (ELCC load-relief / ELCC load-addition adjustment).
    # Clamp post-relief Demand to >= 0: ENS upper bound is set to Demand;
    Demand = max(0.0, Demand + _constant_load_adjustment_RA(am, n, d, h, t, y; constant_load, constant_load_bus_idx, constant_load_distribution, constant_load_distribution_profile, system_peak_scale))
    
    ens_ndht = variable(am, nw, decomp_group, :ens_ndhty, (n,d,h,t,y))
    demand_ndhty = variable(am, nw, decomp_group, :demand_ndhty, (n,d,h,t,y))

    sum_lfl_ldht = JuMP.AffExpr(0.0)
    for lfl in get(inc.lfl_at_bus, n, _EMPTY_INT_VEC)
        JuMP.add_to_expression!(sum_lfl_ldht, 1, variable(am, nw, decomp_group, :lfl_ldhty, (lfl,d,h,t,y)))
    end
    
    # Mirrors the non-RA twin (constraint_PTDF_LoadBalance_ndhty_real): delivered demand
    # plus unmet load (ens) equals nominal Demand plus large flexible load (sum_lfl_ldht).
    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($n,$d,$h,$t,$y)"))) end
    expr = JuMP.@expression(JuMP_model,
        demand_ndhty + ens_ndht - sum_lfl_ldht - Demand
        )
    constraint = JuMP.@constraint(JuMP_model,
        expr == 0
        )
    JuMP.set_name(constraint, string(const_name, "_($n,$d,$h,$t,$y)"))

    JuMP.set_upper_bound(ens_ndht, Demand)
end

function constraint_LoadBalance_ndhty_real_RA(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, n::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, constant_load=0.0, constant_load_bus_idx=nothing, constant_load_distribution="single_bus", constant_load_distribution_profile=Dict{Tuple{Int, Int, Int}, Dict{Int, Float64}}(), system_peak_scale=1.0)

    ids_local_i = am.ref[:nw][0][:bus][n]["aggregation_info"]["new_local_gen_idx"]
    inc = gtep_bus_incidence(am)   # per-bus branch/LFL incidence, precomputed once (mirrors OP twin)

    Demand = _base_ra_bus_demand(am, n, d, h, t, y; system_peak_scale)

    # add constant load (ELCC load-relief / ELCC load-addition adjustment).
    # Clamp post-relief Demand to >= 0: ENS upper bound is set to Demand;
    Demand = max(0.0, Demand + _constant_load_adjustment_RA(am, n, d, h, t, y; constant_load, constant_load_bus_idx, constant_load_distribution, constant_load_distribution_profile, system_peak_scale))
    
    sum_g_idht = JuMP.AffExpr(0.0)
    for i in ids_local_i
        JuMP.add_to_expression!(sum_g_idht, 1, variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y)))
        if parameter(am, 0, :gen_index, "UNIT_CATEGORY", i) == "STORAGE"
            JuMP.add_to_expression!(sum_g_idht, -1, variable(am, nw, decomp_group, :chg_idhty, (i,d,h,t,y)))
        end
    end

    sum_lfl_ldht = JuMP.AffExpr(0.0)
    for lfl in get(inc.lfl_at_bus, n, _EMPTY_INT_VEC)
        JuMP.add_to_expression!(sum_lfl_ldht, 1, variable(am, nw, decomp_group, :lfl_ldhty, (lfl,d,h,t,y)))
    end

    ens_ndht = variable(am, nw, decomp_group, :ens_ndhty, (n,d,h,t,y))

    sum_f_ktd_to_node = JuMP.AffExpr(0.0)    # injection to the node n
    sum_f_ktd_from_node = JuMP.AffExpr(0.0)  # withdrawals from the node n
    for k in get(inc.to_branches, n, _EMPTY_INT_VEC)
        JuMP.add_to_expression!(sum_f_ktd_to_node, 1, variable(am, nw, decomp_group, :f_kdhty, (k,d,h,t,y)))
    end
    for k in get(inc.from_branches, n, _EMPTY_INT_VEC)
        JuMP.add_to_expression!(sum_f_ktd_from_node, 1, variable(am, nw, decomp_group, :f_kdhty, (k,d,h,t,y)))
    end

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($n,$d,$h,$t,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        - (sum_f_ktd_from_node - sum_f_ktd_to_node) + sum_g_idht + ens_ndht - sum_lfl_ldht - Demand
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
    JuMP.set_name(constraint, string(const_name, "_($n,$d,$h,$t,$y)"))    

    JuMP.set_upper_bound(ens_ndht, Demand)
end

function constraint_dc_power_flow_max_kdhty_RA(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, k::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false)


    rate_a = parameter(am, 0, :branch, "rate_a", k) 

    f_kdhty = variable(am, nw, decomp_group, :f_kdhty, (k, d, h, t, y))

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($k,$d,$h,$t,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        f_kdhty - rate_a
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    JuMP.set_name(constraint, string(const_name, "_($k,$d,$h,$t,$y)"))    

end

function constraint_dc_power_flow_min_kdhty_RA(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, k::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false)


    rate_a = parameter(am, 0, :branch, "rate_a", k) 

    f_kdhty = variable(am, nw, decomp_group, :f_kdhty, (k, d, h, t, y))

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($k,$d,$h,$t,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        f_kdhty + rate_a 
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr >= 0
        )    
    JuMP.set_name(constraint, string(const_name, "_($k,$d,$h,$t,$y)"))    

end

function constraint_dc_power_flow_kdhty_RA(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, k::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, constant_load=0.0, system_peak_scale=1.0)

    ptdf_k = parameter(am, 0, :branch, k, "ptdf")

    f_ptdf_kdhty = JuMP.AffExpr(0.0)
    for n in get_index(am, :bus, 0)

        if abs(ptdf_k[n]) >= am.setting["Simulation Setting"]["PTDF_threshold_value"]
        
            p_inj_ndhty = variable(am, nw, decomp_group, :p_inj_ndhty, (n,d,h,t,y))
                
            expr = JuMP.@expression(JuMP_model,  
                ptdf_k[n]*(p_inj_ndhty)
            )  

            JuMP.add_to_expression!(f_ptdf_kdhty, expr)

        end
    end

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($k,$d,$h,$t,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        f_ptdf_kdhty - variable(am, nw, decomp_group, :f_kdhty, (k,d,h,t,y))
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
    JuMP.set_name(constraint, string(const_name, "_($k,$d,$h,$t,$y)"))    
end

function constraint_sum_p_injdhty_RA(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, constant_load=0.0, system_peak_scale=1.0)

    total_p_inj_dhty = JuMP.AffExpr(0.0)

    for n in get_index(am, :bus, 0)
        JuMP.add_to_expression!(total_p_inj_dhty, 1,  variable(am, nw, decomp_group, :p_inj_ndhty, (n,d,h,t,y)))
    end

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($d,$h,$t,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        total_p_inj_dhty 
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
    JuMP.set_name(constraint, string(const_name, "_($d,$h,$t,$y)"))    
end

function constraint_post_cont_es_soc_balance_perfect_foresignt_first_hour_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int, ref_ED_solution, daygroup_first_hour_flag, RA_setting; nw::Int=am.cnw, update::Bool=false, prior_state_solution=nothing, es_avail::Real=1.0)

    idx_string = string("(", i, ", ", d, ", ", h, ", ", t, ", ", y, ")")
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)
    CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")
    BATEFF = sqrt(parameter(am, bus_idx, :gen_bus, tech_idx, "BATEFF"))
    card_T = length(am.setting["run_T"])

    # Hybrid storage: include g_G_ES_idhty from ref_ED (LP omits it as variable).
    hybrid_type = am.ref[:nw][0][:gen_index][i]["hybrid_type"]
    is_hybrid_ES = hybrid_type == "ES"
    hyb_idx = is_hybrid_ES ? am.ref[:nw][0][:gen_index][i]["hybrid_main_gen_idx"] : 0

    # Rolling-horizon: read prior_soc from h-1 commit if present, else from ref.
    has_prior_state = (prior_state_solution !== nothing)

    prior_soc_idht = 0.0
    if daygroup_first_hour_flag == true  # this is the first hour of the day group

        idx_string = string("(", i, ", ", d, ", ", 1, ", ", t, ", ", y, ")")
        prior_g_idht = ref_ED_solution["dispatch"][idx_string]["g_idhty"]
        prior_chg_idht = ref_ED_solution["storage"][idx_string]["chg_idhty"]
        prior_g_G_ES_idht = is_hybrid_ES ? ref_ED_solution["hybrid"][string("(", hyb_idx, ", ", d, ", ", 1, ", ", t, ", ", y, ")")]["g_G_ES_idhty"] : 0.0
        prior_soc_idht = ref_ED_solution["storage"][idx_string]["soc_idhty"] - (BATEFF * (prior_chg_idht + prior_g_G_ES_idht) - (1/BATEFF) * prior_g_idht)

    else
        if (h==1) && (t==1) # hour is 1, but not the first hour of the day group, then link to the last hour of the previous day
            prior_idx_string = string("(", i, ", ", d-1, ", ", 24, ", ", 1, ", ", y, ")")
            if has_prior_state && haskey(prior_state_solution, "storage") && haskey(prior_state_solution["storage"], prior_idx_string)
                prior_soc_idht = prior_state_solution["storage"][prior_idx_string]["soc_idhty"]
            else
                prior_soc_idht = ref_ED_solution["storage"][prior_idx_string]["soc_idhty"]
            end
        elseif (h != 1) && (t == 1) # not the first hour of the day group, and not the first hour, then link to the previous hour
            prior_idx_string = string("(", i, ", ", d, ", ", h-1, ", ", t, ", ", y, ")")
            if has_prior_state && haskey(prior_state_solution, "storage") && haskey(prior_state_solution["storage"], prior_idx_string)
                prior_soc_idht = prior_state_solution["storage"][prior_idx_string]["soc_idhty"]
            else
                prior_soc_idht = ref_ED_solution["storage"][prior_idx_string]["soc_idhty"]
            end
        end
    end

    # Adjust prior_soc_idht using the predefined scaling factor
    if haskey(RA_setting, "prior_soc_adjust")
        prior_soc_idht = min(prior_soc_idht * RA_setting["prior_soc_adjust"], parameter(am, bus_idx, :gen_bus, tech_idx, "ES_MWh"))
    end

    soc_idhty = variable(am, nw, decomp_group, :soc_idhty, (i,d,h,t,y))
    chg_idhty = variable(am, nw, decomp_group, :chg_idhty, (i,d,h,t,y))
    g_idhty = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))

    # Gate onsite renewable->battery charging by the asset's availability: when the hybrid is on
    # outage (es_avail=0) the injection is zero so SOC holds constant instead of overflowing its cap.
    ref_g_G_ES_idhty_now = (is_hybrid_ES ? ref_ED_solution["hybrid"][string("(", hyb_idx, ", ", d, ", ", h, ", ", t, ", ", y, ")")]["g_G_ES_idhty"] : 0.0) * es_avail


    expr = JuMP.@expression(JuMP_model,
        soc_idhty - prior_soc_idht - (BATEFF * (chg_idhty + ref_g_G_ES_idhty_now) - (1/BATEFF) * g_idhty)/card_T
        )
    constraint = JuMP.@constraint(JuMP_model,
        expr == 0
        )
    JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)"))

end

function constraint_post_cont_es_soc_balance_perfect_foresignt_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int, ref_ED_solution, daygroup_first_hour_flag; nw::Int=am.cnw, update::Bool=false, es_avail::Real=1.0)

    idx_string = string("(", i, ", ", d, ", ", h, ", ", t, ", ", y, ")")
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)
    CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")
    BATEFF = sqrt(parameter(am, bus_idx, :gen_bus, tech_idx, "BATEFF"))
    card_T = length(am.setting["run_T"])

    # Hybrid storage: include g_G_ES_idhty from ref_ED (LP omits it as variable).
    hybrid_type = am.ref[:nw][0][:gen_index][i]["hybrid_type"]
    is_hybrid_ES = hybrid_type == "ES"
    hyb_idx = is_hybrid_ES ? am.ref[:nw][0][:gen_index][i]["hybrid_main_gen_idx"] : 0

    soc_idhty = variable(am, nw, decomp_group, :soc_idhty, (i,d,h,t,y))
    chg_idhty = variable(am, nw, decomp_group, :chg_idhty, (i,d,h,t,y))
    if (h==1) && (t==1)
        if daygroup_first_hour_flag == true
            prior_soc_idhty = ref_ED_solution["storage"][string("(", i, ", ", d, ", ", 1, ", ", t, ", ", y, ")")]["soc_idhty"]
        else
            prior_soc_idhty = variable(am, nw, decomp_group, :soc_idhty, (i,d-1,24,1,y))
        end
    else
        prior_soc_idhty = variable(am, nw, decomp_group, :soc_idhty, (i,d,h-1,t,y))
    end
    g_idhty = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))

    # Gate onsite renewable->battery charging by the asset's availability (see first-hour builder).
    ref_g_G_ES_idhty_now = (is_hybrid_ES ? ref_ED_solution["hybrid"][string("(", hyb_idx, ", ", d, ", ", h, ", ", t, ", ", y, ")")]["g_G_ES_idhty"] : 0.0) * es_avail

    expr = JuMP.@expression(JuMP_model,
        soc_idhty - prior_soc_idhty - (BATEFF * (chg_idhty + ref_g_G_ES_idhty_now) - (1/BATEFF) * g_idhty)/card_T
        )
    constraint = JuMP.@constraint(JuMP_model,
        expr == 0
        )
    JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)"))

end


function constraint_cont_deploy_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int, ref_ED_solution, event_first_hour_flag; nw::Int=am.cnw, update::Bool=false, prior_state_solution=nothing)

    idx_string = string("(", i, ", ", d, ", ", h, ", ", t, ", ", y, ")")
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)
    Ramp30_Rate = min(1.0, 30 * parameter(am, bus_idx, :gen_bus, tech_idx, "Ramp"))   # 30 mins ramp rate
    Ramp60_Rate = min(1.0, 60 * parameter(am, bus_idx, :gen_bus, tech_idx, "Ramp"))   # 60 mins ramp rate
    CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")
    unit_group = parameter(am, 0, :gen_index, "UNIT_GROUP", i)

    u_G_iy = parameter(am, bus_idx, :gen_bus, tech_idx, "EXUNITS")

    # VRE shape (available shapes: wind_ons, wind_ofs, pv, rtpv, hydro, csp)
    shape = 1.0
    profile_type = parameter(am, bus_idx, :gen_bus, tech_idx, "Profile_Type")
    type_key = profile_type == "NA" ? "NA" : profile_type * "_shape"

    if type_key != "NA"
        shape = get_vre_zdt_shape(am.ref[:nw][0], y, bus_idx, d, h, type_key)
    end

    flexibility_percent = 0.0
    if (am.setting["Simulation Configuration"]["Hydro_Flexibility_Flag"] == true) && (parameter(am, bus_idx, :gen_bus, tech_idx, "Hydro_Flag") in ("ROR", "IMPOUNDMENT"))
        flexibility_percent = am.setting["Simulation Configuration"]["Hydro_Flexibility_Percent"]
    end

    # Rolling-horizon: anchor ramp on committed prior (60-min), not ref at h.
    use_committed_prior = false
    committed_prior_idx_string = ""
    if event_first_hour_flag && prior_state_solution !== nothing && haskey(prior_state_solution, "dispatch")
        if (h == 1) && (t == 1) && (d > 1)
            committed_prior_idx_string = string("(", i, ", ", d-1, ", ", 24, ", ", 1, ", ", y, ")")
        elseif (h != 1) && (t == 1)
            committed_prior_idx_string = string("(", i, ", ", d, ", ", h-1, ", ", t, ", ", y, ")")
        end
        if committed_prior_idx_string != "" && haskey(prior_state_solution["dispatch"], committed_prior_idx_string)
            use_committed_prior = true
        end
    end

    # 30-min ramp for the true first redispatch hour (mid-hour contingency); 60-min
    # for subsequent hours or once a committed prior exists (full hour since commit).
    ramp_limit = Ramp30_Rate * CAP * u_G_iy
    inter_hour_ramp_limit = min(Ramp60_Rate * CAP * u_G_iy, CAP * u_G_iy)
    if event_first_hour_flag == false || use_committed_prior
        ramp_limit = inter_hour_ramp_limit
    end

    ref_g_idht = ref_ED_solution["dispatch"][idx_string]["g_idhty"]

    g_idhty = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))

    prior_g_idht = 0
    if event_first_hour_flag == false
        # Cross-day boundary: at (d>1, h=1, t=1) the prior hour is (d-1, h=24, t=1)
        if (h==1) && (t==1) && (d > 1)
            prior_g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d-1,24,1,y))
        elseif (h != 1) && (t == 1)
            prior_g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h-1,t,y))
        end
    elseif use_committed_prior
        # Anchor ramp on the committed dispatch at h-1, not on ref at h.
        prior_g_idht = prior_state_solution["dispatch"][committed_prior_idx_string]["g_idhty"]
    end


    expr = JuMP.@expression(JuMP_model,
        (g_idhty) - CAP*u_G_iy*shape*(1 + flexibility_percent)
        )
    constraint = JuMP.@constraint(JuMP_model,
        expr <= 0
        )
    JuMP.set_name(constraint, string(const_name, "_gen_upper_bound_($i,$d,$h,$t,$y)"))

    if event_first_hour_flag == true && !use_committed_prior
        expr = JuMP.@expression(JuMP_model, (g_idhty - ref_g_idht) - ramp_limit)
        constraint = JuMP.@constraint(JuMP_model, expr <= 0)
        JuMP.set_name(constraint, string(const_name, "_first_hour_up_ramp_($i,$d,$h,$t,$y)"))

        expr = JuMP.@expression(JuMP_model, (g_idhty - ref_g_idht) + ramp_limit)
        constraint = JuMP.@constraint(JuMP_model, expr >= 0)
        JuMP.set_name(constraint, string(const_name, "_first_hour_down_ramp_($i,$d,$h,$t,$y)"))
    else

        expr = JuMP.@expression(JuMP_model,  
            (g_idhty - prior_g_idht) - inter_hour_ramp_limit
            )   
        constraint = JuMP.@constraint(JuMP_model, 
            expr <= 0
            )    
        JuMP.set_name(constraint, string(const_name, "_inter_hour_up_ramp_($i,$d,$h,$t,$y)"))    

        expr = JuMP.@expression(JuMP_model,  
            (g_idhty - prior_g_idht) + inter_hour_ramp_limit
            )   
        constraint = JuMP.@constraint(JuMP_model, 
            expr >= 0
            )    
        JuMP.set_name(constraint, string(const_name, "_inter_hour_down_ramp_($i,$d,$h,$t,$y)"))    
    end
    
end

function constraint_redispatch_abs_dev_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int, ref_ED_solution; nw::Int=am.cnw, update::Bool=false)

    idx_string = string("(", i, ", ", d, ", ", h, ", ", t, ", ", y, ")")
    ref_g_idht = ref_ED_solution["dispatch"][idx_string]["g_idhty"]
    g_idhty = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
    redispatch_abs_dev_idhty = variable(am, nw, decomp_group, :redispatch_abs_dev_idhty, (i,d,h,t,y))

    if update
        JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_abs_dev_pos_($i,$d,$h,$t,$y)")))
        JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_abs_dev_neg_($i,$d,$h,$t,$y)")))
    end

    expr = JuMP.@expression(JuMP_model, redispatch_abs_dev_idhty - (g_idhty - ref_g_idht))
    constraint = JuMP.@constraint(JuMP_model, expr >= 0)
    JuMP.set_name(constraint, string(const_name, "_abs_dev_pos_($i,$d,$h,$t,$y)"))

    expr = JuMP.@expression(JuMP_model, redispatch_abs_dev_idhty + (g_idhty - ref_g_idht))
    constraint = JuMP.@constraint(JuMP_model, expr >= 0)
    JuMP.set_name(constraint, string(const_name, "_abs_dev_neg_($i,$d,$h,$t,$y)"))

end

function constraint_cont_deploy_snapshot_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int, ref_ED_solution, risk_ED_solution, event_first_hour_flag; nw::Int=am.cnw, update::Bool=false)
    
    idx_string = string("(", i, ", ", d, ", ", h, ", ", t, ", ", y, ")")
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
    Ramp30_Rate = min(1.0, 30 * parameter(am, bus_idx, :gen_bus, tech_idx, "Ramp"))   # 30 mins ramp rate
    Ramp60_Rate = min(1.0, 60 * parameter(am, bus_idx, :gen_bus, tech_idx, "Ramp"))   # 60 mins ramp rate
    CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")
    unit_group = parameter(am, 0, :gen_index, "UNIT_GROUP", i)

    u_G_iy = parameter(am, bus_idx, :gen_bus, tech_idx, "EXUNITS")

    # VRE shape (available shapes: wind_ons, wind_ofs, pv, rtpv, hydro, csp)
    shape = 1.0
    profile_type = parameter(am, bus_idx, :gen_bus, tech_idx, "Profile_Type")
    type_key = profile_type == "NA" ? "NA" : profile_type * "_shape"

    if type_key != "NA"
        shape = get_vre_zdt_shape(am.ref[:nw][0], y, bus_idx, d, h, type_key)
    end

    flexibility_percent = 0.0
    if (am.setting["Simulation Configuration"]["Hydro_Flexibility_Flag"] == true) && (parameter(am, bus_idx, :gen_bus, tech_idx, "Hydro_Flag") in ("ROR", "IMPOUNDMENT"))
        flexibility_percent = am.setting["Simulation Configuration"]["Hydro_Flexibility_Percent"]
    end

    # Define ramp limit (first hour = 30-min ramp, subsequent hours = 60-min ramp)
    ramp_limit = Ramp30_Rate * CAP * u_G_iy
    if event_first_hour_flag == false
        ramp_limit = min(Ramp60_Rate * CAP * u_G_iy, CAP * u_G_iy)
        inter_hour_ramp_limit = min(Ramp60_Rate * CAP * u_G_iy, CAP * u_G_iy)
    end

    g_idhty = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))

    # Reference and prior solved generation for ramp linkage.
    ref_g_idht = ref_ED_solution["dispatch"][idx_string]["g_idhty"]
    prior_g_idht = 0.0
    if event_first_hour_flag == false
        if (h==1) && (t==1)
            prior_idx_string = string("(", i, ", ", d-1, ", ", 24, ", ", 1, ", ", y, ")")
            if haskey(risk_ED_solution, "dispatch") && haskey(risk_ED_solution["dispatch"], prior_idx_string)
                prior_g_idht = risk_ED_solution["dispatch"][prior_idx_string]["g_idhty"]
            else
                prior_g_idht = ref_ED_solution["dispatch"][prior_idx_string]["g_idhty"]
            end
        elseif (h != 1) && (t == 1)
            prior_idx_string = string("(", i, ", ", d, ", ", h-1, ", ", t, ", ", y, ")")
            if haskey(risk_ED_solution, "dispatch") && haskey(risk_ED_solution["dispatch"], prior_idx_string)
                prior_g_idht = risk_ED_solution["dispatch"][prior_idx_string]["g_idhty"]
            else
                prior_g_idht = ref_ED_solution["dispatch"][prior_idx_string]["g_idhty"]
            end
        end
    end

    expr = JuMP.@expression(JuMP_model, 
        g_idhty - CAP*u_G_iy*shape*(1 + flexibility_percent)
        )
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )
    JuMP.set_name(constraint, string(const_name, "_gen_upper_bound_($i,$d,$h,$t,$y)"))

    if event_first_hour_flag == true
        expr = JuMP.@expression(JuMP_model, 
            (g_idhty - ref_g_idht) - ramp_limit
            )
        constraint = JuMP.@constraint(JuMP_model, 
            expr <= 0
            )
        JuMP.set_name(constraint, string(const_name, "_first_hour_up_ramp_($i,$d,$h,$t,$y)"))

        expr = JuMP.@expression(JuMP_model, 
            (g_idhty - ref_g_idht) + ramp_limit
            )
        constraint = JuMP.@constraint(JuMP_model, 
            expr >= 0
            )
        JuMP.set_name(constraint, string(const_name, "_first_hour_down_ramp_($i,$d,$h,$t,$y)"))
    else
        expr = JuMP.@expression(JuMP_model, 
            (g_idhty - prior_g_idht) - inter_hour_ramp_limit
            )
        constraint = JuMP.@constraint(JuMP_model, 
            expr <= 0
            )
        JuMP.set_name(constraint, string(const_name, "_inter_hour_up_ramp_($i,$d,$h,$t,$y)"))

        expr = JuMP.@expression(JuMP_model, 
            (g_idhty - prior_g_idht) + inter_hour_ramp_limit
            )
        constraint = JuMP.@constraint(JuMP_model, 
            expr >= 0
            )
        JuMP.set_name(constraint, string(const_name, "_inter_hour_down_ramp_($i,$d,$h,$t,$y)"))
    end

end

function constraint_cont_deploy_ES_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int, ref_ED_solution, event_first_hour_flag, RA_setting; nw::Int=am.cnw, update::Bool=false)

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)
    CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")
    u_G_iy = parameter(am, bus_idx, :gen_bus, tech_idx, "EXUNITS")

    idx_string = string("(", i, ", ", d, ", ", h, ", ", t, ", ", y, ")")
    chg_idhty = variable(am, nw, decomp_group, :chg_idhty, (i,d,h,t,y))
    g_idhty   = variable(am, nw, decomp_group, :g_idhty,   (i,d,h,t,y))

    ref_chg_idht = ref_ED_solution["storage"][idx_string]["chg_idhty"]
    ref_g_idht   = ref_ED_solution["dispatch"][idx_string]["g_idhty"]

    # Charging gating
    if RA_setting["Post_Contingency_Charge_method"] == "Not Allowed"

        expr = JuMP.@expression(JuMP_model,
            chg_idhty
            )
        constraint = JuMP.@constraint(JuMP_model,
            expr == 0
            )
        JuMP.set_name(constraint, string(const_name, "_sto_post_charge_up_bound_($i,$d,$h,$t,$y)"))

    elseif RA_setting["Post_Contingency_Charge_method"] == "Same as Reference"
        expr = JuMP.@expression(JuMP_model,
            (chg_idhty - ref_chg_idht)
            )
        constraint = JuMP.@constraint(JuMP_model,
            expr == 0
            )
        JuMP.set_name(constraint, string(const_name, "_sto_post_charge_up_bound_($i,$d,$h,$t,$y)"))

    else # "Optimal" or unrecognised value — allow up to nameplate
        expr = JuMP.@expression(JuMP_model,
            (chg_idhty - u_G_iy * CAP)
            )
        constraint = JuMP.@constraint(JuMP_model,
            expr <= 0
            )
        JuMP.set_name(constraint, string(const_name, "_sto_post_charge_up_bound_($i,$d,$h,$t,$y)"))
    end

    # Discharge gating — symmetric to charge gating above
    if RA_setting["Post_Contingency_Discharge_method"] == "Not Allowed"
        expr = JuMP.@expression(JuMP_model, g_idhty)
        constraint = JuMP.@constraint(JuMP_model, expr == 0)
        JuMP.set_name(constraint, string(const_name, "_sto_post_discharge_up_bound_($i,$d,$h,$t,$y)"))

    elseif RA_setting["Post_Contingency_Discharge_method"] == "Same as Reference"
        expr = JuMP.@expression(JuMP_model, g_idhty - ref_g_idht)
        constraint = JuMP.@constraint(JuMP_model, expr == 0)
        JuMP.set_name(constraint, string(const_name, "_sto_post_discharge_up_bound_($i,$d,$h,$t,$y)"))

    # else "Optimal": no additional bound beyond existing capacity constraints
    end
end


function constraint_cont_deploy_ES_limit_snapshot_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int, risk_flag, ref_ED_solution; nw::Int=am.cnw, update::Bool=false)

    unit_category = parameter(am, 0, :gen_index, "UNIT_CATEGORY", i)
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)   
    big_M = 100000

    if unit_category in ["STORAGE"] && (risk_flag == 1) # if unit is online

        # sum of ENS over all buses; only needed inside this branch, so build it here
        sum_ens = JuMP.AffExpr(0.0)
        for n in get_index(am, :bus, 0)
            JuMP.add_to_expression!(sum_ens, 1, variable(am, nw, decomp_group, :ens_ndhty, (n,d,h,t,y)))
        end

        idx_string = string("(", i, ", ", d, ", ", h, ", ", t, ", ", y, ")")
        g_idhty = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
        chg_idhty = variable(am, nw, decomp_group, :chg_idhty, (i,d,h,t,y))
        ref_g_idht = ref_ED_solution["dispatch"][idx_string]["g_idhty"]
        ref_chg_idht = ref_ED_solution["storage"][idx_string]["chg_idhty"]

        # No ENS => sum_ens = 0 => force storage dispatch/charge to reference values.
        expr = JuMP.@expression(JuMP_model, (g_idhty - ref_g_idht) - big_M * sum_ens)
        constraint = JuMP.@constraint(JuMP_model, expr <= 0)
        JuMP.set_name(constraint, string(const_name, "_sto_post_discharge_dev_up_($i,$d,$h,$t,$y)"))

        expr = JuMP.@expression(JuMP_model, -(g_idhty - ref_g_idht) - big_M * sum_ens)
        constraint = JuMP.@constraint(JuMP_model, expr <= 0)
        JuMP.set_name(constraint, string(const_name, "_sto_post_discharge_dev_dn_($i,$d,$h,$t,$y)"))

        expr = JuMP.@expression(JuMP_model, (chg_idhty - ref_chg_idht) - big_M * sum_ens)
        constraint = JuMP.@constraint(JuMP_model, expr <= 0)
        JuMP.set_name(constraint, string(const_name, "_sto_post_charge_dev_up_($i,$d,$h,$t,$y)"))

        expr = JuMP.@expression(JuMP_model, -(chg_idhty - ref_chg_idht) - big_M * sum_ens)
        constraint = JuMP.@constraint(JuMP_model, expr <= 0)
        JuMP.set_name(constraint, string(const_name, "_sto_post_charge_dev_dn_($i,$d,$h,$t,$y)"))
    end
   
end

function constraint_cont_deploy_ES_snapshot_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int, ref_ED_solution, risk_ED_solution, RA_setting, event_first_hour_flag; nw::Int=am.cnw, update::Bool=false)

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)
    CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")
    idx_string = string("(", i, ", ", d, ", ", h, ", ", t, ", ", y, ")")
    u_G_iy = parameter(am, bus_idx, :gen_bus, tech_idx, "EXUNITS")

    chg_idhty = variable(am, nw, decomp_group, :chg_idhty, (i,d,h,t,y))
    g_idhty   = variable(am, nw, decomp_group, :g_idhty,   (i,d,h,t,y))
    ref_chg_idht = ref_ED_solution["storage"][idx_string]["chg_idhty"]
    ref_g_idht   = ref_ED_solution["dispatch"][idx_string]["g_idhty"]

    # Charging gating
    if RA_setting["Post_Contingency_Charge_method"] == "Not Allowed"

        expr = JuMP.@expression(JuMP_model, chg_idhty)
        constraint = JuMP.@constraint(JuMP_model, expr == 0)
        JuMP.set_name(constraint, string(const_name, "_sto_post_charge_up_bound_($i,$d,$h,$t,$y)"))

    elseif RA_setting["Post_Contingency_Charge_method"] == "Same as Reference"
        expr = JuMP.@expression(JuMP_model, (chg_idhty - ref_chg_idht))
        constraint = JuMP.@constraint(JuMP_model, expr == 0)
        JuMP.set_name(constraint, string(const_name, "_sto_post_charge_up_bound_($i,$d,$h,$t,$y)"))

    else # "Optimal" or unrecognised value — allow up to nameplate
        expr = JuMP.@expression(JuMP_model, (chg_idhty - u_G_iy * CAP))
        constraint = JuMP.@constraint(JuMP_model, expr <= 0)
        JuMP.set_name(constraint, string(const_name, "_sto_post_charge_up_bound_($i,$d,$h,$t,$y)"))
    end

    # Discharge gating — symmetric to charge gating above
    if RA_setting["Post_Contingency_Discharge_method"] == "Not Allowed"
        expr = JuMP.@expression(JuMP_model, g_idhty)
        constraint = JuMP.@constraint(JuMP_model, expr == 0)
        JuMP.set_name(constraint, string(const_name, "_sto_post_discharge_up_bound_($i,$d,$h,$t,$y)"))

    elseif RA_setting["Post_Contingency_Discharge_method"] == "Same as Reference"
        expr = JuMP.@expression(JuMP_model, g_idhty - ref_g_idht)
        constraint = JuMP.@constraint(JuMP_model, expr == 0)
        JuMP.set_name(constraint, string(const_name, "_sto_post_discharge_up_bound_($i,$d,$h,$t,$y)"))

    # else "Optimal": no additional bound beyond existing capacity constraints
    end

end

function constraint_total_cont_deploy_snapshot_dhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, ids_n, ids_i, ids_sto, d::Int, h::Int, t::Int, y::Int, gen_risk_matrix, event_idx, ref_ED_solution; nw::Int=am.cnw, update::Bool=false)

    total_gen_outage = 0.0
    for i in ids_i
        if gen_risk_matrix[event_idx, i] == 0
            total_gen_outage += ref_ED_solution["dispatch"][string("(", i, ", ", d, ", ", h, ", ", t, ", ", y, ")")]["g_idhty"] 
            if i in ids_sto
                total_gen_outage -= ref_ED_solution["storage"][string("(", i, ", ", d, ", ", h, ", ", t, ", ", y, ")")]["chg_idhty"] 
            end
        end
    end

    total_prior_gen = 0.0
    for i in ids_i
        if gen_risk_matrix[event_idx, i] == 1
            total_prior_gen += ref_ED_solution["dispatch"][string("(", i, ", ", d, ", ", h, ", ", t, ", ", y, ")")]["g_idhty"] 
            if i in ids_sto
                total_prior_gen -= ref_ED_solution["storage"][string("(", i, ", ", d, ", ", h, ", ", t, ", ", y, ")")]["chg_idhty"] 
            end
        end
    end

    sum_redispatch_gen = JuMP.AffExpr(0.0)    
    for i in ids_i
        if gen_risk_matrix[event_idx, i] == 1
            JuMP.add_to_expression!(sum_redispatch_gen, 1.0, variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y)))

            if i in ids_sto
                JuMP.add_to_expression!(sum_redispatch_gen, -1, variable(am, nw, decomp_group, :chg_idhty, (i,d,h,t,y)))
            end
        end
    end

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($d,$h,$t,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        (sum_redispatch_gen - total_prior_gen) - total_gen_outage
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    JuMP.set_name(constraint, string(const_name, "_($d,$h,$t,$y)"))    
    
end

function constraint_cont_deploy_dhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, ids_n, ids_i, ids_sto, d::Int, h::Int, t::Int, y::Int, gen_risk_matrix, risk_hour, ref_ED_solution; nw::Int=am.cnw, update::Bool=false)

    ids_sto_set = Set(ids_sto)
    total_gen_outage = 0.0 # total gen outage MW
    total_prior_gen = 0.0 # total generations of available units in the reference solution (i.e., before redispatch)

    for i in ids_i
        string_idx = string("(", i, ", ", d, ", ", h, ", ", t, ", ", y, ")")
        if gen_risk_matrix[risk_hour, i] == 0
            total_gen_outage += ref_ED_solution["dispatch"][string_idx]["g_idhty"] 
            if i in ids_sto_set
                total_gen_outage -= ref_ED_solution["storage"][string_idx]["chg_idhty"] 
            end
        else
            total_prior_gen += ref_ED_solution["dispatch"][string_idx]["g_idhty"] 
            if i in ids_sto_set
                total_prior_gen -= ref_ED_solution["storage"][string_idx]["chg_idhty"] 
            end
        end
    end

    sum_redispatch_gen = JuMP.AffExpr(0.0)    
    for i in ids_i
        if gen_risk_matrix[risk_hour, i] == 1
            JuMP.add_to_expression!(sum_redispatch_gen, 1, variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y)))

            if i in ids_sto_set
                JuMP.add_to_expression!(sum_redispatch_gen, -1, variable(am, nw, decomp_group, :chg_idhty, (i,d,h,t,y)))
            end

        end
    end

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($d,$h,$t,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        (sum_redispatch_gen - total_prior_gen) - total_gen_outage
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    JuMP.set_name(constraint, string(const_name, "_($d,$h,$t,$y)"))    
    
end

function constraint_ES_SOC_Balance_Inter_Hour_NoReserve_Cont_Snapshot_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int, ref_ED_solution, risk_ED_solution, daygroup_first_hour_flag, RA_setting; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=true)

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)
    BATEFF = sqrt(parameter(am, bus_idx, :gen_bus, tech_idx, "BATEFF"))
    card_T = length(am.setting["run_T"])

    # Hybrid storage: include g_G_ES_idhty from ref_ED (snapshot LP omits it as variable).
    hybrid_type = am.ref[:nw][0][:gen_index][i]["hybrid_type"]
    is_hybrid_ES = hybrid_type == "ES"
    hyb_idx = is_hybrid_ES ? am.ref[:nw][0][:gen_index][i]["hybrid_main_gen_idx"] : 0

    prior_soc_idht = 0.0
    if daygroup_first_hour_flag == true  # this is the first hour of the day group

        idx_string = string("(", i, ", ", d, ", ", 1, ", ", t, ", ", y, ")")
        prior_g_idht = ref_ED_solution["dispatch"][idx_string]["g_idhty"]
        prior_chg_idht = ref_ED_solution["storage"][idx_string]["chg_idhty"]
        prior_g_G_ES_idht = is_hybrid_ES ? ref_ED_solution["hybrid"][string("(", hyb_idx, ", ", d, ", ", 1, ", ", t, ", ", y, ")")]["g_G_ES_idhty"] : 0.0
        prior_soc_idht = ref_ED_solution["storage"][idx_string]["soc_idhty"] - (BATEFF * (prior_chg_idht + prior_g_G_ES_idht) - (1/BATEFF) * prior_g_idht)

        # Adjust prior_soc_idht using the predefined scaling factor
        if haskey(RA_setting, "prior_soc_adjust")
            prior_soc_idht = min(prior_soc_idht * RA_setting["prior_soc_adjust"], parameter(am, bus_idx, :gen_bus, tech_idx, "ES_MWh"))
        end

    else
        if (h==1) && (t==1) # hour is 1, but not the first hour of the day group, then link to the last hour of the previous day
            prior_idx_string = string("(", i, ", ", d-1, ", ", 24, ", ", 1, ", ", y, ")")
            if haskey(risk_ED_solution, "storage") && haskey(risk_ED_solution["storage"], prior_idx_string)
                prior_soc_idht = risk_ED_solution["storage"][prior_idx_string]["soc_idhty"]
            else
                prior_soc_idht = ref_ED_solution["storage"][prior_idx_string]["soc_idhty"]
            end
        elseif (h != 1) && (t == 1) # not the first hour of the day group, and not the first hour, then link to the previous hour
            prior_idx_string = string("(", i, ", ", d, ", ", h-1, ", ", t, ", ", y, ")")
            if haskey(risk_ED_solution, "storage") && haskey(risk_ED_solution["storage"], prior_idx_string)
                prior_soc_idht = risk_ED_solution["storage"][prior_idx_string]["soc_idhty"]
            else
                prior_soc_idht = ref_ED_solution["storage"][prior_idx_string]["soc_idhty"]
            end
        end
    end

    soc_idhty = variable(am, nw, decomp_group, :soc_idhty, (i,d,h,t,y))
    chg_idhty = variable(am, nw, decomp_group, :chg_idhty, (i,d,h,t,y))
    g_idhty = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))

    ref_g_G_ES_idhty_now = is_hybrid_ES ? ref_ED_solution["hybrid"][string("(", hyb_idx, ", ", d, ", ", h, ", ", t, ", ", y, ")")]["g_G_ES_idhty"] : 0.0

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end
    expr = JuMP.@expression(JuMP_model,
        soc_idhty - prior_soc_idht - (BATEFF * (chg_idhty + ref_g_G_ES_idhty_now) - (1/BATEFF) * g_idhty)/card_T
        )
    constraint = JuMP.@constraint(JuMP_model,
        expr == 0
        )
    JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)"))

end

function constraint_VREBalance_Budget_NonReserve_RA_snapshot_iy(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, day_group_id::Int, y::Int, ids_d, ids_h, ids_t, simulated_event_hours, risk_ED_solution, type::String; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=true)

    
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
    CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")
    Timeseries_Tag = parameter(am, bus_idx, :gen_bus, tech_idx, "Timeseries_Tag")

    # water budget of day_group_id in bus_idx (MWh/MW)
    budget = am.ref[:nw][0][:hydro_budget][bus_idx][day_group_id]

    sum_generation = JuMP.AffExpr(0.0)
    for d in ids_d
        for h in ids_h
            for t in ids_t
                JuMP.add_to_expression!(sum_generation, 1.0, variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y)))
            end
        end
    end
    u_i = parameter(am, bus_idx, :gen_bus, tech_idx, "EXUNITS")

    prior_generation = 0.0
    for simulated_hours_set in simulated_event_hours
        event_day = simulated_hours_set[1]
        event_hour = simulated_hours_set[2]
        event_time = simulated_hours_set[3]
        idx_string = string("(", i, ", ", event_day, ", ", event_hour, ", ", event_time, ", ", y, ")")
        if haskey(risk_ED_solution, "dispatch") && haskey(risk_ED_solution["dispatch"], idx_string)
            prior_generation += risk_ED_solution["dispatch"][idx_string]["g_idhty"]
        else
            prior_generation += ref_ED_solution["dispatch"][idx_string]["g_idhty"]
        end

    end

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        sum_generation + prior_generation - budget * u_i * CAP
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    JuMP.set_name(constraint, string(const_name, "_($i,$y)"))   

end

function constraint_ES_SOC_Min_NonReserve_idhty_RA(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=true)

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
    STOMIN = parameter(am, bus_idx, :gen_bus, tech_idx, "STOMIN")
    BATEFF = sqrt(parameter(am, bus_idx, :gen_bus, tech_idx, "BATEFF"))

    soc_idhty = variable(am, nw, decomp_group, :soc_idhty, (i,d,h,t,y))
    
    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        soc_idhty - STOMIN
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr >= 0
        )    
    JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)"))

end


function constraint_TherMaxDispatch_NoReserve_RA_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=true)

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)    

    CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")
    
    g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))

    u_i = parameter(am, bus_idx, :gen_bus, tech_idx, "EXUNITS")


    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        g_idht - u_i * CAP
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)"))  end   
end

function constraint_RampUp_InterTemporal_Hour_idhy_RA(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int, day_group; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=true, apply_PMAX::Bool=true)

    HFREQ = am.setting["Simulation Configuration"]["HFREQ"]

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
    CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")
    PMAX = parameter(am, bus_idx, :gen_bus, tech_idx, "PMAX")   # %
    if apply_PMAX == false
        PMAX = 1.0
    end
    Ramp_Rate = min(1.0, 60 * parameter(am, bus_idx, :gen_bus, tech_idx, "Ramp"))   # hourly Ramp Rate (%)

    u_G_iy = parameter(am, bus_idx, :gen_bus, tech_idx, "EXUNITS")

    if (h == 1)
        if am.ref[:nw][0][:repdays][d]["Day"] != am.ref[:nw][0][:repday_groups][day_group]["Start_Day_Id"]
            last_hour = last(am.setting["run_H"])
            
            g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
            prior_g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d-1,last_hour,t,y))
            
            if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$y)"))) end    
            expr = JuMP.@expression(JuMP_model,  
                g_idht - prior_g_idht - (u_G_iy * Ramp_Rate * CAP * PMAX)
                )   
            constraint = JuMP.@constraint(JuMP_model, 
                expr <= 0
                )    
            if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$y)")) end    
        
        elseif am.ref[:nw][0][:repdays][d]["Day"] == am.ref[:nw][0][:repday_groups][day_group]["Start_Day_Id"]
            num_days_in_group = length(am.ref[:nw][0][:repday_groups][day_group]["Day_Idx_List"])
            end_day_id = am.ref[:nw][0][:repday_groups][day_group]["Day_Idx_List"][num_days_in_group]
            last_hour = last(am.setting["run_H"])
            
            g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
            prior_g_idht = variable(am, nw, decomp_group, :g_idhty, (i,end_day_id,last_hour,t,y))
            
            if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$y)"))) end    
            expr = JuMP.@expression(JuMP_model,  
                g_idht - prior_g_idht - (u_G_iy * Ramp_Rate * CAP * PMAX)
                )   
            constraint = JuMP.@constraint(JuMP_model, 
                expr <= 0
                )    
            if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$y)")) end    
        
        end
    elseif (h != 1) 
        HFREQ = am.setting["Simulation Configuration"]["HFREQ"]
        
        g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
        prior_g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h-HFREQ,t,y))
        
        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$y)"))) end    
        expr = JuMP.@expression(JuMP_model,  
            g_idht - prior_g_idht - (u_G_iy * Ramp_Rate * CAP * PMAX)
            )   
        constraint = JuMP.@constraint(JuMP_model, 
            expr <= 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$y)")) end    

    end
end

function constraint_RampDn_InterTemporal_Hour_idhy_RA(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int, day_group; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=true, apply_PMAX::Bool=true)

    HFREQ = am.setting["Simulation Configuration"]["HFREQ"]

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
    CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")
    PMAX = parameter(am, bus_idx, :gen_bus, tech_idx, "PMAX")   # %
    if apply_PMAX == false
        PMAX = 1.0
    end
    Ramp_Rate = min(1.0, 60 * parameter(am, bus_idx, :gen_bus, tech_idx, "Ramp"))   # hourly Ramp Rate (%)

    u_G_iy = parameter(am, bus_idx, :gen_bus, tech_idx, "EXUNITS")

    if (h == 1)
        if am.ref[:nw][0][:repdays][d]["Day"] != am.ref[:nw][0][:repday_groups][day_group]["Start_Day_Id"]
            last_hour = last(am.setting["run_H"])
            
            g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
            prior_g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d-1,last_hour,t,y))
            
            if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$y)"))) end    
            expr = JuMP.@expression(JuMP_model,  
                - g_idht + prior_g_idht - (u_G_iy * Ramp_Rate * CAP * PMAX)
                )   
            constraint = JuMP.@constraint(JuMP_model, 
                expr <= 0
                )    
            if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$y)")) end    
        
        elseif am.ref[:nw][0][:repdays][d]["Day"] == am.ref[:nw][0][:repday_groups][day_group]["Start_Day_Id"]
            num_days_in_group = length(am.ref[:nw][0][:repday_groups][day_group]["Day_Idx_List"])
            end_day_id = am.ref[:nw][0][:repday_groups][day_group]["Day_Idx_List"][num_days_in_group]
            last_hour = last(am.setting["run_H"])
            
            g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
            prior_g_idht = variable(am, nw, decomp_group, :g_idhty, (i,end_day_id,last_hour,t,y))
            
            if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$y)"))) end    
            expr = JuMP.@expression(JuMP_model,  
                - g_idht + prior_g_idht - (u_G_iy * Ramp_Rate * CAP * PMAX)
                )   
            constraint = JuMP.@constraint(JuMP_model, 
                expr <= 0
                )    
            if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$y)")) end    
        
        end
    elseif (h != 1) 
        HFREQ = am.setting["Simulation Configuration"]["HFREQ"]
        
        g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
        prior_g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h-HFREQ,t,y))
        
        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$y)"))) end    
        expr = JuMP.@expression(JuMP_model,  
            - g_idht + prior_g_idht - (u_G_iy * Ramp_Rate * CAP * PMAX)
            )   
        constraint = JuMP.@constraint(JuMP_model, 
            expr <= 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$y)")) end    

    end
end

function constraint_ES_Sto_UC_Limit_idhty_RA(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=true)

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
    
    if parameter(am, bus_idx, :gen_bus, tech_idx, "Storage Commitment") == true

        u_i = parameter(am, bus_idx, :gen_bus, tech_idx, "EXUNITS")
        sto_c_idhty = variable(am, nw, decomp_group, :sto_c_idhty, (i,d,h,t,y))
        
        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
        expr = JuMP.@expression(JuMP_model,  
            sto_c_idhty - u_i 
            )   
        constraint = JuMP.@constraint(JuMP_model, 
            expr <= 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end    

    end

end

function constraint_ES_AET_OP_RA_y(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, y::Int, ids_d; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=true, reserve_flag::Bool=true)

    ids_dhty = [(d,h,t,y) for d in ids_d for h in am.setting["run_H"] for t in am.setting["run_T"]]

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
    
    AET = parameter(am, bus_idx, :gen_bus, tech_idx, "AET")
    BATEFF = sqrt(parameter(am, bus_idx, :gen_bus, tech_idx, "BATEFF"))

    total_num_days = 0.0
    for d in ids_d
        total_num_days += parameter(am, 0, :repdays, "NumDays", d)
    end
    AET_scale = total_num_days / 365    # this will be 1 if PH is not used.

    if AET > 0
    
        sum_AET = JuMP.AffExpr(0.0)
        for (d, h, t, y) in ids_dhty
            
        
            JuMP.add_to_expression!(sum_AET, parameter(am, 0, :repdays, "NumDays", d), variable(am, nw, decomp_group, :chg_idhty, (i,d,h,t,y)))
            JuMP.add_to_expression!(sum_AET, parameter(am, 0, :repdays, "NumDays", d), variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y)))
            
        end

        u_i = parameter(am, bus_idx, :gen_bus, tech_idx, "EXUNITS")

        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$y)"))) end    
        expr = JuMP.@expression(JuMP_model,  
            (1.0 / sqrt(2 * AET * AET_scale)) * sum_AET - sqrt(2 * AET * AET_scale) * u_i
            )   
        constraint = JuMP.@constraint(JuMP_model, 
            expr <= 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$y)")) end    

    end
end

function constraint_ES_SOC_Neutral_idy(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, y::Int, start_day_idx::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  

    last_subhour = last(am.setting["run_T"])
    last_hour = last(am.setting["run_H"])
    
    BATEFF = sqrt(parameter(am, bus_idx, :gen_bus, tech_idx, "BATEFF"))
    card_T = length(am.setting["run_T"])

    soc_end = variable(am, nw, decomp_group, :soc_idhty, (i,d,last_hour,last_subhour,y))
    soc_beginning = variable(am, nw, decomp_group, :soc_idhty, (i,start_day_idx,1,1,y))
    chg_beginning = variable(am, nw, decomp_group, :chg_idhty, (i,start_day_idx,1,1,y))
    g_beginning = variable(am, nw, decomp_group, :g_idhty, (i,start_day_idx,1,1,y))
    
    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        soc_end - soc_beginning + (BATEFF * chg_beginning - (1/BATEFF) * g_beginning)/card_T
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$y)")) end    

end

function constraint_ES_Charge_Max_NonReserve_idhty_RA(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=true)

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
    CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")
    
    chg_idht = variable(am, nw, decomp_group, :chg_idhty, (i,d,h,t,y))
    u_i = parameter(am, bus_idx, :gen_bus, tech_idx, "EXUNITS")
        
    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        chg_idht - u_i * CAP
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end    

end

function constraint_ES_SOC_Max_NonReserve_idhty_RA(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=true)

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)
    BATEFF = sqrt(sqrt(parameter(am, bus_idx, :gen_bus, tech_idx, "BATEFF")))

    soc_idhty = variable(am, nw, decomp_group, :soc_idhty, (i,d,h,t,y))
    u_ESE_iy = parameter(am, bus_idx, :gen_bus, tech_idx, "ES_MWh")

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end
    expr = JuMP.@expression(JuMP_model,
        soc_idhty - u_ESE_iy
        )
    constraint = JuMP.@constraint(JuMP_model,
        expr <= 0
        )
    JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)"))

end

function constraint_VREBalance_Fixed_Profile_NoReserve_idhty_RA(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int, type::String; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=true, apply_PMAX::Bool=true)

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
    CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")
    Timeseries_Tag = parameter(am, bus_idx, :gen_bus, tech_idx, "Timeseries_Tag")
    PMAX = parameter(am, bus_idx, :gen_bus, tech_idx, "PMAX")
    if apply_PMAX == false
        PMAX = 1.0
    end    
    bus_id = parameter(am, 0, :bus, "bus_i", bus_idx)

    operational_option = parameter(am, bus_idx, :gen_bus, tech_idx, "Dispatch")

    # VRE shape (available shapes: wind_ons, wind_ofs, pv, rtpv, hydro, csp)
    shape = 0.0
    profile_type = parameter(am, bus_idx, :gen_bus, tech_idx, "Profile_Type")
    if Timeseries_Tag == "LOCAL"
        type_key = profile_type * "_shape"
        shape = get_vre_zdt_shape(am.ref[:nw][0], y, bus_idx, d, h, type_key)
    else
        if profile_type in ("wind_ons", "wind_ofs", "csp")
            shape = parameter(am, 0, :planning_stages, "repdays", "data", profile_type, d, h, t, y)[Timeseries_Tag]
        end
    end

    g_idht = variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y))
    curt_idht = variable(am, nw, decomp_group, :curt_idhty, (i,d,h,t,y))
    u_i = parameter(am, bus_idx, :gen_bus, tech_idx, "EXUNITS")

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        g_idht + curt_idht - shape * u_i * CAP * PMAX
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr == 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)"))  end   

end

function constraint_VREBalance_Budget_NonReserve_iy_RA(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, day_group_id::Int, y::Int, type::String; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=true)

    ids_ht = [(h,t) for h in am.setting["run_H"] for t in am.setting["run_T"]]
    
    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  
    CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")

    # water budget of day_group_id in bus_idx (MWh/MW)
    budget = am.ref[:nw][0][:hydro_budget][bus_idx][day_group_id]

    sum_generation = JuMP.AffExpr(0.0)
    for d in am.ref[:nw][0][:repday_groups][day_group_id]["Day_Idx_List"]
        for (h, t) in ids_ht
            JuMP.add_to_expression!(sum_generation, 1.0, variable(am, nw, decomp_group, :g_idhty, (i,d,h,t,y)))
        end
    end

    u_i = parameter(am, bus_idx, :gen_bus, tech_idx, "EXUNITS")

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$y)"))) end    
    expr = JuMP.@expression(JuMP_model,  
        sum_generation - budget * u_i * CAP
        )   
    constraint = JuMP.@constraint(JuMP_model, 
        expr <= 0
        )    
    if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$y)"))  end   

end

function constraint_ES_Charge_Max_Sto_UC_NonReserve_idhty(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, i::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false, const_name_flag::Bool=false)

    bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
    tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)  

    if parameter(am, bus_idx, :gen_bus, tech_idx, "Storage Commitment") == true
        
        CAP = parameter(am, bus_idx, :gen_bus, tech_idx, "CAP")
        MAXINVEST = parameter(am, bus_idx, :gen_bus, tech_idx, "MAXINVEST")
        Big_M = max(MAXINVEST, 100)
        
        chg_idht = variable(am, nw, decomp_group, :chg_idhty, (i,d,h,t,y))
            
        sto_c_idhty = variable(am, nw, decomp_group, :sto_c_idhty, (i,d,h,t,y))
    

        if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($i,$d,$h,$t,$y)"))) end    
        expr = JuMP.@expression(JuMP_model,  
            chg_idht - (1 - sto_c_idhty) * CAP * Big_M
            )   
        constraint = JuMP.@constraint(JuMP_model, 
            expr <= 0
            )    
        if const_name_flag JuMP.set_name(constraint, string(const_name, "_($i,$d,$h,$t,$y)")) end    
    end

end

function constraint_b_theta_power_flow_kdhty_RA(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, const_name::String, decomp_group::Int, k::Int, d::Int, h::Int, t::Int, y::Int; nw::Int=am.cnw, update::Bool=false)

    t_bus = parameter(am, 0, :branch, k, "t_bus")
    f_bus = parameter(am, 0, :branch, k, "f_bus")

    b = 1.0 / parameter(am, 0, :branch, k, "br_x_pu")
    # geometric-mean row scaling keeps the ~1e5 susceptance off the matrix ceiling; abs() handles
    # negative-reactance branches (sign stays on the angle term via b/sb). Exact row-scaling, solution-neutral.
    sb = sqrt(abs(b))

    f_kdhty = variable(am, nw, decomp_group, :f_kdhty, (k,d,h,t,y))
    t_bus_angle = variable(am, nw, decomp_group, :bus_angle_ndhty, (t_bus,d,h,t,y))
    f_bus_angle = variable(am, nw, decomp_group, :bus_angle_ndhty, (f_bus,d,h,t,y))

    if update JuMP.delete(JuMP_model, JuMP.constraint_by_name(JuMP_model, string(const_name, "_($k,$d,$h,$t,$y)"))) end
    expr = JuMP.@expression(JuMP_model,
        (1.0 / sb) * f_kdhty - (b / sb) * (f_bus_angle - t_bus_angle)
        )
    constraint = JuMP.@constraint(JuMP_model,
        expr == 0
        )
    JuMP.set_name(constraint, string(const_name, "_($k,$d,$h,$t,$y)"))
end
