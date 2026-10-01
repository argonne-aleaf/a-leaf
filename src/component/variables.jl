# Copyright © 2026, UChicago Argonne, LLC. All Rights Reserved.
# SPDX-License-Identifier: BSD-3-Clause  (see LICENSE)

# Shared JuMP variable builders reused across ALEAF models.
# Keep naming aligned with variable registration keys.


function variable_chg_idhty_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, category::Symbol, variable_name::String, decomp_group, ids_i, ids_d, ids_h, ids_t, ids_y; nw::Int=am.cnw, bounded::Bool=true, report::Bool=true, fix_flag::Bool = false, result::Dict{String,<:Any} = Dict{String,Any}())
    
    ids = [(i, d, h, t, y) for (i) in ids_i for (d) in ids_d for h in ids_h for t in ids_t for y in ids_y]

    var = am.var[:nw][nw][decomp_group][Symbol(variable_name)] = JuMP.@variable(JuMP_model,
        [idx in ids],
        base_name=variable_name,
        integer=false,
        binary=false
    )

    if bounded
        for v in var
            JuMP.set_lower_bound(v, 0)
        end
    end

    if fix_flag == true
        for idx in ids
            fix_value = result["solution"]["storage"][string(idx)]["chg_idhty"]
            JuMP.fix(var[idx], fix_value; force = true)
        end
    end

    report && add_sol_component(am, nw, category, Symbol(variable_name), ids, var, decomp_group)    
end

function variable_soc_idhty_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, category::Symbol, variable_name::String, decomp_group, ids_i, ids_d, ids_h, ids_t, ids_y; nw::Int=am.cnw, bounded::Bool=true, report::Bool=true, fix_flag::Bool = false, result::Dict{String,<:Any} = Dict{String,Any}())
    
    ids = [(i, d, h, t, y) for (i) in ids_i for (d) in ids_d for h in ids_h for t in ids_t for y in ids_y]

    var = am.var[:nw][nw][decomp_group][Symbol(variable_name)] = JuMP.@variable(JuMP_model,
        [idx in ids],
        base_name=variable_name,
        integer=false,
        binary=false
    )

    if bounded
        for v in var
            JuMP.set_lower_bound(v, 0)
        end
    end

    if fix_flag == true
        for idx in ids
            fix_value = result["solution"]["storage"][string(idx)]["soc_idhty"]
            JuMP.fix(var[idx], fix_value; force = true)
        end
    end

    report && add_sol_component(am, nw, category, Symbol(variable_name), ids, var, decomp_group)    
end

function variable_rns_flex_dn_zdhty_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, category::Symbol, variable_name::String, decomp_group, ids_z, ids_d, ids_h, ids_t, ids_y; nw::Int=am.cnw, bounded::Bool=true, report::Bool=true, fix_flag::Bool = false, result::Dict{String,<:Any} = Dict{String,Any}())
    
    ids = [(z, d, h, t, y) for z in ids_z for d in ids_d for h in ids_h for t in ids_t for y in ids_y]
        
    var = am.var[:nw][nw][decomp_group][Symbol(variable_name)] = JuMP.@variable(JuMP_model,
        [idx in ids],
        base_name=variable_name,
        integer=false,
        binary=false
    )

    if bounded
        for v in var
            JuMP.set_lower_bound(v, 0)
        end
    else
        for idx in ids
            JuMP.fix(var[idx], 0.0; force = true)
        end
    end

    if fix_flag == true
        for idx in ids
            fix_value = result["solution"]["scarcity"][string(idx)]["rns_flex_dn_zdhty"]
            JuMP.fix(var[idx], fix_value; force = true)
        end
    end

    report && add_sol_component(am, nw, category, Symbol(variable_name), ids, var, decomp_group)    
end

function variable_rns_flex_up_zdhty_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, category::Symbol, variable_name::String, decomp_group, ids_z, ids_d, ids_h, ids_t, ids_y; nw::Int=am.cnw, bounded::Bool=true, report::Bool=true, fix_flag::Bool = false, result::Dict{String,<:Any} = Dict{String,Any}())
    
    ids = [(z, d, h, t, y) for z in ids_z for d in ids_d for h in ids_h for t in ids_t for y in ids_y]
        
    var = am.var[:nw][nw][decomp_group][Symbol(variable_name)] = JuMP.@variable(JuMP_model,
        [idx in ids],
        base_name=variable_name,
        integer=false,
        binary=false
    )

    if bounded
        for v in var
            JuMP.set_lower_bound(v, 0)
        end
    else
        for idx in ids
            JuMP.fix(var[idx], 0.0; force = true)
        end
    end

    if fix_flag == true
        for idx in ids
            fix_value = result["solution"]["scarcity"][string(idx)]["rns_flex_up_zdhty"]
            JuMP.fix(var[idx], fix_value; force = true)
        end
    end

    report && add_sol_component(am, nw, category, Symbol(variable_name), ids, var, decomp_group)    
end

function variable_rns_nonspin_zdhty_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, category::Symbol, variable_name::String, decomp_group, ids_z, ids_d, ids_h, ids_t, ids_y; nw::Int=am.cnw, bounded::Bool=true, report::Bool=true, fix_flag::Bool = false, result::Dict{String,<:Any} = Dict{String,Any}())
    
    ids = [(z, d, h, t, y) for z in ids_z for d in ids_d for h in ids_h for t in ids_t for y in ids_y]
        
    var = am.var[:nw][nw][decomp_group][Symbol(variable_name)] = JuMP.@variable(JuMP_model,
        [idx in ids],
        base_name=variable_name,
        integer=false,
        binary=false
    )

    if bounded
        for v in var
            JuMP.set_lower_bound(v, 0)
        end
    else
        for idx in ids
            JuMP.fix(var[idx], 0.0; force = true)
        end
    end

    if fix_flag == true
        for idx in ids
            fix_value = result["solution"]["scarcity"][string(idx)]["rns_nonspin_zdhty"]
            JuMP.fix(var[idx], fix_value; force = true)
        end
    end

    report && add_sol_component(am, nw, category, Symbol(variable_name), ids, var, decomp_group)   
end

function variable_rns_reg_dn_zdhty_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, category::Symbol, variable_name::String, decomp_group, ids_z, ids_d, ids_h, ids_t, ids_y; nw::Int=am.cnw, bounded::Bool=true, report::Bool=true, fix_flag::Bool = false, result::Dict{String,<:Any} = Dict{String,Any}())
    
    ids = [(z, d, h, t, y) for z in ids_z for d in ids_d for h in ids_h for t in ids_t for y in ids_y]
        
    var = am.var[:nw][nw][decomp_group][Symbol(variable_name)] = JuMP.@variable(JuMP_model,
        [idx in ids],
        base_name=variable_name,
        integer=false,
        binary=false
    )

    if bounded
        for v in var
            JuMP.set_lower_bound(v, 0)
        end
    else
        for idx in ids
            JuMP.fix(var[idx], 0.0; force = true)
        end
    end

    if fix_flag == true
        for idx in ids
            fix_value = result["solution"]["scarcity"][string(idx)]["rns_reg_dn_zdhty"]
            JuMP.fix(var[idx], fix_value; force = true)
        end
    end

    report && add_sol_component(am, nw, category, Symbol(variable_name), ids, var, decomp_group)   
end

function variable_rns_reg_up_zdhty_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, category::Symbol, variable_name::String, decomp_group, ids_z, ids_d, ids_h, ids_t, ids_y; nw::Int=am.cnw, bounded::Bool=true, report::Bool=true, fix_flag::Bool = false, result::Dict{String,<:Any} = Dict{String,Any}())
    
    ids = [(z, d, h, t, y) for z in ids_z for d in ids_d for h in ids_h for t in ids_t for y in ids_y]
        
    var = am.var[:nw][nw][decomp_group][Symbol(variable_name)] = JuMP.@variable(JuMP_model,
        [idx in ids],
        base_name=variable_name,
        integer=false,
        binary=false
    )

    if bounded
        for v in var
            JuMP.set_lower_bound(v, 0)
        end
    else
        for idx in ids
            JuMP.fix(var[idx], 0.0; force = true)
        end
    end

    if fix_flag == true
        for idx in ids
            fix_value = result["solution"]["scarcity"][string(idx)]["rns_reg_up_zdhty"]
            JuMP.fix(var[idx], fix_value; force = true)
        end
    end

    report && add_sol_component(am, nw, category, Symbol(variable_name), ids, var, decomp_group)   
end

function variable_rns_cont_zdhty_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, category::Symbol, variable_name::String, decomp_group, ids_z, ids_d, ids_h, ids_t, ids_y; nw::Int=am.cnw, bounded::Bool=true, report::Bool=true, fix_flag::Bool = false, result::Dict{String,<:Any} = Dict{String,Any}())
    
    ids = [(z, d, h, t, y) for z in ids_z for d in ids_d for h in ids_h for t in ids_t for y in ids_y]
        
    var = am.var[:nw][nw][decomp_group][Symbol(variable_name)] = JuMP.@variable(JuMP_model,
        [idx in ids],
        base_name=variable_name,
        integer=false,
        binary=false
    )

    if bounded
        for v in var
            JuMP.set_lower_bound(v, 0)
        end
    else
        for idx in ids
            JuMP.fix(var[idx], 0.0; force = true)
        end
    end

    if fix_flag == true
        for idx in ids
            fix_value = result["solution"]["scarcity"][string(idx)]["rns_cont_zdhty"]
            JuMP.fix(var[idx], fix_value; force = true)
        end
    end

    report && add_sol_component(am, nw, category, Symbol(variable_name), ids, var, decomp_group)   
end

function variable_max_ENS_MWh_y_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, category::Symbol, variable_name::String, decomp_group, ids_y; nw::Int=am.cnw, bounded::Bool=true, report::Bool=true, fix_flag::Bool = false, result::Dict{String,<:Any} = Dict{String,Any}())
    
    var = am.var[:nw][nw][decomp_group][Symbol(variable_name)] = JuMP.@variable(JuMP_model,
        [idx in ids_y],
        base_name=variable_name,
        integer=false,
        binary=false
    )

    if bounded
        for idx in ids_y
            JuMP.set_lower_bound(var[idx], 0)
        end
    end

    report && add_sol_component(am, nw, category, Symbol(variable_name), ids_y, var, decomp_group) 
end

function variable_slack_RPS_ny_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, category::Symbol, variable_name::String, decomp_group, ids_n, ids_y; nw::Int=am.cnw, bounded::Bool=true, report::Bool=true, fix_flag::Bool = false, result::Dict{String,<:Any} = Dict{String,Any}())
    
    ids = [(n, y) for n in ids_n for y in ids_y]
        
    var = am.var[:nw][nw][decomp_group][Symbol(variable_name)] = JuMP.@variable(JuMP_model,
        [idx in ids],
        base_name=variable_name,
        integer=false,
        binary=false
    )

    if bounded
        for v in var
            JuMP.set_lower_bound(v, 0)
        end
    end

    report && add_sol_component(am, nw, category, Symbol(variable_name), ids, var, decomp_group) 
end

function variable_slack_CEG_ndy_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, category::Symbol, variable_name::String, decomp_group, ids_n, ids_d, ids_y; nw::Int=am.cnw, bounded::Bool=true, report::Bool=true, fix_flag::Bool = false, result::Dict{String,<:Any} = Dict{String,Any}())
    
    ids = [(n, d, y) for n in ids_n for d in ids_d for y in ids_y]
        
    var = am.var[:nw][nw][decomp_group][Symbol(variable_name)] = JuMP.@variable(JuMP_model,
        [idx in ids],
        base_name=variable_name,
        integer=false,
        binary=false
    )

    if bounded
        for v in var
            JuMP.set_lower_bound(v, 0)
        end
    end

    report && add_sol_component(am, nw, category, Symbol(variable_name), ids, var, decomp_group) 
end

function variable_slack_CEG_ny_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, category::Symbol, variable_name::String, decomp_group, ids_n, ids_y; nw::Int=am.cnw, bounded::Bool=true, report::Bool=true, fix_flag::Bool = false, result::Dict{String,<:Any} = Dict{String,Any}())
    
    ids = [(n, y) for n in ids_n for y in ids_y]
        
    var = am.var[:nw][nw][decomp_group][Symbol(variable_name)] = JuMP.@variable(JuMP_model,
        [idx in ids],
        base_name=variable_name,
        integer=false,
        binary=false
    )

    if bounded
        for v in var
            JuMP.set_lower_bound(v, 0)
        end
    end

    report && add_sol_component(am, nw, category, Symbol(variable_name), ids, var, decomp_group) 
end

function variable_p_inj_ndhty_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, category::Symbol, variable_name::String, decomp_group, ids_n, ids_d, ids_h, ids_t, ids_y; nw::Int=am.cnw, bounded::Bool=true, report::Bool=true, fix_flag::Bool = false, result::Dict{String,<:Any} = Dict{String,Any}())
    
    ids = [(n, d, h, t, y) for n in ids_n for d in ids_d for h in ids_h for t in ids_t for y in ids_y]
        
    var = am.var[:nw][nw][decomp_group][Symbol(variable_name)] = JuMP.@variable(JuMP_model,
        [idx in ids],
        base_name=variable_name,
        integer=false,
        binary=false
    )


    if fix_flag == true
        for idx in ids
            fix_value = result["solution"]["powerflow"][string(idx)]["p_inj_ndhty"]
            JuMP.fix(var[idx], fix_value; force = true)
        end
    end

    report && add_sol_component(am, nw, category, Symbol(variable_name), ids, var, decomp_group) 
end

function variable_flow_relax_fdhty_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, category::Symbol, variable_name::String, decomp_group, ids_k, ids_d, ids_h, ids_t, ids_y; nw::Int=am.cnw, bounded::Bool=true, report::Bool=true, fix_flag::Bool = false, result::Dict{String,<:Any} = Dict{String,Any}())
    
    ids = [(k, d, h, t, y) for k in ids_k for d in ids_d for h in ids_h for t in ids_t for y in ids_y]
        
    var = am.var[:nw][nw][decomp_group][Symbol(variable_name)] = JuMP.@variable(JuMP_model,
        [idx in ids],
        base_name=variable_name,
        integer=false,
        binary=false
    )

    if bounded
        for v in var
            JuMP.set_lower_bound(v, 0)
        end
    end

    report && add_sol_component(am, nw, category, Symbol(variable_name), ids, var, decomp_group) 
end

function variable_demand_ndhty_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, category::Symbol, variable_name::String, decomp_group, ids_n, ids_d, ids_h, ids_t, ids_y; nw::Int=am.cnw, bounded::Bool=true, report::Bool=true, fix_flag::Bool = false, result::Dict{String,<:Any} = Dict{String,Any}())
    
    ids = [(n, d, h, t, y) for n in ids_n for d in ids_d for h in ids_h for t in ids_t for y in ids_y]
        
    var = am.var[:nw][nw][decomp_group][Symbol(variable_name)] = JuMP.@variable(JuMP_model,
        [idx in ids],
        base_name=variable_name,
        integer=false,
        binary=false
    )

    if bounded
        for v in var
            JuMP.set_lower_bound(v, 0)
        end
    end

    if fix_flag == true
        for idx in ids
            fix_value = result["solution"]["powerflow"][string(idx)]["demand_ndhty"]
            JuMP.fix(var[idx], fix_value; force = true)
        end
    end

    report && add_sol_component(am, nw, category, Symbol(variable_name), ids, var, decomp_group) 
end

function variable_ens_ndhty_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, category::Symbol, variable_name::String, decomp_group, ids_n, ids_d, ids_h, ids_t, ids_y; nw::Int=am.cnw, bounded::Bool=true, report::Bool=true, fix_flag::Bool = false, result::Dict{String,<:Any} = Dict{String,Any}())
    
    ids = [(n, d, h, t, y) for n in ids_n for d in ids_d for h in ids_h for t in ids_t for y in ids_y]
        
    var = am.var[:nw][nw][decomp_group][Symbol(variable_name)] = JuMP.@variable(JuMP_model,
        [idx in ids],
        base_name=variable_name,
        integer=false,
        binary=false
    )

    if bounded
        for v in var
            JuMP.set_lower_bound(v, 0)
        end
    else
        for idx in ids
            JuMP.fix(var[idx], 0.0; force = true)
        end
    end

    if fix_flag == true
        for idx in ids
            fix_value = result["solution"]["scarcity"][string(idx)]["ens_ndhty"]
            JuMP.fix(var[idx], fix_value; force = true)
        end
    end

    report && add_sol_component(am, nw, category, Symbol(variable_name), ids, var, decomp_group) 
end

function variable_curt_idhty_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, category::Symbol, variable_name::String, decomp_group, ids_i, ids_d, ids_h, ids_t, ids_y; nw::Int=am.cnw, bounded::Bool=true, report::Bool=true, fix_flag::Bool = false, result::Dict{String,<:Any} = Dict{String,Any}())
    
    ids_i_curt = []

    for i in ids_i
        bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
        tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)

        if parameter(am, bus_idx, :gen_bus, tech_idx, "VRE_Flag") == true
            if (parameter(am, bus_idx, :gen_bus, tech_idx, "FUEL_LIMIT") == "Fixed Profile")
                push!(ids_i_curt, i)
            end
        end
    end
    
    
    ids = [(i, d, h, t, y) for (i) in ids_i_curt for (d) in ids_d for h in ids_h for t in ids_t for y in ids_y]

    var = am.var[:nw][nw][decomp_group][Symbol(variable_name)] = JuMP.@variable(JuMP_model,
        [idx in ids],
        base_name=variable_name,
        integer=false,
        binary=false
    )

    if bounded
        for v in var
            JuMP.set_lower_bound(v, 0)
        end
    end

    if fix_flag == true
        for idx in ids
            fix_value = result["solution"]["dispatch"][string(idx)]["curt_idhty"]
            JuMP.fix(var[idx], fix_value; force = true)
        end
    end

    report && add_sol_component(am, nw, category, Symbol(variable_name), ids, var, decomp_group)    
end

function variable_spin_idhty_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, category::Symbol, variable_name::String, decomp_group, ids_i, ids_d, ids_h, ids_t, ids_y; nw::Int=am.cnw, bounded::Bool=true, report::Bool=true, fix_flag::Bool = false, result::Dict{String,<:Any} = Dict{String,Any}())
    
    ids = [(i, d, h, t, y) for (i) in ids_i for (d) in ids_d for h in ids_h for t in ids_t for y in ids_y]

    var = am.var[:nw][nw][decomp_group][Symbol(variable_name)] = JuMP.@variable(JuMP_model,
        [idx in ids],
        base_name=variable_name,
        integer=false,
        binary=false
    )

    if bounded
        for v in var
            JuMP.set_lower_bound(v, 0)
        end
    end

    if fix_flag == true
        for idx in ids
            fix_value = result["solution"]["reserve"][string(idx)]["spin_idhty"]
            JuMP.fix(var[idx], fix_value; force = true)
        end
    end

    report && add_sol_component(am, nw, category, Symbol(variable_name), ids, var, decomp_group)    
end

function variable_flex_dn_idhty_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, category::Symbol, variable_name::String, decomp_group, ids_i, ids_d, ids_h, ids_t, ids_y; nw::Int=am.cnw, bounded::Bool=true, report::Bool=true, fix_flag::Bool = false, result::Dict{String,<:Any} = Dict{String,Any}())
    
    ids = [(i, d, h, t, y) for (i) in ids_i for (d) in ids_d for h in ids_h for t in ids_t for y in ids_y]

    var = am.var[:nw][nw][decomp_group][Symbol(variable_name)] = JuMP.@variable(JuMP_model,
        [idx in ids],
        base_name=variable_name,
        integer=false,
        binary=false
    )

    if bounded
        for v in var
            JuMP.set_lower_bound(v, 0)
        end
    end

    if fix_flag == true
        for idx in ids
            fix_value = result["solution"]["reserve"][string(idx)]["flex_dn_idhty"]
            JuMP.fix(var[idx], fix_value; force = true)
        end
    end

    report && add_sol_component(am, nw, category, Symbol(variable_name), ids, var, decomp_group)    
end

function variable_flex_up_idhty_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, category::Symbol, variable_name::String, decomp_group, ids_i, ids_d, ids_h, ids_t, ids_y; nw::Int=am.cnw, bounded::Bool=true, report::Bool=true, fix_flag::Bool = false, result::Dict{String,<:Any} = Dict{String,Any}())
    
    ids = [(i, d, h, t, y) for (i) in ids_i for (d) in ids_d for h in ids_h for t in ids_t for y in ids_y]

    var = am.var[:nw][nw][decomp_group][Symbol(variable_name)] = JuMP.@variable(JuMP_model,
        [idx in ids],
        base_name=variable_name,
        integer=false,
        binary=false
    )

    if bounded
        for v in var
            JuMP.set_lower_bound(v, 0)
        end
    end

    if fix_flag == true
        for idx in ids
            fix_value = result["solution"]["reserve"][string(idx)]["flex_up_idhty"]
            JuMP.fix(var[idx], fix_value; force = true)
        end
    end

    report && add_sol_component(am, nw, category, Symbol(variable_name), ids, var, decomp_group)    
end

function variable_reg_dn_idhty_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, category::Symbol, variable_name::String, decomp_group, ids_i, ids_d, ids_h, ids_t, ids_y; nw::Int=am.cnw, bounded::Bool=true, report::Bool=true, fix_flag::Bool = false, result::Dict{String,<:Any} = Dict{String,Any}())
    
    ids = [(i, d, h, t, y) for (i) in ids_i for (d) in ids_d for h in ids_h for t in ids_t for y in ids_y]

    var = am.var[:nw][nw][decomp_group][Symbol(variable_name)] = JuMP.@variable(JuMP_model,
        [idx in ids],
        base_name=variable_name,
        integer=false,
        binary=false
    )

    if bounded
        for v in var
            JuMP.set_lower_bound(v, 0)
        end
    end

    if fix_flag == true
        for idx in ids
            fix_value = result["solution"]["reserve"][string(idx)]["reg_dn_idhty"]
            JuMP.fix(var[idx], fix_value; force = true)
        end
    end

    report && add_sol_component(am, nw, category, Symbol(variable_name), ids, var, decomp_group)    
end

function variable_reg_up_idhty_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, category::Symbol, variable_name::String, decomp_group, ids_i, ids_d, ids_h, ids_t, ids_y; nw::Int=am.cnw, bounded::Bool=true, report::Bool=true, fix_flag::Bool = false, result::Dict{String,<:Any} = Dict{String,Any}())
    
    ids = [(i, d, h, t, y) for (i) in ids_i for (d) in ids_d for h in ids_h for t in ids_t for y in ids_y]

    var = am.var[:nw][nw][decomp_group][Symbol(variable_name)] = JuMP.@variable(JuMP_model,
        [idx in ids],
        base_name=variable_name,
        integer=false,
        binary=false
    )

    if bounded
        for v in var
            JuMP.set_lower_bound(v, 0)
        end
    end

    if fix_flag == true
        for idx in ids
            fix_value = result["solution"]["reserve"][string(idx)]["reg_up_idhty"]
            JuMP.fix(var[idx], fix_value; force = true)
        end
    end

    report && add_sol_component(am, nw, category, Symbol(variable_name), ids, var, decomp_group)    
end

function variable_g_idhty_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, category::Symbol, variable_name::String, decomp_group, ids_i, ids_d, ids_h, ids_t, ids_y; nw::Int=am.cnw, bounded::Bool=true, report::Bool=true, fix_flag::Bool = false, result::Dict{String,<:Any} = Dict{String,Any}())
    
    ids = [(i, d, h, t, y) for (i) in ids_i for (d) in ids_d for h in ids_h for t in ids_t for y in ids_y]

    var = am.var[:nw][nw][decomp_group][Symbol(variable_name)] = JuMP.@variable(JuMP_model,
        [idx in ids],
        base_name=variable_name,
        integer=false,
        binary=false
    )

    if bounded
        for v in var
            JuMP.set_lower_bound(v, 0)
        end
    end

    if fix_flag == true
        for idx in ids
            fix_value = result["solution"]["dispatch"][string(idx)]["g_idhty"]
            JuMP.fix(var[idx], fix_value; force = true)
        end
    end

    report && add_sol_component(am, nw, category, Symbol(variable_name), ids, var, decomp_group)    
end

function variable_su_idhty_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, category::Symbol, variable_name::String, decomp_group, ids_i, ids_d, ids_h, ids_t, ids_y; nw::Int=am.cnw, bounded::Bool=true, report::Bool=true, fix_flag::Bool = false, result::Dict{String,<:Any} = Dict{String,Any}())
    
    ids = [(i, d, h, t, y) for (i) in ids_i for (d) in ids_d for h in ids_h for t in ids_t for y in ids_y]

    var = am.var[:nw][nw][decomp_group][Symbol(variable_name)] = JuMP.@variable(JuMP_model,
        [idx in ids],
        base_name=variable_name,
        integer=false,
        binary=false
    )

    if bounded
        for v in var
            JuMP.set_lower_bound(v, 0)
        end
    end

    if fix_flag == true
        for idx in ids
            fix_value = result["solution"]["commitment"][string(idx)]["su_idhty"]
            JuMP.fix(var[idx], fix_value; force = true)
        end
    end

    report && add_sol_component(am, nw, category, Symbol(variable_name), ids, var, decomp_group)    
end

function variable_sto_c_idhty_integer(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, category::Symbol, variable_name::String, decomp_group, ids_i, ids_d, ids_h, ids_t, ids_y; nw::Int=am.cnw, bounded::Bool=true, report::Bool=true, fix_flag::Bool = false, result::Dict{String,<:Any} = Dict{String,Any}())
    
    ids_i_commit = []

    for i in ids_i
        bus_idx = am.ref[:nw][0][:gen_index][i]["bus_idx"]
        tech_idx = am.ref[:nw][0][:gen_index][i]["genco_tech_id"]
        if parameter(am, 0, :gen_index, "UNIT_CATEGORY", i) == "STORAGE"
            if parameter(am, bus_idx, :gen_bus, tech_idx, "Storage Commitment") == true
                push!(ids_i_commit, i)
            end
        end
    end
    
    ids = [(i, d, h, t, y) for i in ids_i_commit for (d) in ids_d for h in ids_h for t in ids_t for y in ids_y]

    if fix_flag == true # c_idhty is non-integer when we fix it.
        var = am.var[:nw][nw][decomp_group][Symbol(variable_name)] = JuMP.@variable(JuMP_model,
            [idx in ids],
            base_name=variable_name,
            integer=false,
            binary=false
        )
    else
        var = am.var[:nw][nw][decomp_group][Symbol(variable_name)] = JuMP.@variable(JuMP_model,
            [idx in ids],
            base_name=variable_name,
            integer=false,
            binary=true
        )
    end

    if bounded
        for idx in ids
            JuMP.set_lower_bound(var[idx], 0)
            JuMP.set_upper_bound(var[idx], 1)
        end
    end

    if fix_flag == true
        for idx in ids
            fix_value = round(Int, result["solution"]["storage_commitment"][string(idx)]["sto_c_idhty"])
            JuMP.fix(var[idx], fix_value; force = true)
        end
    end

    report && add_sol_component(am, nw, category, Symbol(variable_name), ids, var, decomp_group)    
end

function variable_water_use_jidhty_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, category::Symbol, variable_name::String, decomp_group, ids_d, ids_h, ids_t, ids_y; nw::Int=am.cnw, bounded::Bool=true, report::Bool=true, fix_flag::Bool = false, result::Dict{String,<:Any} = Dict{String,Any}())
    
    # Define generators that needs water level variables
    ids_i = am.ref[:nw][0][:water_management]["gen_index"]
    ids_j = am.ref[:nw][0][:water_management]["segment_index"]

    ids = [(i, j, d, h, t, y) for i in ids_i for j in ids_j for (d) in ids_d for h in ids_h for t in ids_t for y in ids_y]

    var = am.var[:nw][nw][decomp_group][Symbol(variable_name)] = JuMP.@variable(JuMP_model,
            [idx in ids],
            base_name=variable_name,
            integer=false,
            binary=false
        )

    if bounded
        for v in var
            JuMP.set_lower_bound(v, 0)
        end
    end

    if fix_flag == true
        for idx in ids
            fix_value = result["solution"]["water_management"][string(idx)]["water_use_ijdhty"]
            JuMP.fix(var[idx], fix_value; force = true)
        end
    end

    report && add_sol_component(am, nw, category, Symbol(variable_name), ids, var, decomp_group)    
end

function variable_c_idhty_integer(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, category::Symbol, variable_name::String, decomp_group, ids_i, ids_d, ids_h, ids_t, ids_y; nw::Int=am.cnw, bounded::Bool=true, report::Bool=true, fix_flag::Bool = false, result::Dict{String,<:Any} = Dict{String,Any}())
    
    
    ids = [(i, d, h, t, y) for (i) in ids_i for (d) in ids_d for h in ids_h for t in ids_t for y in ids_y]

    if fix_flag == true # c_idhty is non-integer when we fix it.
        var = am.var[:nw][nw][decomp_group][Symbol(variable_name)] = JuMP.@variable(JuMP_model,
            [idx in ids],
            base_name=variable_name,
            integer=false,
            binary=false
        )
    else
        var = am.var[:nw][nw][decomp_group][Symbol(variable_name)] = JuMP.@variable(JuMP_model,
            [idx in ids],
            base_name=variable_name,
            integer=true,
            binary=false
        )
    end

    if bounded
        for v in var
            JuMP.set_lower_bound(v, 0)
        end
    end

    if fix_flag == true
        for idx in ids
            fix_value = result["solution"]["commitment"][string(idx)]["c_idhty"]
            if fix_value > 0.0
                JuMP.fix(var[idx], fix_value; force = true)
            else
                JuMP.fix(var[idx], 0.0; force = true)
            end
        end
    end

    report && add_sol_component(am, nw, category, Symbol(variable_name), ids, var, decomp_group)    
end

function variable_u_ESE_iy_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, category::Symbol, variable_name::String, decomp_group, ids_i, ids_y; nw::Int=am.cnw, bounded::Bool=true, report::Bool=true, fix_flag::Bool = false, result::Dict{String,<:Any} = Dict{String,Any}(), PH_val::Bool=false)
    
    ids = [(i, y) for (i) in ids_i for y in ids_y]
    integer_flag_PH = false

    var = am.var[:nw][nw][decomp_group][Symbol(variable_name)] = JuMP.@variable(JuMP_model,
        [(i, y) in ids],
        base_name = variable_name,
        integer = false,
        binary = false
    )

    if bounded
        for idx in ids

            JuMP.set_lower_bound(var[idx], 0)
        end
    end

    if fix_flag == true
        for idx in ids

            if am.setting["Simulation Configuration"]["Run_expansion_flag"] == true
                pu_power_base = am.setting["Simulation Setting"]["per_unit_base_value"]
                fix_value = result["solution"]["expansion"][string(idx)]["u_ESE_iy"] / pu_power_base    # the u_ESE_iy value in the expansion solution is not per-unit when we run OP after EXP reporting
             
                JuMP.fix(var[idx], fix_value; force = true)
            
            elseif am.setting["Simulation Configuration"]["Use_predefined_expansion_data_for_OP_flag"] == true
                fix_value = result["solution"]["expansion"][string(idx)]["u_ESE_iy"]                    # the u_ESE_iy value in the pre-defined expansion solution is per-unit
             
                JuMP.fix(var[idx], fix_value; force = true)
            end
            
        end
    end

    report && add_sol_component(am, nw, category, Symbol(variable_name), ids, var, decomp_group)    
end


function variable_u_T_ky_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, category::Symbol, variable_name::String, decomp_group, ids_k, ids_y; nw::Int=am.cnw, bounded::Bool=true, report::Bool=true, fix_flag::Bool = false, result::Dict{String,<:Any} = Dict{String,Any}(), PH_val::Bool=false)
    
    ids = [(k, y) for (k) in ids_k for y in ids_y]
    integer_flag_PH = false
    
    var = am.var[:nw][nw][decomp_group][Symbol(variable_name)] = JuMP.@variable(JuMP_model,
        [(k, y) in ids],
        base_name = variable_name,
        integer = false,
        binary = false
    )

    if bounded
        for v in var
            JuMP.set_lower_bound(v, 0)
        end
    end

    if fix_flag == true
        for idx in ids
            fix_value = result["solution"]["expansion"][string(idx)]["u_T_ky"]
             
            JuMP.fix(var[idx], fix_value; force = true)
        end
    end

    report && add_sol_component(am, nw, category, Symbol(variable_name), ids, var, decomp_group)    
end


function variable_u_G_iy_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, category::Symbol, variable_name::String, decomp_group, ids_i, ids_y; nw::Int=am.cnw, bounded::Bool=true, report::Bool=true, fix_flag::Bool = false, result::Dict{String,<:Any} = Dict{String,Any}(), PH_val::Bool=false)
    
    ids = [(i, y) for (i) in ids_i for y in ids_y]
    integer_flag_PH = false
    
    var = am.var[:nw][nw][decomp_group][Symbol(variable_name)] = JuMP.@variable(JuMP_model,
        [(i, y) in ids],
        base_name = variable_name,
        integer = false,
        binary = false
    )

    if bounded
        for v in var
            JuMP.set_lower_bound(v, 0)
        end
    end

    if fix_flag == true
        for idx in ids
            fix_value = result["solution"]["expansion"][string(idx)]["u_G_iy"]
             
            JuMP.fix(var[idx], fix_value; force = true)
        end
    end

    report && add_sol_component(am, nw, category, Symbol(variable_name), ids, var, decomp_group)    
end

function variable_u_newESH_iy_integer_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, category::Symbol, variable_name::String, decomp_group, ids_i, ids_y; nw::Int = am.cnw, bounded::Bool = true, report::Bool = true, fix_flag::Bool = false, result::Dict{String,<:Any} = Dict{String,Any}(), PH_val::Bool=false)

    ids = [(i, y) for (i) in ids_i for y in ids_y]

    var = am.var[:nw][nw][decomp_group][Symbol(variable_name)] = JuMP.@variable(JuMP_model,
        [(i, y) in ids],
        base_name = variable_name,
        integer = false,
        binary = false
    )
    
    if bounded
        for idx in ids

            (i, y) = idx

            bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
            tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)

            MAXINVEST = parameter(am, bus_idx, :gen_bus, tech_idx, "MAXINVEST")

            lower_limit = 0
            upper_limit = parameter(am, bus_idx, :gen_bus, tech_idx, "STOHR_MAX") * MAXINVEST
            
            JuMP.set_lower_bound(var[idx], lower_limit)

            JuMP.set_upper_bound(var[idx], upper_limit)
        end
    end

    integer_flag_PH = false
    for idx in ids

        (i, y) = idx

        bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
        tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)
        integer_flag = parameter(am, bus_idx, :gen_bus, tech_idx, "Integrality")

        if (integer_flag == true) & (fix_flag == false)
            integer_flag_PH = true
        end
    end

    if fix_flag == true
        for idx in ids
            fix_value = result["solution"]["expansion"][string(idx)]["u_new_ESH_iy"]
            if fix_value > 0.0
                JuMP.fix(var[idx], fix_value; force = true)
            else
                JuMP.fix(var[idx], 0.0; force = true)
            end
        end
    end

    report && add_sol_component(am, nw, category, Symbol(variable_name), ids, var, decomp_group)    
end

function variable_u_retG_iy_integer_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, category::Symbol, variable_name::String, decomp_group, ids_i, ids_y; nw::Int = am.cnw, bounded::Bool = true, report::Bool = true, fix_flag::Bool = false, result::Dict{String,<:Any} = Dict{String,Any}(), PH_val::Bool=false)

    ids = [(i, y) for (i) in ids_i for y in ids_y]

    var = am.var[:nw][nw][decomp_group][Symbol(variable_name)] = JuMP.@variable(JuMP_model,
        [(i, y) in ids],
        base_name = variable_name,
        integer = false,
        binary = false
    )

    if bounded
        for idx in ids
            
            (i, y) = idx

            bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
            tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)
            ret_flag = parameter(am, bus_idx, :gen_bus, tech_idx, "RET_FLAG")
            EXUNITS = parameter(am, bus_idx, :gen_bus, tech_idx, "EXUNITS")
            
            if y == -1
                JuMP.fix(var[idx], 0.0; force = true)
            else
                JuMP.set_lower_bound(var[idx], 0.0)
            end
        end
    end

    integer_flag_PH = false
    for idx in ids

        (i, y) = idx

        bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
        tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)
        integer_flag = parameter(am, bus_idx, :gen_bus, tech_idx, "Integrality")

        if (integer_flag == true) & (fix_flag == false)
            JuMP.set_integer(var[idx])
            integer_flag_PH = true
        end
    end

    if fix_flag == true
        for idx in ids
            fix_value = result["solution"]["expansion"][string(idx)]["u_ret_G_iy"]
            if fix_value > 0.0
                JuMP.fix(var[idx], fix_value; force = true)
            else
                JuMP.fix(var[idx], 0.0; force = true)
            end
        end
    end

    report && add_sol_component(am, nw, category, Symbol(variable_name), ids, var, decomp_group)    
end

function variable_u_newG_iy_integer_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, category::Symbol, variable_name::String, decomp_group, ids_i, ids_y; nw::Int = am.cnw, bounded::Bool = true, report::Bool = true, fix_flag::Bool = false, result::Dict{String,<:Any} = Dict{String,Any}(), PH_val::Bool=false, ref_expansion_result::Dict{String,<:Any} = Dict{String,Any}())

    ids = [(i, y) for (i) in ids_i for y in ids_y]

    var = am.var[:nw][nw][decomp_group][Symbol(variable_name)] = JuMP.@variable(JuMP_model,
            [(i, y) in ids],
            base_name = variable_name,
            integer = false,
            binary = false
        )

    MIP_relaxed_solution_bounds_step = am.setting["Simulation Setting"]["MIP_relaxed_solution_bounds_step_value"]

    integer_flag_PH = false
    if bounded
        if isempty(ref_expansion_result)
            for idx in ids
            
                (i, y) = idx
    
                bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
                tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)
                invest_flag = parameter(am, bus_idx, :gen_bus, tech_idx, "INVEST_FLAG")
                integer_flag = parameter(am, bus_idx, :gen_bus, tech_idx, "Integrality")
    
                if y == -1
                    JuMP.fix(var[idx], 0.0; force = true)
                else
                    if invest_flag == true
    
                        JuMP.set_lower_bound(var[idx], parameter(am, bus_idx, :gen_bus, tech_idx, "MININVEST"))
                        JuMP.set_upper_bound(var[idx], parameter(am, bus_idx, :gen_bus, tech_idx, "MAXINVEST"))

                        if (integer_flag == true) & (fix_flag == false)
                            JuMP.set_integer(var[idx])
                            integer_flag_PH = true
                        end
    
                    else
                        JuMP.fix(var[idx], 0.0; force = true)
                    end
                end
            end

        else

            for idx in ids
            
                (i, y) = idx
    
                bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
                tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)
                invest_flag = parameter(am, bus_idx, :gen_bus, tech_idx, "INVEST_FLAG")
                integer_flag = parameter(am, bus_idx, :gen_bus, tech_idx, "Integrality")
    
                if y == -1
                    JuMP.fix(var[idx], 0.0; force = true)
                else

                    if invest_flag == true
    
                        ref_u_newG_iy_solution = ref_expansion_result[string("(", i, ", ", y, ")")]["u_new_G_iy"]

                        upper_bound = 0.0
                        lower_bound = 0.0
                        warm_start_value = 0.0
                        if ref_u_newG_iy_solution > 0.0001
                            upper_bound = round(ref_u_newG_iy_solution, RoundUp) + MIP_relaxed_solution_bounds_step
                            lower_bound = min(0.0, round(ref_u_newG_iy_solution, RoundDown) - MIP_relaxed_solution_bounds_step)
                            warm_start_value = round(ref_u_newG_iy_solution)
                        end

                        if upper_bound > 0
                            JuMP.set_lower_bound(var[idx], lower_bound)
                            JuMP.set_upper_bound(var[idx], upper_bound)
                            JuMP.set_start_value(var[idx], warm_start_value)
                        else
                            JuMP.fix(var[idx], 0.0; force = true)
                            JuMP.set_start_value(var[idx], 0.0)
                        end

                        if upper_bound > 0
                            if (integer_flag == true) & (fix_flag == false)
                                JuMP.set_integer(var[idx])
                                integer_flag_PH = true
                            end
                        end
    
                    else
                        JuMP.fix(var[idx], 0.0; force = true)
                    end
                end
            end


        end


        
    end

    
    for idx in ids

        (i, y) = idx

        bus_idx = parameter(am, 0, :gen_index, "bus_idx", i)
        tech_idx = parameter(am, 0, :gen_index, "genco_tech_id", i)
        integer_flag = parameter(am, bus_idx, :gen_bus, tech_idx, "Integrality")

        if (integer_flag == true) & (fix_flag == false)
            JuMP.set_integer(var[idx])
            integer_flag_PH = true
        end
    end

    if fix_flag == true
        for idx in ids
            fix_value = result["solution"]["expansion"][string(idx)]["u_new_G_iy"]
            JuMP.fix(var[idx], fix_value; force = true)
        end
    end

    report && add_sol_component(am, nw, category, Symbol(variable_name), ids, var, decomp_group)
end

function variable_u_newT_ky_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, category::Symbol, variable_name::String, decomp_group, ids_k, ids_y; nw::Int = am.cnw, bounded::Bool = true, report::Bool = true, fix_flag::Bool = false, result::Dict{String,<:Any} = Dict{String,Any}(), PH_val::Bool=false, cycling_flag::Bool=false)
    
    ids_ky = [(k, y) for (k) in ids_k for y in ids_y]
    integer_flag_PH = false

    var = am.var[:nw][nw][decomp_group][Symbol(variable_name)] = JuMP.@variable(JuMP_model,
        [(k, y) in ids_ky],
        base_name = variable_name,
        integer = false,
        binary = false
    )
        
    if bounded
        for idx in ids_ky
            if am.setting["Simulation Configuration"]["transmission_expansion_flag"] == true
                JuMP.set_lower_bound(var[idx], 0)
                JuMP.set_upper_bound(var[idx], am.setting["Planning Design"]["transmission_expansion_limit_value"] - 1.0)
            else
                JuMP.fix(var[idx], 0.0; force=true)
            end
        end
    end

    report && add_sol_component(am, nw, category, Symbol(variable_name), ids_ky, var, decomp_group)    
end

function variable_f_kdhty_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, category::Symbol, variable_name::String, decomp_group, ids_k, ids_d, ids_h, ids_t, ids_y; nw::Int=am.cnw, bounded::Bool=false, report::Bool=true, fix_flag::Bool = false, result::Dict{String,<:Any} = Dict{String,Any}())
    
    ids = [(k, d, h, t, y) for (k) in ids_k for (d) in ids_d for (h) in ids_h for (t) in ids_t for y in ids_y]
        
    var = am.var[:nw][nw][decomp_group][Symbol(variable_name)] = JuMP.@variable(JuMP_model,
        [idx in ids],
        base_name=variable_name,
        integer=false,
        binary=false
    )

    if bounded
        for v in var
            JuMP.set_lower_bound(v, 0)
        end
    end

    if fix_flag == true
        for idx in ids
            fix_value = result["solution"]["expansion"][string(idx)]["f_kdhty"]
            JuMP.fix(var[idx], fix_value; force = true)
        end
    end

    report && add_sol_component(am, nw, category, Symbol(variable_name), ids, var, decomp_group)    
end

function variable_f_exp_kdhty_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, category::Symbol, variable_name::String, decomp_group, ids_k, ids_d, ids_h, ids_t, ids_y; nw::Int=am.cnw, report::Bool=true)

    # Expansion-increment flow for the enhanced-hybrid B-theta formulation. Free (signed:
    # counterflow is allowed within the unbuilt headroom); bounded only by the capacity and
    # coupling rows. Registered ONLY for eligible AC corridors (see hybrid_exp_branch).
    ids = [(k, d, h, t, y) for (k) in ids_k for (d) in ids_d for (h) in ids_h for (t) in ids_t for y in ids_y]

    var = am.var[:nw][nw][decomp_group][Symbol(variable_name)] = JuMP.@variable(JuMP_model,
        [idx in ids],
        base_name=variable_name,
        integer=false,
        binary=false
    )

    report && add_sol_component(am, nw, category, Symbol(variable_name), ids, var, decomp_group)
end

function variable_nonspin_idhty_real_EXP(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, category::Symbol, variable_name::String, decomp_group, ids_i, ids_d, ids_h, ids_t, ids_y; nw::Int=am.cnw, bounded::Bool=true, report::Bool=true, fix_flag::Bool = false, result::Dict{String,<:Any} = Dict{String,Any}())
    
    ids = [(i, d, h, t, y) for (i) in ids_i for (d) in ids_d for h in ids_h for t in ids_t for y in ids_y]

    var = am.var[:nw][nw][decomp_group][Symbol(variable_name)] = JuMP.@variable(JuMP_model,
        [idx in ids],
        base_name=variable_name,
        integer=false,
        binary=false
    )

    if bounded
        for idx in ids
            JuMP.set_lower_bound(var[idx], 0)


            # Set to zero if UC is not modeled
            if am.setting["Simulation Configuration"]["Dispatch_Mode_in_EXP"] == "Economic Dispatch"
                JuMP.fix(var[idx], 0.0; force = true)
            end
        end
    end

    if fix_flag == true
        for idx in ids
            fix_value = result["solution"]["reserve"][string(idx)]["nonspin_idhty"]
            JuMP.fix(var[idx], fix_value; force = true)
        end
    end

    report && add_sol_component(am, nw, category, Symbol(variable_name), ids, var, decomp_group)    
end

function variable_nonspin_idhty_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, category::Symbol, variable_name::String, decomp_group, ids_i, ids_d, ids_h, ids_t, ids_y; nw::Int=am.cnw, bounded::Bool=true, report::Bool=true, fix_flag::Bool = false, result::Dict{String,<:Any} = Dict{String,Any}())
    
    ids = [(i, d, h, t, y) for (i) in ids_i for (d) in ids_d for h in ids_h for t in ids_t for y in ids_y]

    var = am.var[:nw][nw][decomp_group][Symbol(variable_name)] = JuMP.@variable(JuMP_model,
        [idx in ids],
        base_name=variable_name,
        integer=false,
        binary=false
    )

    if bounded
        for idx in ids
            JuMP.set_lower_bound(var[idx], 0)


            # Set to zero if UC is not modeled
            if am.setting["Simulation Configuration"]["Dispatch_Mode_in_OP"] == "Economic Dispatch"
                JuMP.fix(var[idx], 0.0; force = true)
            end
        end
    end

    if fix_flag == true
        for idx in ids
            fix_value = result["solution"]["reserve"][string(idx)]["nonspin_idhty"]
            JuMP.fix(var[idx], fix_value; force = true)
        end
    end

    report && add_sol_component(am, nw, category, Symbol(variable_name), ids, var, decomp_group)    
end

function variable_bus_angle_ndhty_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, category::Symbol, variable_name::String, decomp_group, ids_n, ids_d, ids_h, ids_t, ids_y; nw::Int=am.cnw, bounded::Bool=true, report::Bool=true, fix_flag::Bool = false, result::Dict{String,<:Any} = Dict{String,Any}(), ref_bus_ids::Union{Nothing,Set{Int}}=nothing)

    ids = [(n, d, h, t, y) for n in ids_n for d in ids_d for h in ids_h for t in ids_t for y in ids_y]

    var = am.var[:nw][nw][decomp_group][Symbol(variable_name)] = JuMP.@variable(JuMP_model,
        [idx in ids],
        base_name=variable_name,
        integer=false,
        binary=false
    )

    # Anchor one reference bus per synchronous island (angle=0). Default (no ref_bus_ids)
    # keeps the legacy single global reference at bus 1.
    if bounded
        anchor = ref_bus_ids === nothing ? Set([1]) : ref_bus_ids
        for idx in ids

            (n, d, h, t, y) = idx
            if n in anchor
                JuMP.fix(var[idx], 0.0; force = true)
            end
        end
    end

    if fix_flag == true
        for idx in ids
            fix_value = result["solution"]["powerflow"][string(idx)]["bus_angle_ndhty"]
            JuMP.fix(var[idx], fix_value; force = true)
        end
    end

    report && add_sol_component(am, nw, category, Symbol(variable_name), ids, var, decomp_group) 
end

function variable_lfl_ind_lsdhty_binary(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, category::Symbol, variable_name::String, decomp_group, ids_lfl, ids_d, ids_h, ids_t, ids_y; nw::Int=am.cnw, bounded_lower::Bool=false, lower_bound=0, bounded_upper::Bool=false, upper_bound=0, report::Bool=true, fix_flag::Bool = false, result::Dict{String,<:Any} = Dict{String,Any}())
    

    ids = [] 
    for lfl in ids_lfl 
        segment_vector = [s for s in 1:am.ref[:nw][0][:demand][lfl]["Num_DR_Segments"]]
        for s in segment_vector
            if s > 1    # skip the first segment indicator
                append!(ids, [(lfl, s, d, h, t, y) for d in ids_d for h in ids_h for t in ids_t for y in ids_y])
            end
        end
    end

    var = am.var[:nw][nw][decomp_group][Symbol(variable_name)] = JuMP.@variable(JuMP_model,
        [idx in ids],
        base_name=variable_name
    )

    for idx in ids
        if am.ref[:nw][0][:demand][idx[1]]["Integer_Flag"] == true
            JuMP.set_binary(var[idx])
        end
    end
    
    if bounded_lower
        for idx in ids
            JuMP.set_lower_bound(var[idx], lower_bound)
        end
    end

    if bounded_upper
        for idx in ids
            JuMP.set_upper_bound(var[idx], upper_bound)
        end
    end

    if fix_flag == true
        for idx in ids
            fix_value = result["solution"]["demand"][string(idx)]["lfl_ind_lsdhty"]
            JuMP.fix(var[idx], fix_value; force = true)
        end
    end


    report && add_sol_component(am, nw, category, Symbol(variable_name), ids, var, decomp_group)   
end

function variable_lfl_seg_lsdhty_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, category::Symbol, variable_name::String, decomp_group, ids_lfl, ids_d, ids_h, ids_t, ids_y; nw::Int=am.cnw, bounded_lower::Bool=false, lower_bound=0, bounded_upper::Bool=false, upper_bound=0, report::Bool=true, fix_flag::Bool = false, result::Dict{String,<:Any} = Dict{String,Any}())
    

    ids = [] 
    for lfl in ids_lfl 
        segment_vector = [s for s in 1:am.ref[:nw][0][:demand][lfl]["Num_DR_Segments"]]
        for s in segment_vector
            append!(ids, [(lfl, s, d, h, t, y) for d in ids_d for h in ids_h for t in ids_t for y in ids_y])
        end
    end

    var = am.var[:nw][nw][decomp_group][Symbol(variable_name)] = JuMP.@variable(JuMP_model,
        [idx in ids],
        base_name=variable_name,
        integer=false,
        binary=false
    )

    if bounded_lower
        for idx in ids
            JuMP.set_lower_bound(var[idx], lower_bound)
        end
    end

    if bounded_upper
        for idx in ids
            JuMP.set_upper_bound(var[idx], upper_bound)
        end
    end

    if fix_flag == true
        for idx in ids
            fix_value = result["solution"]["demand"][string(idx)]["lfl_seg_lsdhty"]
            JuMP.fix(var[idx], fix_value; force = true)
        end
    end


    report && add_sol_component(am, nw, category, Symbol(variable_name), ids, var, decomp_group)   
end

function define_variable_idhty_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, category::Symbol, variable_name::String, decomp_group, ids_i, ids_d, ids_h, ids_t, ids_y; nw::Int=am.cnw, bounded_lower::Bool=false, lower_bound=0, bounded_upper::Bool=false, upper_bound=0, report::Bool=true)
    
    ids = [(i, d, h, t, y) for (i) in ids_i for d in ids_d for h in ids_h for t in ids_t for y in ids_y]

    var = am.var[:nw][nw][decomp_group][Symbol(variable_name)] = JuMP.@variable(JuMP_model,
        [idx in ids],
        base_name=variable_name,
        integer=false,
        binary=false
    )

    if bounded_lower
        for idx in ids
            JuMP.set_lower_bound(var[idx], lower_bound)
        end
    end

    if bounded_upper
        for idx in ids
            JuMP.set_upper_bound(var[idx], upper_bound)
        end
    end

    report && add_sol_component(am, nw, category, Symbol(variable_name), ids, var, decomp_group)      
end


function variable_cont_deploy_idhty_real(JuMP_model::JuMP.AbstractModel, am::Abstract_ALEAF_Model, category::Symbol, variable_name::String, decomp_group, ids_i, ids_d, ids_h, ids_t, ids_y; nw::Int=am.cnw, bounded::Bool=true, report::Bool=true, fix_flag::Bool = false, result::Dict{String,<:Any} = Dict{String,Any}())
    
    ids = [(i, d, h, t, y) for (i) in ids_i for (d) in ids_d for h in ids_h for t in ids_t for y in ids_y]

    var = am.var[:nw][nw][decomp_group][Symbol(variable_name)] = JuMP.@variable(JuMP_model,
        [idx in ids],
        base_name=variable_name,
        integer=false,
        binary=false
    )

    if bounded
        for idx in ids

            (i, d, h, t, y) = idx

            JuMP.set_lower_bound(var[idx], 0)

            CAP = parameter(am, 0, :gen_index, "CAP", i)
            JuMP.set_upper_bound(var[idx], CAP)
        end
    end

    if fix_flag == true
        for idx in ids
            fix_value = result["solution"]["dispatch"][string(idx)]["g_idhty"]
            JuMP.fix(var[idx], fix_value; force = true)
        end
    end

    report && add_sol_component(am, nw, category, Symbol(variable_name), ids, var, decomp_group)    
end


