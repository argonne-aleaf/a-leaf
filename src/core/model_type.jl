# Copyright © 2026, UChicago Argonne, LLC. All Rights Reserved.
# SPDX-License-Identifier: BSD-3-Clause  (see LICENSE)

# ALEAF model type hierarchy and shared model-instance structs.

"Jonghwan Kwon; Argonne National Laboratory; kwonj@anl.gov"

##### Top Level Abstract Types #####

"Root of the ALEAF model formulation type hierarchy"
abstract type Abstract_ALEAF_Model end

"individual model type"
abstract type Abstract_LC_GTEP_Model <: Abstract_ALEAF_Model end

"constructor for ALEAF_Model_Structure_with_PH"
mutable struct ALEAF_Model_Structure_PH <: Abstract_ALEAF_Model
    model::Dict{Symbol,<:Any}
    model_type::String

    data::Dict{String,<:Any}
    setting::Dict{String,<:Any}
    solution::Dict{String,<:Any}

    ref::Dict{Symbol,<:Any}
    var::Dict{Symbol,<:Any}
    con::Dict{Symbol,<:Any}

    sol::Dict{Symbol,<:Any}

    cnw::Int

    PH::Dict{Symbol,<:Any}
end

"constructor for ALEAF_Model_Structure_for_RA"
mutable struct ALEAF_Model_Structure_RA <: Abstract_ALEAF_Model
    model::Dict{Symbol,<:Any}
    model_type::String

    setting::Dict{String,<:Any}
    solution::Dict{String,<:Any}

    ref::Dict{Symbol,<:Any}
    var::Dict{Symbol,<:Any}
    con::Dict{Symbol,<:Any}

    sol::Dict{Symbol,<:Any}

    cnw::Int
end
1
