# Copyright © 2026, UChicago Argonne, LLC. All Rights Reserved.
# SPDX-License-Identifier: BSD-3-Clause  (see LICENSE)

# Root ALEAF module. Include order below is intentional — do not reorder.

"Jonghwan Kwon; Argonne National Laboratory; kwonj@anl.gov"

module ALEAF

using JuMP
using MathOptInterface
using LinearAlgebra
const _MOI = MathOptInterface

# Solvers
using HiGHS

# Optional CPLEX (commercial licence) via ext/ALEAFCplexExt.jl: the extension activates when the CPLEX package is loaded
# (`using CPLEX`, before `ALEAF.run_ALEAF`). Default solver is HiGHS; these fallbacks are used when CPLEX is not loaded.
_cplex_missing() = error("solver_name = CPLEX needs the CPLEX package loaded: `using CPLEX` " *
    "(requires IBM CPLEX Studio and CPLEX_STUDIO_BINARIES). execute_ALEAF.jl loads it automatically once the package is " *
    "installed (`Pkg.add(\"CPLEX\")`). Or set solver_name = HiGHS.")
cplex_optimizer_type(v::Val) = _cplex_missing()
cplex_direct_model(v::Val; pass_names::Bool=false) = _cplex_missing()
cplex_silence_error_channel!(v::Val, model) = nothing
cplex_loaded(v::Val) = false

# Optional GPU solvers (cuOpt, MadNLP/ExaModels) via ext/ package extensions;
# load the matching package(s) before `using ALEAF`.
gpu_optimizer(v::Val) = error("Selected solver needs its GPU package loaded first: " *
    "`using cuOpt` for cuOpt, or `using MadNLP, MadNLPGPU, CUDA, CUDSS, ExaModels` for MadNLP, then `using ALEAF`.")
# 2-arg form (MadNLP has constructor-argument options); fallback errors if package absent.
gpu_optimizer(v::Val, solver_setting) = gpu_optimizer(v)

gpu_setting_optimizer(v::Val) = gpu_optimizer(v)

# Device memory a finished GPU solve still holds; the ext for the loaded solver returns it.
# No-op for CPU solvers and whenever the GPU package is absent.
gpu_reclaim!(v::Val) = nothing

# (free, total) device bytes, or nothing when no GPU solver is loaded.
gpu_memory_info(v::Val) = nothing

# Load order matters: core runtime/utilities, then shared components, then models.
include("io/logging.jl")
include("core/model_type.jl")
include("io/data_io.jl")
include("core/base_functions.jl")
include("network/ptdf_reduction.jl")
include("network/generate_network.jl")
include("io/export_network.jl")
include("runtime/model_handler.jl")
include("solvers/solve_and_collect.jl")

# Shared components used by multiple models
include("component/variables.jl")
include("component/constraints.jl")
    
# Model workflows
include("model/LCO_GTEP.jl")
include("model/RA.jl")

# Supporting utilities
include("util/outage_scenarios_model.jl")
include("util/scenario_reduction_GTEP.jl")

end


