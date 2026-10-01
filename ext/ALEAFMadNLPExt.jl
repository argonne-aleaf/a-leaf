# Copyright © 2026, UChicago Argonne, LLC. All Rights Reserved.
# SPDX-License-Identifier: BSD-3-Clause  (see LICENSE)

module ALEAFMadNLPExt

using ALEAF
using ExaModels, MadNLP, MadNLPGPU, CUDA, CUDSS

# cuOpt's options are MOI attributes (set after the solver is attached); MadNLP's
# options are ExaModels.Optimizer *constructor* arguments, so they cannot go through
# the same attribute mechanism. The solve path reads them from the "MadNLP Setting"
# sheet here, at construction time. See the GPU Solvers docs page (configuration/GPU_Solvers.md).

# kkt_system is a Julia type, so the sheet stores its name and we map it here.
const MADNLP_KKT_SYSTEMS = Dict(
    "SparseCondensed" => MadNLP.SparseCondensedKKTSystem,   # GPU default (pairs with cuDSS)
    "DenseCondensed"  => MadNLP.DenseCondensedKKTSystem,
    "Sparse"          => MadNLP.SparseKKTSystem,
    "ScaledSparse"    => MadNLP.ScaledSparseKKTSystem,
    "SparseUnreduced" => MadNLP.SparseUnreducedKKTSystem,
    "Dense"           => MadNLP.DenseKKTSystem,
)

# Read one MadNLP option from the Solver Setting dict by its Parameter name.
# Uses the sheet Value only when that row's Flag is true; otherwise the default.
function _madnlp_opt(solver_setting, name, default)
    for (_, row) in solver_setting
        if row isa AbstractDict && get(row, "Parameter", nothing) == name
            return get(row, "Flag", false) == true ? row["Value"] : default
        end
    end
    return default
end

# Solve path (solve_and_collect.jl): GPU condensed-KKT MadNLP via ExaModels,
# with kkt_system and tol taken from the MadNLP Setting sheet.
function ALEAF.gpu_optimizer(::Val{:MadNLP}, solver_setting)
    kkt_name = _madnlp_opt(solver_setting, "kkt_system", "SparseCondensed")
    haskey(MADNLP_KKT_SYSTEMS, kkt_name) ||
        error("Unknown MadNLP kkt_system \"$kkt_name\". Valid: $(sort(collect(keys(MADNLP_KKT_SYSTEMS))))")
    kkt_type = MADNLP_KKT_SYSTEMS[kkt_name]
    tol = _madnlp_opt(solver_setting, "tol", 1e-6)
    return () -> ExaModels.Optimizer(madnlp, CUDABackend(); kkt_system = kkt_type, tol = tol)
end

# Setting path (LCO_GTEP.jl): default MadNLP optimizer; overridden by the solve path above.
ALEAF.gpu_setting_optimizer(::Val{:MadNLP}) = () -> ExaModels.MadNLPOptimizer()

end
