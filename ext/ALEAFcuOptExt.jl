# Copyright © 2026, UChicago Argonne, LLC. All Rights Reserved.
# SPDX-License-Identifier: BSD-3-Clause  (see LICENSE)

module ALEAFcuOptExt

using ALEAF
using cuOpt
import MathOptInterface as MOI

ALEAF.gpu_optimizer(::Val{:cuOpt}) = cuOpt.Optimizer

# Found through whichever CUDA the session already loaded, so the ext trigger stays `cuOpt` alone
# and a missing binding degrades to a no-op rather than an error mid-run.
function _loaded_cuda()
    isdefined(cuOpt, :CUDA) && return getfield(cuOpt, :CUDA)
    for (pkg, mod) in Base.loaded_modules
        pkg.name == "CUDA" && return mod
    end
    return nothing
end

function ALEAF.gpu_reclaim!(::Val{:cuOpt})
    cuda = _loaded_cuda()
    (cuda === nothing || !isdefined(cuda, :reclaim)) && return nothing
    return Base.invokelatest(getfield(cuda, :reclaim))
end

# nvidia-smi is the fallback because cuOpt.jl links libcuopt directly; CUDA.jl need not be installed,
# and without it neither the reclaim nor a device query is reachable through Julia.
function _nvidia_smi_memory()
    out = try
        read(`nvidia-smi --query-gpu=memory.free,memory.total --format=csv,noheader,nounits`, String)
    catch
        return nothing
    end
    fields = split(strip(first(split(strip(out), '\n'))), ',')
    length(fields) == 2 || return nothing
    free_mib = tryparse(Float64, strip(fields[1]))
    total_mib = tryparse(Float64, strip(fields[2]))
    (free_mib === nothing || total_mib === nothing) && return nothing
    return (free_mib * 2^20, total_mib * 2^20)
end

function ALEAF.gpu_memory_info(::Val{:cuOpt})
    cuda = _loaded_cuda()
    if cuda !== nothing && isdefined(cuda, :memory_info)
        info = try Base.invokelatest(getfield(cuda, :memory_info)) catch; nothing end
        info === nothing || return info
    end
    return _nvidia_smi_memory()
end

# cuOpt.jl computes LP duals (Optimizer.dual_solution) but ships no MOI.ConstraintDual
# getter, so JuMP.dual() throws and ALEAF's collectors record 0.0. Surface the stored
# duals; sign convention mirrors HiGHS.jl, the wrapper cuOpt.jl was forked from.
_cuopt_sense_mult(m::cuOpt.Optimizer) = m.objective_sense == MOI.MAX_SENSE ? -1.0 : 1.0

# PDLP returns an interior (non-basic) point, so clamp inequality duals to the MOI cone
# sign; equality/interval duals are free-signed (this is HiGHS.jl's no-basis branch).
_cuopt_cone_dual(d::Float64, ::Type{<:MOI.LessThan}) = min(d, 0.0)
_cuopt_cone_dual(d::Float64, ::Type{<:MOI.GreaterThan}) = max(d, 0.0)
_cuopt_cone_dual(d::Float64, ::Type) = d

# cuOpt.jl applies no MOI sign convention to its raw row duals; the HiGHS-fork convention
# is replicated above/below. The one bit unverifiable without a GPU run: if a Swing
# calibration vs a CPLEX/HiGHS reference shows duals uniformly flipped, set this to -1.0.
const _RAW_DUAL_SIGN = 1.0

function MOI.get(
    model::cuOpt.Optimizer,
    attr::MOI.ConstraintDual,
    ci::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},S},
) where {S}
    attr.result_index == 1 || throw(MOI.ResultIndexBoundsError(attr, 1))
    info = model.affine_constraint_info[cuOpt._ConstraintKey(ci.value)]
    raw = _RAW_DUAL_SIGN * model.dual_solution[info.row + 1]   # info.row is 0-based
    return _cuopt_cone_dual(_cuopt_sense_mult(model) * raw, S)
end

# The getter reads cuOpt internals (dual_solution / affine_constraint_info /
# _ConstraintInfo.row); warn loudly if a cuOpt upgrade moves them so duals don't
# silently regress to 0 again. Keep the cuOpt compat bound in Project.toml in sync.
function __init__()
    ok = isdefined(cuOpt, :_ConstraintKey) && isdefined(cuOpt, :_ConstraintInfo) &&
        all(f -> f in fieldnames(cuOpt.Optimizer),
            (:dual_solution, :affine_constraint_info, :objective_sense)) &&
        :row in fieldnames(cuOpt._ConstraintInfo)
    ok || @warn "ALEAFcuOptExt: cuOpt internals changed; the ConstraintDual glue may be " *
                "broken (duals could regress to 0). Re-verify against jump-dev/cuOpt.jl."
end

end
