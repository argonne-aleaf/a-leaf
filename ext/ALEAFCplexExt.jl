# Copyright © 2026, UChicago Argonne, LLC. All Rights Reserved.
# SPDX-License-Identifier: BSD-3-Clause  (see LICENSE)

module ALEAFCplexExt

using ALEAF
using CPLEX
using JuMP
import MathOptInterface as MOI

# Optional CPLEX support (commercial licence). The extension activates when the CPLEX package is loaded
# (`using CPLEX`), then set solver_name = CPLEX in the setting workbook. The default solver is HiGHS.

ALEAF.cplex_loaded(::Val{:CPLEX}) = true
ALEAF.cplex_optimizer_type(::Val{:CPLEX}) = CPLEX.Optimizer

# A JuMP direct-mode model backed by a fresh CPLEX optimizer (CPLEX misbehaves in auto mode for the
# big expansion/operation models, so the model builders use direct mode when the CPLEX sheet asks for it).
function ALEAF.cplex_direct_model(::Val{:CPLEX}; pass_names::Bool=false)
    cplex_model = CPLEX.Optimizer()
    pass_names && (cplex_model.pass_names = true)
    return JuMP.direct_model(cplex_model)
end

# Reach the underlying CPLEX.Optimizer of a JuMP model, or nothing if not CPLEX-backed.
# Handles both solver modes ALEAF uses:
#   * direct mode:  JuMP.backend(model) is the CPLEX.Optimizer itself.
#   * caching mode (solver_direct_mode_flag=False): the backend is a CachingOptimizer whose inner optimizer is a
#     LazyBridgeOptimizer{CPLEX.Optimizer}. The CPLEX.Optimizer (and its env) already exists after set_optimizer,
#     even while the CachingOptimizer state is EMPTY_OPTIMIZER (before optimize!), and is the same instance MOI copies
#     into and solves, so channel changes made now persist to the solve.
function _cplex_backend(model::JuMP.AbstractModel)
    CPX = CPLEX.Optimizer
    b = try JuMP.backend(model) catch; return nothing end
    b isa CPX && return b
    # Peel a LazyBridgeOptimizer / CachingOptimizer wrapper via its `.optimizer`/`.model` fields
    # (works regardless of attach state, unlike attached_optimizer).
    cur = b
    for _ in 1:6
        cur isa CPX && return cur
        if cur isa MOI.Utilities.CachingOptimizer
            cur = cur.optimizer === nothing ? break : cur.optimizer
        elseif hasproperty(cur, :model)
            cur = getfield(cur, :model)
        elseif hasproperty(cur, :optimizer)
            cur = getfield(cur, :optimizer)
        else
            break
        end
    end
    inner = try MOI.Utilities.attached_optimizer(b) catch; nothing end
    (inner isa CPX) ? inner : nothing
end

# Silence ONLY CPLEX's error message channel (cpxerror) before optimize!.
#
# CPLEX.jl unconditionally calls CPXgetnummipstarts inside MOI.optimize! (MOI_wrapper.jl), which returns
# CPXERR_NOT_MIP (3003) on any continuous (LP) problem and makes CPLEX print a one-line "CPLEX Error 3003: Not a
# mixed-integer problem" to whatever channels are wired to the screen. That message is harmless (the solve still
# returns OPTIMAL) but noisy, and it is emitted by CPLEX itself, not by any parameter ALEAF sets. It cannot be
# suppressed via a solver attribute without CPX_PARAM_SCRIND, which would also hide the barrier/solver log.
# Disconnecting only the cpxerror channel drops the 3003 while leaving the results/log/warning channels (the barrier
# iteration log) intact. Genuine solve failures are still detected via return codes / termination status, so nothing
# is masked.
function ALEAF.cplex_silence_error_channel!(::Val{:CPLEX}, model::JuMP.AbstractModel)
    opt = _cplex_backend(model)
    opt === nothing && return
    try
        env = opt.env
        res = Ref{Ptr{Cvoid}}(); warn = Ref{Ptr{Cvoid}}()
        err = Ref{Ptr{Cvoid}}(); logc = Ref{Ptr{Cvoid}}()
        CPLEX.CPXgetchannels(env, res, warn, err, logc)
        CPLEX.CPXdisconnectchannel(env, err[])
    catch
        # Best-effort only; never let console hygiene break a solve.
    end
    return
end

end # module
