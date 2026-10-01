# Copyright © 2026, UChicago Argonne, LLC. All Rights Reserved.
# SPDX-License-Identifier: BSD-3-Clause  (see LICENSE)

# ============================================================================================
# ptdf_reduction.jl — live fine->zonal (region) network reduction for aggregated AC corridors.
#
# Corridor susceptances are estimated from synthetic operating snapshots solved on the FINE network
# (reactances + bus peak loads + plant capacities only; no dispatch, no time series), then refined so
# the REDUCED network reproduces the fine network's cross-border flows. Lines kept at native resolution
# (mixed-resolution runs) are fixed edges of the reduced network.
#
# CORE       : reduce_network_by_snapshots -> corridor br_x_pu for the aggregated network
# PTDF utils : reduced_network_ptdf / attach_branch_ptdf! -> PTDF of the (final) region network
#              fine_to_zonal_ptdf                          -> zonal PTDF of the FINE network
#
# All identifiers are Strings. Each connected component of a region graph is grounded at its
# lowest-index node; exported quantities are gauge-invariant for balanced injections.
# ============================================================================================

using LinearAlgebra
using SparseArrays

Base.@kwdef struct NetworkReductionOpts
    n_snapshots::Int = 250          # training snapshots
    n_holdout::Int = 50             # held-out snapshots (diagnostics and the per-island safety check)
    snapshot_chunk::Int = 50        # fine DC solves per block (bounds memory on large islands)
    refine_iters::Int = 30
    refine_tol::Float64 = 1e-6      # stop when max |Δ log b| < tol
    seed::UInt64 = 0x00a1eaf5eed00001
    load_range::Tuple{Float64,Float64} = (0.5, 1.0)     # load level (x bus peak load)
    wind_range::Tuple{Float64,Float64} = (0.0, 0.95)    # wind output (x capacity)
    solar_range::Tuple{Float64,Float64} = (0.0, 0.85)   # solar output (x capacity)
    avail_range::Tuple{Float64,Float64} = (0.3, 1.0)    # per-bus dispatchable availability
    bnd::Float64 = 50.0             # b kept in [b0/bnd, b0*bnd] (numerical guard)
    x_sigdigits::Int = 10           # rounding makes x identical across BLAS thread counts
    min_island_regions::Int = 2
    init::Symbol = :stage_a         # :stage_a (closed-form start) or :parallel (diagnostics only)
end

# a single aggregated AC corridor: region pair (ca<=cb) + its crossing fine legs (fine_f, fine_t, b=1/|x|)
struct PTDFCorridor
    ca::String
    cb::String
    members::Vector{Tuple{String,String,Float64}}
end

# a line kept at native resolution (both ends at the finest level): fixed edge of the reduced network
struct KeptLine
    key::String
    fine_f::String
    fine_t::String
    b::Float64
end

struct NetworkReductionResult
    x_by_pair::Dict{Tuple{String,String},Float64}
    n_islands::Int
    n_corridors::Int
    n_kept::Int
    n_fallback::Int           # corridors whose slope stayed non-positive (previous b kept)
    n_clamped::Int            # corridors at the [b0/bnd, b0*bnd] bound
    n_island_fallback::Int    # islands reverted to parallel-combine by the held-out safety check
    err_parallel::Float64     # held-out cross-border flow error, parallel-combine x
    err_final::Float64        # held-out cross-border flow error, estimated x
    err_kept::Float64         # held-out flow error on kept lines
    leakage::Float64          # share of zone injection carried by fine lines that are not reduced edges
    coverage::Float64         # share of injection data (load + capacity) on reduced-region buses
    max_components::Int       # largest number of connected components in a reduced island graph
    checksum::UInt64
    elapsed::Float64
end

function _uf_find!(parent::Dict{String,String}, x::String)
    while parent[x] != x
        parent[x] = parent[parent[x]]
        x = parent[x]
    end
    return x
end
function _uf_union!(parent::Dict{String,String}, a::String, b::String)
    ra = _uf_find!(parent, a); rb = _uf_find!(parent, b)
    ra == rb && return
    ra < rb ? (parent[rb] = ra) : (parent[ra] = rb)
end

# counter-based generator keyed by bus id: the same draws in every process, platform and Julia version,
# independent of bus ordering (parallel workers rebuild the network; recipients re-aggregate locally)
@inline function _splitmix64(x::UInt64)
    z = x + 0x9e3779b97f4a7c15
    z = (z ⊻ (z >> 30)) * 0xbf58476d1ce4e5b9
    z = (z ⊻ (z >> 27)) * 0x94d049bb133111eb
    return z ⊻ (z >> 31)
end
function _fnv1a64(s::AbstractString)
    h = 0xcbf29ce484222325
    for c in codeunits(s)
        h = (h ⊻ UInt64(c)) * 0x00000100000001b3
    end
    return h
end
@inline _u01(seed::UInt64, s::Int, key::UInt64) =
    Float64(_splitmix64(seed ⊻ _splitmix64(UInt64(s) ⊻ _splitmix64(key))) >> 11) * 2.0^-53
@inline _draw(r::Tuple{Float64,Float64}, u::Float64) = r[1] + (r[2] - r[1]) * u

# MW value from an xlsx cell (Int / Float / numeric String); missing, Bool and text -> 0
function _mw(v)
    (v === missing || v === nothing || v isa Bool) && return 0.0
    v isa Real && return isfinite(Float64(v)) ? Float64(v) : 0.0
    x = tryparse(Float64, strip(string(v)))
    return (x === nothing || !isfinite(x)) ? 0.0 : x
end

# unit class from UNIT_CATEGORY: 2 = wind, 3 = solar, 4 = dispatchable, 0 = storage (no net injection)
function _unit_class(cat::AbstractString)
    c = uppercase(strip(cat))
    c in ("WIND_ONS", "WIND_OFS") && return 2
    c in ("PV", "CSP") && return 3
    c == "STORAGE" && return 0
    return 4
end

"""
    fine_injection_data(load, plants) -> Dict(bus => (load, wind cap, solar cap, dispatchable cap))

`load`: fine bus -> peak MW load. `plants`: (fine bus, UNIT_CATEGORY, CAP MW) per unit. Storage is excluded.
"""
function fine_injection_data(load::AbstractDict, plants::AbstractVector)
    acc = Dict{String,Vector{Float64}}()
    for (b, v) in load
        a = get!(acc, string(b), zeros(4)); a[1] += max(_mw(v), 0.0)
    end
    for (b, cat, cap) in plants
        k = _unit_class(string(cat)); k == 0 && continue
        a = get!(acc, string(b), zeros(4)); a[k] += max(_mw(cap), 0.0)
    end
    return Dict{String,NTuple{4,Float64}}(k => (v[1], v[2], v[3], v[4]) for (k, v) in acc)
end

# connected components of a region graph; component ids follow the lowest node index
function _components(nreg::Int, pairs::AbstractVector{Tuple{Int,Int}})
    parent = collect(1:nreg)
    find(x) = (while parent[x] != x; parent[x] = parent[parent[x]]; x = parent[x]; end; x)
    for (i, j) in pairs
        ri = find(i); rj = find(j)
        ri == rj && continue
        ri < rj ? (parent[rj] = ri) : (parent[ri] = rj)
    end
    comp = zeros(Int, nreg); id = Dict{Int,Int}()
    for i in 1:nreg
        comp[i] = get!(id, find(i), length(id) + 1)
    end
    return comp, length(id)
end

# angles of a region network for injections R (nodes x columns). Each connected component is grounded
# at its lowest-index node; with `balance`, a component's net injection (flow on lines that are not
# edges of this network) is removed uniformly first so every block stays consistent.
function _grounded_solve(pairs::AbstractVector{Tuple{Int,Int}}, bvec::AbstractVector{Float64}, nreg::Int,
                         R::AbstractMatrix{Float64}; balance::Bool = true)
    comp, ncomp = _components(nreg, pairs)
    isref = falses(nreg); seen = falses(ncomp)
    for i in 1:nreg
        seen[comp[i]] || (isref[i] = true; seen[comp[i]] = true)
    end
    Rb = Matrix{Float64}(R)
    if balance && ncomp > 1
        for c in 1:ncomp
            mem = findall(==(c), comp)
            @views Rb[mem, :] .-= sum(Rb[mem, :], dims = 1) ./ length(mem)
        end
    end
    Θ = zeros(nreg, size(R, 2))
    kc = findall(!, isref)
    isempty(kc) && return Θ, ncomp
    I = Int[]; J = Int[]; V = Float64[]
    for (e, (i, j)) in enumerate(pairs)
        b = bvec[e]
        push!(I, i, j, i, j); push!(J, i, j, j, i); push!(V, b, b, -b, -b)
    end
    Lk = sparse(I, J, V, nreg, nreg)[kc, kc]
    fac = try
        cholesky(Symmetric(Lk))
    catch
        lu(Lk)
    end
    Θ[kc, :] = fac \ Rb[kc, :]
    return Θ, ncomp
end

"""
    reduced_network_ptdf(pairs, bvec, nnode) -> Matrix (edges x nnode)

PTDF of a network given as an edge list `pairs[e] = (i, j)` with susceptances `bvec[e]` (parallel edges
allowed). Column r = edge flows for 1 MW injected at node r and withdrawn at the grounded (lowest-index)
node of r's connected component, whose own column is 0. Rows are oriented i -> j.
"""
function reduced_network_ptdf(pairs::AbstractVector{Tuple{Int,Int}}, bvec::AbstractVector{Float64}, nnode::Int)
    Θ, _ = _grounded_solve(pairs, bvec, nnode, Matrix{Float64}(I, nnode, nnode); balance = false)
    F = zeros(length(pairs), nnode)
    for (e, (i, j)) in enumerate(pairs)
        @views F[e, :] .= bvec[e] .* (Θ[i, :] .- Θ[j, :])
    end
    return F
end

"""
    attach_branch_ptdf!(branches, nbus; max_entries) -> number of connected components

Set `branch["ptdf"]` (length-`nbus` vector indexed by model bus key) for every branch from the final
branch list (f_bus -> t_bus, br_x_pu), so PTDF mode is identical to B-theta on the same network. DC ties
get zero rows (not supported in PTDF mode).
"""
function attach_branch_ptdf!(branches::AbstractDict, nbus::Int; max_entries::Int = 50_000_000)
    ks = sort!(collect(keys(branches)), by = k -> parse(Int, string(k)))
    ac = [k for k in ks if get(branches[k], "dc_line", false) != true]
    length(ac) * nbus > max_entries &&
        error("PTDF mode: $(length(ac)) AC branches x $nbus buses exceeds the dense PTDF limit ($max_entries); use B-theta at this resolution")
    pairs = Tuple{Int,Int}[(Int(branches[k]["f_bus"]), Int(branches[k]["t_bus"])) for k in ac]
    bvec = Float64[1.0 / abs(Float64(branches[k]["br_x_pu"])) for k in ac]
    F = reduced_network_ptdf(pairs, bvec, nbus)
    for (e, k) in enumerate(ac); branches[k]["ptdf"] = F[e, :]; end
    for k in ks
        get(branches[k], "dc_line", false) == true && (branches[k]["ptdf"] = zeros(nbus))
    end
    return _components(nbus, pairs)[2]
end

struct _Island
    buses::Vector{String}
    nidx::Dict{String,Int}
    B::SparseMatrixCSC{Float64,Int}
    keep::Vector{Int}
    regions::Vector{String}          # regions touched by a corridor or kept line
    ridx::Dict{String,Int}
    region_of_bus::Vector{Int}       # 0 = passive bus (its region is not a reduced node)
    corr::Vector{PTDFCorridor}
    kept::Vector{KeptLine}
end

function _partition_islands(nodal_ac::Vector{Tuple{String,String,Float64}}, bus2region::Dict{String,String},
                            corridors::Vector{PTDFCorridor}, kept::Vector{KeptLine}, min_regions::Int)
    busset = Set{String}()
    for (f, t, _x) in nodal_ac
        (haskey(bus2region, f) && haskey(bus2region, t)) || continue
        push!(busset, f); push!(busset, t)
    end
    parent = Dict{String,String}(b => b for b in busset)
    for (f, t, _x) in nodal_ac
        (haskey(parent, f) && haskey(parent, t)) || continue
        _uf_union!(parent, f, t)
    end
    root_of = Dict{String,String}(b => _uf_find!(parent, b) for b in busset)
    isl_buses = Dict{String,Vector{String}}()
    for b in busset; push!(get!(isl_buses, root_of[b], String[]), b); end
    # a corridor goes to the island holding most of its susceptance (legs can straddle islands)
    isl_corr = Dict{String,Vector{PTDFCorridor}}()
    for cr in sort(corridors, by = c -> (c.ca, c.cb))
        wt = Dict{String,Float64}()
        for (f, _t, bs) in cr.members
            r = get(root_of, f, ""); r == "" && continue
            wt[r] = get(wt, r, 0.0) + bs
        end
        isempty(wt) && continue
        r = first(sort(collect(wt), by = x -> (-x[2], x[1])))[1]
        push!(get!(isl_corr, r, PTDFCorridor[]), cr)
    end
    isl_kept = Dict{String,Vector{KeptLine}}()
    for kl in sort(kept, by = k -> k.key)
        r = get(root_of, kl.fine_f, "")
        (r == "" || get(root_of, kl.fine_t, "") != r) && continue
        bus2region[kl.fine_f] == bus2region[kl.fine_t] && continue
        push!(get!(isl_kept, r, KeptLine[]), kl)
    end

    islands = _Island[]
    for root in sort!(collect(keys(isl_buses)))
        corr = get(isl_corr, root, PTDFCorridor[]); kl = get(isl_kept, root, KeptLine[])
        (isempty(corr) && isempty(kl)) && continue
        buses = sort!(isl_buses[root])
        # reduced nodes = regions touched by a corridor or kept line; other buses stay passive (e.g. the
        # unmatched placeholder "0")
        rset = Set{String}()
        for cr in corr; push!(rset, cr.ca); push!(rset, cr.cb); end
        for k in kl; push!(rset, bus2region[k.fine_f]); push!(rset, bus2region[k.fine_t]); end
        present = Set{String}(bus2region[b] for b in buses)
        regions = sort!([r for r in rset if r in present])
        length(regions) < min_regions && continue
        ridx = Dict{String,Int}(r => i for (i, r) in enumerate(regions))
        # corridors whose regions both sit in this island (a region can span islands)
        corr = [cr for cr in corr if haskey(ridx, cr.ca) && haskey(ridx, cr.cb)]
        nidx = Dict{String,Int}(b => i for (i, b) in enumerate(buses)); n = length(buses)
        I = Int[]; J = Int[]; V = Float64[]
        for (f, t, x) in nodal_ac
            (get(root_of, f, "") == root && get(root_of, t, "") == root) || continue
            i = nidx[f]; j = nidx[t]; b = 1.0 / abs(x)
            push!(I, i, j, i, j); push!(J, i, j, j, i); push!(V, b, b, -b, -b)
        end
        B = sparse(I, J, V, n, n)
        keep = collect(2:n)                    # fine gauge: the lexicographically-min bus
        region_of_bus = [get(ridx, bus2region[b], 0) for b in buses]
        push!(islands, _Island(buses, nidx, B, keep, regions, ridx, region_of_bus, corr, kl))
    end
    return islands
end

"""
    fine_to_zonal_ptdf(nodal_ac, bus2region, corridors; bus_weight=nothing) -> Dict(pair => Dict(region => PTDF))

Zonal PTDF of the FINE network: 1 MW spread over a region's buses in proportion to `bus_weight`
(uniform when `nothing` or when a region's weights sum to 0), withdrawn the same way from the island's
lowest region (whose column is 0). Rows are the summed flows over each corridor's legs, oriented ca -> cb.
"""
function fine_to_zonal_ptdf(nodal_ac::Vector{Tuple{String,String,Float64}}, bus2region::Dict{String,String},
                            corridors::Vector{PTDFCorridor}; bus_weight::Union{Nothing,AbstractDict}=nothing,
                            min_island_regions::Int = 2)
    out = Dict{Tuple{String,String},Dict{String,Float64}}()
    for isl in _partition_islands(nodal_ac, bus2region, corridors, KeptLine[], min_island_regions)
        n = length(isl.buses); nreg = length(isl.regions)
        w = [bus_weight === nothing ? 1.0 : max(_mw(get(bus_weight, b, 0.0)), 0.0) for b in isl.buses]
        W = _region_weights(isl, w)
        RHS = zeros(n, nreg)
        for j in 1:n
            r = isl.region_of_bus[j]; r == 0 && continue
            RHS[j, r] += W[j]
            RHS[j, :] .-= r == 1 ? W[j] : 0.0
        end
        RHS[:, 1] .= 0.0
        Θ = zeros(n, nreg)
        Θ[isl.keep, :] = _factor(isl.B[isl.keep, isl.keep]) \ RHS[isl.keep, :]
        for cr in isl.corr
            row = zeros(nreg)
            for (f, t, bs) in cr.members
                (haskey(isl.nidx, f) && haskey(isl.nidx, t)) || continue
                s = bus2region[f] == cr.ca ? bs : -bs
                @views row .+= s .* (Θ[isl.nidx[f], :] .- Θ[isl.nidx[t], :])
            end
            out[(cr.ca, cr.cb)] = Dict(isl.regions[r] => row[r] for r in 1:nreg)
        end
    end
    return out
end

# SPD sparse factorization (Cholesky; LU if the matrix is not numerically positive definite)
_factor(A::SparseMatrixCSC{Float64,Int}) = try
    cholesky(Symmetric(A))
catch
    lu(A)
end

# per-bus share of its region (weights w; uniform when a region's weights sum to 0); 0 for passive buses
function _region_weights(isl::_Island, w::Vector{Float64})
    nreg = length(isl.regions); wsum = zeros(nreg); cnt = zeros(nreg)
    for (j, r) in enumerate(isl.region_of_bus)
        r == 0 && continue
        wsum[r] += w[j]; cnt[r] += 1
    end
    return [r == 0 ? 0.0 : (wsum[r] > 0 ? w[j] / wsum[r] : 1.0 / cnt[r]) for (j, r) in enumerate(isl.region_of_bus)]
end

# synthetic injections (buses x snapshots) for snapshot ids `ids`: load at `load_range` of peak, wind /
# solar at a uniform output level, dispatchable capacity at a per-bus availability covers the residual in
# proportion to available capacity; each column is balanced
function _snapshot_block(ids, keys_b::Vector{UInt64}, load::Vector{Float64}, capW::Vector{Float64},
                         capS::Vector{Float64}, capD::Vector{Float64}, opts::NetworkReductionOpts)
    n = length(load); P = zeros(n, length(ids))
    Ltot = sum(load); Wtot = sum(capW); Stot = sum(capS)
    for (c, s) in enumerate(ids)
        lf = _draw(opts.load_range, _u01(opts.seed, s, UInt64(1)))
        cw = _draw(opts.wind_range, _u01(opts.seed, s, UInt64(2)))
        cs = _draw(opts.solar_range, _u01(opts.seed, s, UInt64(3)))
        D = lf * Ltot; Vre = cw * Wtot + cs * Stot
        vscale = (Vre > D && Vre > 0) ? D / Vre : 1.0
        resid = D - Vre * vscale
        A = 0.0
        @inbounds for j in 1:n
            if capD[j] > 0
                a = capD[j] * _draw(opts.avail_range, _u01(opts.seed, s, keys_b[j]))
                P[j, c] = a; A += a
            end
        end
        dscale = A > 0 ? resid / A : 0.0
        lscale = lf
        if A <= 0 && resid > 0                 # nothing dispatchable: shrink load to the available VRE
            lscale = Ltot > 0 ? (Vre * vscale) / Ltot : 0.0
        end
        @inbounds for j in 1:n
            P[j, c] = vscale * (cw * capW[j] + cs * capS[j]) + dscale * P[j, c] - lscale * load[j]
        end
    end
    return P
end

# squared flow error of edges `pairs` (susceptances b) at angles Θ against reference flows f
function _edge_sq_err(pairs::AbstractVector{Tuple{Int,Int}}, b::AbstractVector{Float64},
                      Θ::AbstractMatrix{Float64}, f::AbstractMatrix{Float64})
    se = 0.0
    for (e, (i, j)) in enumerate(pairs)
        @inbounds for c in 1:size(f, 2)
            se += (b[e] * (Θ[i, c] - Θ[j, c]) - f[e, c])^2
        end
    end
    return se
end

"""
    reduce_network_by_snapshots(nodal_ac, bus2region, corridors, kept, bus_inj; opts, verbose)

Estimate one susceptance per aggregated corridor so the reduced (region) network reproduces the fine
network's cross-border flows over synthetic snapshots built from `bus_inj` (see `fine_injection_data`).

- `nodal_ac`: fine AC lines (fine_f, fine_t, x); `bus2region`: fine bus -> region id.
- `corridors`: aggregated corridors (estimated); `kept`: native-resolution lines (fixed edges).

Stage A: b = Σ f·Δθ / Σ Δθ² with zone angles = weighted mean of member-bus angles (weights = load +
capacity). Refinement: the same slope against the REDUCED network's own angles, damped b ← √(b·b_new).
The candidate with the lowest training flow error is kept; an island falls back to parallel-combine if
its held-out error is not better.
"""
function reduce_network_by_snapshots(nodal_ac::Vector{Tuple{String,String,Float64}},
                                     bus2region::Dict{String,String},
                                     corridors::Vector{PTDFCorridor},
                                     kept::Vector{KeptLine},
                                     bus_inj::Dict{String,NTuple{4,Float64}};
                                     opts::NetworkReductionOpts = NetworkReductionOpts(),
                                     verbose::Bool = false)
    t_start = time()
    x_by_pair = Dict{Tuple{String,String},Float64}()
    n_isl = 0; n_corr = 0; n_kept = 0; n_fb = 0; n_clamp = 0; n_isl_fb = 0; max_comp = 0
    se_p = 0.0; se_f = 0.0; sr = 0.0; se_k = 0.0; sr_k = 0.0; lk_num = 0.0; lk_den = 0.0
    inj_total = sum(v[1] + v[2] + v[3] + v[4] for (b, v) in bus_inj if haskey(bus2region, b); init = 0.0)
    inj_used = 0.0
    xr(b) = round(1.0 / b, sigdigits = opts.x_sigdigits)

    for isl in _partition_islands(nodal_ac, bus2region, corridors, kept, opts.min_island_regions)
        n_isl += 1
        corr = isl.corr; kl = isl.kept; ncorr = length(corr); nk = length(kl)
        nreg = length(isl.regions); n = length(isl.buses)
        n_corr += ncorr; n_kept += nk
        b0 = [clamp(sum(m[3] for m in cr.members), 1e-6, 1e6) for cr in corr]

        local Fk
        try
            Fk = _factor(isl.B[isl.keep, isl.keep])
        catch err
            verbose && @warn "network reduction: island factorization failed ($err) -> parallel-combine"
            for (e, cr) in enumerate(corr); x_by_pair[(cr.ca, cr.cb)] = xr(b0[e]); end
            n_isl_fb += 1
            continue
        end

        # injection data, masked to buses of reduced regions (passive buses inject nothing)
        load = zeros(n); capW = zeros(n); capS = zeros(n); capD = zeros(n)
        for (j, b) in enumerate(isl.buses)
            isl.region_of_bus[j] == 0 && continue
            v = get(bus_inj, b, nothing); v === nothing && continue
            load[j] = max(v[1], 0.0); capW[j] = max(v[2], 0.0); capS[j] = max(v[3], 0.0); capD[j] = max(v[4], 0.0)
        end
        inj_used += sum(load) + sum(capW) + sum(capS) + sum(capD)
        wz = _region_weights(isl, load .+ capW .+ capS .+ capD)
        keys_b = UInt64[_fnv1a64(b) for b in isl.buses]

        pairs  = Tuple{Int,Int}[(isl.ridx[cr.ca], isl.ridx[cr.cb]) for cr in corr]
        kpairs = Tuple{Int,Int}[(isl.ridx[bus2region[k.fine_f]], isl.ridx[bus2region[k.fine_t]]) for k in kl]
        kb = Float64[k.b for k in kl]
        allpairs = vcat(pairs, kpairs)
        max_comp = max(max_comp, _components(nreg, allpairs)[2])

        S = opts.n_snapshots; H = opts.n_holdout
        ftr = zeros(ncorr, S); Ztr = zeros(nreg, S); num = zeros(ncorr); den = zeros(ncorr)
        fho = zeros(ncorr, H); Zho = zeros(nreg, H); kho = zeros(nk, H)

        # fine DC power flow over snapshot blocks: cross-border flows, zone injections, stage-A sums
        function run_block!(ids, cols, fout, Zout, kout, stage_a::Bool)
            P = _snapshot_block(ids, keys_b, load, capW, capS, capD, opts)
            Θ = zeros(n, length(ids))
            Θ[isl.keep, :] = Fk \ P[isl.keep, :]
            for (e, cr) in enumerate(corr)
                for (f, t, bs) in cr.members
                    (haskey(isl.nidx, f) && haskey(isl.nidx, t)) || continue
                    s = bus2region[f] == cr.ca ? bs : -bs
                    i = isl.nidx[f]; j = isl.nidx[t]
                    @inbounds for (c, col) in enumerate(cols); fout[e, col] += s * (Θ[i, c] - Θ[j, c]); end
                end
            end
            if kout !== nothing
                for (e, k) in enumerate(kl)
                    i = isl.nidx[k.fine_f]; j = isl.nidx[k.fine_t]
                    @inbounds for (c, col) in enumerate(cols); kout[e, col] = k.b * (Θ[i, c] - Θ[j, c]); end
                end
            end
            θz = zeros(nreg, length(ids))
            @inbounds for j in 1:n
                r = isl.region_of_bus[j]; r == 0 && continue
                for (c, col) in enumerate(cols)
                    Zout[r, col] += P[j, c]
                    θz[r, c] += wz[j] * Θ[j, c]
                end
            end
            if stage_a
                for (e, (ia, ib)) in enumerate(pairs)
                    @inbounds for (c, col) in enumerate(cols)
                        d = θz[ia, c] - θz[ib, c]
                        num[e] += fout[e, col] * d; den[e] += d * d
                    end
                end
            end
            return nothing
        end
        ch = max(opts.snapshot_chunk, 1)
        for lo in 1:ch:S
            run_block!(lo:min(lo + ch - 1, S), lo:min(lo + ch - 1, S), ftr, Ztr, nothing, true)
        end
        for lo in 1:ch:H
            hi = min(lo + ch - 1, H)
            run_block!((S + lo):(S + hi), lo:hi, fho, Zho, kho, false)
        end

        # stage A -> refinement; keep the candidate with the lowest training flow error
        lo_b = b0 ./ opts.bnd; hi_b = b0 .* opts.bnd
        b = copy(b0); nonpos = 0
        if ncorr > 0 && opts.init == :stage_a
            dmax = maximum(den)
            for e in 1:ncorr
                if num[e] > 0 && den[e] > 1e-12 * dmax && isfinite(num[e] / den[e])
                    b[e] = clamp(num[e] / den[e], lo_b[e], hi_b[e])
                else
                    nonpos += 1
                end
            end
        end
        best_b = copy(b0); best_err = Inf
        if ncorr > 0 && S > 0
            best_err = _edge_sq_err(pairs, b0, first(_grounded_solve(allpairs, vcat(b0, kb), nreg, Ztr)), ftr)
            for it in 0:opts.refine_iters
                Θr, _ = _grounded_solve(allpairs, vcat(b, kb), nreg, Ztr)
                e_tr = _edge_sq_err(pairs, b, Θr, ftr)
                e_tr < best_err && (best_err = e_tr; best_b = copy(b))
                it == opts.refine_iters && break
                maxd = 0.0; nonpos = 0; bnext = similar(b)
                for (e, (i, j)) in enumerate(pairs)
                    nu = 0.0; de = 0.0
                    @inbounds for c in 1:S
                        d = Θr[i, c] - Θr[j, c]
                        nu += ftr[e, c] * d; de += d * d
                    end
                    bn = b[e]
                    if nu > 0 && de > 0 && isfinite(nu / de)
                        bn = nu / de
                    else
                        nonpos += 1
                    end
                    v = clamp(sqrt(b[e] * bn), lo_b[e], hi_b[e])
                    maxd = max(maxd, abs(log(v) - log(b[e])))
                    bnext[e] = v
                end
                b = bnext
                if maxd < opts.refine_tol          # converged: score the last iterate, then stop
                    e_tr = _edge_sq_err(pairs, b, first(_grounded_solve(allpairs, vcat(b, kb), nreg, Ztr)), ftr)
                    e_tr < best_err && (best_err = e_tr; best_b = copy(b))
                    break
                end
            end
        end
        b = best_b

        # held-out check: an island keeps parallel-combine unless the estimate is better on unseen snapshots
        if H > 0 && (ncorr > 0 || nk > 0)
            Θh, _ = _grounded_solve(allpairs, vcat(b, kb), nreg, Zho)
            if ncorr > 0
                Θp, _ = _grounded_solve(allpairs, vcat(b0, kb), nreg, Zho)
                ep = _edge_sq_err(pairs, b0, Θp, fho); ef = _edge_sq_err(pairs, b, Θh, fho)
                if ef > ep
                    verbose && @warn "network reduction: island held-out error not improved -> parallel-combine"
                    b = copy(b0); Θh = Θp; ef = ep; n_isl_fb += 1
                end
                se_p += ep; se_f += ef; sr += sum(abs2, fho)
            end
            se_k += _edge_sq_err(kpairs, kb, Θh, kho); sr_k += sum(abs2, kho)
        end
        # leakage: zone injection not carried by the reduced edges' fine flows
        if H > 0
            res = copy(Zho)
            for (e, (i, j)) in enumerate(pairs), c in 1:H
                res[i, c] -= fho[e, c]; res[j, c] += fho[e, c]
            end
            for (e, (i, j)) in enumerate(kpairs), c in 1:H
                res[i, c] -= kho[e, c]; res[j, c] += kho[e, c]
            end
            lk_num += sum(abs2, res); lk_den += sum(abs2, Zho)
        end

        n_fb += nonpos
        for (e, cr) in enumerate(corr)
            (b[e] <= lo_b[e] * (1 + 1e-12) || b[e] >= hi_b[e] * (1 - 1e-12)) && (n_clamp += 1)
            x_by_pair[(cr.ca, cr.cb)] = xr(b[e])
        end
        verbose && @info "network reduction island: buses=$n regions=$nreg corridors=$ncorr kept=$nk non-positive=$nonpos"
    end

    ck = _fnv1a64(join(["$(k[1])|$(k[2])|$(v)" for (k, v) in sort!(collect(x_by_pair), by = first)], ";"))
    rel(a, b) = b > 0 ? sqrt(a / b) : NaN
    return NetworkReductionResult(x_by_pair, n_isl, n_corr, n_kept, n_fb, n_clamp, n_isl_fb,
                                  rel(se_p, sr), rel(se_f, sr), rel(se_k, sr_k), rel(lk_num, lk_den),
                                  inj_total > 0 ? inj_used / inj_total : NaN, max_comp, ck, time() - t_start)
end
