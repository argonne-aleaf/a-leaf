# Copyright © 2026, UChicago Argonne, LLC. All Rights Reserved.
# SPDX-License-Identifier: BSD-3-Clause  (see LICENSE)

# Aggregated network DB export: the live-aggregated network (e.g. County built from a private nodal DB) becomes
# the finest level of a standalone DB, so the nodal source never has to be shared.

# aggregate_branches_using_database! converts branch Length km -> miles with this factor; the export inverts it.
const _EXPORT_MILES_PER_KM = 0.621371
# network/config tables are read with first_row_value = 2 (row 1 holds optional category labels)
const _EXPORT_HEADER_ROW = 2
const _EXPORT_BRANCH_COLUMNS = ("UID", "f_bus", "t_bus", "rate_a", "rate_c", "Length", "model_flag", "expansion_flag",
                                "max_rate_a", "br_x_pu", "br_r_pu", "transformer", "dc_line", "Notes")
const _EXPORT_BUS_COLUMNS = ("bus_i", "bus_name", "MW load", "TimeSeriesStatus", "Notes", "full  bus name", "Longitude", "Latitude")

const _ExportCells = Dict{Tuple{Int,Int},Any}


# Non-empty cells keyed by (row, col); rows are streamed so a stale <dimension> tag cannot clip data.
function _export_read_cells(ws::XLSX.Worksheet)
    cells = _ExportCells()
    for sr in XLSX.eachrow(ws)
        isempty(sr) && continue
        r = XLSX.row_number(sr)
        c1, c2 = XLSX.column_bounds(sr)
        for c in c1:c2
            v = sr[c]
            ismissing(v) || (cells[(r, c)] = v)
        end
    end
    return cells
end

function _export_read_workbook(path::String)
    return XLSX.openxlsx(path) do xf
        [name => _export_read_cells(xf[name]) for name in XLSX.sheetnames(xf)]
    end
end

_export_header(cells::AbstractDict, row::Int) =
    Dict{String,Int}(string(v) => c for ((r, c), v) in cells if r == row && v isa AbstractString)

_export_data_rows(cells::AbstractDict, header_row::Int) = sort!(unique(r for (r, _) in keys(cells) if r > header_row))

_export_isnum(v) = v isa Real && !(v isa Bool) && isfinite(v)

function _export_cell_value(v)
    v isa Bool && return v
    v isa Integer && return Int(v)
    if v isa Real
        isfinite(v) || error("Aggregated network export: cannot write non-finite value $(v) to xlsx.")
        return Float64(v)
    end
    v isa AbstractString && return String(v)
    v isa XLSX.CellValueType && return v
    error("Aggregated network export: unsupported cell value $(repr(v)) of type $(typeof(v)).")
end

function _export_write_workbook(path::String, sheets::Vector{Pair{String,_ExportCells}})
    xf = XLSX.openxlsx(path, mode = "w")
    for (i, (name, cells)) in enumerate(sheets)
        ws = i == 1 ? xf[1] : XLSX.addsheet!(xf, name)
        i == 1 && XLSX.rename!(ws, name)
        for (r, c) in sort!(collect(keys(cells)))
            v = cells[(r, c)]
            (ismissing(v) || v === nothing || (v isa AbstractString && isempty(v))) && continue
            ws[r, c] = _export_cell_value(v)
        end
    end
    # saved only once every sheet is built: the do-block form of openxlsx writes the file even when it throws
    XLSX.writexlsx(path, xf, overwrite = true)
    return path
end


# network_data carries no coordinates, so each bus gets the load-weighted centroid of its member fine buses.
function _export_bus_cells(src::_ExportCells, network_data::Dict{String,<:Any})
    hr = _EXPORT_HEADER_ROW
    h = _export_header(src, hr)
    for n in ("bus_i", "MW load", "Longitude", "Latitude")
        haskey(h, n) || error("Aggregated network export: the source bus sheet has no '$n' column.")
    end
    extra = setdiff(keys(h), _EXPORT_BUS_COLUMNS)
    isempty(extra) || @aleaf_warn "Aggregated network export: bus columns $(collect(extra)) have no aggregated value and are left empty."

    fine = Dict{String,Tuple{Float64,Any,Any}}()
    for r in _export_data_rows(src, hr)
        id = get(src, (r, h["bus_i"]), missing)
        ismissing(id) && continue
        load = max(_mw(get(src, (r, h["MW load"]), 0.0)), 0.0)
        fine[string(id)] = (load, get(src, (r, h["Longitude"]), missing), get(src, (r, h["Latitude"]), missing))
    end

    out = _ExportCells(k => v for (k, v) in src if k[1] <= hr)
    n_without_coords = 0
    for (i, k) in enumerate(sort!(collect(keys(network_data["bus"])), by = x -> parse(Int, x)))
        b = network_data["bus"][k]
        pts = [fine[m] for m in string.(b["aggregation_info"]["aggregated_regions_bus_i"])
               if haskey(fine, m) && _export_isnum(fine[m][2]) && _export_isnum(fine[m][3])]
        lon = lat = missing
        if isempty(pts)
            n_without_coords += 1
        else
            w = [p[1] for p in pts]
            sum(w) > 0 || (w = ones(length(pts)))   # no member load to weight by
            lon = round(sum(w .* [Float64(p[2]) for p in pts]) / sum(w), digits = 6)
            lat = round(sum(w .* [Float64(p[3]) for p in pts]) / sum(w), digits = 6)
        end
        # region ids keep the config's cell type: the model matches bus_i to sub_area_mapping by raw value
        id = b["bus_i"]
        vals = Dict{String,Any}("bus_i" => id, "bus_name" => id, "MW load" => b["MW load"], "TimeSeriesStatus" => true,
                                "full  bus name" => id, "Longitude" => lon, "Latitude" => lat)
        for (n, c) in h
            v = get(vals, n, missing)
            ismissing(v) || (out[(hr + i, c)] = v)
        end
    end
    return out, n_without_coords
end


function _export_branch_cells(src::_ExportCells, network_data::Dict{String,<:Any})
    hr = _EXPORT_HEADER_ROW
    h = _export_header(src, hr)
    absent = [n for n in _EXPORT_BRANCH_COLUMNS if !haskey(h, n)]
    isempty(absent) || error("Aggregated network export: the source branch sheet lacks columns $(absent).")
    extra = setdiff(keys(h), _EXPORT_BRANCH_COLUMNS)
    isempty(extra) || @aleaf_warn "Aggregated network export: branch columns $(collect(extra)) have no aggregated value and are left empty."

    bus = network_data["bus"]
    out = _ExportCells(k => v for (k, v) in src if k[1] <= hr)
    for (i, k) in enumerate(sort!(collect(keys(network_data["branch"])), by = x -> parse(Int, x)))
        br = network_data["branch"][k]
        f = bus[string(br["f_bus"])]["bus_i"]
        t = bus[string(br["t_bus"])]["bus_i"]
        (isequal(f, get(br, "f_bus_id", f)) && isequal(t, get(br, "t_bus_id", t))) ||
            error("Aggregated network export: branch $k endpoint ids disagree with its bus keys.")
        dc = get(br, "dc_line", false) == true
        n = length(get(br, "merged_line_list", []))
        note = get(br, "nodal_kept", false) == true ? "kept" : string(dc ? "dc" : "agg", " n=", n)
        vals = Dict{String,Any}(
            "UID" => i, "f_bus" => f, "t_bus" => t, "rate_a" => br["rate_a"], "rate_c" => br["rate_a"],
            "Length" => br["length"] / _EXPORT_MILES_PER_KM, "model_flag" => br["model_flag"] == true,
            "expansion_flag" => br["expansion_flag"] == true, "max_rate_a" => br["max_rate_a"], "br_x_pu" => br["br_x_pu"],
            "br_r_pu" => 0, "transformer" => false, "dc_line" => dc, "Notes" => note)
        for (name, c) in h
            v = get(vals, name, missing)
            ismissing(v) || (out[(hr + i, c)] = v)
        end
    end
    return out
end


# Rows follow the source sheet because the model re-keys entries in Dict order; bus_ID/bus_name become the region id.
function _export_entity_cells(src::_ExportCells, entries::AbstractDict, sheet::String)
    hr = _EXPORT_HEADER_ROW
    hdr = sort!([(c, string(v)) for ((r, c), v) in src if r == hr && v isa AbstractString])
    any(n -> n == "bus_ID", last.(hdr)) || error("Aggregated network export: the source '$sheet' sheet has no bus_ID column.")
    keys_e = sort!(collect(keys(entries)), by = x -> parse(Int, x))
    # match on the columns the reader actually loaded (it stops at a blank header cell)
    sig_cols = [(c, n) for (c, n) in hdr if all(k -> haskey(entries[k], n), keys_e)]
    by_sig = Dict{Any,Vector{Int}}()
    src_rows = _export_data_rows(src, hr)
    for r in src_rows
        push!(get!(by_sig, Tuple(get(src, (r, c), missing) for (c, _) in sig_cols), Int[]), r)
    end
    order = Tuple{Int,Int,String}[]
    n_unmatched = 0
    for k in keys_e
        rows = get(by_sig, Tuple(get(entries[k], n, missing) for (_, n) in sig_cols), Int[])
        isempty(rows) && (n_unmatched += 1)
        push!(order, (isempty(rows) ? typemax(Int) : popfirst!(rows), parse(Int, k), k))
    end
    sort!(order)
    n_unmatched > 0 && @aleaf_warn "Aggregated network export: $n_unmatched '$sheet' entries match no source row; written after the matched rows."

    bool_cols = Set(c for ((r, c), v) in src if r > hr && v isa Bool)
    out = _ExportCells(k => v for (k, v) in src if k[1] <= hr)
    for (i, (_, _, k)) in enumerate(order)
        e = entries[k]
        for (c, n) in hdr
            v = n in ("bus_ID", "bus_name") ? e["bus_i"] : get(e, n, missing)
            if c in bool_cols && v isa AbstractString && uppercase(strip(v)) in ("TRUE", "FALSE")
                v = uppercase(strip(v)) == "TRUE"
            end
            ismissing(v) || (out[(hr + i, c)] = v)
        end
    end
    return out, (exported = length(order), source_rows = length(src_rows), unmatched = n_unmatched)
end


# Drops column `name` (header on `header_row`) and compacts all-empty rows, which would end the reader's table early.
function _export_drop_column(cells::_ExportCells, header_row::Int, name::String; dedup::Bool)
    h = _export_header(cells, header_row)
    dc = get(h, name, nothing)
    dc === nothing && error("Aggregated network export: no '$name' column to drop on config row $header_row.")
    ncol = maximum(c for (_, c) in keys(cells))
    out = _ExportCells(((r, c > dc ? c - 1 : c) => v) for ((r, c), v) in cells if r <= header_row && c != dc)
    seen = Set{Any}()
    nr = header_row
    data_rows = _export_data_rows(cells, header_row)
    for r in data_rows
        t = Tuple(get(cells, (r, c), missing) for c in 1:ncol if c != dc)
        all(ismissing, t) && continue
        if dedup
            t in seen && continue
            push!(seen, t)
        end
        nr += 1
        for (j, v) in enumerate(t)
            ismissing(v) || (out[(nr, j)] = v)
        end
    end
    return out, (rows_in = length(data_rows), rows_out = nr - header_row)
end


# Level 1 (nodal) is removed so the exported buses are the finest level; every level number shifts down by one.
function _export_config_sheets(config_file::String)
    hr = _EXPORT_HEADER_ROW
    src = _export_read_workbook(config_file)
    cells = Dict(src)
    for s in ("network_resolution_level", "Network Setting", "sub_area_mapping", "sub_area_list")
        haskey(cells, s) || error("Aggregated network export: config workbook has no '$s' sheet.")
    end

    nrl = cells["network_resolution_level"]
    h = _export_header(nrl, hr)
    (haskey(h, "Network_Resolution_Level") && haskey(h, "Network_Resolution_ID")) ||
        error("Aggregated network export: network_resolution_level needs Network_Resolution_Level/Network_Resolution_ID headers.")
    lc, ic = h["Network_Resolution_Level"], h["Network_Resolution_ID"]
    rows = _export_data_rows(nrl, hr)
    # the model keys levels by row order, so rows must be levels 1..n in order
    (length(rows) >= 2 && all(j -> isequal(get(nrl, (rows[j], lc), missing), j), eachindex(rows))) ||
        error("Aggregated network export: network_resolution_level rows must be levels 1..n in order (found at least one gap).")
    removed = string(nrl[(rows[1], ic)])
    finest = string(nrl[(rows[2], ic)])
    new_nrl = _ExportCells(k => v for (k, v) in nrl if k[1] <= hr)
    for (j, r) in enumerate(rows[2:end]), ((rr, c), v) in nrl
        rr == r && (new_nrl[(hr + j, c)] = c == lc ? v - 1 : v)
    end

    ns = copy(cells["Network Setting"])
    level_changes = Pair{String,Any}[]
    for ((r, c), key) in cells["Network Setting"]
        (c == 1 && key isa AbstractString) || continue
        v = get(ns, (r, 2), missing)
        if endswith(key, "_level")
            if _export_isnum(v) && isinteger(v)
                v > 1 || error("Aggregated network export: Network Setting '$key' = $v is the '$removed' level, which the export removes.")
                ns[(r, 2)] = v - one(v)
                push!(level_changes, key => (v => v - one(v)))
            else
                @aleaf_warn "Aggregated network export: Network Setting '$key' = $(repr(v)) is not a level number; copied unchanged."
            end
        elseif endswith(key, "_type") && isequal(v, removed)
            error("Aggregated network export: Network Setting '$key' uses the '$removed' level, which the export removes.")
        end
    end

    level_of(n) = (m = match(r"^Network Data Level (\d+)$", n); m === nothing ? nothing : parse(Int, m[1]))
    present = Set(k for k in level_of.(first.(src)) if k !== nothing)
    absent = [k for k in 2:length(rows) if !(k in present)]
    isempty(absent) || error("Aggregated network export: config lacks 'Network Data Level' sheet(s) $(absent).")

    out = Pair{String,_ExportCells}[]
    mapping_info = list_info = nothing
    for (name, sheet_cells) in src
        k = level_of(name)
        if k !== nothing
            k == 1 || push!(out, string("Network Data Level ", k - 1) => sheet_cells)
        elseif name == "network_resolution_level"
            push!(out, name => new_nrl)
        elseif name == "Network Setting"
            push!(out, name => ns)
        elseif name == "sub_area_mapping"
            c, mapping_info = _export_drop_column(sheet_cells, hr, removed; dedup = true)
            push!(out, name => c)
        elseif name == "sub_area_list"
            c, list_info = _export_drop_column(sheet_cells, hr, removed; dedup = false)
            push!(out, name => c)
        else
            push!(out, name => sheet_cells)
        end
    end
    return out, (removed_level = removed, finest_level = finest, n_levels = length(rows) - 1, level_changes = level_changes,
                 sub_area_mapping = mapping_info, sub_area_list = list_info)
end


"""
    export_aggregated_network(network_data, ALEAF_setting, out_dir, out_name; network_id="Base", config_id="Base")

Writes `network_data` (after `get_network_data!`) as `out_dir/network_<out_name>_<network_id>.xlsx` and the source
config minus its finest level as `out_dir/config_<out_name>_<config_id>.xlsx`; returns a summary NamedTuple.
"""
function export_aggregated_network(network_data::Dict{String,<:Any}, ALEAF_setting::Dict{String,<:Any},
                                   out_dir::String, out_name::String; network_id = "Base", config_id = "Base")

    src_net = ALEAF_setting["network_data_file_location"]
    src_cfg = get(ALEAF_setting, "network_config_file_location", nothing)
    src_cfg === nothing && error("Aggregated network export needs a network config workbook (Network_Configuration_File_ID).")
    net_out = joinpath(out_dir, string("network_", out_name, "_", network_id, ".xlsx"))
    cfg_out = joinpath(out_dir, string("config_", out_name, "_", config_id, ".xlsx"))
    sources = abspath.((src_net, src_cfg))
    (abspath(net_out) in sources || abspath(cfg_out) in sources) &&
        error("Aggregated network export: output would overwrite its own source workbook.")

    config_sheets, cfg_info = _export_config_sheets(src_cfg)

    # exported bus ids become the finest level, so every bus must sit exactly at the old level 2
    for b in values(network_data["bus"])
        rc = b["region_config"]
        res = rc[rc["region_lookup_type"]]
        res == cfg_info.finest_level || error("Aggregated network export: bus $(b["bus_i"]) is aggregated at '$res'; " *
            "the export needs every bus at '$(cfg_info.finest_level)' (level 2 of the source config).")
    end
    bus_ids = [string(b["bus_i"]) for b in values(network_data["bus"])]
    allunique(bus_ids) || error("Aggregated network export: aggregated bus ids are not unique.")

    src_sheets = _export_read_workbook(src_net)
    sheet_names = first.(src_sheets)
    for s in ("bus", "branch", "plant", "hybrid", "demand")
        s in sheet_names || error("Aggregated network export: source network workbook has no '$s' sheet.")
    end

    net_sheets = Pair{String,_ExportCells}[]
    entity_info = Dict{String,Any}()
    n_without_coords = 0
    for (name, cells) in src_sheets
        if name == "bus"
            c, n_without_coords = _export_bus_cells(cells, network_data)
            push!(net_sheets, name => c)
        elseif name == "branch"
            push!(net_sheets, name => _export_branch_cells(cells, network_data))
        elseif name in ("plant", "hybrid", "demand")
            c, entity_info[name] = _export_entity_cells(cells, network_data[name], name)
            push!(net_sheets, name => c)
        else
            push!(net_sheets, name => cells)
        end
    end

    mkpath(out_dir)
    _export_write_workbook(net_out, net_sheets)
    _export_write_workbook(cfg_out, config_sheets)

    return (network_file = net_out, config_file = cfg_out, source_network_file = src_net, source_config_file = src_cfg,
            buses = length(network_data["bus"]), branches = length(network_data["branch"]),
            plant = entity_info["plant"], hybrid = entity_info["hybrid"], demand = entity_info["demand"],
            buses_without_coords = n_without_coords, config = cfg_info)
end
