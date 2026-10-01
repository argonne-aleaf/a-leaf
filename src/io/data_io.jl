# Copyright © 2026, UChicago Argonne, LLC. All Rights Reserved.
# SPDX-License-Identifier: BSD-3-Clause  (see LICENSE)

# Shared data I/O helpers. Keep parsing backward compatible with input templates.

using XLSX
using DataFrames


function check_and_create_path(Report_path::String)
    if !ispath(Report_path)
        mkpath(Report_path)
    end
end


function read_xlsx_return_dict_string_any_for_sim_config(file_location::String, tap_name::String; first_row_value=1)

    data = DataFrame(XLSX.readtable(file_location, tap_name, first_row=first_row_value))

    simulation_columns = names(data)[3:end]

    sim_dict = Dict{String, Any}()
    for sim_col in simulation_columns
        sim_dict[sim_col] = Dict{String, Any}()
        for i in 1:nrow(data)
            sim_dict[sim_col][data.Setting[i]] = data[i, sim_col]
        end
    end

    # remove ID if Run_Flag == false 
    for sim_col in simulation_columns 
        if sim_dict[sim_col]["Run_Flag"] == false 
            delete!(sim_dict, sim_col)
        end
    end

    return sim_dict

end


# Shared body so the string-path and open-handle methods build byte-identical dicts.
function _build_dict_string_any_from_dataframe(data::DataFrame)
    li = Dict{String, Any}()
    for row in eachrow(data)
        di = Dict{String, Any}()
        for name in names(row)
            di[string(name)] = row[name]
        end
        li[string(DataFrames.row(row))] = di
    end
    return li
end


function read_xlsx_return_dict_string_any(file_location::String, tap_name::String; first_row_value=1)

    data = DataFrame(XLSX.readtable(file_location, tap_name, first_row=first_row_value))
    return _build_dict_string_any_from_dataframe(data)

end


# Read from an already-open workbook handle to avoid re-unzipping/parsing per sheet.
function read_xlsx_return_dict_string_any(xf::XLSX.XLSXFile, tap_name::String; first_row_value=1)

    data = DataFrame(XLSX.gettable(xf[tap_name]; first_row=first_row_value))
    return _build_dict_string_any_from_dataframe(data)

end


# Shared body: string-path and open-handle methods build byte-identical dicts.
function _build_setting_dict_string_any_from_dataframe(data::DataFrame)
    dic = Dict{String, Any}()
    for row in eachrow(data)
        dic[values(row)[1]] = values(row)[2]
    end
    return dic
end


function read_setting_xlsx_return_dict_string_any(file_location::String, tap_name::String; first_row_value=1)

    data = DataFrame(XLSX.readtable(file_location, tap_name, first_row=first_row_value))
    return _build_setting_dict_string_any_from_dataframe(data)

end


# Read from an already-open workbook handle to avoid re-unzipping/parsing per sheet.
function read_setting_xlsx_return_dict_string_any(xf::XLSX.XLSXFile, tap_name::String; first_row_value=1)

    data = DataFrame(XLSX.gettable(xf[tap_name]; first_row=first_row_value))
    return _build_setting_dict_string_any_from_dataframe(data)

end


function resolve_network_input_files!(ALEAF_setting::Dict{String,<:Any}, case_id)

    data_location = joinpath(pwd(), "data", ALEAF_setting["Simulation Setting"]["test_system_name"])
    ALEAF_setting["data_location"] = data_location
    ALEAF_setting["data_location_timeseries"] = joinpath(data_location, "timeseries_data_files")

    simulation_config = ALEAF_setting["Simulation Configuration"][string(case_id)]
    network_data_id = simulation_config["Network_Data_File_ID"]
    test_system_name = get(
        ALEAF_setting["Simulation Setting"],
        "original_test_system_name",
        ALEAF_setting["Simulation Setting"]["test_system_name"],
    )
    # The split-input workflow requires a base network workbook and may include a separate config workbook.
    network_file_location = joinpath(data_location, string("network_", test_system_name, "_", network_data_id, ".xlsx"))
    if !isfile(network_file_location)
        error("Failed to find network data workbook $(basename(network_file_location)) in $(data_location).")
    end

    network_config_file_location = ""
    network_config_id = get(simulation_config, "Network_Configuration_File_ID", "")
    if !ismissing(network_config_id)
        network_config_id = strip(string(network_config_id))
        if !isempty(network_config_id) && uppercase(network_config_id) != "NA"
            network_config_file_location = joinpath(data_location, string("config_", test_system_name, "_", network_config_id, ".xlsx"))
            if !isfile(network_config_file_location)
                error("Failed to find network config workbook $(basename(network_config_file_location)) in $(data_location).")
            end
        end
    end

    ALEAF_setting["network_data_file_location"] = network_file_location
    ALEAF_setting["network_config_file_location"] = isempty(network_config_file_location) ? nothing : network_config_file_location
    ALEAF_setting["network_config_available"] = !isempty(network_config_file_location)

    # Keep the legacy setting key aligned with the new base network workbook.
    ALEAF_setting["test_system_file_location"] = network_file_location

    return nothing

end


function load_file_path_sheet!(ALEAF_setting::Dict{String,<:Any}, file_location::String)

    ALEAF_setting["File Path"] = Dict{String, Any}()
    # File-path references now live with the base network workbook regardless of config availability.
    data = DataFrame(XLSX.readtable(file_location, "File Path"))
    for row in eachrow(data)
        ALEAF_setting["File Path"][values(row)[1]] = values(row)[2]
    end

    return nothing

end


function read_csv_return_dict_string_any_per_round(file_location::String, round_idx)

    data = CSV.read(file_location, DataFrame)

    data = data[data.round .== round_idx, Not(:round)]

    li = Dict{String, Any}()
    for row in eachrow(data)
        di = Dict{String, Any}()
        for name in names(row)
            di[string(name)] = row[name]
        end
        li[string(DataFrames.row(row))] = di
    end

    return li

end


function read_csv_return_dict_string_any(file_location::String)

    data = CSV.read(file_location, DataFrame)

    li = Dict{String, Any}()
    for row in eachrow(data)
        di = Dict{String, Any}()
        for name in names(row)
            di[string(name)] = row[name]
        end
        li[string(DataFrames.row(row))] = di
    end

    return li

end


function convert_dataFrame_to_dict_string_any(input_dataframe)

    li = Dict{String, Any}()
    for row in eachrow(input_dataframe)
        di = Dict{String, Any}()
        for name in names(row)
            di[string(name)] = row[name]
        end
        li[string(DataFrames.row(row))] = di
    end

    return li

end



function define_output_path(ALEAF_setting::Dict{String,<:Any}, case_id)

    output_path = joinpath(pwd(), "output")
    ALEAF.check_and_create_path(output_path)

    if haskey(ALEAF_setting["Simulation Setting"], "original_test_system_name")
        output_path = joinpath(output_path, ALEAF_setting["Simulation Setting"]["original_test_system_name"])
    else
        output_path = joinpath(output_path, ALEAF_setting["Simulation Setting"]["test_system_name"])
    end
    ALEAF.check_and_create_path(output_path)

    SNAME = ALEAF_setting["Simulation Configuration"][string(case_id)]["Case_ID"]
    output_path = joinpath(output_path, string("case_id_", case_id, "_", SNAME))
    ALEAF.check_and_create_path(output_path)

    # Keep trailing separator for legacy string(output_path, file_name) call sites.
    return string(output_path, "/")

end



