# Copyright © 2026, UChicago Argonne, LLC. All Rights Reserved.
# SPDX-License-Identifier: BSD-3-Clause  (see LICENSE)

# Activate the project environment and load A-LEAF as a package
using Pkg
Pkg.activate(".")
# Optional solver: load CPLEX when it is installed in this environment (needed only for solver_name = CPLEX; the
# default solver HiGHS needs nothing).
Base.find_package("CPLEX") === nothing || @eval using CPLEX
using ALEAF

#=
# run_ALEAF

To run the ALEAF process, call the "run_ALEAF" function with the name of the master setting file
(an .xlsx workbook in the setting/ folder, without the extension).

-	master_setting_file_name::String: (Required) The name of the setting file to be used in the ALEAF process.

Edit the `setting` line below, or override it without editing this file:

	julia --project=. execute_ALEAF.jl ALEAF_Simulation_Setting_NorthAmerica
	ALEAF_SETTING=ALEAF_Simulation_Setting_NorthAmerica julia --project=. execute_ALEAF.jl
=#

# <-- EDIT: default setting file
setting = "ALEAF_Simulation_Setting_NorthAmerica"

# Command-line argument or ALEAF_SETTING environment variable takes precedence.
setting = !isempty(ARGS) ? ARGS[1] : get(ENV, "ALEAF_SETTING", setting)
println("A-LEAF setting file: $setting")

ALEAF.run_ALEAF(; master_setting_file_name = setting)
