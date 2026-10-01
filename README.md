<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="data/common/aleaf_logo_dark.svg">
    <img src="data/common/aleaf_logo.svg" alt="A-LEAF: Argonne Large-Scale Electricity Analysis Framework" width="480">
  </picture>
</p>

# Argonne Large-Scale Electricity Analysis Framework (A-LEAF)

A-LEAF is an integrated power system simulation framework for long-term capacity expansion, operations, and reliability assessment workflows. It provides three linked model families: generation and transmission expansion planning, system operation (unit commitment and economic dispatch), and reliability and resource adequacy assessment.

- ANL page: https://www.anl.gov/esia/a-leaf
- Documentation: https://argonne-aleaf.github.io/a-leaf-docs/

## What A-LEAF does

**Capabilities.** A-LEAF combines advanced optimization (least-cost planning and operations, sub-hourly dispatch, multiday representative periods, and joint generation, transmission, and storage expansion), a detailed U.S. grid database, climate and weather data, wholesale market design, policies and regulations, multi-sector interdependency, and reliability and resource adequacy assessment.

<p align="center">
  <a href="data/common/aleaf_capabilities.svg"><img src="data/common/aleaf_capabilities.svg" alt="Power system modeling with A-LEAF: advanced optimization, detailed U.S. grid database, climate and weather data, wholesale market design, policies and regulations, multi-sector interdependency, and reliability and resource adequacy assessment." width="100%"></a>
</p>

**Applications.** Typical studies include electricity system modernization, wholesale market analysis, long-term planning, short-term operations and reliability assessment, technoeconomic valuation, and extreme weather resilience.

<p align="center">
  <a href="data/common/aleaf_applications.svg"><img src="data/common/aleaf_applications.svg" alt="A-LEAF applications: electricity system modernization, wholesale market analysis, long-term planning, short-term operations and reliability assessment, technoeconomic valuation, and extreme weather resilience." width="100%"></a>
</p>

The figures show the full A-LEAF framework. Some items, such as the county-level Texas system, coupling with the TIMES energy systems model, and weather years derived from climate models, use data or models that are not part of this public release. The bundled example is the North America database in `data/NorthAmerica/`. See the [documentation](https://argonne-aleaf.github.io/a-leaf-docs/) for details.

## Requirements
- macOS or Linux (Windows is untested)
- Julia `1.12` or newer (tested with `1.12.5`)
- Solver: **HiGHS** is installed with the Julia dependencies and is the default. IBM ILOG CPLEX (license required) and GPU solvers are optional; see [Solvers](#solvers)
- This repository, which includes the bundled input data (about 200 MB in `data/`) and the `setting/` control workbook
- Optional: an NVIDIA GPU for the GPU solvers

## Installation (from scratch)

**1. Install Julia.** The simplest way is [juliaup](https://github.com/JuliaLang/juliaup):
```bash
curl -fsSL https://install.julialang.org | sh
juliaup add 1.12
juliaup default 1.12
```

**2. Clone the repository:**
```bash
git clone https://github.com/argonne-aleaf/a-leaf.git a-leaf
cd a-leaf
```

**3. Install the Julia dependencies:**
```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```
This installs the exact package versions recorded in `Manifest.toml` and takes several minutes the first time
(it downloads and precompiles about 130 packages).

**4. Check the installation:**
```bash
julia --project=. -e 'using ALEAF; println("LOAD OK")'
```

### Troubleshooting

| Symptom | Fix |
|---|---|
| `setting/<name>.xlsx not found` | Run from the repository root; A-LEAF reads `setting/` and `data/` relative to the current directory. |
| Very slow first run | The first `using ALEAF` precompiles the project (and again after you edit files under `src/`). Later runs are fast. |

## Quick Start

### Run a case
A run is controlled by one Excel workbook in `setting/`. From the repository root, pass its file name without the `.xlsx` extension:

```bash
julia --project=. execute_ALEAF.jl ALEAF_Simulation_Setting_NorthAmerica
```

Or edit the `setting` line at the top of `execute_ALEAF.jl` and run `julia --project=. execute_ALEAF.jl`. You can also set the `ALEAF_SETTING` environment variable. A-LEAF is loaded as a Julia package from this clone (`using ALEAF`); after you edit files under `src/`, the next launch precompiles the change automatically.

The equivalent call from a Julia REPL:

```julia
using Pkg
Pkg.activate(".")
using ALEAF

ALEAF.run_ALEAF(; master_setting_file_name = "ALEAF_Simulation_Setting_NorthAmerica")
```

`run_ALEAF` has one required argument, `master_setting_file_name::String`, the setting workbook name. Results are written under `output/`.

Bundled setting workbooks:

| Setting workbook (`setting/`) | Data (`data/`) |
|---|---|
| `ALEAF_Simulation_Setting_NorthAmerica` | `NorthAmerica/` |

The bundled workbook enables one case, `Test_EXP` (expansion planning for the Texas system at balancing-authority resolution), as a first run. Full-size cases are large and can run for hours; check the `Simulation Configuration` sheet for which cases have `Run_Flag` enabled. For what each setting means, see the [documentation](https://argonne-aleaf.github.io/a-leaf-docs/).

### Alternative: install as a Julia package
The same prerequisites apply.
```julia
using Pkg
Pkg.add(url="https://github.com/argonne-aleaf/a-leaf.git")
```
Then call `using ALEAF; ALEAF.run_ALEAF(; master_setting_file_name = "...")`. This uses the installed copy of the package, not a clone, so edits to a clone's `src/` are not used. Run from a directory that contains `setting/` and `data/`. To update, run `Pkg.update()` then `Pkg.precompile()`.

## Solvers
- Default: `HiGHS` (installed automatically; set `solver_name = HiGHS` in the `Simulation Setting` sheet)
- Optional, needs a license: `CPLEX`, loaded as a package extension (see below)
- Optional GPU backends, off by default: `cuOpt` (PDLP, LP only) and `MadNLP`/`ExaModels`, loaded as package extensions
- Optimization stack: `JuMP`, `MathOptInterface`

See `Project.toml` for the full dependency list.

### Optional: CPLEX
CPLEX is not installed by default and is not needed to run A-LEAF. If you have an IBM ILOG CPLEX Studio license (supported versions 12.10, 20.1 and 22.1.x; recommended 22.1.x):
1. Set `CPLEX_STUDIO_BINARIES` to the folder that contains the CPLEX binaries, then add and build the package. When working from a clone, run `julia --project=. -e 'using Pkg; Pkg.add("CPLEX")'` (this edits the clone's `Project.toml`/`Manifest.toml`; keep that change local and do not commit it). In package mode (`Pkg.add(url=...)`), run `Pkg.add("CPLEX")` in your own environment instead.
2. Load it. With `execute_ALEAF.jl` nothing needs editing: it loads CPLEX automatically when the package is installed. In your own script or the Julia REPL, add `using CPLEX` (before or after `using ALEAF`, but before `ALEAF.run_ALEAF`).
3. Set `solver_name = CPLEX` in the `Simulation Setting` sheet; the solver parameters are read from the `CPLEX Setting` sheet.

If `solver_name = CPLEX` is selected without `using CPLEX`, A-LEAF stops with a message saying so.

### GPU solvers (optional)
Set `solver_name` in the setting workbook (`Simulation Setting` sheet, cell `B11`) to `cuOpt` or `MadNLP`. The GPU packages are not part of the default environment: install `cuOpt` (or the MadNLP stack) into a Julia environment together with A-LEAF and load it with `using ALEAF`, so the matching extension activates. cuOpt is LP-only, so integer build variables must be relaxed. cuOpt and MadNLP cannot be loaded in the same Julia session. Setup, options and caveats: [GPU Solvers](https://argonne-aleaf.github.io/a-leaf-docs/configuration/GPU_Solvers/).

## Environment Variables
| Variable | Used by | Purpose |
|---|---|---|
| `ALEAF_SETTING` | `execute_ALEAF.jl` | Setting workbook name (no `.xlsx`); a command-line argument takes precedence, the default in the file is used last. |
| `ALEAF_CASE_ID` | `run_ALEAF` | Optional: run only this `Run_Flag` case (one case per Julia process, useful for GPU runs). |
| `CPLEX_STUDIO_BINARIES` | `Pkg.build("CPLEX")` | Path to the CPLEX binaries (optional CPLEX support only). |

## Repository Layout
- `setting/`: simulation setting workbooks
- `data/`: network, technology and time-series data
- `output/`: run results
- `src/`: source code
  - `core`, `component`: base types, variables and constraints
  - `io`: data I/O and logging
  - `network`: network generation and PTDF reduction
  - `model`: `LCO_GTEP` (expansion and operation), `RA` (resource adequacy)
  - `runtime`: run orchestration and model dispatch
  - `solvers`: solve and result collection
  - `util`: scenario reduction and outage scenarios
- `ext/`: optional solver extensions (CPLEX, cuOpt, MadNLP)
- `execute_ALEAF.jl`: entry script

## Citing A-LEAF
If you use A-LEAF in your work, please cite it. GitHub shows the citation under "Cite this repository" in the sidebar, generated from [CITATION.cff](CITATION.cff). A BibTeX entry:

```bibtex
@software{aleaf,
  title   = {Argonne Large-Scale Electricity Analysis Framework (A-LEAF)},
  author  = {Kwon, Jonghwan and Levin, Todd and Mann, Neal},
  year    = {2026},
  version = {1.0.0},
  url     = {https://github.com/argonne-aleaf/a-leaf},
  license = {BSD-3-Clause}
}
```

## Contact
Questions, bug reports, and feature requests are welcome as [GitHub issues](https://github.com/argonne-aleaf/a-leaf/issues). For other inquiries, contact Jonghwan Kwon ([kwonj@anl.gov](mailto:kwonj@anl.gov)) or Todd Levin ([tlevin@anl.gov](mailto:tlevin@anl.gov)).

## License
A-LEAF, including the input data distributed under `data/`, is released under the BSD 3-Clause license of UChicago Argonne, LLC. See [LICENSE](LICENSE).

## Data Sources
Bundled inputs draw on public sources, including:
- the [Annual Technology Baseline (ATB) 2024](https://atb.nlr.gov/electricity/2024/data) cost and performance projections from the National Laboratory of the Rockies (NLR, formerly NREL). `data/common/ATB_2024.csv` is a subset of the ATB data download: only the columns and parameters A-LEAF reads (CAPEX, fixed and variable O&M, FCR, WACC) are kept, for all technologies, cases and scenarios. To use other ATB parameters or a different ATB release, download the full data from the link above and replace the file.
- data from the U.S. Energy Information Administration (EIA), such as Annual Energy Outlook (AEO) fuel-price projections (`data/NorthAmerica/timeseries_data_files/Fuel/`).
