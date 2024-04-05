# EPI50 analysis

An analytical framework to estimate the global public health impact of 50 years of the Expanded Programme on Immunization (EPI). This framework estimates deaths averted, years of live saved, and years of full health gained for 14 pathogens within the EPI portfolio. Exclusions include SARS-CoV-2/COVID-19 and HPV. The analysis is performed at country level for 194 UN member states for the timeframe mid-1974 t0 mid-2024, with results summarised at the global and WHO-region level. 

## The repository

This open-source repository contains all code needed to fully reproduce the analysis described in the peer-reviewed publication:

[Contribution of vaccination to improved child survival: quantifying 50 years of the Expanded Programme on Immunization](https://www.sciencedirect.com/science/article/pii/S1755436521000785)

This repository is primarily written in R, and is stable for R versions 4.3.0 and 4.3.2. It may be possible to run alternative versions of R, but these have not been tested and verified. Configuration files are primarily written in YAML markup. The use of YAML files requires no additional software beyond the R packages on which this repositiory depends.

All package dependencies will be automatically installed (if necessary) and loaded when the pipeline is run. See `dependencies.R` in the main directory for a list of all package dependencies. All of the data used to generate this analysis is publicly available. In some cases, this data is stored in the `/input/` directory, whilst in other cases data is pulled from internet. As such, a live internet connection is required to run this analysis.

No special cluster computing resources are required to run this analysis. Local parallelisation is used throughout the pipeline to improve runtime, however this parallelisation is only available on UNIX operating systems. Running this analysis on Windows operating systems is still possible, yet will be slower. Expect the full analysis (exluding any package installations) to take between 45 to 60 minutes on a UNIX machine and between 60 and 120 minutes on a Windows machine, dependant on machine specifications.

We invite any potential collaborators interested in expanding upon this analysis to fork this repository. Please submit a pull request if you feel the parent repository would benefit from any changes to code or configuration files.

## Directory structure

All code required to run the pipeline sits within the main repository directory (called `epi50` by default). Initially, three sub-directories can also be found within this main directory: 
 1. a `/config/` directory containing a series of configuration files written in YAML.
 2. an `/input/` directory containing the input data required to run the pipeline.
 3. an `/extern/` directory containing processed results from the measles and polio models used in this analysis.
A fourth directory, `/output/`, is automatically created when the pipeline is first launched. All final intermediary and final results are stored in this output folder. The figures presented in the peer-reviewed publication are stored within the `/output/figures/` sub-directory. 

## Running the pipeline

All pipeline modules are launched from `launch.R`. Preferred usage is to 'source' this file (without 'echo' is ideal). When sourced, the current working directory is automatically reset the EPI50 repository. Alternative UNIX command line usage is to `cd` to the EPI50 repository then call `sh launch.sh`.

The pipeline consists of a series of modules that, in general, should be run consecutively. The module/s to be run can be set in `launch.R` (line 20). There are 8 modules in the entire pipeline, each described below. By default, all 8 modules will be run, and in the process producing all outputs described in the associated manuscript.

### Step 1: Generate a model calibration
- Set `do_step = 1` in `launch.R` and 'source' to calibrate the model
- There are three broad options for model calibration available (by setting `calibration_type` in your analysis file):
  1. Calibrate the model to a user-defined effective reproduction number at the start of the simulation period (`calibration_type: "r_user"`, the default)
  2. Calibrate the model to an effective reproduction number that is calculated from setting-specific epidemiological data (`calibration_type: "r_data"`)
  3. Calibrate the model to a set of setting-specific epidemiological metrics (such as confirmed cases, hospitalisations, and deaths) over time (`calibration_type: "epi_data"`)

### Step 2: Run scenarios
- Set `do_step = 2` in `launch.R` and 'source' to simulate the baseline and any user-defined alternative scenarios detailed in the analysis file

## Authors

#### Development and maintenance:
* Andrew J. Shattock (shattocka@who.int)

#### Contributors:
* Helen C. Johnson
* So Yoon Sim
* Austin Carter
