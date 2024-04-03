# EPI50 analysis

An analytical framework to estimate the global public health impact of 50 years of the Expanded Programme on Immunization (EPI). This framework estimates deaths averted, years of live saved, and years of full health gained for 14 pathogens within the EPI portfolio. Exclusions include SARS-CoV-2/COVID-19 and HPV. The analysis is performed at country level for 194 UN member states for the timeframe mid-1974 t0 mid-2024, with results summarised at the global and WHO-region level. 

## The repository

This open-source repository contains all code needed to fully reproduce the analysis described in the peer-reviewed publication:

[Contribution of vaccination to improved child survival: quantifying 50 years of the Expanded Programme on Immunization](https://www.sciencedirect.com/science/article/pii/S1755436521000785)

This repository is primarily written in R, and is stable for R version 4.3.0. It may be possible to run alternative versions of R, but these have not been tested and verified. Configuration files are primarily written in YAML markup. The use of YAML files requires no additional software beyond the R packages on which this repositiory depends.

All package dependencies will be automatically installed (if necessary) and loaded when the pipeline is run. See `dependencies.R` in the main directory for a list of all package dependencies. All of the data used to generate this analysis is publicly available. In some cases, this data is stored in the `/input/` directory, whilst in other cases data is pulled from internet. As such, a live internet connection is required to run this analysis.

No special cluster computing resources are required to run this analysis. Local parallelisation is used throughout the pipeline to improve runtime, however this parallelisation is only available on UNIX operating systems. Running this analysis on Windows operating systems is still possible, yet will be slower. Expect the full analysis (exluding any package installations) to take between 45 to 60 minutes on a UNIX machine and between 60 and 120 minutes on a Windows machine, dependant on machine specifications.

We invite any potential collaborators interested in expanding upon this analysis to fork this repository. Please submit a pull request if you feel the parent repository would benefit from any changes to code or configuration files.

## The pipeline

### Getting started
- All pipeline modules are launched from `launch.R`
  - Preferred usage is to 'source' this file (without 'echo' is ideal)
- Alternative UNIX command line usage is to `cd` to the code directory then call `sh launch.sh`
  
### Directory structure
- This is all taken care of when `launch.R` (or `bash_launch.sh`) is called
- The model will automatically re-point the current working directory
- The model will also automatically create the necessary output directory structure

#### Development and maintenance:
* Andrew J. Shattock (shattocka@who.int)

#### Contributors:
* Helen C. Johnson
* So Yoon Sim
* Austin Carter
