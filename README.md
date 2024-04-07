# EPI50 analysis

An analytical framework to estimate the global public health impact of 50 years of the Expanded Programme on Immunization (EPI). This framework estimates deaths averted, years of live saved, and years of full health gained for 14 pathogens within the EPI portfolio. Exclusions include SARS-CoV-2/COVID-19 and HPV. The analysis is performed at country level for 194 UN member states for the timeframe mid-1974 to mid-2024, with results summarised at the global and WHO-region level. This analysis builds upon infectious disease modelling estimates produced by the Vaccine Impact Modelling Consortium (VIMC) and the Global Burden of Disease study (GBD).

## The repository

This open-source repository contains all code needed to fully reproduce the analysis presented in the peer-reviewed publication:

[Contribution of vaccination to improved child survival: quantifying 50 years of the Expanded Programme on Immunization](https://www.sciencedirect.com/science/article/pii/S1755436521000785)

This repository is primarily written in R, and is stable for R versions 4.3.0 and 4.3.2. It may be possible to run alternative versions of R, but these have not been tested and verified. Configuration files are primarily written in YAML markup. The use of YAML files requires no additional software beyond the R packages on which this repositiory depends.

We invite any potential collaborators interested in expanding upon this analysis to fork this repository. Please submit a pull request if you feel the parent repository would benefit from any changes to code or configuration files.

#### Package and data dependencies

All R package dependencies will be automatically installed (if necessary) and loaded the first time the pipeline is run. If many packages need to be installed, the installation process could take up to one hour. See `dependencies.R` in the main directory for a full list of all R package dependencies, all of which are available on [CRAN](https://cran.r-project.org/). All of the data used to generate this analysis is publicly available. In some cases, this data is stored in the `/input/` directory, whilst in other cases data is pulled from the internet. As such, a live internet connection is required to run this analysis.

#### Computing resource requirements

No special cluster computing resources are required to run this analysis. Local parallelisation is used throughout the pipeline to improve runtime, however this parallelisation is only available on UNIX operating systems. Running this analysis on Windows operating systems is still possible, yet will be slower. Expect the full analysis (exluding any package installations) to take between 45 to 60 minutes on a UNIX machine and between 60 and 120 minutes on a Windows machine, dependant on machine specifications. R console output will keep the user updated with runtime progress.

## Directory structure

All code required to run the pipeline sits within the main repository directory (called `epi50` by default). Initially, three sub-directories exist within this main directory:

 1. A `/config/` directory containing a series of configuration files written in YAML.
 2. An `/input/` directory containing the input data required to run the pipeline.
 3. An `/extern/` directory containing processed results from the measles and polio models used in this analysis.

A fourth directory, `/output/`, is automatically created when the pipeline is first launched. All intermediary and final results are stored in this output folder. The figures presented in the peer-reviewed publication are stored within the `/output/figures/manuscript/` sub-directory. 

## Analysis configuration

The files contained in the `/config/` directory, all written in YAML mark up, configure the analysis to be run. These configuration files set key options such as diseases and vaccines to be modelled, countries to be run, and covariates data to be used for regression modelling. These configuration files also define data dictionaries for converting between various data sources and EPI variable naming conventions. Configured data sets include the WHO Immunization Information System (WIISE), World Population Prospects (WPP), Global Burden of Disease (GBD), and GapMinder. Several general algorithm options and model assumptions are also defined in `options.R`.

As standard, the configuration files in this repository will fully reproduce the analysis presented in the corresponding publication.

## Running the pipeline

The analysis pipeline consists of eight *modules*, indentified by numbers 1 to 8, each described below. In general, these modules should be run consecutively from 1 to 8. The module/s to be run are defined by the `run_module` variable in line 20 of `launch.R`. Use a single value to run/re-run a specific module, or use a vector to run/re-run a subset of modules or all modules. 

```{r}
# All of the following are are valid syntax for defining modules to be run (line 20, launch.R)
run_modules = 1
run_modules = c(1, 2, 4)
run_modules = 2 : 5
run_modules = 1 : 8
```

Modules are launched by running `launch.R`. Preferred usage is to 'source' this file (without 'echo' is ideal). When sourced, the current working directory is automatically reset the EPI50 repository. Alternative UNIX command line usage is to `cd` to the EPI50 repository then call `sh launch.sh`. By default all 8 modules will be run, in the process producing all outputs presented in the corresponding publication.

#### Module 1: Prepare inputs
The `run_prepare` module loads and formats modelling estimates from VIMC and GBD, vaccine coverage data from WHO (WIISE and SIA databases), demography data from WPP, and pubic health covariate data from GapMinder. Further, this module interprets all configuration files, using data dictionaries to convert variables to EPI50 naming conventions where appropriate.

#### Module 2: External models
The `run_external` module loads and formats outcomes from transmission models simulated outside of VIMC scope. Both measles and polio are considered 'externally modelled' for the purpose of this analysis. Only minor formatting jobs are applied to externally modelled diseases, as this analysis attempts to use the outcomes in their purest form.

Note that it is possible to simulate the [DynaMICE measles model](https://pubmed.ncbi.nlm.nih.gov/37474227/) directly from within this module. However this requires the user to clone another [github repository](https://github.com/ashattock/dynamice) and also have access to a computing cluster that uses the SLURM queueing system. 

#### Module 3: Static models
The `run_static` module produces vaccine impact estimates for diseases outside of VIMC scope and which haven't been 'externally modelled'. In this analysis, static modelling is used for diphtheria, tetanus, pertussis, and tuberculosis. Broadly, vaccine impact is calculated using a combination of GBD burden estimates, vaccine efficacy profiles, and vaccine coverage. A full description of the static modelling approach used in this analysis is provided in the Supplementary Material of the publication. 

#### Module 4: Regression (geographical-imputation)
The `run_regression("impute")` module uses time-series regression models to produce vaccine impact estimates for countries outside of VIMC scope, for diseases within VIMC scope. A series of models that vary in terms of predictor covariates are evaluated for each vaccine in each country, with the most parsimonious model selected. For countries without impact estimates, such estimates are imputed using the most commonly selected model in each WHO region.

#### Module 5: Impact functions (temporal-extrapolation)
The `run_impact` module seeks to determine the relationship between vaccine coverage and vaccine impact for each vaccine in each country. This is acheived by fitting four pre-defined statistical functions (straight line, exponential growth, logarithmic growth, sigmoidal growth) in cumulative space, and selecting the most parsimonious function.

#### Module 6: Historical evaluation
The `run_history` module compiles the final vaccine impact results by evaluating impact functions for each vaccine in each country using all coverage values. Temporal extrapolation takes effect at this stage, where vaccine coverage values that do not have corresponding impact estimates are evaluated.

#### Module 7: Regression (impact inference)
The `run_regression("infer")` module reapplies the regression models used for geographical imputation, but for all vaccines (including statically modelled and externally modelled diseases). The purpose of this additional regression step is to infer the predictors of vaccine impact for all diseases. Note that outcomes from this module were not presented in the final publication.

#### Module 8: Produce results
The `run_results` module produces all result and diagnostic figures presented in the publication.

## Authors

#### Development and maintenance:
* Andrew J. Shattock (shattocka@who.int)

#### Contributors:
* Helen C. Johnson
* Austin Carter
* So Yoon Sim
