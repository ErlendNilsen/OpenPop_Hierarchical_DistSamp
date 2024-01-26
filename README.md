
[![License: AGPL
v3](https://img.shields.io/badge/License-AGPL_v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)

# An Integrated Open Population Distance Sampling Model

## What is in this repository?
This repository contains code for a workflow analysing line transect data
using an integrated distance sampling model (IDSM). 
The model utilizes age-structured survey data and auxiliary
data from marked individuals to jointly estimate changes in population density and
temporal variation in underlying demographic rates (recruitment rate and
survival probability). It is a multi-area model, meaning it simultaneously
models processes across a defined number of areas, and shares information
both across space and time. 

The IDSM workflow is set up as a targets pipeline (https://books.ropensci.org/targets/)
and was specifically written for data collected through the [Norwegian
monitoring program for tetraonid
birds](https://honsefugl.nina.no/Innsyn/en) (mainly Willow Ptarmigan
*Lagopus lagopus*), but can be used for other systems that collect
age-structured distance sampling data.

The model itself is written and implemented in NIMBLE 
(see [here](https://r-nimble.org/) for more information about NIMBLE for R). 
As of v.2.0, the NIMBLE model code is written by a function called "writeModelCode.R"
in the "R" folder. 
Previous versions relied on NIMBLE code found in the "NIMBLE code" folder. 

Additional R functions used for downloading and wrangling the data, simulating data, 
preparing data in correct format, setting up and running the model, as well 
as post-hoc analyses and plotting and quality control of results is contained in 
the R folder. Refer to each function's roxygen documentation for detailed documentation. 

The complete workflows for analysing simulated data, real data from the 
Lierne area only, and data from all areas for which public data on willow
ptarmigan in Norway are available are provided in the following master
scripts: 

- "Analysis_SimData.R" for a single simulated dataset
- "Analysis_SimData_Replicates.R" for multiple simulated datasets including replicate runs
- "Analysis_RealData_LierneVest.R" for real data from Lierna area only
- "Analysis_RealData.R" for real data from all areas
- "_targets.R" for real data from all area, organized in a targets pipeline.

Note that only the multi-area workflows ("Analysis_RealData.R" and "_targets.R") have been
further developed in the the v2.0 release, i.e. use a newer and reparameterized version
of the model (see "R/writeModelCode.R"). 
Workflows for Lierne only and for simulated data are fully functional, but
use the original parameterisation of the model (see "NIMBLE Code/"). For these latter 
workflows, we therefore recommend using v1.3 of the code. 
  
## Additional dependencies
There are a three dependencies that need to be manually installed to
run the workflow. 
First, you need to install NIMBLE (follow instructions given here: [https://r-nimble.org/download](https://r-nimble.org/download)). 
Second, the analysis uses code from the nimbleDistance package (https://github.com/scrogster/nimbleDistance). 
to estimate the half normal detection distribution. 
Third, we use the LivingNorwayR package (https://livingnorway.github.io/LivingNorwayR/) to download and
wrangle the line transect distance sampling survey data.

Finally, running the workflow requires access to additional data (radio-telemetry data on ptarmigan, 
rodent occupancy data, and shapefiles for municipalities in Norway). 
These can be downloaded from OSF: [https://osf.io/7326r/](https://osf.io/7326r/).

## Citation
Releases of this repository are archived and issued DOIs via Zenodo: https://zenodo.org/records/10462269
The citation of the latest version is: 

ChloeRNater, ErlendNilsen, christofferhohi, Matthew Grainger, & Bernardo Brandão Niebuhr. (2024). 
ErlendNilsen/OpenPop_Integrated_DistSamp: Ptarmigan IDSM v2.0 (v.2.0). 
Zenodo. https://doi.org/10.5281/zenodo.10462269

## Contact
Erlend Nilsen: erlend.nilsen@nina.no

Chloé Nater: chloe.nater@nina.no
