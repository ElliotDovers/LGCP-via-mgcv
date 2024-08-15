# LGCP-via-mgcv
This repo accompanies the manuscript "Fitting Log-Gaussian Cox Processes Using Generalized Additive Model Software" published in The American Statistician:\
[https://www.tandfonline.com/doi/full/10.1080/00031305.2024.2316725](https://www.tandfonline.com/doi/full/10.1080/00031305.2024.2316725)\
The code can be used to replicate both the simulations and real data analysis therein.

## LGCP via mgcv - Supplementary.pdf contains appendices:
- S1 contains additional simulation settings and results.
- S2 contains a vignette for analyzing the gorillas nesting data presented in the manuscript - with code to obtain the freely available dataset.

## LGCP via mgcv - Code.zip is directory of all code used for simulations and analyses presented:
This includes the following folders:

### SIMULATIONS
This folder contains the R scripts and other files used to run the simulation study.
NOTE: Simulations were run in parallel on a High Performance Computing Cluster.
Please run "sim_analysis_iterated.R" for a demonstration of the simulation process and analysis (this has been set to "job=1" for such a demonstration)

### APPLICATION
This folder contains the R script to analyze the gorillas nesting data.

see additional README files within subfolders for further details.
