Rtools needs to be installed and properly linked to run Nimble
The folder is built as follows:

- \code includes all code needed to run the analysis
- Run a basic mcmc including diagnostics on all or parts of the data with run_nimble_model_base
- create output and plots with create_plots
- out_ext_functions_base.R is a set of functions needed for create_plots.R
- run_nimble_model_rX.R are full runs of different chain lengths that I am currently running in parallel because it is not quite clear how many iterations are needed for convergence

- model output is saved in \samples as part of the run_nimble_model_base with naming depending on model specification and time
- plots are saved in a subfolder of \plots named using the naming of the saved samples
- data needed for the project are in \data
- all other explanations are included as comments. Let me know if something is unclear.