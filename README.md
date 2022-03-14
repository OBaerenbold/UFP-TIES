# Context

This is the code used to build and run the Nimble model for Source Apportionment, as
discussed in:

Baerenbold, O. *et al* (2022). 
A dependent Bayesian Dirichlet Process model for source apportionment of particle number size distribution.
Special Issue on Environmental Data Science, Environmetrics. 

All plotting routines and several sample data sets are also included here. The full output of the 
Nimble run as summarized and discussed in the paper is available on [Zenodo](https://doi.org/10.5281/zenodo.6352716).

# Installation and Setup

If running on Windows, you will need **RTools** installed and functioning properly
in order to build, link and run **Nimble**. You can install Nimble from CRAN
using

    install.packages("nimble")

and it has two dependencies: **coda** and **igraph**. **igraph** has several odd dependencies,
so we recommend installing the binary ersion for Windows rather than compiling it. Do compile
**nimble** - if both **code** and **igraph** are installed, it should work reasonably easily
with **Rtools**. These steps have been tested on both R 3.6.1 and 4.1.2 on Windows 7 and 10.

# Running a Test Version

All data needed for the models to load and run is contained in ./data/

In the ./code directory are a number of R scripts. The **run_nimble_model_noar_test.R**
script is a test routine which you can run which should ensure the directories are set
up correctly, the data is accessible, and the Nimble framework is working. This produces
a minimal example of the source apportionment model implementation, suitable for testing.
The full model can then be run as in the next. 

**Note**: the full script takes quite a while to run. The output of 
**run_nimble_model_noar_full.R** is what is contained in the Zenodo archive linked above.
You can download this file (900+MB), place it in the ./samples directory, and then 
the model diagnostics code at the end of the full.R file will work, and you can also
regenerate all of the plots presented in the paper.

The **out_ext_functions_base.R** file contains a set of utility functions used in the 
**create_plots.R** script.

All model outputs are saved in a ./samples directory in the parent directory of ./code
(the parent assumed to be the working directory), and plots are saved in ./plots 
in similar fashion. 

