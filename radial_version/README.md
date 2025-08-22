## Fitting code - Radial Version

These scripts are designed to work together and be used on the CANFAR science platform to conduct visibility modelling of galaxies. This version of the tool utilizes galario's sampleProfile and chi2Profile functions, and therefore only models axisymmetric profiles defined along a radial vector. Currently, only the gaussian ring and jinc function profiles have been implemented, but it is relatively straightforward to add in other profiles by following the formatting used for the existing two. 

This tool is designed to run fittings, but also to explore how the initial guess provided to the sampler affects the convergence of the sampling in the case of fitting visibilities from galactic observations. For each fitting profile specified a 'column' is created, and the first 'row' of the explored grid uses the provided initial guess for parameters. In each subsequent 'row', each initial parameter position in the starting vector is perturbed by a random percentage within $\pm 5\%$ of the explorable range of that parameter. The size of this perturbation can be adjusted on line 88 of the `launch_headless_sessions.py` file. _If only a single fit is needed, then the number of rows can be set to one_. 

A number of fitting parameters are common between fittings in each row of the explored parameter grid, and these are set in the `global_fitting_params.txt` file. The parameters are labelled and most have a short explanation to help in correctly setting them. 

In order to interact with CANFAR and manage the headless sessions running these fits, the skaha package is used. This requires that the first script, `launch_headless_sessions.py`, be run in an environment with python $\geq 3.9$, and the skaha, numpy, and pandas packages installed. 

In order to generate visibilities from a model of sky brightness, the galario package is used. This requires that the fitting script, `run_fittings.py`, be run in an environment with python $\leq 3.8$, and the galario, corner, numpy, pandas, emcee, scipy, astropy, uvplot, and matplotlib packages installed. 

After these two environments have been configured, **please adjust the paths in the following lines**:

**launch_headless_sessions.py:** lines 36, 56, 65

**launch_fittings.sh:** line 2 

**global_fitting_params.txt:** lines 1, 13

Because this tool runs each fitting in a headless session on CANFAR, the instruction to begin the fitting must be passed in the `cmd` argument as a shell script. The `launch_fittings.sh` script acts as an interface between the session launching and parameter assignment that happens in the environment with python $\geq 3.9$, and the fitting procedures which are completed in the environment with python $=3.8$.

The file `visualizeModel.py` contains several functions used to provide a quick visualization of the quality of the fit, which are called in `run_fittings.py`. If a FITS image of the data you are modelling is provided, then a figure comparing the image to the model is created. A plot of Re(V) and Im(V) with uvdistance is also generated using the uvplot package. An example of this comparison plot is provided. 

![alt text](https://github.com/emcarver/ALMA-Galaxy-Visibility-Modelling/edit/main/radial_version/example_plot.png "Example of the image to model comparison plot")

Each of these files also contain comments, and `launch_headless_sessions.py` has a header that goes into detail to the flow of the tool, similar to what I have written here. 

***IN PROGRESS. UNFINISHED. KEEP WRITING.*** 
