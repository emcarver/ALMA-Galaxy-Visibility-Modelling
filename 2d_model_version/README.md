## Fitting code - Two-Dimensional Model Version

These scripts are designed to work together and be used on the CANFAR science platform to conduct visibility modelling of galaxies. The general workflow is to generate a model and generate visibilities using the uv sampling of the observed data, then find $\chi^2$ of this model. This is used to guide an MCMC sampler, which conducts a fitting by finding the most likely value for each model parameter. After this, a corner plot showing the probability distribution of fitting parameters is generated.

This version of the tool utilizes galario's sampleImage and chi2Image functions, and therefore can work with any two-dimensional model. Keep in mind that this model has to be repeatedly generated for the purpose of fitting though, so analytic models will be fastest. Currently, the gaussian ring model as well as a gaussian with one or two additional gaussian components ("blobs") have been implemented as models, but it is relatively straightforward to add in other models by following the formatting used for the existing three. 

Overall, this tool is incomplete and does not run as expected, but it is provided for reference. The galario sampleImage and chi2Image functions are more computationally expensive than their profile counterparts, and the model itself is simply a larger array. As such, these fits took a very long time to run. Part of the issue was in the fact that arguments to the lnpostfn function have to be passed between processes, and when this data is large, there is a slowdown. The ugly fix to this is to implement global parameters and use the initializer in the Pool, which allowed for an almost 10x increase in speed on my personal laptop. However, the fitting models are not implemented quite correctly (specifically when it comes to the intensity scaling) so this code is incomplete. I also implemented all the model making with numpy, but I did not explore any other options for model generation that might be simpler or faster. 

This tool is designed to run fittings, but also to explore how the initial guess provided to the sampler affects the convergence of the sampling in the case of fitting visibilities from galactic observations. For each fitting profile specified a 'column' is created, and the first 'row' of the explored grid uses the provided initial guess for parameters. In each subsequent 'row', each initial parameter position in the starting vector is perturbed by a random percentage within $\pm 5$% of the explorable range of that parameter. The size of this perturbation can be adjusted on line 90 of the `launch_headless_fittings.py` file. If only a single fit is needed, then the number of rows can be set to one. 

A number of fitting parameters are common between fittings in each row of the explored parameter grid, and these are set in the `global_fitting_params.txt` file. The parameters are labelled within the file.

In order to interact with CANFAR and manage the headless sessions running these fits, the skaha package is used. This requires that the first script, `launch_headless_fittings.py`, be run in an environment with python $\geq 3.9$, and the skaha, numpy, and pandas packages installed. 

In order to generate visibilities from a model of sky brightness, the galario package is used. This requires that the fitting script, `run_fittings.py`, be run in an environment with python $\leq 3.8$, and the galario, corner, numpy, pandas, emcee, scipy, astropy, uvplot, and matplotlib packages installed. 

After these two environments have been configured, **please adjust the paths in the following lines**:

- launch_headless_fittings.py: lines 27, 58, and 67

- launch_fittings.sh: line 2 

- global_fitting_params.txt: line 1

Because this tool runs each fitting in a headless session on CANFAR, the instruction to begin the fitting must be passed in the `cmd` argument as a shell script. The `launch_fittings.sh` script acts as an interface between the session launching and parameter assignment that happens in the environment with python $\geq 3.9$, and the fitting procedures which are completed in the environment with python $\leq 3.8$.
