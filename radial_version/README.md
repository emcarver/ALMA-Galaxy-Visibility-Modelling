## Fitting code - Radial Version

These scripts are designed to work together and be used on the CANFAR science platform to conduct visibility modelling of galaxies. This version of the tool utilizes galario's sampleProfile and chi2Profile functions, and therefore only models axisymmetric profiles defined along a radial vector. Currently, only the gaussian ring and jinc function profiles have been implemented, but it is relatively straightforward to add in other profiles. 

This tool is designed to run fittings, but also to explore how the initial guess provided to the sampler affects the convergence of the sampling in the case of fitting visibilities from galactic observations. For each fitting profile specified a 'column' is created, and the first 'row' of the explored grid uses the provided initial guess for parameters. In each subsequent 'row', each initial parameter position in the starting vector is perturbed by a random percentage within $\pm 5\%$ of the explorable range of that parameter. 

In order to interact with CANFAR and manage the headless sessions running these fits, the skaha package is used. This requires that the first script, `launch_headless_sessions.py`, be run in an environment with python $\geq$ 3.9, and the skaha, numpy, and pandas packages installed. 

In order to generate visibilities from a model of sky brightness, the galario package is used. This requires that the fitting script, `run_fittings.py`, be run in an environment with python $\leq$ 3.8, and the galario, corner, numpy, pandas, emcee, scipy, astropy, uvplot, and matplotlib packages installed. 

After these two environments have been configured, please adjust the paths in the following places:

**launch_headless_sessions.py:** lines 36, 56, 65

**launch_fittings.sh:** line 2 

**global_fitting_params.txt:** lines 1, 13

***IN PROGRESS. UNFINISHED. KEEP WRITING.***
