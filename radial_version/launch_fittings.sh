#!/bin/sh
/arc/home/ecarver/.conda/envs/galario_casa/bin/python /arc/projects/uvdisk_fit/radial_fit_code/run_fittings.py $1 $2 $3 $4 $5
#This should follow the format:
#/complete/path/to/env /complete/path/to/run_fittings.py $1 $2 $3 $4 $5
#Please provide the path to an environment py<=3.8, with galario, casa, corner, numpy, pandas, emcee, scipy, astropy, uvplot, and matplotlib installed.