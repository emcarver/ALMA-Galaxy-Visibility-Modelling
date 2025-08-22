#!/bin/sh
/arc/home/ecarver/.conda/envs/galario_casa/bin/python /arc/projects/uvdisk_fit/approved_code_only/run_fittings.py $1 $2 $3 $4 $5

#Please put the paths to a python environment that is py3.8 or lower, that has the following packages (plus their dependencies) installed:
#corner
#numpy
#pandas
#emcee
#galario