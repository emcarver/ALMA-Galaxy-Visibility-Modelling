###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%
#In its current form, this script is structured so that for each row of each fit assigned in model_fits, the initial guess parameter is perturbed by a random amount in \pm 5% 
#of the total value range. By changing this initial guess parameter, my hope is to find a better understanding of the errors on these fittings, and the fragility of these convergences.

#In order to launch multiple headless sessions in a nice manner, you must use skaha, which requires python>=3.9. Galario, used for the fitting, requires python<=3.8. 
#Therefore, a series of .py and .sh scripts must be used to go about launching multiple headless sessions in parallel to use in fitting. 
#launch_headless_sessions.py -> this file dictates the number of fittings to run, the fittypes, the allowed parameter ranges, and interacts with skaha.Session. 
#                               In the session.create call, the file launch_fittings.sh is passed as the command
#launch_fittings.sh ----------> this file allows the fittings to be run, and specifically in an environment where galario (and other fitting dependencies) are installed. 
#                               One of the arguments passed in this file is run_fittings.py
#run_fittings.py -------------> this file actually implements all the model equations, passes them into galario, and runs the MCMC sampling to fit the models to the observed visibilities.
#                               This file utilizes multiprocessing, so that the fits can converge more quickly. In combination with the fact that several fits can be performed 
#                               at once by utilizing the structure in launch_headless_sess.py, this should allow the parameter space to be explored more efficiently. 
#                               One of the arguments passed to this script is global_params.txt
#                               #Aside on multiprocessing: I have found that passing NO arguments to the lnpostfn function in the sampler results in a MAJOR speedup (~10x on my personal computer), but this means using a ton of global variables which is not fantastic. If you need a major speedup, feel free to try it out. 
#global_params.txt -----------> this file contains a handful of parameters global to all fits run, such as the visibilities table and observing frequency.

#Again, please note that this should be run in an environment py>=3.9, with skaha, numpy, and pandas installed. 
#In launch_fittings.sh, please provide the path to an environment py<=3.8, with galario, casa, corner, numpy, pandas, emcee, and matplotlib installed.
###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%

import os
import numpy as np
import pandas as pd
#from skaha.session import Session
import subprocess

#Some intermediate and product files will have to be generated. Define a directory to put them in.
productdir='/Users/school/Desktop/Research/2025/final_code_radial/testproducts/'

#Name the types of fits you want to perform + iterate over. 
model_fits = ['gaussring'] #'gaussring', 'jinc'

#Define the explorable ranges for all the parameters in each model. These are set somewhat narrow. Adjust as needed?
#Note that passing an array in a command line argument is difficult, and a 2d array is just awful. 
#These will get converted to csv files and those files will be read in as arguments instead. They will be (re-)written in the first row of the grid (in case bounds change between runs).
model_fits_ranges = {'gaussring': [[-5, 15], [0, 10],[0, 20], [0, 90], [0, 180],[-5, 5], [-5, 5]],
                     'jinc': [[0, 15], [0, 20],[-1, 3], [0, 90], [0, 180],[-2, 2], [-2, 2]]}

#Define the initial guess for fitting to begin at for each model. 
model_fits_initialguesses = {'gaussring': [6, 1, 7, 60, 15, 0, 0],
                             'jinc': [4, 5, 0.5, 25, 18, 0, 0]}

#Define the output product name. Suffixes such as _rowx_corner.png will be added
model_fits_outputs = {'gaussring': productdir+'gaussring_testtwo',
                     'jinc': productdir+'jinc'}

#Path to the param file
global_paramfile='/Users/school/Desktop/Research/2025/final_code_radial/global_fitting_params.txt'

#choose number of times to perturb the initial guess
gridrows=1

#set the session computing parameters
cores=8
mem=6
image='images.canfar.net/skaha/astroml:latest'
cmd = '/Users/school/Desktop/Research/2025/final_code_radial/launch_fittings.sh'

for fit in model_fits:
    for row in range(gridrows):
        if row==0:
            #Save the array arguments to text files for easier passing. Do this only once for each fit type. 
            pd.DataFrame(np.array(model_fits_ranges[fit])).to_csv(productdir+fit+'_param_ranges.csv')
            pd.DataFrame(np.array(model_fits_initialguesses[fit])).to_csv(productdir+fit+'_param_guess_row'+str(row)+'.csv')
            arglist = [global_paramfile, fit, productdir+fit+'_param_ranges.csv', productdir+fit+'_param_guess_row'+str(row)+'.csv', model_fits_outputs[fit]+'_row'+str(row)]
            args = ' '.join(arglist)
            subprocess.call(['sh', cmd, args])
            '''
            session = Session()
            session_id = session.create(name=fit+'-row'+str(row),
                                        image=image,
                                        cores=cores,
                                        ram=mem,
                                        kind='headless',
                                        cmd=cmd,
                                        args=args,
                                        env={'sessiontype':'headless'})
            print("Session ID: {}".format(session_id))
            '''
        else:
            
            prange = np.array(model_fits_ranges[fit])
            #In order to evaluate the reliance of the fit on the initial guess, I will adjust each entry in the parameter grid by a random percentage of the parameter range (pm 5%)
            adjust = [np.random.uniform(-0.05, 0.05, size=1)[0]*(prange[i,1] - prange[i,0]) for i in range(prange.shape[0])]
            guess = np.array(model_fits_initialguesses[fit]) + adjust
            pd.DataFrame(guess).to_csv(productdir+fit+'_param_guess_row'+str(row)+'.csv')
            #in testing, I have found that sometimes the order of the sessions get mixed up. Failsafe check if the ranges csv exists yet?
            if not os.path.exists(productdir+fit+'_param_ranges.csv'):
                pd.DataFrame(np.array(model_fits_ranges[fit])).to_csv(productdir+fit+'_param_ranges.csv')
            arglist = [global_paramfile, fit, productdir+fit+'_param_ranges.csv', productdir+fit+'_param_guess_row'+str(row)+'.csv', model_fits_outputs[fit]+'_row'+str(row)]
            args = ' '.join(arglist)
            subprocess.call(['sh', cmd, args])
            '''
            session = Session()
            session_id = session.create(name=fit+'-row'+str(row),
                                        image=image,
                                        cores=cores,
                                        ram=mem,
                                        kind='headless',
                                        cmd=cmd,
                                        args=args,
                                        env={'sessiontype':'headless'})
            print("Session ID: {}".format(session_id))
            '''

#Notes to self:
#monitor the headless session events and logs to look for successful launches.
# https://ws-uv.canfar.net/skaha/v0/session/[id]?view=events
# https://ws-uv.canfar.net/skaha/v0/session/[id]?view=logs