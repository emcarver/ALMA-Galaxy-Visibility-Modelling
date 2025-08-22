###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%
#In its current form, this script is structured so that for each row of each fit assigned in model_fits, the initial guess parameter is perturbed by a random amount in \pm 5% 
#of the total value range. By changing this initial guess parameter, my hope is to find a better understanding of the errors on these fittings, and the fragility of these convergences.

#In order to launch multiple headless sessions in a nice manner, you must use skaha, which requires python>=3.9. Galario, used for the fitting, requires python<=3.8. 
#Therefore, a series of .py and .sh scripts must be used to go about launching multiple headless sessions in parallel to use in fitting. 
#launch_headless_sess.py -----> this file dictates the number of fittings to run, the fittypes, the allowed parameter ranges, and interacts with skaha.Session. 
#                               In the session.create call, the file launch_fittings.sh is passed as the command
#launch_fittings.sh ----------> this file allows the fittings to be run, and specifically in an environment where galario (and other fitting dependencies) are installed. 
#                               One of the arguments passed in this file is run_fittings.py
#run_fittings.py -------------> this file actually implements all the model equations, passes them into galario, and runs the MCMC sampling to fit the models to the observed visibilities.
#                               This file (attempts) to utilize multiprocessing, so that the fits can converge more quickly. In combination with the fact that several fits can be performed 
#                               at once by utilizing the structure in launch_headless_sess.py, this should allow the parameter space to be explored more efficiently. 
#                               One of the arguments passed to this script is global_params.txt
#global_params.txt -----------> this file contains a handful of parameters global to all fits run, such as the visibilities table and observing frequency.

#Again, please note that this should be run in an environment py>=3.9, with skaha, numpy, pandas, and uvplot installed. 
#In launch_fittings.sh, please provide the path to an environment py<=3.8, with galario, casa, corner, numpy, pandas, schwimmbad, emcee, and matplotlib installed.
###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%

import os
import numpy as np
import pandas as pd
from skaha.session import Session

#Some intermediate and product files will have to be generated. Define a directory to put them in.
productdir='/arc/projects/uvdisk_fit/'

#Name the types of fits you want to perform + iterate over. 
model_fits = ['noblob', 'oneblob'] #'noblob', 'oneblob', 'twoblobs'

#Define the explorable ranges for all the parameters in each model. These are set somewhat narrow. Adjust as needed?
#Note that passing an array in a command line argument is difficult, and a 2d array is just awful. Especially with variable size between fit models.
#These will get converted to txt files and those files will be read in as arguments instead. They will be (re-)written in the first row of the grid (in case bounds change between runs).
#note that according to the galario docs, the intensities in a 2d model should be in units of Jy/pix. But in testing, this doesn't make a lot of sense with what I've seen? 
#Ranges defined here are quite wide but this is something I wasn't able to go and investigate more deeply.
model_fits_ranges = {'noblob': [[-8, 8], [0, 5], [0, 30], [0, 180], [0, 90], [-1, 1], [-1, 1]],
                     'oneblob': [[-8, 8], [0, 5], [0, 30], [0, 180], [0, 90], 
                                 [-8, 8], [0, 5], [0, 3], [-30, 30], [-30, 30], [0, 180], [-1, 1], [-1,1]],
                     'twoblobs': [[-8, 8], [0, 5], [0, 30], [0, 180], [0, 90], 
                                  [-8, 8], [0, 3], [0, 3], [-30, 30], [-30, 30], [0, 180], 
                                  [-8, 8], [0, 3], [0, 3], [-30, 30], [-30, 30], [0, 180], [-1, 1], [-1,1]]}

#Define the initial guess for fitting to begin at for each model. 
model_fits_initialguesses = {'noblob': [-4, 1, 7, 13, 60, 0, 0],
                             'oneblob': [-4, 1, 7, 13, 60, 
                                         -4, 1, 1, 0.25, 6.5, 0, 0, 0],
                             'twoblobs': [-2, 1, 7, 13, 60, 
                                          -2, 1, 1, 0.25, 6.5, 0, 
                                          -2, 1, 1, 0, -7, 0, 0, 0]}

#Define the output product name. Suffixes such as _rowx_corner.png will be added
model_fits_outputs = {'noblob': productdir+'noblobs',
                    'oneblob': productdir+'oneblob',
                    'twoblobs': productdir+'twoblobs'}

#Path to the param file
global_paramfile='/arc/projects/uvdisk_fit/approved_code_only/global_fitting_params.txt'

#choose number of times to perturb the initial guess
gridrows=2

#set the session computing parameters
cores=16
mem=8
image='images.canfar.net/skaha/astroml:latest'
cmd = '/arc/projects/uvdisk_fit/approved_code_only/launch_fittings.sh'

for fit in model_fits:
    for row in range(gridrows):
        if row==0:
            #Save the array arguments to text files for easier passing. Do this only once for each fit type. 
            pd.DataFrame(np.array(model_fits_ranges[fit])).to_csv(productdir+fit+'_param_ranges.csv')
            pd.DataFrame(np.array(model_fits_initialguesses[fit])).to_csv(productdir+fit+'_param_guess_row'+str(row)+'.csv')
            arglist = [global_paramfile, fit, productdir+fit+'_param_ranges.csv', productdir+fit+'_param_guess_row'+str(row)+'.csv', model_fits_outputs[fit]+'_row'+str(row)]
            args = ' '.join(arglist)
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
        else:
            prange = np.array(model_fits_ranges[fit])
            #In order to evaluate the reliance of the fit on the initial guess, I will adjust each entry in the parameter grid by a random percentage of the available range (pm 5%)
            adjust = [np.random.uniform(-0.05, 0.05, size=1)[0]*(prange[i,1] - prange[i,0]) for i in range(prange.shape[0])]
            guess = np.array(model_fits_initialguesses[fit]) + adjust
            pd.DataFrame(guess).to_csv(productdir+fit+'_param_guess_row'+str(row)+'.csv')
            #in testing, I have found that 
            if not os.path.exists(productdir+fit+'_param_ranges.csv'):
                pd.DataFrame(np.array(model_fits_ranges[fit])).to_csv(productdir+fit+'_param_ranges.csv')
            arglist = [global_paramfile, fit, productdir+fit+'_param_ranges.csv', productdir+fit+'_param_guess_row'+str(row)+'.csv', model_fits_outputs[fit]+'_row'+str(row)]
            args = ' '.join(arglist)
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

#Notes to self:
#monitor the headless session events and logs to look for successful launches.
# https://ws-uv.canfar.net/skaha/v0/session/[id]?view=events
# https://ws-uv.canfar.net/skaha/v0/session/[id]?view=logs