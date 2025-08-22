#Note that in the emcee docs, it is reminded that any arguments passed to the lnpostfn function must be picklable, but this obviously
#can be quite demanding when the passed data is large. In order to fix issues of slowdown (by about a factor of 10!!!), I have switched the 
#code to work with GLOBALS so that there is less to pass between. This is super ugly and I hate it but it works. 
###%%%###%%%###%%%###%%%###%%%###

###%%%###%%%### Package Imports ###%%%###%%%###
import os
import sys
import ast
import time
import corner
import logging
import numpy as np
import pandas as pd
from multiprocessing import Pool
from emcee import EnsembleSampler as es
from galario.double import chi2Image, get_image_size
###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###

###%%%###%%%### Necessary to stop numpy from tripping over itself? ###%%%###%%%###
os.environ["OMP_NUM_THREADS"] = "1"
###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###

###%%%###%%%### Fitting profiles ###%%%###%%%### 
def gauss_ring(peak_ring, sigma_ring, rad_ring, radius):
    return 10**peak_ring * np.exp((-1/2)*((radius-rad_ring)/sigma_ring)**2)

def gauss2d(peak, sigmax, axisratio, x, y):
    sigmay = sigmax * axisratio
    return 10**peak * np.exp((-1/2) * ((x/sigmax)**2 + (y/sigmay)**2))

def add_ring(xx, yy, peak, sigma, rad, pa, inc):
    #xx and yy should be meshgrid objects? That represent the plane on which the model will be evaluated
    #apply PA and inc rotations for the ring
    xpa = xx*np.cos(pa) + yy*np.sin(pa)
    ypa = -xx*np.sin(pa) + yy*np.cos(pa)
    xinc = xpa/np.cos(inc)

    #make a radius vector, evaluate the ring model
    radius_vec = np.hypot(xinc, ypa)
    ring_model = gauss_ring(peak, sigma, rad, radius_vec)
    return ring_model

def add_blob(xx, yy, peak, sigma, axratio, pa, offRA, offDec):
    #xx and yy should be meshgrid objects? That represent the plane on which the model will be evaluated
    #Implement the offset (with respect to the centre of the ring)
    xoff = xx - offRA
    yoff = yy - offDec

    #apply PA rotation for the blob
    xpa = xoff*np.cos(pa) + yoff*np.sin(pa)
    ypa = -xoff*np.sin(pa) + yoff*np.cos(pa)

    #evaluate the model
    blob_model = gauss2d(peak, sigma, axratio, xpa, ypa)
    return blob_model

def gaussring_noblob(peak_r, sigma_r, rad_r, pa_r, inc_r, npix, pixscale):
    #Initialize the image plane
    image_size = npix * pixscale
    x = np.linspace(-image_size/2, image_size/2, npix)
    y = np.linspace(-image_size/2, image_size/2, npix)
    xx, yy = np.meshgrid(x, y)

    ring_model = add_ring(xx, yy, peak_r, sigma_r, rad_r, pa_r, inc_r)
    return ring_model

def gaussring_oneblob(peak_r, sigma_r, rad_r, pa_r, inc_r, peak_b, sigma_b, axrat_b, offRA_b, offDec_b, pa_b, npix, pixscale):    
    #Initialize the image plane
    image_size = npix * pixscale
    x = np.linspace(-image_size/2, image_size/2, npix)
    y = np.linspace(-image_size/2, image_size/2, npix)
    xx, yy = np.meshgrid(x, y)

    ring_model = add_ring(xx, yy, peak_r, sigma_r, rad_r, pa_r, inc_r)
    blob_model = add_blob(xx, yy, peak_b, sigma_b, axrat_b, pa_b, offRA_b, offDec_b)

    return ring_model + blob_model

def gaussring_twoblob(peak_r, sigma_r, rad_r, pa_r, inc_r, peak_ba, sigma_ba, axrat_ba, offRA_ba, offDec_ba, pa_ba, peak_bb, sigma_bb, axrat_bb, offRA_bb, offDec_bb, pa_bb, npix, pixscale):
    #Initialize the image plane
    image_size = npix * pixscale
    x = np.linspace(-image_size/2, image_size/2, npix)
    y = np.linspace(-image_size/2, image_size/2, npix)
    xx, yy = np.meshgrid(x, y)

    ring_model = add_ring(xx, yy, peak_r, sigma_r, rad_r, pa_r, inc_r)

    blobA = add_blob(xx, yy, peak_ba, sigma_ba, axrat_ba, pa_ba, offRA_ba, offDec_ba)
    blobB = add_blob(xx, yy, peak_bb, sigma_bb, axrat_bb, pa_bb, offRA_bb, offDec_bb)

    return ring_model + blobA + blobB

###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###

'''###%%%###%%%### Fitting profiles ###%%%###%%%###  #alternate versions where the unit conversions are performed in the MODEL (good??) (not tested)
def gauss_ring(peak_ring, sigma_ring, rad_ring, radius):
    return 10**peak_ring * np.exp((-1/2)*((radius-rad_ring)/sigma_ring)**2)

def gauss2d(peak, sigmax, axisratio, x, y):
    sigmay = sigmax * axisratio
    return 10**peak * np.exp((-1/2) * ((x/sigmax)**2 + (y/sigmay)**2))

def add_ring(xx, yy, peak, sigma, rad, pa, inc):
    #xx and yy should be meshgrid objects? That represent the plane on which the model will be evaluated
    #apply PA and inc rotations for the ring
    xpa = xx*np.cos(pa) + yy*np.sin(pa)
    ypa = -xx*np.sin(pa) + yy*np.cos(pa)
    xinc = xpa/np.cos(inc)

    #make a radius vector, evaluate the ring model
    radius_vec = np.hypot(xinc, ypa)
    ring_model = gauss_ring(peak, sigma, rad, radius_vec)
    return ring_model

def add_blob(xx, yy, peak, sigma, axratio, pa, offRA, offDec):
    #xx and yy should be meshgrid objects? That represent the plane on which the model will be evaluated
    #Implement the offset (with respect to the centre of the ring)
    xoff = xx - offRA
    yoff = yy - offDec

    #apply PA rotation for the blob
    xpa = xoff*np.cos(pa) + yoff*np.sin(pa)
    ypa = -xoff*np.sin(pa) + yoff*np.cos(pa)

    #evaluate the model
    blob_model = gauss2d(peak, sigma, axratio, xpa, ypa)
    return blob_model

def gaussring_noblob(peak_r, sigma_r, rad_r, pa_r, inc_r, npix, pixscale):
    arcsecs = [np.deg2rad(item/3600) for item in [sigma_r, rad_r]]
    sigma_r, rad_r = arcsecs
    degrees = [np.deg2rad(item) for item in [pa_r, inc_r]]
    pa_r, inc_r = degrees
    
    #Initialize the image plane
    image_size = npix * pixscale
    x = np.linspace(-image_size/2, image_size/2, npix)
    y = np.linspace(-image_size/2, image_size/2, npix)
    xx, yy = np.meshgrid(x, y)

    ring_model = add_ring(xx, yy, peak_r, sigma_r, rad_r, pa_r, inc_r)
    return ring_model

def gaussring_oneblob(peak_r, sigma_r, rad_r, pa_r, inc_r, peak_b, sigma_b, axrat_b, offRA_b, offDec_b, pa_b, npix, pixscale):
    #perform unit conversions. 
    arcsecs = [np.deg2rad(item/3600) for item in [sigma_r, rad_r, sigma_b, offRA_b, offDec_b]]
    sigma_r, rad_r, sigma_b, offRA_b, offDec_b = arcsecs
    degrees = [np.deg2rad(item) for item in [pa_r, inc_r, pa_b]]
    pa_r, inc_r, pa_b = degrees
    
    #Initialize the image plane
    image_size = npix * pixscale
    x = np.linspace(-image_size/2, image_size/2, npix)
    y = np.linspace(-image_size/2, image_size/2, npix)
    xx, yy = np.meshgrid(x, y)

    ring_model = add_ring(xx, yy, peak_r, sigma_r, rad_r, pa_r, inc_r)
    blob_model = add_blob(xx, yy, peak_b, sigma_b, axrat_b, pa_b, offRA_b, offDec_b)

    return ring_model + blob_model

def gaussring_twoblob(peak_r, sigma_r, rad_r, pa_r, inc_r, peak_ba, sigma_ba, axrat_ba, offRA_ba, offDec_ba, pa_ba, peak_bb, sigma_bb, axrat_bb, offRA_bb, offDec_bb, pa_bb, npix, pixscale):
    #perform unit conversions. 
    arcsecs = [np.deg2rad(item/3600) for item in [sigma_r, rad_r, sigma_ba, offRA_ba, offDec_ba, sigma_bb, offRA_bb, offDec_bb]]
    sigma_r, rad_r, sigma_ba, offRA_ba, offDec_ba, sigma_bb, offRA_bb, offDec_bb= arcsecs
    degrees = [np.deg2rad(item) for item in [pa_r, inc_r, pa_ba, pa_bb]]
    pa_r, inc_r, pa_ba, pa_bb = degrees
    
    #Initialize the image plane
    image_size = npix * pixscale
    x = np.linspace(-image_size/2, image_size/2, npix)
    y = np.linspace(-image_size/2, image_size/2, npix)
    xx, yy = np.meshgrid(x, y)

    ring_model = add_ring(xx, yy, peak_r, sigma_r, rad_r, pa_r, inc_r)

    blobA = add_blob(xx, yy, peak_ba, sigma_ba, axrat_ba, pa_ba, offRA_ba, offDec_ba)
    blobB = add_blob(xx, yy, peak_bb, sigma_bb, axrat_bb, pa_bb, offRA_bb, offDec_bb)

    return ring_model + blobA + blobB

###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###'''

#Note that the order of the function definitions is strange since this is the version of the code working with globals that I needed to just run

###%%%###%%%### Emcee evaluation functions ###%%%###%%%###
def lnpriorfn(p, par_ranges):
    for i in range(len(p)):
        if isinstance(par_ranges[i][0], str) or isinstance(par_ranges[i][1], str):
            logging.error("Bad param_ranges at index {i}: {par_ranges[i]}")
            raise ValueError("Non-numeric bounds detected in param_ranges")
        if p[i] < par_ranges[i][0] or p[i] > par_ranges[i][1]:
            return -np.inf
    return 0.0
###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###

###%%%###%%%### Functions to do all the work ###%%%###%%%###
def read_param_file(file):
    params={}
    with open(file, 'r') as f:
        for line in f:
            line = line.split("#", 1)[0].strip() #This will skip in-line comments
            if "=" in line:
                pkey, pvalue = line.split("=", maxsplit=1)
                pkey = pkey.strip()
                pvalue = ast.literal_eval(pvalue.strip())
                if isinstance(pvalue, list):
                    pvalue = np.array(pvalue, dtype=float)
                params[pkey] = pvalue
    return params

def initialize_data(datatable, frequency):
    u, v, re, im, w = np.require(np.loadtxt(datatable, unpack=True), requirements='C')
    wavelength = 299792458/frequency
    u /= wavelength
    v /= wavelength

    nx, dx = get_image_size(u, v)
    return  nx, dx, u, v, re, im, w

def save_results(fittype, sampler, pos, outname):
    np.savetxt(outname+'_pos.txt', pos, fmt='%10.6e', delimiter='\t')
    #I previously had been saving ALL the outputs, as I thought I would use them all. I did not. Uncomment and add prob, state to the args if you want to save these 
    #np.savetxt(outname+'_prob.txt', prob, fmt='%10.6e', delimiter='\t')
    #The state is a listbut each of the entries is a different type (string, numpy array) so I think I have to manually write to a file if we want to use it again later.
    #with open(outname+'_state.txt', 'w') as newfile:
    #    for entry in state:
    #        newfile.write(f"{entry}\n")
    
    ndim = pos.shape[1]
    try:
        samples_shape = sampler.chain.shape
        samples = sampler.chain.reshape((-1, ndim))
        np.savetxt(outname+'_chain.txt', samples, fmt='%10.6e', delimiter='\t', header=str(samples_shape))
    except:
        logging.warning('Saving the chain failed. Falling back to the documentation style of doing things.')
    logging.info('Successfully saved the sampler output variables to text files for later use.')

    samples = sampler.chain[:, -1000:, :].reshape((-1, ndim))

    if fittype=='noblob':
        label=["Peak", r"$\sigma$", r"R$_ring$", "PA", "Inc",r"$\Delta$RA", r"$\Delta$Dec"]
    elif fittype=='oneblob':
        label=["Peak R", r"$\sigma_R$", r"rad$_R$", "PA R", "Inc R", 
               "Peak B1", r"$\sigma_B1$", "Ax Ratio B1", r"$\Delta$RA B1", r"$\Delta$Dec B1", "PA B1", r"$\Delta$RA Net", r"$\Delta$Dec Net"]
    elif fittype=='twoblobs':
        label=["Peak R", r"$\sigma_R$", r"rad$_R$", "PA R", "Inc R", 
               "Peak B1", r"$\sigma_B1$", "Ax Ratio B1", r"$\Delta$RA B1", r"$\Delta$Dec B1", "PA B1", 
               "Peak B2", r"$\sigma_B2$", "Ax Ratio B2", r"$\Delta$RA B2", r"$\Delta$Dec B2", "PA B2" ,r"$\Delta$RA Net", r"$\Delta$Dec Net"]
    else:
        logging.warning('Please choose a valid fitting model, or add a new one into the code.')
        return None

    fig = corner.corner(samples, labels=label,
                    show_titles=True, quantiles=[0.16, 0.50, 0.84],
                    label_kwargs={'labelpad':20, 'fontsize':0}, fontsize=8)
    fig.savefig(outname+'_corner.png')
    logging.info('Successfully saved the corner plot')

def total_initializer(paramfile, ranges, guesses):
    #Declare that I'm going to make a bunch of global variables
    global nxy_global, dxy_global, u_global, v_global, Re_global, Im_global, w_global, ranges_global, paramdict_global, ndim_global, guesses_global

    #Do the parts of the code that I used to have at the top of main, but now I need the globals to exist before I can declare lnpostfn
    paramdict_global = read_param_file(paramfile)
    nxy_global, dxy_global, u_global, v_global, Re_global, Im_global, w_global = initialize_data(paramdict_global['visFile'], paramdict_global['obsFreq'])
    ranges_global = pd.read_csv(ranges, header=0, index_col=0).values
    guesses_global = pd.read_csv(guesses, header=0, index_col=0).values[:,0]
    ndim_global = len(ranges_global)

#Now that I can't have any arguments passed to lnpostfn, I must deal with each case with different profiles OUTSIDE of lnpostfn. Therefore all the definitions. 
def lnpostfn_noblob(params):
    lnprior = lnpriorfn(params, ranges_global)  # apply prior
    if not np.isfinite(lnprior):
        return -np.inf

    a, b, c, d, e, f, g= params

    arcsecs = [np.deg2rad(item/3600) for item in [b, c]]
    b, c = arcsecs
    degrees = [np.deg2rad(item) for item in [d, e]]
    d, e = degrees

    model_img = gaussring_noblob(a, b, c, d, e, nxy_global, dxy_global) 
    chi2 = chi2Image(model_img, dxy_global, u_global, v_global, Re_global, Im_global, w_global, dRA=f, dDec=g, origin='lower')

    return -0.5 * chi2 #+ lnprior

def lnpostfn_oneblob(params):
    lnprior = lnpriorfn(params, ranges_global)  # apply prior
    if not np.isfinite(lnprior):
        return -np.inf
    a, b, c, d, e, f, g, h, i, k, l, m, n= params
    #perform unit conversions. 
    arcsecs = [np.deg2rad(item/3600) for item in [b, c, g, i, k]]
    b, c, g, i, k = arcsecs
    degrees = [np.deg2rad(item) for item in [d, e, l]]
    d, e, l = degrees

    model_img = gaussring_oneblob(a, b, c, d, e, f, g, h, i, k, l, nxy_global, dxy_global)
    chi2 = chi2Image(model_img, dxy_global, u_global, v_global, Re_global, Im_global, w_global, dRA=m, dDec=n, origin='lower')
    return -0.5 * chi2 #+ lnprior

def lnpostfn_twoblobs(params):
    lnprior = lnpriorfn(params, ranges_global)  # apply prior
    if not np.isfinite(lnprior):
        return -np.inf
    a, b, c, d, e, f, g, h, i, k, l, m, n, o, p, q, s, t, z= params

    #perform unit conversions. 
    arcsecs = [np.deg2rad(item/3600) for item in [b, c, g, i, k, n, p, q]]
    b, c, g, i, k, n, p, q= arcsecs
    degrees = [np.deg2rad(item) for item in [d, e, l, s]]
    d, e, l, s = degrees


    model_img = gaussring_twoblob(a, b, c, d, e, f, g, h, i, k, l, m, n, o, p, q, s, nxy_global, dxy_global)
    chi2 = chi2Image(model_img, dxy_global, u_global, v_global, Re_global, Im_global, w_global, dRA=t, dDec=z, origin='lower')
    return -0.5 * chi2 #+ lnprior

def main(paramfile, fitType, ranges, guesses, productname): 
    #Set up the logger
    logging.basicConfig(filename=productname+'.log', filemode='w', level=logging.INFO)
    #the initializer function must be called outside the Pool, as well as in the initializer argument of the Pool. 
    total_initializer(paramfile, ranges, guesses)

    with Pool(initializer=total_initializer, initargs=(paramfile, productname, ranges, guesses)) as pool:
        if fitType=='noblob':
            sampler=es(paramdict_global['numWalkers'], ndim_global, lnpostfn_noblob, pool=pool)
        elif fitType=='oneblob':
            sampler=es(paramdict_global['numWalkers'], ndim_global, lnpostfn_oneblob, pool=pool)
        elif fitType=='twoblobs':
            sampler=es(paramdict_global['numWalkers'], ndim_global, lnpostfn_twoblobs, pool=pool)
        pos = [guesses_global + 1e-5*np.random.randn(ndim_global) for i in range(paramdict_global['numWalkers'])]
        logging.info('Beginning the burn-in...')
        bi_start = time.time()
        pos0, _, _ = sampler.run_mcmc(pos, paramdict_global['burnTime'])
        sampler.reset()
        bi_end=time.time()
        logging.info("Duration of the burn-in is {0:.1f} seconds".format(bi_end-bi_start))
        logging.info('Beginning the fitting run...')
        fit_start=time.time()
        pos1, _, _ = sampler.run_mcmc(pos0, paramdict_global['runTime']) #put in prob1, state1 if you want to save them and use them for anything else?
        fit_end=time.time()
        logging.info("Duration of the fitting run is {0:.1f} seconds".format(fit_end-fit_start))

    save_results(fitType, sampler, pos1, productname)

###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###

###%%%###%%%### Call main, wrap in name=main to prevent recursion ###%%%###%%%###
if __name__=='__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])