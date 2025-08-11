#This script is designed to take in the output products of my fitting script chain, and use it to generate model visibilities, residual visibilities,
#and a CLEAN image of the residual visibilities.

#requires that you input the same ms file that you generated the uvtable from (split with keepFlags=False, uniform size, etc.)
#how can I call the same logger that I set up in the other python file? just call separately but don't overwrite??

###%%%###%%%###%%%###%%%###%%%###

###%%%###%%%### Package Imports ###%%%###%%%###
import os
import sys
import ast
import logging
import numpy as np
from scipy import ndimage
from astropy.io import fits
from astropy.wcs import WCS
from scipy.special import j1
from astropy import units as u
from uvplot import UVTable, COLUMNS_V0
from astropy.coordinates import SkyCoord
from astropy.nddata.utils import Cutout2D
from galario.double import sampleProfile, get_image_size

#I like when my matplotlib figures have a certain style so I set that here. 
from matplotlib import pyplot as plt
from matplotlib import rc
rc('text', usetex=True)
font = {'family' : 'serif',
        'weight' : 'bold',
        'size'   : '14'}
rc('font', **font)
###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###

###%%%###%%%### Fitting profiles ###%%%###%%%###  
def radial_gaussian_ring(peak, sigma, ring_rad, start, step, numsteps):
    radius = np.linspace(start, start + numsteps*step, numsteps)
    return peak*np.exp((-1/2)*((radius-ring_rad)/sigma)**2)

def radial_jinc(normalization, width, offset, start, step, numsteps):
    radius = np.linspace(start, start + numsteps*step, numsteps)
    return normalization*j1(width*radius)/radius + offset

def gaussian_ring(peak, sigma, ring_rad, radius):
    return 10**peak*np.exp((-1/2)*((radius-ring_rad)/sigma)**2)

def jinc(normalization, width, offset, radius):
    return 10**normalization*j1(width*radius)/radius + offset

def make_2d(paramdict, posfile, fitType):
    #The outputs of a model fitting are parameters that describe the profile, in physical units such as Jy/sr, arcseconds, etc. 
    #The model should be defined in the same units as the image is in. This means that I will do everything in arcseconds/pix. 
    target_hdr = fits.open(paramdict['dataFITSFile'])[0].header
    numpix = target_hdr['NAXIS2']
    pixscale = target_hdr['CDELT2']*3600 #put into arcseconds to match the model

    image_size = numpix * pixscale
    x = np.linspace(-image_size/2, image_size/2, numpix)
    y = np.linspace(-image_size/2, image_size/2, numpix)
    xx, yy = np.meshgrid(x, y) #note that this is now defined in arcseconds/pixel. 

    #Find the fitted values from the walkers end position
    endpos = np.loadtxt(posfile)
    ndim = endpos.shape[1]
    bestparams = [np.percentile(endpos[:,i], 50) for i in range(ndim)]

    if fitType=='gaussring':
        peak, sigma, ringrad, inc, pa, dRA, dDec = bestparams
        #The fitted parameters will come out in the "nice" units of degrees. Need to convert them to pixels for use in the model image creation. 
        degrees2rad = [np.deg2rad(item) for item in [inc, pa]]
        inc, pa = degrees2rad
        degrees2pix = [item*3600*pixscale for item in [dRA, dDec]]
        dRA, dDec = degrees2pix

        #Apply the dRA and dDec transforms. 
        xx_shifted = xx + dRA #positive dRA = move left
        yy_shifted = yy - dDec #positive dDEc = move up 

        #Apply the PA and inc rotations. 
        xxpa = xx_shifted*np.cos(pa) + yy_shifted*np.sin(pa)
        yypa = -xx_shifted*np.sin(pa) + yy_shifted*np.cos(pa)
        xxinc = xxpa/np.cos(inc)

        #Make the grid into a radius vector that we can then evaluate the model over. 
        radiusvec = np.hypot(xxinc, yypa)
        model_jysr = gaussian_ring(peak, sigma, ringrad, radiusvec)
    elif fitType=='jinc':
        normalization, width, offset, inc, pa, dRA, dDec = bestparams
        degrees = [np.deg2rad(item) for item in [inclination, posangle, dRA, dDec]]
        inclination, posangle, dRA, dDec = degrees

    #Apply the dRA and dDec transforms. 
        xx_shifted = xx + dRA #positive dRA = move left
        yy_shifted = yy - dDec #positive dDEc = move up 

        #Apply the PA and inc rotations. 
        xxpa = xx_shifted*np.cos(pa) + yy_shifted*np.sin(pa)
        yypa = -xx_shifted*np.sin(pa) + yy_shifted*np.cos(pa)
        xxinc = xxpa/np.cos(inc)

        #Make the grid into a radius vector that we can then evaluate the model over. 
        radiusvec = np.hypot(xxinc, yypa)
        model_jysr = jinc(normalization, width, offset, radiusvec)
    else:
        logging.warning('Please choose a valid fitting model, or add a new one into the code.')
        return None
    return model_jysr
###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###
###%%%###%%%### Dealing with visibilities function definitions ###%%%###%%%###
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
###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###

###%%%###%%%### Dealing with plotting images functions ###%%%###%%%###
#Because the model was created on a grid that is identical to that on which the fits image already exists, we don't have to worry about changing the image or pixel size. 
#The image is in units of Jy/bm though, which has to be converted over to Jy/sr to match the intensity scale of the model (and provide a unit better made for comparison)
def jybm_to_jysr(infile):
    hdr = fits.open(infile)[0].header
    bmaj = hdr['BMAJ']
    bmin = hdr['BMIN']

    omega_bm_deg2 = (np.pi/(2*np.log(2))) * bmaj*bmin
    omega_bm_sr = omega_bm_deg2 / (180/np.pi)**2
    return omega_bm_sr #divide the Jy/bm measurement by this factor

def img_prepper(fitsimg, center, box):
    im = fits.open(fitsimg)[0]
    im_wcs = WCS(im.header, naxis=2)
    conv = jybm_to_jysr(infile=fitsimg)
    im_plot = Cutout2D(im.data/conv, center, box*u.arcsecond, wcs=im_wcs)

    return im_plot.data, im_wcs

def total_prepper(fitsimg, model, center, box):
    if center.__class__.__name__!='SkyCoord': #needed this earlier when i thought I could pass this in as a SkyCoord object, but now I do the conversion manually. 
        logging.error('The image center parameter must be an astropy SkyCoord object. Please check this and run again.')
        return None
    img_final, img_wcs = img_prepper(fitsimg, center, box)
    model_final = Cutout2D(model, center, box, wcs=img_wcs).data
    return img_final, img_wcs, model_final

def plot_panel_pretty(ax, plotdata, label, vmin=None, vmax=None, cmap='inferno', text_color='w'):
    im = ax.imshow(plotdata, vmin=vmin, vmax=vmax, origin='lower', cmap=cmap)
    ax.text(5, 5, label, color=text_color, fontsize=16)
    #cbar = plt.colorbar(mappable=im, ax=ax, orientation='vertical', location='right', pad=0.05, shrink=0.8, aspect=15)
    return im#, cbar

def plot_img_model_residimg(paramdict, model_array, productprefix):
    center = SkyCoord(paramdict['plotCenterRA'], paramdict['plotCenterDec'], frame='icrs')
    data_img, data_wcs, model_plot = total_prepper(paramdict['dataFITSFile'], model_array, center, [paramdict['plotBoxSide']*u.arcsecond, paramdict['plotBoxSide']*u.arcsecond])

    fig = plt.figure(figsize=(12,5))

    ax1 = plt.subplot(131, projection=data_wcs)
    im1 = plot_panel_pretty(ax1, data_img, 'CLEAN Image', np.percentile(data_img, 1), np.percentile(data_img, 99.95))
    ax2 = plt.subplot(132, projection=data_wcs)
    plot_panel_pretty(ax2, model_plot, 'Model Sky Intensity', np.percentile(data_img, 1), np.percentile(data_img, 99.95))

    axes=[ax1, ax2]
    for ax in axes:
        if ax != ax1:
            ax.set_ylabel(r'', size=0)
            ax.coords[1].set_ticklabel(size=0)
        else:
            ax.set_ylabel(r"Declination (J2000)", size=22,labelpad=1)
        if ax != ax2:
            ax.set_xlabel(r'', size=0)
        else:
            ax.set_xlabel(r"Right Ascension (J2000)", size=22)

    cbar = fig.colorbar(im1, ax=fig.axes, orientation='horizontal', anchor=(0, 8.5),aspect=30)
    cbar.ax.set_title(paramdict['plotName']+' (Jy/sr)',fontsize=26,ha='center')
    cbar.ax.xaxis.set_ticks_position('top')
    cbar.ax.tick_params(direction='in',length=5,bottom=True,top=True)
    cbar.ax.xaxis.set_tick_params(labelsize=18)

    #plt.suptitle(paramdict['plotName'], fontsize=22)
    fig.subplots_adjust(hspace=0.1,wspace=0.05)
    plt.savefig(fname=productprefix+'_'+paramdict['plotName'].replace(' ', '_')+'_comparison.pdf', bbox_inches='tight')
###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###

###%%%###%%%### Define functions that take care of the model visibilities ###%%%###%%%### 
def sample_fitmodel(fittype, paramdict, posfile):
    #use the mean position value for each param at the end of the run to determine the 'most likely value' for that fit parameter
    endpos = np.loadtxt(posfile)
    ndim = endpos.shape[1]
    avgvalues = [np.percentile(endpos[:,i], 50) for i in range(ndim)]
    dRa=avgvalues[3]
    dDec=avgvalues[4]
    PA=avgvalues[5]
    inc=avgvalues[6]

    nxy, dxy, u, v, _, _, _ = initialize_data(paramdict['visFile'], paramdict['obsFreq'])
    

    if fittype=='gaussring':
        fitmodel = radial_gaussian_ring(avgvalues[0], avgvalues[1], avgvalues[2], paramdict['radiusStart'], paramdict['radiusStep'], paramdict['radiusNumSteps'])
        model_vis = np.array(sampleProfile(fitmodel, paramdict['radiusStart'], paramdict['radiusStep'], paramdict['radiusNumSteps'], nxy, dxy, u, v, dRA=avgvalues[5], dDec=avgvalues[4], PA=avgvalues[5], inc=avgvalues[6]), dtype=np.complex256) #specifying datatype here because I suspecting some rounding errors earlier and implemented this. it wasn't the fix, but it also isn't broken
    elif fittype=='jinc':
        fitmodel = radial_jinc(avgvalues[0], avgvalues[1], avgvalues[2], paramdict['radiusStart'], paramdict['radiusStep'], paramdict['radiusNumSteps'])
        model_vis = np.array(sampleProfile(fitmodel, paramdict['radiusStart'], paramdict['radiusStep'], paramdict['radiusNumSteps'], nxy, dxy, u, v, dRA=avgvalues[3], dDec=avgvalues[4], PA=avgvalues[5], inc=avgvalues[6]), dtype=np.complex256) #specifying datatype here because I suspecting some rounding errors earlier and implemented this. it wasn't the fix, but it also isn't broken
    else:
        logging.warning('Please choose a valid fitting model, or add a new one into the code.')
        return None
    return model_vis, PA, inc, dRa, dDec

def modelVdata_uvplot(paramdict, modelvis, pa, inc, dra, ddec, productprefix):
    _, _, u, v, datRe, datIm, w = initialize_data(paramdict['visFile'], paramdict['obsFreq'])
    modelvisRe = np.real(modelvis)
    modelvisIm = np.imag(modelvis)
    model_uvtable = UVTable(uvtable=[u, v, modelvisRe, modelvisIm, w], columns=COLUMNS_V0)
    data_uvtable = UVTable(uvtable=[u, v, datRe, datIm, w], columns=COLUMNS_V0)

    for uvtable in [model_uvtable, data_uvtable]:
        uvtable.apply_phase(dRA=dra, dDec=ddec)
        uvtable.deproject(inc=inc, PA=pa)

    #Now have the visibilities, u and v in units of wavelengths. Deprojected (I think that's a fair thing to do for uvdist/amp. Don't use that version when/if I do amp/u only plots)
    
    #Need to set bins such that I don't have a ridiculous number of points and plotting fails. Determine the upper bound and then choose bin sizes from there?
    uvdist = np.hypot(datRe, datIm)
    maxdist = np.max(uvdist)

    uvbin_size = paramdict['uvbinSize']
    if uvbin_size==0:
        if maxdist>=1e6:
            uvbin_size=50e3
        elif maxdist>=1e5:
            uvbin_size=25e3
        else:
            uvbin_size=10e3
    
    axes = data_uvtable.plot(label='Data', uvbin_size=uvbin_size, color='black')
    model_uvtable.plot(label='Model', uvbin_size=uvbin_size, axes=axes, yerr=False, linestyle='-', color='red')
    axes[0].figure.savefig(productprefix+'_'+paramdict['plotName']+'_uvplot.pdf')
###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###

###%%%###%%%### Define main function ###%%%###%%%### 
def visualize_main(paramsin, fitType, positionFile, productname):
    #Set the logger, but use the same filename as used in the main fitting script, in append mode. 
    logging.basicConfig(filename=productname+'.log', filemode='a', level=logging.INFO)

    #Evaluate the best fit model
    logging.info('Beginning to make the 2D model array...')
    model_array = make_2d(paramsin, positionFile, fitType)
    logging.info('Successfully generated the 2D model array.')

    #Create visualizations of the model compared to the image to help understand the fit
    logging.info('Beginning to make the comparison plot...')
    plot_img_model_residimg(paramsin, model_array, productname)
    logging.info('Successfully generated the comparison plot.')

    #Take care of generating the model visibilities, making a uvplot
    #THIS IS BROKEN THIS IS BROKEN DO NOT USE IT RIGHT NOW
    #fakevis, pa, inc, dra, ddec = sample_fitmodel(fitType, paramsin, positionFile)
    #modelVdata_uvplot(paramsin, fakevis, pa, inc, dra, ddec, productname)

###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###%%%###

###%%%###%%%### Call main, wrap to prevent issues with being called by run_fittings.py ###%%%###%%%### 
#This technically can be called as a script with command line arguments, or imported into the runfittings file now???? Decide which is better? Separate for ease?
if __name__=='__main__':
    visualize_main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])