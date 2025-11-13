'''
GOAL: CONVERT JM'S SGA2025 BITMASK TO A FORMAT COMPATIBLE WITH GALFIT.
    - need *ellipse-griz.fits for each galaxy
    - convert to bool mask
    - save as *-r-mask.fits
    -reproject onto WISE pixel scale
    -save *-wise-mask.fits
'''

import os
import sys
import numpy as np

from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from reproject import reproject_interp

homedir = os.getenv("HOME")


#for reference...just in case.
pixel_scale = {'FUV':1.5,'NUV':1.5,'g':0.262,'r':0.262,'z':0.262,'W1':2.75,'W2':2.75,'W3':2.75,'W4':2.75}
psf_oversampling = {'FUV':1,'NUV':1,'g':1,'r':1,'g':1,'W1':1,'W2':1,'W3':1,'W4':1}
mag_zeropoint = {'FUV':22.5,'NUV':22.5,'g':22.5,'r':22.5,'g':22.5,'W1':22.5,'W2':22.5,'W3':22.5,'W4':22.5}
image_resolution = {'FUV':6,'NUV':6,'g':1.5,'r':1.5,'z':1.5,'W1':6.1,'W2':6.4,'W3':6.5,'W4':12}
minmag2fit = {'FUV':10,'NUV':10,'g':17,'r':17,'z':17,'W1':10,'W2':10,'W3':10,'W4':10}


def remove_galaxy_SGA2020(maskfile):
    
    hdu = fits.open(maskfile)
    
    #remove 4096 (2**12, corresponds to pixel touching SGA galaxy (https://www.legacysurvey.org/dr10/bitmasks/)
    #change 4096 pixels to 0 so that they are ignored.
    hdu[0].data[hdu[0].data == 4096] = 0

    # write out updated mask
    hdu.writeto(maskfile,overwrite=True)

    
def remove_galaxy_SGA2025(maskfile):
    '''
    Open maskfile! For SGA2025, should be *-ellipse-griz.fits
    
    Note that the maskbits are:
        - 2^0 --> bright star
        - 2^1 --> Gaia star
        - 2^2 --> extended source (galaxy)
        - 2^3 --> SGA source
    
    We want to keep 2^0 and 2^1 pixel masks, as well as 2^2 pixel masks. 
    We want to REMOVE 2^3 and 2^2+2^3 pixel masks, as these are the 
        galaxies on which we will run GALFIT modeling.
    '''
    
    #read the maskfile
    hdu = fits.open(maskfile)
    
    #the maskbits are the second extension of the ellipse-griz fits.
    mask = hdu[1]
    
    #create a value flag
    source_flag = (mask.data == 2**2) | (mask.data == (2**2 + 2**3))
    
    #omit from the mask (i.e., set equal to zero)
    mask.data[source_flag] = 0
    
    mask.writeto(maskfile,overwrite=True)   #this will overwrite our copy of the *ellipse-griz.fits file

    
def reproject_mask(maskfile, reffile):
    '''
    maskfile: the mask to reproject
    reffile: the file with header information for reprojection! 
    
    e.g., if I want to reproject r-band onto WISE pixel scale, 
        -  maskfile should be r-band mask
        -  reffile should be WISE image
    
    This function will write the wisemask to the same directory as where reproject_mask() is run
    
    Output: output filename for the reprojected mask
    '''
    
    #remove galaxy pixels from mask...will save.
    remove_galaxy_SGA2025(maskfile)
    
    with fits.open(maskfile) as hmask, fits.open(reffile) as href:

        #reproject using HDU
        wisemask, footprint = reproject_interp(hmask[1], href[0].header)
        
        outname = maskfile.replace('r-mask', 'wise-mask')
        
        fits.writeto(outname, wisemask, href[0].header, overwrite=True)
    
    return outname      