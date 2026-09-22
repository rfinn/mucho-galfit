'''
GOAL: TRIM SGA2025 PSFs FROM 325x325 PIXELS TO 31x31 PIXELS
    - need W1, W3 PSF images in FITS format
    - x axis range will be 0 to 30, y axis range 0 to 30
        - center of psf should lie at (15, 15) -- must be at ~center!
Steps:
    - load .fz PSF image from parent directory
    - use astropy to clip to 
        - x1 = (im_length/2) - 15, 
        - x2 = (im_length/2) + 15,
        - y1 = (im_height/2) - 15,
        - y2 = (im_height/2) + 15.
    - for a 325 x 325 image, that is (148 to 178), (148 to 178)
'''

from astropy.io import fits
from astropy.nddata import Cutout2D
import numpy as np
import glob
import os

def load_fz(im_path):
    '''
    note: image data should be in hdu=1?
    '''
    return fits.open(im_path)


def get_psf_data(hdu):
    return hdu[1].data

def get_psf_header(hdu):
    return hdu[1].header


def trim_psf_data(data):
    '''
    we want 31x31 image, centered on the psf
    '''
    
    #center of psf is the current center of the image
    position = (len(data)/2, len(data)/2)
    
    #size of image centered on position
    size=31
    
    #trim...
    cutout = Cutout2D(data, position=position, size=size)
    
    return cutout


def save_trimmed_psf(cropped_data, header, output_destination):
    
    hdu = fits.PrimaryHDU(data=cropped_data, header=header)
    hdu.writeto(output_destination, overwrite=True)
    

def crop_save_psf(path_to_psf, psf_image_name, output_dir):    
    
    full_path = os.path.join(path_to_psf,psf_image_name)
    
    try:
        hdu = load_fz(full_path)
        data = get_psf_data(hdu)
        header = get_psf_header(hdu)
    except:
        print(f'{full_path} does not exist!')
        return
    
    output_path = os.path.join(output_dir,psf_image_name.replace('.fz',''))
    
    #for the moment, ignore cropping instructions for r-band PSF. just save.
    if ('W1' not in psf_image_name) and ('W3' not in psf_image_name):
        save_trimmed_psf(data, header, output_path)
        print('PSF saved!')
        return
    
    cropped_data = trim_psf_data(data)
    
    save_trimmed_psf(cropped_data, header, output_path)
    print('PSF trimmed and saved!')
    return