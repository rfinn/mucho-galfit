#!/usr/bin/env python

'''
GOAL:
* set up directories for running galfit on WISESize sample

PROCEDURE:
* read .fits tables

* for every galaxy, create a directory

* unpack and move band images to this directory

* generate unWISE psfs, move to this directory

USAGE:
* on draco, move to /mnt/astrophysics/mg-output-wisesize/

* call as:
    python /mnt/astrophysics/wisesize/github/wisesize/mucho-galfit/setup_galfit.py

    * will create directories of galaxies. each directory will contain the image, noise, unWISE psfs

'''

import os
import sys
import glob
from astropy.io import fits
import numpy as np
from astropy.table import Table

#import create objid function
from merge_ns_catalogs import create_OBJIDs

#importing mask util functions...
sys.path.insert(0,'utils')
from convert_mask import reproject_mask
from galaxies_in_fov import get_galaxies_in_fov


##########################################################################     
### FUNCTIONS
##########################################################################     

#functions to change .fits.fz to .fits
    
def funpack_image(input_, output):
    hdu_list = fits.open(input_)
    #find first HDU with data
    for h in hdu_list:
        if h.data is not None:
            fits.writeto(output, h.data, header=h.header, overwrite=True)
            hdu_list.close()
            return
    print(f"Warning: no data found in {input_}")
    hdu_list.close()

def funpack_all(start_dir, output_dir):
    for filename in os.listdir(start_dir):
        if '.fz' in filename:
            funpack_image(start_dir+filename, output_dir+filename.replace('.fz',''))


#prepare the masks!
def move_masks(start_dir, output_dir, wise_image_file):
    '''
    Run AFTER the *maskbits.fits.fz is converted to .fits and relocated to the OBJIDxxxxx directory
    '''
    
    for filename in os.listdir(output_dir):
        
        if 'maskbits' in filename:

            rmask_file = output_dir+filename.replace('maskbits','image-r-mask') #replace 'maskbits' with 'image-r-mask' in mask filename (the image- is for cadence purposes)
            os.system(f'cp {filename} {rmask_file}') #make copy with new filename!
            
            #takes r-band mask (maskfile), converts to wise mask (reffile header) with the name outname
            #ALSO removes 4096 bitmask -- the SGA galaxy! -- from the mask.
                        
            reproject_mask(rmask_file, wise_image_file)
            
            return
    
    
#convert invvar image to noise
def convert_invvar_noise(invvar_image, noise_image):
    import warnings
    
    # read in invvar image
    # print('invvar image = ',invvar_image, os.path.basename(invvar_image))
    hdu = fits.open(invvar_image)
    data = hdu[0].data
    header = hdu[0].header
    hdu.close()
    
    warnings.simplefilter("ignore", RuntimeWarning)
    
    # operate on pixels to take sqrt(1/x)
    noise_data = np.sqrt(1/data)
    
    # check for bad values, nan, inf
    # set bad pixels to very high value, like 1e6
    noise_data = np.nan_to_num(noise_data,nan=1.e6)

    # write out as noise image
    fits.writeto(noise_image,noise_data,header=header,overwrite=True)

    
def radec_to_groupname(ra, dec, prefix=''):
    
    #36-arcsec precision (0.01 degrees)
    ra1, dec1 = np.atleast_1d(ra), np.atleast_1d(dec)
    
    group_name = '{}{:05d}{}{:04d}'.format(
            prefix, int(100*ra1.item()), 'm' if dec1 < 0 else 'p',
            int(100*np.abs(dec1).item()))

    return group_name  


def get_wise_psfs(param_dict, path_to_image_dir):

    main_dir = param_dict['main_dir']
    outdir = main_dir+param_dict['path_to_images']
    
    #load functions from pull_unwise_psfs (read tile table, get galaxy image's coadd_id, pull associated psf for W1-4
    #and save to path_to_image_dir
    sys.path.append(main_dir+'github/wisesize/unwise_PSFs/')
    from pull_unwise_psfs import read_tiles, get_coadd_id, pull_unwise_psf
    
    #read in tiles. the rest will follow at the bottom of the forthcoming loop
    tile_path = main_dir+param_dict['tile_path']   #contains COADD IDs and the RA+DEC of tile centers
    tile_table = read_tiles(tile_path)  
    
    #directory of (primary) galaxy
    path_to_image_dir = outdir+obj_id+'/'

    #get coadd id of (primary) galaxy image
    coadd_id = get_coadd_id(tile_path, path_to_image_dir, tile_table=tile_table)

    #pulls psf for W1-4 bands
    for band in range(1,5):
        pull_unwise_psf(path_to_image_dir, coadd_id, band)
    
    
#path_to_repos e.g., /mnt/astrophysics/wisesize/
def get_images(objid,ra,dec,output_loc,data_root_dir):
    ###############################################################################
    ### GET IMAGES
    ###############################################################################

    print(objid)
    
    #output_loc is the directory holding the individual galaxy output directories (which GALFIT will be pulling from!)
    #e.g., /mnt/astrophysics/wisesize/mg_output_wisesize/OBJ10000/
    output_dir = os.path.join(output_loc,objid+'/')
    if not os.path.exists(output_dir):
        print("making the output directory ",output_dir)
        os.mkdir(output_dir)

    #data_root_dir is where JM's input_ images are initially stored
    if not os.path.exists(data_root_dir):
        print(f"could not find data_root_dir {data_root_dir} - exiting")
        sys.exit()
    
    #just pulling the RA directory name...extracts integer from ra, then puts in xxx format (three integer places)
    #ra_slice = f'{np.trunc(ra):03.0f}'
    ra_slice = f'{int(ra):03d}'
        
    if dec>32.:   #if DEC>32 degrees, then galaxy is in "north" catalog. else, south catalog.
        data_dir = f'{data_root_dir}dr11-north/{ra_slice}/'
    if dec<32.:
        data_dir = f'{data_root_dir}dr11-south/{ra_slice}/'
        
    if not os.path.exists(data_dir):
        print(f"could not find data_dir {data_dir} - exiting")
        sys.exit()

    group_name = radec_to_groupname(ra, dec, prefix='')
    
    data_dir = data_dir+group_name+'/'
    
    funpack_all(data_dir, output_dir)
    
    #masks! rename maskbits to rband mask, remove 4096 (SGA galaxy) mask; create WISE mask
    wise_image = glob.glob(f'{output_dir}*-image-W3.fits')[0]   #will output the image path+filename
    move_masks(data_dir, output_dir, wise_image) 
    
    #move Legacy Survey Viewer JPG image (if it exists)
    try:
        ls_im = glob.glob(f'{data_dir}*image.jpg')[0]
        os.system(f'cp {ls_im} {output_dir}')
    except:
        print(f'LS Viewer image not found in {data_dir}. Skipping.')
    
    #define invvar image names; if the std does not exist, then convert invvar to std and save to output_dir
    for bandpass in ['g','r','z','W1','W2','W3','W4']:
        invvar_image = f'SGA2025_{group_name}-invvar-{bandpass}.fits'

        # check if noise image exists in output_dir, if not make it from invvar 
        sigma_image = invvar_image.replace('invvar','std')
        if not os.path.exists(output_dir+sigma_image):
            convert_invvar_noise(os.path.join(output_dir,invvar_image),os.path.join(output_dir,sigma_image))


    ###############################################################################
    ### END GET IMAGES
    ###############################################################################


##########################################################################     
### END FUNCTIONS
##########################################################################     

if __name__ == '__main__':
    
    
    ######################
    ### Parameter File ###
    ######################
    
    param_file = '/mnt/astrophysics/wisesize/github/mucho-galfit/paramfile.txt'
        
    #create dictionary with keyword and values from param textfile
    param_dict={}
    with open(param_file) as f:
        for line in f:
            try:
                key = line.split()[0]
                val = line.split()[1]
                param_dict[key] = val
            except:
                continue
    
    main_dir = param_dict['main_dir']
    path_to_pyscripts = main_dir+param_dict['path_to_scripts']
    
    outdir = main_dir+param_dict['path_to_images']
    data_root_dir = param_dict['data_root_dir']
    
    main_catalog_path = param_dict['main_catalog']
    
    objid_col = param_dict['objid_col']

    group_name_col = param_dict['group_name_col']
    objname_col = param_dict['objname_col']

    maintab = Table.read(param_dict['main_catalog'])

    primary_group_col = param_dict['primary_group_col']
    
    ###########################################
    # Check if main catalog has OBJID column. #
    # If not, create one (and save result)! #
    ###########################################
    
    if objid_col not in maintab.columns:
        maintab = create_OBJIDs(maintab)
        maintab.write(main_catalog_path, overwrite=True)
    
    ############################
    # Isolate Primary Galaxies #
    ############################
    
    #trim to only include primary galaxies; if no such flag exists, assume all galaxies are primary.
    try:
        primary_flag = maintab[primary_group_col]
    except:
        primary_flag = np.ones(len(maintab),dtype=bool)   #all true

    maintab = maintab[primary_flag]
        
    #############################
    # Create Primary Galaxy Dirs #
    ##############################
    
    #check that outdir (where the individual primary galaxy directories will live) exists! if not, create it.
    if os.path.exists(outdir):
        os.chdir(outdir)
    else:
        os.system(f'mkdir {outdir}')
    
    # for each primary galaxy, create a directory
    for i in range(len(maintab)):
        
        obj_id = maintab[objid_col][i]
        ra = maintab['RA'][i]
        dec = maintab['DEC'][i]
        objname = maintab[objname_col][i]
        #group_name = etab[group_name_col][i] # this is either the objname, or objname_GROUP for groups
        
        # if etab[primary_group_col][i] & (etab[group_mult_col][i] > 0): # make directory for primary targets
        # galpath = outdir+etab[objid_col][i]
        path_to_image_dir = outdir+obj_id+'/'
        
        # make directory if it doesn't already exist
        if not os.path.exists(path_to_image_dir):
            os.mkdir(path_to_image_dir)
        os.chdir(path_to_image_dir)

        #copy images
        get_images(obj_id,ra,dec,outdir,data_root_dir)
        
        #get galaxes in FOV, save to galsFOV.txt in path_to_image_dir
        get_galaxies_in_fov(maintab, path_to_image_dir)
        
        ############
        ### PSFs ###
        ############
        
        #get_wise_psfs(param_dict, path_to_image_dir)

        # for testing
        #if i == 1:
        #    os.chdir(outdir)
        #    sys.exit()
        
    os.chdir(outdir)   #return to the main output directory
    