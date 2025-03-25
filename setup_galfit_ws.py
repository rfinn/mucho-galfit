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
    python /mnt/astrophysics/kconger_wisesize/github/wisesize/mucho-galfit/setup_galfit.py

    * will create directories of galaxies. each directory will contain the image, noise, unWISE psfs

'''

import os
import sys
from astropy.io import fits
import numpy as np
from astropy.table import Table

##########################################################################     
### FUNCTIONS
##########################################################################     

#functions to change .fits.fz to .fits
def funpack_image_cfitsio(input,output):
    command = 'funpack -O {} {}'.format(output,input)
    print(command)
    os.system(command)
def funpack_image(input,output,nhdu=1):
    hdu = fits.open(input)
    print('input file = ',input)
    fits.writeto(output,data=hdu[nhdu].data, header=hdu[nhdu].header, overwrite=True)
    hdu.close()
    #print('finished unpacking image')

#unpack composite images into their constituent wavelength bands
#save results in the same spot as where JM images are stored? (path_to_im=data_root_dir)
def extract_bands(path_to_im,im_name=None,grz=False,WISE=False):
    if grz:
        ims=fits.open(path_to_im+im_name)
        header_im = ims[0].header
        header_invvar = ims[1].header
        g_im, r_im, z_im = ims[0].data[0], ims[0].data[1], ims[0].data[2]
        g_inv, r_inv, z_inv = ims[1].data[0], ims[1].data[1], ims[1].data[2]
        ims.close()
        
        grz_im_names = [im_name.replace('.fits','-im-g.fits'),im_name.replace('.fits','-im-r.fits'),
                        im_name.replace('.fits','-im-z.fits'),im_name.replace('.fits','-invvar-g.fits'),
                        im_name.replace('.fits','-invvar-r.fits'),im_name.replace('.fits','-invvar-z.fits')]
        grz_ims = [g_im,r_im,z_im,g_inv,r_inv,z_inv]
        
        for i,filename in enumerate(grz_im_names):
            header = header_im if i<=2 else header_invvar
            fits.writeto(path_to_im+filename,grz_ims[i],header=header,overwrite=True)
        
    if WISE:
        ims=fits.open(path_to_im+im_name)
        header_im = ims[0].header
        header_invvar = ims[1].header
        w1_im, w2_im, w3_im, w4_im = ims[0].data[0], ims[0].data[1], ims[0].data[2], ims[0].data[3]
        w1_inv, w2_inv, w3_inv, w4_inv = ims[1].data[0], ims[1].data[1], ims[1].data[2], ims[1].data[3]
        ims.close()
        
        wise_im_names = [im_name.replace('-unwise.fits','-im-W1.fits'),im_name.replace('-unwise.fits','-im-W2.fits'),
                        im_name.replace('-unwise.fits','-im-W3.fits'),im_name.replace('-unwise.fits','-im-W4.fits'),
                        im_name.replace('-unwise.fits','-invvar-W1.fits'),im_name.replace('-unwise.fits','-invvar-W2.fits'),
                        im_name.replace('-unwise.fits','-invvar-W3.fits'),im_name.replace('-unwise.fits','-invvar-W4.fits')]
        wise_ims = [w1_im,w2_im,w3_im,w4_im,w1_inv,w2_inv,w3_inv,w4_inv]
        
        for i,filename in enumerate(wise_im_names):
            header = header_im if i<=2 else header_invvar
            fits.writeto(path_to_im+filename,wise_ims[i],header=header,overwrite=True)
    
#convert invvar image to noise
def convert_invvar_noise(invvar_image, noise_image):
    # read in invvar image
    print('invvar image = ',invvar_image, os.path.basename(invvar_image))
    hdu = fits.open(invvar_image)
    data = hdu[0].data
    header = hdu[0].header
    hdu.close()
    # operate on pixels to take sqrt(1/x)
    noise_data = np.sqrt(1/data)
    
    # check for bad values, nan, inf
    # set bad pixels to very high value, like 1e6
    noise_data = np.nan_to_num(noise_data,nan=1.e6)

    # write out as noise image
    fits.writeto(noise_image,noise_data,header=header,overwrite=True)
    
#path_to_repos e.g., /mnt/astrophysics/kconger_wisesize/
def get_images(objid,ra,dec,output_loc,data_root_dir):
    ###############################################################################
    ### GET IMAGES
    ###############################################################################

    #output_loc is the directory holding the individual galaxy output directories (which GALFIT will be pulling from!)
    #e.g., /mnt/astrophysics/kconger_wisesize/mg_output_wisesize/OBJ1000/
    output_dir = os.path.join(output_loc,objid+'/')
    if not os.path.exists(output_dir):
        print("making the output directory ",output_dir)
        os.mkdir(output_dir)

    #data_root_dir is where JM's input images are initially stored
    if not os.path.exists(data_root_dir):
        print(f"could not find data_root_dir - exiting")
        sys.exit()

    ra_val = str(int(ra)) if len(str(int(ra)))==3 else '0'+str(int(ra))
    dec_val = str(int(dec)) if len(str(int(dec)))==2 else '0'+str(int(dec))
    
    if dec>32.:   #if DEC>32 degrees, then galaxy is in "north" catalog
        data_dir = f'{data_root_dir}/dr9-north/native/{ra_val}/'
    if dec<32.:
        data_dir = f'{data_root_dir}/dr9-south/native/{ra_val}/'
    
    
    if not os.path.exists(data_dir):
        print(f"could not find data_dir - exiting")
        sys.exit()

    print("source directory for JM images = ",data_dir)
    
    #np.modf()[0] isolates the decimals
    #str() converts to string
    #+'000' ensures that there are at least 4 decimal places (including the required 1 from np.modf()
    #[1:6] isolates '.xxxx', where xxxx are the 4 decimal places 
    ra_string = ra_val + (str(np.modf(ra)[0])+'000')[1:6]    
    dec_string = dec_val + (str(np.modf(dec)[0])+'000')[1:6]
    
    im_name_grz = f'SGA2025_J{ra_string}+{dec_string}.fits'
    extract_bands(data_dir,im_name=im_name_grz,grz=True)
    
    im_name_wise = im_name_grz.replace('.fits','-unwise.fits')
    extract_bands(data_dir,im_name=im_name_wise,WISE=True)
    
    for bandpass in ['r','g','z','W1','W2','W3','W4']:
        image = f'{im_name_grz.replace('.fits','')}-im-{bandpass}.fits'
        invvar_image = f'{im_name.replace('.fits','')}-invvar-{bandpass}.fits'
        #psf_image = f'{group_name}-custom-psf-{bandpass}.fits.fz'

        os.system(f'cp {image} {output_dir}')
        os.system(f'cp {invvar_image} {output_dir}')
        
        # created images
        sigma_image = invvar_image.replace('invvar','std')

        # check if noise image exists, if not make it from invvar    
        if not os.path.exists(sigma_image):
            convert_invvar_noise(os.path.join(output_dir,invvar_image),os.path.join(output_dir,sigma_image))
        

    ###############################################################################
    ### END GET IMAGES
    ###############################################################################


##########################################################################     
### END FUNCTIONS
##########################################################################     

if __name__ == '__main__':
    
    param_file = '/mnt/astrophysics/kconger_wisesize/github/mucho-galfit/paramfile.txt'
        
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
    
    main_catalog = param_dict['main_catalog']
    #phot_catalog = param_dict['phot_catalog']
    
    objid_col = param_dict['objid_col']
    #primary_group_col = param_dict['primary_group_col']
    #group_mult_col = param_dict['group_mult_col']
    group_name_col = param_dict['group_name_col']
    objname_col = param_dict['objname_col']

    try:
        #etab = Table.read(phot_catalog)
        maintab = Table.read(main_catalog)        
    
    except FileNotFoundError:
        print("ERROR: problem locating catalogs - exiting")
        sys.exit()
    
    #load functions from pull_unwise_psfs (read tile table, get galaxy image's coadd_id, pull associated psf for W1-4
    #and save to path_to_image_dir
    sys.path.append(main_dir+'github/wisesize/unwise_PSFs/')
    from pull_unwise_psfs import read_tiles, get_coadd_id, pull_unwise_psf
        
    #read in tiles. the rest will follow at the bottom of the forthcoming loop
    tile_path = main_dir+param_dict['tile_path']   #contains COADD IDs and the RA+DEC of tile centers
    tile_table = read_tiles(tile_path)    
    
    #check that outdir exists! if not, create it.
    if os.path.exists(outdir):
        os.chdir(outdir)
    else:
        os.system(f'mkdir {outdir}')
    
    # for each galaxy, create a directory (and write sourcelist?)
    for i in range(len(maintab)):
        
        obj_id = maintab[objid_col][i]
        ra = maintab['RA'][i]
        dec = maintab['DEC'][i]
        objname = maintab[objname_col][i]
        #group_name = etab[group_name_col][i] # this is either the objname, or objname_GROUP for groups
        
        #if etab[primary_group_col][i] & (etab[group_mult_col][i] > 0): # make directory for primary targets
        #galpath = outdir+etab[objid_col][i]
        path_to_image_dir = outdir+obj_id+'/'
        # make directory if it doesn't already exist
        if not os.path.exists(path_to_image_dir):
            os.mkdir(path_to_image_dir)
        os.chdir(path_to_image_dir)

        #sourcefile = path_to_image_dir+'/{}sourcelist'.format(maintab['VFID'][i])
        #sourcelist = open(sourcefile,'w')
        ## write out one line with OBJID, objname, RA, DEC, wavelength
        #output_string = maintab['OBJID'][i] + ' ' + maintab['objname'][i] + ' ' + str(maintab['RA'][i]) + ' ' + str(maintab['DEC'][i]) + ' ' + str(wavelength) + ' \n'.format()
        #sourcelist.write(output_string)
        #sourcelist.close()

        #copy images
        get_images(obj_id,ra,dec,outdir,data_root_dir)
        
        ############
        ### PSFs ###
        ############
        
        #directory of (primary) galaxy
        path_to_image_dir = outdir+obj_id+'/'
        
        #get coadd id of (primary) galaxy image
        coadd_id = get_coadd_id(tile_path, path_to_image_dir, tile_table=tile_table)
        
        #pulls psf for W1-4 bands
        for band in range(1,5):
            pull_unwise_psf(path_to_image_dir, coadd_id, band)
        
        #try:
        #    #runs code to pull unwise psf images and save to outdir/objid   
        #    os.system(f'python3 pull_unwise_psfs.py -objid {objid}')    
        #except:
        #    os.system(f'python pull_unwise_psfs.py -objid {objid}')
       
        # for testing
        #if i > 1:
        #    os.chdir(outdir)
        #    sys.exit()
        
    os.chdir(outdir)