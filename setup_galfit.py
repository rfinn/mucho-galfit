#!/usr/bin/env python
'''
GOAL:
* set up directories for running galfit on the Virgo sample

PROCEDURE:
* read the vf_v2_legacy_ephot.fits
* for all galaxies with a PRIMARY flag, create a directory

USAGE:

* On the virgo vms, move to /mnt/astrophysics/muchogalfit-output
 
* Call as:

python ~/github/muchogalfit/setup_galfit.py 

-will create directories for each primary galaxy, with each directory containing the image, invvar, and psf

* then run the parallel processing on grawp

NOTES: 
* 2023-06-20 : RF updating the program to run on the full virgo sample.  The program 
was originally created to run on Kim's sample. 

* 2024-03-05 : RF rewriting to move all the paths
'''
import os
import sys

##########################################################################     
### FUNCTIONS
##########################################################################     


def funpack_image(input,output,nhdu=1):
    command = 'funpack -O {} {}'.format(output,input)
    print(command)
    os.system(command)
def funpack_image_astropy(input,output,nhdu=1):
    from astropy.io import fits
    hdu = fits.open(input)
    print('input file = ',input)
    #print(hdu.info())
    fits.writeto(output,data=hdu[nhdu].data, header=hdu[nhdu].header, overwrite=True)
    hdu.close()
    print('finished unpacking image')

def convert_invvar_noise(invvar_image, noise_image):
    # convert invvar image to noise
    from astropy.io import fits
    import numpy as np
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
    
def get_images(vfid,ra,objname,group_name,aphys='/mnt/astrophysics/rfinn/'):
    ###############################################################################
    ### GET IMAGES
    ###############################################################################
    
    # set up path name for image directory
    # directory where galaxy images are
    output_dir = os.path.join(aphys,'muchogalfit-output',vfid+'/')
    if not os.path.exists(output_dir):
        print("making the output directory ",output_dir)
        os.mkdir(output_dir)
    data_root_dir = '/mnt/virgofilaments-data/' #directory on draco or virgo vms
    if not os.path.exists(data_root_dir):

        data_root_dir = '/mnt/astrophysics/virgofilaments-data/'
        if not os.path.exists(data_root_dir):
            print(f"could not find data_dir - exiting")
            sys.exit()
    
    data_dir = f'/mnt/astrophysics/virgofilaments-data/{int(ra)}/{group_name}/'
    if not os.path.exists(data_dir):
        data_dir = f'/mnt/virgofilaments-data/{int(ra)}/{group_name}/'
        if not os.path.exists(data_dir):
            
            print(f"could not find data_dir - exiting")
            sys.exit()

    data_dir = f'{data_root_dir}/{int(ra)}/{group_name}/'
    print("source directory for JM images = ",data_dir)

    for bandpass in ['r','g','z','W1','W2','W3','W4']:
        image = f'{group_name}-custom-image-{bandpass}.fits.fz'
        invvar_image = f'{group_name}-custom-invvar-{bandpass}.fits.fz'    
        psf_image = f'{group_name}-custom-psf-{bandpass}.fits.fz'
    
        # need to copy image from data directory if it doesn't exist
        if not os.path.exists(image.replace('.fz','')):
            funpack_image(os.path.join(data_dir,image),os.path.join(output_dir,image.replace('.fz','')))
        else:
            print("image is in place - no need to recopy")
        image = image.replace('.fz','')


        if not os.path.exists(invvar_image.replace('.fz','')):
            funpack_image(os.path.join(data_dir,invvar_image),os.path.join(output_dir,invvar_image.replace('.fz','')),nhdu=0)        
        invvar_image = invvar_image.replace('.fz','')

        if not os.path.exists(psf_image.replace('.fz','')):
            funpack_image(os.path.join(data_dir,psf_image),os.path.join(output_dir,psf_image.replace('.fz','')),nhdu=0)        
        psf_image = psf_image.replace('.fz','')

        # created images
        sigma_image = f'{group_name}-custom-std-{bandpass}.fits'

    
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


    from astropy.table import Table

    # define environment variable so funpack can find the correct variables
    os.environ["LD_LIBRARY_PATH"]="/opt/ohpc/pub/compiler/gcc/9.4.0/lib64:/home/siena.edu/rfinn/software/cfitsio-4.2.0/"

    # read in virgo catalog
    tablename = 'vf_v2_legacy_ephot.fits'
    maintablename = 'vf_v2_main.fits'
    # directory structure on grawp
    catalogdir='/mnt/astrophysics/rfinn/catalogs/Virgo/v2/'
    try:
        
        etab = Table.read(catalogdir+tablename)
        maintab = Table.read(catalogdir+maintablename)        
        
        aphys = '/mnt/astrophysics/rfinn/'
    except FileNotFoundError:
        try:
            # test to see if running on Virgo VMS or draco
            catalogdir='/mnt/astrophysics/catalogs/Virgo/v2/'
            etab = Table.read(catalogdir+tablename)
            maintab = Table.read(catalogdir+maintablename)        
            aphys = '/mnt/astrophysics/'        
        except FileNotFoundError: # 
            print("ERROR: problem locating astrophysics drive - exiting")
            sys.exit()




    outdir = os.path.join(aphys,'muchogalfit-output/')
    os.chdir(outdir)
    # for each galaxy, create a directory and write sourcelist

    # could run in parallel by just sending in the line number
    line_number = sys.argv[1]
    for i in [line_number]:
        if etab['GROUP_PRIMARY'][i] & (etab['GROUP_MULT'][i] > 0): # make directory for primary targets
            galpath = outdir+etab['VFID'][i]

            # make directory if it doesn't already exist
            if not os.path.exists(galpath):
                os.mkdir(galpath)
            os.chdir(galpath)
            #sourcefile = galpath+'/{}sourcelist'.format(maintab['VFID'][i])
            #sourcelist = open(sourcefile,'w')
            ## write out one line with VFID, objname, RA, DEC, wavelength
            #output_string = maintab['VFID'][i] + ' ' + maintab['objname'][i] + ' ' + str(maintab['RA'][i]) + ' ' + str(maintab['DEC'][i]) + ' ' + str(wavelength) + ' \n'.format()
            #sourcelist.write(output_string)
            #sourcelist.close()

            # copy images and funpack
            vfid = maintab['VFID'][i]
            ra = maintab['RA'][i]
            dec = maintab['DEC'][i]
            objname = maintab['objname'][i]
            group_name = etab['GROUP_NAME'][i] # this is either the objname, or objname_GROUP for groups

            get_images(vfid,ra,objname,group_name,aphys=aphys)

            os.chdir(outdir)

            # for testing
            #if i > 1:
            #    sys.exit()
