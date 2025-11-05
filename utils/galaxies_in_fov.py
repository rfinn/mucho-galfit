import sys
import os
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.io import fits
import glob


def get_header(path_to_image):
    
    if not os.path.exists(path_to_image):
        print(f'File {path_to_image} does not exist! Unable to create galsFOV.txt')
        sys.exit()
    
    hdu = fits.open(path_to_image)
    header = hdu[0].header
    hdu.close()
    
    return header


def get_image_size(path_to_image):

    if not os.path.exists(path_to_image):
        print(f'File {path_to_image} does not exist! Unable to create galsFOV.txt')
        sys.exit()
    
    data = fits.getdata(path_to_image)
    data_shape = data.shape
    return data_shape


def get_galaxies_in_fov(gtab, output_dir):
        """
        get catalog galaxies in FOV. needed for building galaxy webpages! 
        gtab --> main galaxy table
        output_dir --> path to primary galaxy directory where galsFOV.txt will be saved
        """
        
        # create a SkyCoord object from RA and DEC of galaxies
        galcoord = SkyCoord(gtab['RA'],gtab['DEC'],frame='icrs',unit='deg')
        print(output_dir)
        #grab some filename in output_dir with 'image' and 'fits' in the name
        im_path = glob.glob(f'{output_dir}*-image-W3.fits')[0]
        
        #grab the image header information
        im_header = get_header(im_path)
        
        #set up image wcs
        image_wcs = WCS(im_header)

        # get the size of the image
        xmax, ymax = get_image_size(im_path)
    
        # find galaxies on cutout
        x,y = image_wcs.world_to_pixel(galcoord)

        # create flag to save galaxies on the image
        flag = (x > 0) & (x < xmax) & (y>0) & (y < ymax)        
        galids = gtab['OBJID'][flag]
        x,y = x[flag],y[flag]
        # write out file containing OBJID, x, y
        ofilename = f'galsFOV.txt'
        outfile = open(ofilename,'w')
        for i in range(len(galids)):
            outfile.write(f'{galids[i]}, {x[i]:.2f}, {y[i]:.2f} \n')
        outfile.close()