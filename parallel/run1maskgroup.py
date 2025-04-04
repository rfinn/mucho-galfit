#!/usr/bin/env python

'''
GOAL: 
- this creates a mask for a group image
- the mask is built from the r-band image
- the r-band mask is then reprojected for use with the unWISE images

INPUT:
- directory of the primary group number comes in as sys.argv[1]

PROCEDURE:
- load galaxy catalog, which must have RA and DEC columns 
  - vf_v2_main.fits for Virgo
  - wisesize_main.fits for WISEsize
- get galaxies in FOV
- create a mask
- remove mask for galaxies in FOV
  - for virgo, we use JM's ephot catalog to remove any masked pixels within galaxy's
    ellipse, but we keep the masks for gaia stars.
  - for WISEsize, as of april 04, 2025, we don't have photometry from JM.  
    So we can just remove the SExtractor sources that are at each galaxy position.  
    This won't be perfect, but it will allow us to move forward until SGA2025 is 
    released.

USAGE:
This should be run from a directory that contains
- image
- noise image

To run:
Run on r-band first!!!

python ~/github/virgowise/run1galfitgroup.py galname 

- galname (should be the VFID for Virgo and) should correspond to the 
  subdirectories in 
  - Virgo: /mnt/astrophysics/rfinn/muchogalfit-output (grawp directory)
    - each subdirectory has an input file that is created by setup_galfit.py
  - WISEsize: 
    - /mnt/astrophysics/wisesize/mg_output_wisesize (draco directory)
    - /mnt/astrophysics/rfinn/wisesize/mg_output_wisesize (grawp directory)
    - each directory has grz, W1-W4 images; all invvar images; all std images; 
      wise psfs
    - Kim has script setup_galfit_ws.py that creates these directories.  
      in wisesize branch of mucho-galfit

'''
import os
import sys
#import glob
import numpy as np

from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from reproject import reproject_interp
from scipy.stats import scoreatpercentile


homedir = os.getenv("HOME")
# add in masking from halphagui
sys.path.append(homedir+'/github/halphagui/')
from maskwrapper import buildmask
import imutils
from mask1galaxy import get_galaxy_params
#import reproject_mask

### DICTIONARIES

pixel_scale = {'FUV':1.5,'NUV':1.5,'g':0.262,'r':0.262,'z':0.262,'W1':2.75,'W2':2.75,'W3':2.75,'W4':2.75}
psf_oversampling = {'FUV':1,'NUV':1,'g':1,'r':1,'g':1,'W1':1,'W2':1,'W3':1,'W4':1}
mag_zeropoint = {'FUV':22.5,'NUV':22.5,'g':22.5,'r':22.5,'g':22.5,'W1':22.5,'W2':22.5,'W3':22.5,'W4':22.5}
image_resolution = {'FUV':6,'NUV':6,'g':1.5,'r':1.5,'z':1.5,'W1':6.1,'W2':6.4,'W3':6.5,'W4':12}
minmag2fit = {'FUV':10,'NUV':10,'g':17,'r':17,'z':17,'W1':10,'W2':10,'W3':10,'W4':10}

# this is approx the peak mag of COG_MTOT_ histogram
guess_mag = {'FUV':17.5,'NUV':17.5,'g':15.5,'r':15.5,'z':15.5,'W1':15,'W2':15,'W3':14.5,'W4':14.5}


# set up a dictionary for the radius to use for the first guess of the sersic profile
# a better way is to use a constant angular size, and then we can translate that into pixels using the pixel scale


###########################################
## SET UP CATALOGS
###########################################

# setting up to look for Kim's param file
# get current directory
cdir = os.getcwd()

# first option is for running with wiseSize
if 'mg_output_wisesize' in cdir:
    # open paramfile.txt
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
            
    MAIN_CATALOG = main_catalog
    ELLIPSE_CATALOG = None
    OBJID_COLNAME = 'OBJID'

else:
    catalogdir='/mnt/astrophysics/catalogs/Virgo/v2/'
    # read in virgo catalog
    etablename = 'vf_v2_legacy_ephot.fits'
    maintablename = 'vf_v2_main.fits'

    MAIN_CATALOG = catalogdir+maintablename
    ELLIPSE_CATALOG = catalogdir+etablename
    OBJID_COLNAME = 'VFID'

# directory structure on grawp



class buildgroupmask(buildmask):
    def __init__(self,image):
        self.image_name = image
        self.image, self.imheader = fits.getdata(self.image_name,header = True)
        self.ymax,self.xmax = self.image.shape        
        self.xc = self.xmax/2.
        self.yc = self.ymax/2.
        self.image_wcs = WCS(self.imheader)
        
        # another parameter to initialize
        self.off_center_flag = False
        self.remove_center_object_flag = False

        # to track usr masks - need this for writing functions
        self.deleted_objects = []
        self.usr_mask = np.zeros_like(self.image)

        
        # SE parameters for masking
        self.threshold=0.005
        self.snr=10
        self.snr_analysis=10        
        self.sepath = homedir+'/github/halphagui/astromatic/'
        self.config='default.sex.HDI.mask'

        # set auto to true to run masking in automatic mode
        self.auto = True
        
        # GAIA catalog, with path set for grawp
        # TODO - need to implement gaia query to get the faint stars as well
        self.gaiapath = '/mnt/astrophysics/rfinn/catalogs/gaia-mask-dr9.virgo.fits'
        self.gaia_mask = None
        self.add_gaia_stars = True        

        t = self.image_name.split('.fit')
        self.mask_image=t[0]+'-mask.fits'
        self.mask_inv_image=t[0]+'-inv-mask.fits'
        print('saving mask image as: ',self.mask_image)

        # default se parameters
        self.objra = None
        self.objdec = None
        self.objsma = None
        self.objBA = None
        self.objPA = None
        self.ngrow = 4
        self.sesnr = 10
        self.minarea = 5
        

        # get the pixel scale in the image
        self.pscalex,self.pscaley = self.image_wcs.proj_plane_pixel_scales() # appears to be degrees/pixel
        
        # get image dimensions in deg,deg
        self.dxdeg,self.dydeg = imutils.get_image_size_deg(self.image_name)
        

        # Get coord of image center.  will use when getting gaia stars
        self.racenter,self.deccenter = imutils.get_image_center_deg(self.image_name)                



        # TODO - add ellipse param for each galaxy in FOV
        self.ellipseparams = None

        # other params that need to be initialized        
        self.v1,self.v2=scoreatpercentile(self.image,[.5,99.5])
        self.adjust_mask = True
        self.figure_size = (10,5)
        self.cmap = 'gist_heat_r'
    def build_mask(self):

        ##################################
        ### METHODS FROM MASKWRAPPER.PY
        ##################################
        # SET UP AND RUN SOURCE EXTRACTOR
        # link source extractor files
        self.link_files() # maskwrapper method
        self.runse() # maskwrapper method
        self.get_gaia_stars() # maskwrapper method
        self.add_gaia_masks() # maskwrapper method
        self.grow_mask() # maskwrapper method
        self.grow_mask() # maskwrapper method
        self.grow_mask() # maskwrapper method
        self.grow_mask() # maskwrapper method

        ##################################
        ### METHODS FROM THIS CLASS
        ##################################
        self.get_galaxies_in_fov() # this class, updated to work with any gal catalog with RA and DEC columns

        # the ellipse params are used to clear and SE detection on galaxy
        # while leaving gaia stars
        if ELLIPSE_CATALOG is not None:
            self.get_ellipse_params() # this class, using JM ephot table to get ellipse params for each galaxy

        ##################################
        ### METHODS FROM MASKWRAPPER.PY
        ##################################
        self.remove_center_object() # remove center object flag is false in the init function
        #m.remove_gals(xgals,ygals)
        #self.write_mask()        
        self.show_mask_mpl()
        
    def get_galaxies_in_fov(self):
        """get catalog galaxies in FOV.  updated to work with Virgo or WISEsize """
        from astropy.wcs import WCS
        from astropy.coordinates import SkyCoord
        from astropy.table import Table

        # read in galaxy catalog, where MAIN_CATALOG is defined at the start of this program
        catalog=MAIN_CATALOG
        if os.path.exists(catalog):
            gtab = Table.read(catalog)
        else:
            #try:
            #    # test to see if running on Virgo VMS
            #    catalog='/mnt/astrophysics/catalogs/Virgo/v2/vf_v2_main.fits'
            #    vtab = Table.read(catalog)
            #except FileNotFoundError: # adding case for testing on laptop
            #    catalog=os.getenv("HOME")+'/research/Virgo/tables-north/v2/vf_v2_main.fits'
            #    vtab = Table.read(catalog)
            print(f"WARNING: could not find {MAIN_CATALOG}")
            print("Terminating program...")
            sys.exit()

        # create a SkyCoord object from RA and DEC of virgo galaxies
        galcoord = SkyCoord(vtab['RA'],vtab['DEC'],frame='icrs',unit='deg')

        # set up image wcs

        # get the size of the image
        xmax, ymax = get_image_size(self.image_name)
    
        # find galaxies on cutout
        x,y = self.image_wcs.world_to_pixel(galcoord)

        # create flag to save galaxies on the image
        flag = (x > 0) & (x < xmax) & (y>0) & (y < ymax)        
        galids = gtab[OBJ][flag]
        x,y = x[flag],y[flag]
        # write out file containing VFID, x, y
        ofilename = f'galsFOV.txt'
        outfile = open(ofilename,'w')
        for i in range(len(galids)):
            outfile.write(f'{vfids[i]}, {x[i]:.2f}, {y[i]:.2f} \n')
        outfile.close()
        self.keepflag = flag
        self.galids = galids
        self.xpixel = x
        self.ypixel = y

        
    
    def get_ellipse_params(self):
        """

        call after get_galaxies_in_fov()

        this create lists of ellipse parameters to use with maskwrapper

        """
        # get ellipse params as well
    
        #self.objparams = [self.defcat.cat['RA'][self.igal],self.defcat.cat['DEC'][self.igal],mask_scalefactor*self.radius_arcsec[self.igal],self.BA[self.igal],self.PA[self.igal]+90]
        gRA = []
        gDEC = []
        gRAD = []
        gBA = []
        gPA = []
        for vf in self.galids:
            t = get_galaxy_params(vf)
            gRA.append(t[0])
            gDEC.append(t[1])
            gRAD.append(t[2])
            gBA.append(t[3])
            gPA.append(float(t[4])+90)
        # convert radius to pixels

        self.objBA = gBA
        self.objPA = gPA
        self.objsma = gRAD
        self.objsma_pixels = self.objsma/(self.pscalex.value*3600)
        
                

    def remove_gals(self,xgals,ygals):
        '''
        method to remove SE segmentation id that corresponds to each galaxy position

        '''

        ## TODO - should use ellipse for each galaxy like I do in the regular masking routine
        
        for x,y in zip(xgals,ygals):
            # get mask value at location of galaxy
            maskval = self.maskdat[int(y),int(x)]
            
            print(f"galaxy at {x:.1f},{y:.1f}, value={maskval:.0f}")
            # get pixels where the mask value = value at location of galaxy

            # and 
            # set those pixels to zero
            # adding commit just because
            self.maskdat[self.maskdat == maskval]=0
        self.write_mask()
    def remove_group_gals(self):
        """ copied this from maskwrapper """ 

        # need to replace this with a function that will remove any objects within the specificed central ellipse

        if self.off_center_flag:
            print('setting center object to objid ',self.galaxy_id)
            self.center_object = self.galaxy_id
        else:
            self.center_object = self.read_se_cat()
        if self.center_object is not np.nan:
            self.maskdat[self.maskdat == self.center_object] = 0
        if self.objsma is not None:
            # remove central objects within elliptical aperture
            print("ellipse params in remove_central_object :",self.xpixel,self.ypixel,self.objsma_pixels,self.objBA,self.objPA)
            self.maskdat,self.ellipseparams = remove_central_objects(self.maskdat, sma=self.objsma_pixels, BA=self.objBA, PA=self.objPA, xc=self.xpixel,yc=self.ypixel)
        else:
            print("no ellipse params")
            self.ellipseparams = None
        self.update_mask()
        
### FUNCTIONS




def reproject_mask(maskfile, reffile):
    hmask = fits.open(maskfile)

    href = fits.open(reffile)

    # reproject r-band mask onto W3 header

    wisemask,footprint = reproject_interp(hmask,href[0].header)

    # all wise images have the same pixel scale, so we only need one wise mask
    outname = maskfile.replace('r-mask','wise-mask')
    fits.writeto(outname,wisemask,href[0].header,overwrite=True)
    return outname



def get_image_size(image):
    from astropy.io import fits
    data = fits.getdata(image)
    data_shape = data.shape
    return data_shape

def get_xy_from_wcs(ra,dec,image):
    from astropy.io import fits
    from astropy.wcs import WCS
    from astropy.coordinates import SkyCoord
    header = fits.getheader(image)
    w = WCS(header)
    c = SkyCoord(ra,dec,unit='deg',frame='icrs')
    xobj,yobj = w.world_to_pixel(c)
    return xobj,yobj

def get_galaxies_in_fov(image,bandpass='W3'):
    """get galaxies in FOV using virgo catalog """
    from astropy.wcs import WCS
    from astropy.coordinates import SkyCoord
    from astropy.table import Table

    # read in virgo catalog
    catalog=MAIN_CATALOG
    if os.path.exists(catalog):
        vtab = Table.read(catalog)
    else:
        #try:
        #    # test to see if running on Virgo VMS
        #    catalog='/mnt/astrophysics/catalogs/Virgo/v2/vf_v2_main.fits'
        #    vtab = Table.read(catalog)
        #except FileNotFoundError: # adding case for testing on laptop
        #    catalog=os.getenv("HOME")+'/research/Virgo/tables-north/v2/vf_v2_main.fits'
        #    vtab = Table.read(catalog)
        print(f"WARNING: could not find {MAIN_CATALOG}")
        print("Terminating program...")
        sys.exit()
        
    # create a SkyCoord object from RA and DEC of virgo galaxies
    galcoord = SkyCoord(vtab['RA'],vtab['DEC'],frame='icrs',unit='deg')

    # set up image wcs
    image_wcs = WCS(image)

    # get the size of the image
    xmax, ymax = get_image_size(image)
    
    # find galaxies on cutout
    x,y = image_wcs.world_to_pixel(galcoord)

    # create flag to save galaxies on the image
    flag = (x > 0) & (x < xmax) & (y>0) & (y < ymax)        
    vfids = vtab[OBJID_COLNAME][flag]
    x,y = x[flag],y[flag]
    # write out file containing VFID, x, y
    ofilename = f'galsFOV-{bandpass}.txt'
    outfile = open(ofilename,'w')
    for i in range(len(vfids)):
        outfile.write(f'{vfids[i]}, {x[i]:.2f}, {y[i]:.2f} \n')
    outfile.close()


    # get ellipse params as well

    #self.objparams = [self.defcat.cat['RA'][self.igal],self.defcat.cat['DEC'][self.igal],mask_scalefactor*self.radius_arcsec[self.igal],self.BA[self.igal],self.PA[self.igal]+90]    
    return x,y

def get_rband_mask(image):
    mask_image = get_maskname(image)
    
    #xgals, ygals = get_galaxies_in_fov(image,bandpass='r')
    
    m = buildgroupmask(image)
    m.build_mask()
    os.rename(m.mask_image,mask_image)



def get_group_mask(image,bandpass=None,overwrite=True):
    """
    create the mask for the group image, remove galaxies to be fitted  

    the logic in this function assumes that the r-band mask is made first!!!

    """
    # create the mask
    # remove objects at the positions of the galaxies to be fitted

    # this will return the r-band mask for g,r,z and the 'wise' mask for W1-W4
    mask_image = get_maskname(image)
    
    # if bandpass='r', make the mask
    if (bandpass == 'r'):
        if os.path.exists(mask_image) and not overwrite:
            print("found r-band mask  - no remaking it")
        else:
            print("did not find r-band mask - remaking it!")
            mask_image = get_rband_mask(image)

    else: # reproject r-band mask onto
        
        # look for r-band mask
        if bandpass in ['W1','W2','W3','W4']:
            rmask = mask_image.replace('wise','r')
        print("looking for an r-band masked called ",rmask)
        rimage = image.replace(bandpass,'r')        
        # if r-band mask doesn't exist, build it
        if not os.path.exists(rmask):
            print("did not find r-band mask - remaking it!")
            rmask = get_rband_mask(rimage)
        else:
            print("\tfound it!")
        if bandpass in ['W1','W2','W3','W4']:
            mask_image = reproject_mask(rmask,image)
    

    # return mask_image
    return mask_image

def get_maskname(image):
    """ return the wise mask for wise images and the r mask for legacy images  """
    # check for wise extenstion
    elist = ['W1','W2','W3','W4']
    
    for e in elist:
        if e in image:
            maskname = image.replace(e+".fits","wise-mask.fits")
            return maskname
    # check for legacy
    llist = ['image-g','image-r','image-z']
    for l in llist:
        if l in image:
            maskname = image.replace(l+".fits","image-r-mask.fits")
            return maskname
        

                                                                                                  
if __name__ == '__main__':

    # TODO : if ellip phot file is NONE, remove the SE objects that are at each galaxy position
    
    # run this from /mnt/astrophysics

    # DONE: move galfit output to a new destination
    # /mnt/astrophysics/rfinn/muchogalfit-output

    # TODO : change to run from current directory which should be path_to_images, don't assume or change to a directory

    # TODO : remove any other path dependencies to virgo and fix with statements at the beginning
    topdir = '/mnt/astrophysics/rfinn/muchogalfit-output/'
    try:
        os.chdir(topdir)
    except FileNotFoundError: # assuming that we are running on virgo vms or draco
        topdir = '/mnt/astrophysics/muchogalfit-output/'
        os.chdir(topdir)
    # take as input the galaxy name
    vfid = sys.argv[1]
    bandpass = 'r'
    # move to muchogalfit-output directory
    output_dir = topdir+vfid+'/'
    if not os.path.exists(output_dir):
        print('WARNING: {} does not exist\n Be sure to run setup_galfit.py first')
        os.chdir(topdir)
        sys.exit()
    
    os.chdir(output_dir)

    try:
        
        etab = Table.read(catalogdir+etablename)
        maintab = Table.read(catalogdir+maintablename)        
        
        aphys = '/mnt/astrophysics/rfinn/'
    except FileNotFoundError:
        #try:
        #    # test to see if running on Virgo VMS or draco
        #    catalogdir='/mnt/astrophysics/catalogs/Virgo/v2/'
        #    etab = Table.read(catalogdir+tablename)
        #    maintab = Table.read(catalogdir+maintablename)        
        #    aphys = '/mnt/astrophysics/'        
        #except FileNotFoundError: # 
        #    print("ERROR: problem locating astrophysics drive - exiting")
        #    sys.exit()
        print(f"WARNING: could not find {MAIN_CATALOG}")
        print("Exiting program...")
        sys.exit

    matchflag = etab[OBJID_COLNAME] == vfid
    
    if np.sum(matchflag) < 1:
        print("ERROR: did not find a matching VFID for ",vfid)
    matchindex = np.arange(len(etab))[matchflag][0]
    
        
    for bandpass in ['r','W1']:
        objname = etab['GALAXY'][matchindex]
        # look in vf tables to find if file is group or not
        if etab['GROUP_MULT'][matchindex] > 1:
            image = f'{objname}_GROUP-custom-image-{bandpass}.fits'
            invvar_image = f'{objname}_GROUP-custom-invvar-{bandpass}.fits'    
            psf_image = f'{objname}_GROUP-custom-psf-{bandpass}.fits'
        else:
            image = f'{objname}-custom-image-{bandpass}.fits'
            invvar_image = f'{objname}-custom-invvar-{bandpass}.fits'    
            psf_image = f'{objname}-custom-psf-{bandpass}.fits'
            
        print("image = ",image)
            
        # get mask
        mask_image = get_group_mask(image,bandpass=bandpass,overwrite=True)

    os.chdir(topdir)
