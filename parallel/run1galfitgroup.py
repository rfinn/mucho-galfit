#!/usr/bin/env python

'''
GOAL: 
- this program will run galfit on one of JM's group images
- galfit is run once with generic input parameters, no convolution
- the output is used as input for a second run, w/convolution enabled

INPUT:
- group 

PROCEDURE:
- load virgo catalog (vf_v2_main.fits)
- get virgo galaxies in FOV
- create input file for galfit, adding sersic profile at each location
- run galfit w/out convolution

- parse output from galfit, getting details on all sources
- create new input for galfit, with all sources and convolution
- run galfit with convolution


USAGE:
This should be run from a directory that contains
- galfit input file
- image
- noise image

To run:
Run on r-band first!!!

python ~/github/virgowise/run1galfitgroup.py galname r

- galname should be the VFID and should correspond to the subdirectories in 
/mnt/astrophysics/rfinn/muchogalfit-output (grawp directory)

/mnt/astrophysics/muchogalfit-output (how directory is mounted on virgo vms)

- each subdirectory has an input file that is created by setup_galfit.py


This is called my mksbatchgroups.py for running in parallel.

'''
import os
import sys
#import glob
import numpy as np

from astropy.io import fits
from astropy.wcs import WCS
from reproject import reproject_interp
from astropy.table import Table

homedir = os.getenv("HOME")
# add in masking from halphagui
sys.path.append(homedir+'/github/halphagui/')
sys.path.append(homedir+'/github/mucho-galfit/parallel/')
from maskwrapper import buildmask
#from run1galfit import get_maskname
#import run1maskgroup as mg

#import reproject_mask

### DICTIONARIES

pixel_scale = {'FUV':1.5,'NUV':1.5,'g':0.262,'r':0.262,'z':0.262,'W1':2.75,'W2':2.75,'W3':2.75,'W4':2.75}
psf_oversampling = {'FUV':1,'NUV':1,'g':1,'r':1,'g':1,'W1':1,'W2':1,'W3':1,'W4':1}
mag_zeropoint = {'FUV':22.5,'NUV':22.5,'g':22.5,'r':22.5,'g':22.5,'W1':22.5,'W2':22.5,'W3':22.5,'W4':22.5}
image_resolution = {'FUV':6,'NUV':6,'g':1.5,'r':1.5,'z':1.5,'W1':6.1,'W2':6.4,'W3':6.5,'W4':12}

# check this - what is the min mag we want to fit????
minmag2fit = {'FUV':18,'NUV':18,'g':20,'r':20,'z':20,'W1':17,'W2':17,'W3':17,'W4':17}

# this is approx the peak mag of COG_MTOT_ histogram
guess_mag = {'FUV':15.5,'NUV':15.5,'g':14.5,'r':14.5,'z':14.5,'W1':9,'W2':9,'W3':9,'W4':10}


# set up a dictionary for the radius to use for the first guess of the sersic profile
# a better way is to use a constant angular size, and then we can translate that into pixels using the pixel scale




##########################################################################     
### FUNCTIONS
##########################################################################     

def get_galfit_nsersic(galfit_output_image,ngal=1):

    """
    INPUT:
    * galfit_output_image : will have model fits in header 2
    * ngal : number of galaxies modeled in galfit image.  the default is 1

    RETURN:
    nsersic : sersic index for each galaxy that was fit 

    """
    # read in image header
    hdu = fits.open(galfit_output_image)
    # extension 2 has the model info
    imheader = hdu[2].header
    #print(imheader)
    hdu.close()
    
    nsersic = np.zeros(ngal)
    high_nsersic_flag = np.zeros(ngal,'bool')
    for i in range(ngal):
        # build header keywords
        hkey = f"{i+1}_N"


        
        # query galheader
        t = imheader[hkey]
        # parse output : '*16.1052* +/- *1545.5958*'
        n = float(t.replace('*','').split('+/-')[0])
        if n > 6:
            high_nsersic_flag[i] = True
        nsersic[i] = n

    return nsersic, high_nsersic_flag
    

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
    

def parse_galfit_output(galfit_outimage):
    from astropy.io import fits
    #import numpy as np
    numerical_error_flag=0
    header_keywords=['1_XC','1_YC','1_MAG','1_RE','1_N','1_AR','1_PA','2_SKY','CHI2NU']
    fit_parameters=[]
    working_dir=os.getcwd()+'/'
    image_header = fits.getheader(galfit_outimage,2)
    for hkey in header_keywords:
        s=str(image_header[hkey])
        #print hkey
        if s.find('[') > -1:
            s=s.replace('[','')
            s=s.replace(']','')
            t=s.split('+/-')
            values=(float(t[0]),0.)# fit and error
        else:
            t=s.split('+/-')
            try:
                values=(float(t[0]),float(t[1]))# fit and error
            except ValueError:
                # look for * in the string, which indicates numerical problem
                if t[0].find('*') > -1:
                    numerical_error_flag=1
                    t[0]=t[0].replace('*','')
                    t[1]=t[1].replace('*','')
                    values=(float(t[0]),float(t[1]))# fit and error
            except IndexError: # for CHI2NU
                chi2nu=float(t[0])
                continue
        fit_parameters.append(values)
    fit_parameters.append(numerical_error_flag)
    fit_parameters.append(chi2nu)
    return fit_parameters


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
    catalog='/mnt/astrophysics/rfinn/catalogs/Virgo/v2/vf_v2_main.fits'
    try:
        vtab = Table.read(catalog)
    except FileNotFoundError:
        try:
            # test to see if running on Virgo VMS
            catalog='/mnt/astrophysics/catalogs/Virgo/v2/vf_v2_main.fits'
            vtab = Table.read(catalog)
        except FileNotFoundError: # adding case for testing on laptop
            catalog=os.getenv("HOME")+'/research/Virgo/tables-north/v2/vf_v2_main.fits'
            vtab = Table.read(catalog)
            
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
    vfids = vtab['VFID'][flag]
    x,y = x[flag],y[flag]
    # write out file containing VFID, x, y
    ofilename = f'galsFOV-{bandpass}.txt'
    outfile = open(ofilename,'w')
    for i in range(len(vfids)):
        outfile.write(f'{vfids[i]}, {x[i]:.2f}, {y[i]:.2f} \n')
    outfile.close()
    return x,y,vfids




        
def write_galfit_input(output_dir, image,sigma_image,psf_image,bandpass,xgal=None,ygal=None, mask_image=None, firstpass=True,\
                       rBA=None,fixBA=False,rPA=None,fixPA=False,nsersic_flag=None):
    """
    
    PARAMS:
    * BA : BA value to start with
    * fixBA : hold BA fixed
    * PA : PA value to start with
    * fixPA : hold PA fixed; use this if holding r-band values fixed
    * rBA : a list of BA from r-band fits, will have more than one entry for a group image
    * rOA : a list of PA from r-band fits, will have more than one entry for a group image
    * nsersic_flag : boolean to hold n=6 for particular galaxies
    """
 
    if mask_image is not None:
        maskfound = True
    else:
        maskfound = False
        

    t = image.split('-')
    if firstpass:
        if fixBA:
            output_image = f'{t[0]}-{bandpass}-fixBA-out1.fits'
        else:
            output_image = f'{t[0]}-{bandpass}-out1.fits'
    else:
        if fixBA:
            output_image = f'{t[0]}-{bandpass}-fixBA-out2.fits'
        else:
            output_image = f'{t[0]}-{bandpass}-out2.fits'

    # prepend output directory to all images
    sigma_image = output_dir+'/'+sigma_image

    sigma_image = os.path.basename(sigma_image)

        
    psf_sampling = psf_oversampling[bandpass]
    pscale = pixel_scale[bandpass]
    magzp = mag_zeropoint[bandpass]
    minmag = minmag2fit[bandpass]    
    # TODO: add mask to galfit input
    # have updated mask wrapper in halpha gui

    # TODONE: add mask to galfit input
    # have updated mask wrapper in halpha gui
    #maskfound = False
    #mask_image = get_maskname(image)
    #print()
    #print(f"mask name = {mask_image}")
    #print()
    #print(f"cwd = {os.getcwd()}")
    #print()
    #print(f"listdir = {os.listdir()}")
    #print()
    if os.path.exists(mask_image):
        maskfound = True
        print(f"found mask {mask_image}.  Will implement masking in galfit!")
    else:
        print()
        print(f"no mask found for {image} {mask_image}- will NOT implement masking in galfit :(")
        maskfound = False
        print()
    

    # TODONE: need to get xmaxfit,ymaxfit
    xmax,ymax = get_image_size(image) # am I mixing x and y dimensions here?
    print('image dimensions = ',xmax,ymax)
    xminfit = 1
    xmaxfit = xmax 
    yminfit = 1
    ymaxfit = ymax

    
    #if firstpass:
    #    # read in image header and convert RA, DEC to xpixel,ypixel
    
    convolution_size = int(image_resolution[bandpass]/pixel_scale[bandpass]*20)
    if firstpass:
        # remove an straggler galfit files
        # this ensures that the first pass results are in galfit.01
        os.system('rm galfit.??')
        # DONE: need to get (x,y) center of object



        BA=1
        fitBA = 1
        PA=0
        fitPA = 1
        nsersic=2
        fitn = 1

        # set the magnitude based on the typical galaxy parameters
        mag=guess_mag[bandpass]
        fitmag = 1
        sky = 0
        # set initial guess to 25 arcsec,
        # and translate into pixels based on pixelscale of bandpass 
        rad = 30/pixel_scale[bandpass]
        fitrad = 1
        

        
        outfile = open('galfit.input1','w')
        outfile.write('# IMAGE PARAMETERS\n')
        outfile.write('A) '+image+'              # Input data image (FITS file)\n')
        outfile.write('B) '+output_image+'       # Name for the output image\n')
        outfile.write('C) %s                # Sigma image name (made from data if blank or "none") \n'%(sigma_image))
        if maskfound:
            outfile.write(f'F) {mask_image}     # Bad pixel mask (FITS image or ASCII coord list)\n')
        outfile.write('H) '+str(int(round(xminfit)))+' '+str(int(round(xmaxfit)))+' '+str(int(round(yminfit)))+' '+str(int(round(ymaxfit)))+'     # Image region to fit (xmin xmax ymin ymax)\n')
        outfile.write('J) %5.2f              # Magnitude photometric zeropoint \n'%(magzp))
        outfile.write('K) %6.5f   %6.5f         # Plate scale (dx dy)  [arcsec/pix]\n'%(pscale,pscale))
        outfile.write('O) regular                # Display type (regular, curses, both)\n')
        outfile.write('P) 0                   # Create output image only? (1=yes; 0=optimize) \n')
        outfile.write('S) 0                   # Modify/create objects interactively?\n')
        outfile.write(' \n')
        # write object
        # github test
        objnum = 1
        for xi, yi in zip(xgal,ygal):
            outfile.write(f'# Object number: {objnum} \n')
            outfile.write(' 0) sersic             # Object type \n')
            outfile.write(' 1) %8.1f  %8.1f 1 1  # position x, y        [pixel] \n'%(xi,yi))
            outfile.write(' 3) %5.2f      %i       # total magnitude     \n'%(mag,fitmag))
            outfile.write(' 4) %8.2f       %i       #     R_e              [Pixels] \n'%(rad,fitrad))
            outfile.write(' 5) %5.2f       %i       # Sersic exponent (deVauc=4, expdisk=1)   \n'%(nsersic,int(fitn)))
            if rBA is not None:
                if fixBA:
                    print("Fixing BA")
                    outfile.write(' 9) %5.2f       %i       # axis ratio (b/a)    \n'%(rBA[objnum-1],0))
                else:
                    outfile.write(' 9) %5.2f       %i       # axis ratio (b/a)    \n'%(rBA[objnum-1],1))
            else:
                outfile.write(' 9) %5.2f       %i       # axis ratio (b/a)    \n'%(BA,int(fitBA)))
            if rPA is not None:
                if fixPA:
                    outfile.write('10) %5.2f       %i       # position angle (PA)  [Degrees: Up=0, Left=90] \n'%(rPA[objnum-1],0))
                else:
                    outfile.write('10) %5.2f       %i       # position angle (PA)  [Degrees: Up=0, Left=90] \n'%(rPA[objnum-1],1))
            else:
                outfile.write('10) %5.2f       %i       # position angle (PA)  [Degrees: Up=0, Left=90] \n'%(PA,int(fitPA)))
            outfile.write(" Z) 0                  # Output option (0 = residual, 1 = Don't subtract)  \n")
            objnum += 1

        outfile.write(' \n')
        outfile.write(f'# Object number: {objnum} \n')
        outfile.write(' 0) sky             # Object type \n')
        outfile.write(' 1) %8.1f   1  # sky background at center of fitting region [ADUs] \n'%(sky))
        outfile.write(' 2) 0      0       # dsky/dx (sky gradient in x)    \n')
        outfile.write(' 3) 0      0       # dsky/dy (sky gradient in y) \n')
        outfile.write(" Z) 0                  # Output option (0 = residual, 1 = Don't subtract)  \n")
                  
    
        outfile.close()
    else:
        # read in output from first pass
        # and use results as input into second pass
        
        input = open('galfit.01','r')
        all_lines = input.readlines()
        outlines = []
        

        i = 0
        holdfixed = False
        skyobject = True
        objnum=0
        for i in range(len(all_lines)):
            
            line = all_lines[i]
            if '0) sersic' in line:
                objnum += 1
                if nsersic_flag is None:
                    fixn = False
                else:
                    fixn = nsersic_flag[objnum-1]
                
            if line.startswith('B)'):
                outlines.append(line.replace('out1.fits','out2.fits'))
                holdfixed = False
                
            elif line.startswith('D)'):
                
                outlines.append('D) '+psf_image+'     # Input PSF image and (optional) diffusion kernel\n')
            elif line.startswith('E)'):
                outlines.append('E) %i                   # PSF oversampling factor relative to data\n'%(psf_sampling))
            elif line.startswith('I)'):
                outlines.append('I) '+str(int(round(convolution_size)))+' '+str(int(round(convolution_size)))+'             # Size of convolution box (x y)\n')

            elif line.startswith(' 0)'):
                if line.find('sky') > -1:
                    skyobject = True
                else:
                    skyobject = False
                outlines.append(line)
            
            
                
            elif line.startswith(' 3)'):# check if mag is too faint, then don't fit for new parameters
                t = line.split()
                # check the magnitude returned from no covolution
                #print(f"checking magnitude {float(t[1])} compared to {minmag}",holdfixed)
                if float(t[1]) > minmag and not skyobject:
                    # keep the positions fixed
                    #print(i,len(outlines))
                    outlines[i-2] = outlines[i-2].replace(' 1 1',' 0 0')
                                
                    # set all parameters to fixed
                    # fix the x and y coords
                    print("WARNING: object too faint to fit!!!!!!!!")
                                
                                
                                
                    holdfixed=True
                    outlines.append(line.replace(' 1 ',' 0 '))
                else:
                    holdfixed=False
                    outlines.append(line)
                #print(f"checking magnitude {float(t[1])} compared to {minmag}",holdfixed)                
            elif line.startswith(' 4)'):
                print("result from round 1 = ",line)
                if holdfixed and not skyobject:
                    outlines.append(line.replace(' 1 ',' 0 '))
                else:
                    gparams = line.split()
                    rad = float(gparams[1].replace('*',''))
                    print(f"HEY::::: radius from round one = {rad:.3f}")
                    if rad < 2:
                        print("setting radius to min value")
                        outlines.append(' 4) 10      1       #     R_e              [Pixels] \n')
                    else:
                        outlines.append(line)

                # set the min radius to 5 arcsec/pscale

                 
            elif line.startswith(' 5)'):
                if holdfixed and not skyobject:
                    outlines.append(line.replace(' 1 ',' 0 '))
                else:
                    gparams = line.split()
                    n = float(gparams[1].replace('*',''))
                    #print(f"HEY::::: ba from round one = {rad:.3f}")
                    if fixn:
                        outlines.append(' 5) 6      0       #     Sersic index n \n')
                    elif (n < 0.5) | (n > 5):
                        #print("setting radius to min value")
                        outlines.append(' 5) 2      1       #     Sersic index n \n')
                    else:
                        outlines.append(line)

            elif line.startswith(' 9)'):
                if rBA is not None:
                    if fixBA:
                        print("fixing BA!")                        
                        outlines.append(' 9) %5.2f       %i       # axis ratio (b/a)    \n'%(rBA[objnum-1],0))
                    else:
                        outlines.append(' 9) %5.2f       %i       # axis ratio (b/a)    \n'%(rBA[objnum-1],1))
                        
                elif holdfixed and not skyobject:
                    outlines.append(line.replace(' 1 ',' 0 '))
                else:
                    gparams = line.split()
                    ba = float(gparams[1].replace('*',''))
                    #print(f"HEY::::: ba from round one = {rad:.3f}")
                    if ba < 0.1:
                        print("setting radius to min value")
                        outlines.append(' 9) 1      1       #     Axis Ratio (b/a) \n')
                    else:
                        outlines.append(line)



            elif line.startswith('10)'):
                if rPA is not None:
                    if fixPA:
                        print("fixing PA!")
                        outlines.append('10) %5.2f       %i       # position angle (PA)  [Degrees: Up=0, Left=90] \n'%(rPA[objnum-1],0))
                    else:
                        outlines.append('10) %5.2f       %i       # position angle (PA)  [Degrees: Up=0, Left=90] \n'%(rPA[objnum-1],1))
                
                elif holdfixed and not skyobject:
                    outlines.append(line.replace(' 1 ',' 0 '))

                else:
                    outlines.append(line)

            else:
                outlines.append(line)


        outfile = open('galfit.input2','w')
        for l in outlines:
            outfile.write(l)
        outfile.close()
        input.close()
                                                                                                  
if __name__ == '__main__':


    # import argparse

    # parser = argparse.ArgumentParser(
    #     description='Program to create create input file for galfit and then run galfit. ')
    # parser.add_argument('--fix2rparams',
    #                     dest='fix2rparams',
    #                     default=False,
    #                     action='store_true',
    #                     help='set this to run WISE images using r-band PA and BA')


    # args = parser.parse_args()
    
    # run this from /mnt/astrophysics

    # DONE: move galfit output to a new destination
    # /mnt/astrophysics/rfinn/muchogalfit-output


    # DONE: move galfit output to a new destination
    # /mnt/astrophysics/rfinn/muchogalfit-output
    topdir = '/mnt/astrophysics/rfinn/muchogalfit-output/'
    try:
        os.chdir(topdir)
    except FileNotFoundError: # assuming that we are running on virgo vms or draco
        topdir = '/mnt/astrophysics/muchogalfit-output/'
        os.chdir(topdir)
        
    # take as input the galaxy name
    vfid = sys.argv[1]

    # set bandpass to the wavelength that is passed in from the command line
    bandpass = sys.argv[2]
    
    # move to muchogalfit-output directory
    output_dir = topdir+vfid+'/'
    if not os.path.exists(output_dir):
        print('WARNING: {} does not exist\n Be sure to run setup_galfit.py first')
        os.chdir(topdir)
        sys.exit()
    
    os.chdir(output_dir)

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

    matchindex_primary = etab['VFID'] == vfid

    print(f"number of matches to {vfid} = {np.sum(matchindex_primary)}")
    objname = etab['GALAXY'][matchindex_primary][0]
    print(f"objname = {objname}")
    # look in vf tables to find if file is group or not
    if etab['GROUP_MULT'][matchindex_primary] > 1:
        image = f'{objname}_GROUP-custom-image-{bandpass}.fits'
        invvar_image = f'{objname}_GROUP-custom-invvar-{bandpass}.fits'
        std_image = f'{objname}_GROUP-custom-std-{bandpass}.fits'            
        psf_image = f'{objname}_GROUP-custom-psf-{bandpass}.fits'
    else:
        image = f'{objname}-custom-image-{bandpass}.fits'
        invvar_image = f'{objname}-custom-invvar-{bandpass}.fits'
        std_image = f'{objname}-custom-std-{bandpass}.fits'    
        psf_image = f'{objname}-custom-psf-{bandpass}.fits'
            
    print("image = ",image)


    # get galaxies in the FOV
    # def get_galaxies_in_fov(image,bandpass='W3'):    
    x,y,vfids = get_galaxies_in_fov(image, bandpass=bandpass)
    ngal = len(vfids)
    matchflag = np.zeros(len(etab),'bool')
    # TODO - need to fix this for group pointings
    for v in vfids:
        mindex = etab['VFID'] == v
        matchflag[mindex] = True 

    
    
    if np.sum(matchflag) < 1:
        print("ERROR: did not find a matching VFID for ",vfid)

    # this can now contain multiple ids if this is a group image
    matchindex = np.arange(len(etab))[matchflag]


    ########################################################
    ## USE r-band PARAMETERS
    #  if a 3rd argument is passed into program, then
    #  keep BA and PA fixed at r-band values
    ########################################################

    if len(sys.argv) > 3:
        # TODO - need to fix this to properly account for the group galaxies
        # open the galfit output table from rband
        rgalfit = Table.read(topdir+'/vf_v2_galfit_r.fits')
        # check numerical error flag
        print(rgalfit['Numerical_Error'][matchindex_primary])
        if rgalfit['Numerical_Error'][matchindex_primary][0]:
            print("not using r-band params b/c they are not reliable")
            rPA = None
            rBA = None
            fixBA = False
            fixPA = False
        else:
            # get the PA and BA for this galaxy
            rPA = rgalfit['PA'][matchindex]
            rBA = rgalfit['AR'][matchindex]
            print("r-band axis ratio(s): ",rBA)
            fixBA = True
            fixPA = True
            
        if rgalfit['CNumerical_Error'][matchindex_primary][0]:
            print("not using r-band params b/c they are not reliable")
            rCPA = None
            rCBA = None
            fixCBA = False
            fixCPA = False
            
        else:
            # get the PA and BA for this galaxy
            rCPA = rgalfit['CPA'][matchindex]
            rCBA = rgalfit['CAR'][matchindex]            
            fixCBA = True
            fixCPA = True
    else:
        rPA = None
        rBA = None
        fixBA = False
        fixPA = False
        rCPA = None
        rCBA = None
        fixCBA = False
        fixCPA = False
        


        

    # masks are made first, but run1maskgroup.py
    # so we just need the correct mask name
    mask_image = get_maskname(image)
    
    # TODONE: remove previous galfit files if they exist
    os.system('rm galfit.??')
    os.system('rm galfit.input?')

    
    # TODONE: add code to generate galfit input for first run, no convolution, generic starting point
    #write_galfit_input(output_dir, image, sigma_image,psf_image,bandpass,xgal=None,ygal=None, mask_image=None, firstpass=True):
    #rBA=None,fixBA=False,rPA=None,fixPA=False):    
    write_galfit_input(output_dir, image, std_image, psf_image, bandpass, xgal=x, ygal=y, mask_image=mask_image,rPA=rPA,fixPA=fixPA,rBA=rBA,fixBA=fixBA)
    
    # code to run galfit
    print('running galfit')
    os.system(f"galfit galfit.input1")


    
    # TODONE: read galfit output, and create new input to run with convolution
    write_galfit_input(output_dir, image, std_image, psf_image, bandpass, xgal=x, ygal=y, mask_image=mask_image,firstpass=False,rPA=rCPA,fixPA=fixCPA,rBA=rCBA,fixBA=fixCBA)

    # TODO: make sure I am using the correct PSF images
    print('running galfit second time')
    os.system(f"galfit galfit.input2")


    # TODO : add a third iteration if sersic n > 6
    # rerun galfit and hold n fixed at n=6

    # check index from prior run
    # TODO : check output sersic index.  If it's n > 6, then rerun and hold n=6.

    t = image.split('-')
    if fixBA:
        output_image = f'{t[0]}-{bandpass}-fixBA-out2.fits'
    else:
        output_image = f'{t[0]}-{bandpass}-out2.fits'
    
    nsersic, high_nsersic_flag = get_galfit_nsersic(output_image,ngal=ngal)

    if np.sum(high_nsersic_flag) > 0:
        print("running again and holding the sersic index fixed for {np.sum(high_nsersic_flag)}")
        write_galfit_input(output_dir, image, std_image, psf_image, bandpass, xgal=x, ygal=y, mask_image=mask_image,firstpass=False,rPA=rCPA,fixPA=fixCPA,rBA=rCBA,fixBA=fixCBA,nsersic_flag=high_nsersic_flag)

        print('running galfit third time to hold fixed n=6')
        os.system(f"galfit galfit.input2")

    os.chdir(topdir)
