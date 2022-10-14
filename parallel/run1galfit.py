#!/usr/bin/env python

'''
GOAL: 
- this program will run galfit on one galaxy
- galfit is run once with generic input parameters, no convolution
- the output is used as input for a second run, w/convolution enabled

USAGE:
This should be run from a directory that contains
- galfit input file
- image
- noise image

To run:
python ~/github/virgowise/run1galfit.py galname

- galname should be the VFID and should correspond to the subdirectories in 
/mnt/astrophysics/rfinn/muchogalfit-output (grawp directory)

/mnt/astrophysics/muchogalfit-output (how directory is mounted on virgo vms)

- each subdirectory has an input file that is created by setup_galfit.py

'''
import os
import sys
#import glob
import numpy as np
### DICTIONARIES

pixel_scale = {'FUV':2.5,'NUV':2.5,'g':0.262,'r':0.262,'z':0.262,'W1':2.5,'W2':2.5,'W3':2.5,'W4':2.5}
psf_oversampling = {'FUV':1,'NUV':1,'g':1,'r':1,'g':1,'W1':1,'W2':1,'W3':1,'W4':1}
mag_zeropoint = {'FUV':22.5,'NUV':22.5,'g':22.5,'r':22.5,'g':22.5,'W1':22.5,'W2':22.5,'W3':22.5,'W4':22.5}

# TODO: set up a dictionary for the radius to use for the first guess of the sersic profile

### FUNCTIONS
def funpack_image(input,output,nhdu=1):
    from astropy.io import fits
    hdu = fits.open(input)
    print('input file = ',input)
    print(hdu.info())
    fits.writeto(output,data=hdu[nhdu].data, header=hdu[nhdu].header, overwrite=True)
    hdu.close()
    
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

def convert_invvar_noise(invvar_image, noise_image):
    # convert invvar image to noise
    from astropy.io import fits
    import numpy as np
    # read in invvar image
    #print('invvar image = ',invvar_image)
    hdu = fits.open(invvar_image)
    data = hdu[1].data
    header = hdu[1].header
    hdu.close()
    # operate on pixels to take sqrt(1/x)
    noise_data = np.sqrt(1/data)
    
    # check for bad values, nan, inf
    # set bad pixels to very high value, like 1e6
    noise_data = np.nan_to_num(noise_data,nan=1.e6)

    # write out as noise image
    fits.writeto(noise_image,noise_data,header=header,overwrite=True)

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

def write_galfit_input(galdir, output_dir, ra, dec, bandpass, firstpass=True):
    galname = os.path.basename(galdir)
    image = f'{galname}-custom-image-{bandpass}.fits.fz'
    invvar_image = f'{galname}-custom-invvar-{bandpass}.fits.fz'    
    psf_image = f'{galname}-custom-psf-{bandpass}.fits.fz'

    # created images
    sigma_image = f'{galname}-custom-std-{bandpass}.fits'    
    if firstpass:
        output_image = f'{galname}-{bandpass}-out1.fits'
    else:
        output_image = f'{galname}-{bandpass}-out2.fits'    

    
    if firstpass:

        
        # funpack the .fz images so galfit can read them
        # save them in output_dir
        funpack_image(os.path.join(galdir,image),funpacked_image)
        funpack_image(os.path.join(galdir,psf_image),os.path.join(output_dir,psf_image.replace('.fz','')),nhdu=0)

    # get information from the image, like the image size and position of gal
    image = image.replace('.fz','')
    psf_image = psf_image.replace('.fz','')
                  
    # prepend output directory to all images

    
    # check if noise image exists, if not make it from invvar
    sigma_image = output_dir+'/'+sigma_image
    invvar_image = galdir+'/'+invvar_image    
    if not os.path.exists(sigma_image):
        convert_invvar_noise(invvar_image,sigma_image)
    sigma_image = os.path.basename(sigma_image)


        
    psf_sampling = psf_oversampling[bandpass]
    pscale = pixel_scale[bandpass]
    magzp = mag_zeropoint[bandpass]
    
    # TODO: add mask to galfit input

    # TODO: need to get xmaxfit,ymaxfit

    # TODO: need to get (x,y) center of object

    # make of values for xminfit, etc for now
    # get image size
    
    xmax,ymax = get_image_size(image) # am I mixing x and y dimensions here?
    print('image dimensions = ',xmax,ymax)
    xminfit = 1
    xmaxfit = xmax 
    yminfit = 1
    ymaxfit = ymax
    
    if firstpass:
        # read in image header and convert RA, DEC to xpixel,ypixel
        xobj, yobj = get_xy_from_wcs(ra,dec,image)

    if firstpass:
        BA=1
        fitBA = 1
        PA=0
        fitPA = 1
        nsersic=2
        fitn = 1
        mag=14
        fitmag = 1
        rad = 5
        fitrad = 1
        sky = 0
    else:
        # read in output from first pass run of galfit
        t = parse_galfit_output(output_image.replace('out2','out1'))
        print(t)
        # header_keywords=['1_XC','1_YC','1_MAG','1_RE','1_N','1_AR','1_PA','2_SKY','CHI2NU']
        xc, yc, mag, rad, nsersic, BA, PA, sky = t[0][0],t[1][0],t[2][0],t[3][0],t[4][0],t[5][0],t[6][0],t[7][0]
        fitBA = 1
        fitPA = 1
        fitn = 1
        fitmag = 1
        fitrad = 1
        # get convolution size - set to cutout size?
    # TODO: need to decide on the right size for this
    # Chien had recommended the full image, but we could do something like 8-10x pixel size
    # a smaller convolution size should make galfit run faster
    convolution_size = xobj
    if firstpass:
        outfile = open('galfit.input1','w')
    else:
        outfile = open('galfit.input2','w')        
    outfile.write('# IMAGE PARAMETERS\n')
    outfile.write('A) '+image+'              # Input data image (FITS file)\n')
    outfile.write('B) '+output_image+'       # Name for the output image\n')
    outfile.write('C) %s                # Sigma image name (made from data if blank or "none") \n'%(sigma_image))
    if not firstpass:
        outfile.write('D) '+psf_image+'     # Input PSF image and (optional) diffusion kernel\n')
        outfile.write('E) %i                   # PSF oversampling factor relative to data\n'%(psf_sampling))
    outfile.write('H) '+str(int(round(xminfit)))+' '+str(int(round(xmaxfit)))+' '+str(int(round(yminfit)))+' '+str(int(round(ymaxfit)))+'     # Image region to fit (xmin xmax ymin ymax)\n')
    if not firstpass:
        outfile.write('I) '+str(int(round(convolution_size)))+' '+str(int(round(convolution_size)))+'             # Size of convolution box (x y)\n')
    outfile.write('J) %5.2f              # Magnitude photometric zeropoint \n'%(magzp))
    outfile.write('K) %6.5f   %6.5f         # Plate scale (dx dy)  [arcsec/pix]\n'%(pscale,pscale))
    outfile.write('O) regular                # Display type (regular, curses, both)\n')
    outfile.write('P) 0                   # Create output image only? (1=yes; 0=optimize) \n')
    outfile.write('S) 0                   # Modify/create objects interactively?\n')
    outfile.write(' \n')
    # write object
    outfile.write('# Object number: 1 \n')
    outfile.write(' 0) sersic             # Object type \n')
    outfile.write(' 1) %8.1f  %8.1f 1 1  # position x, y        [pixel] \n'%(xobj,yobj))
    outfile.write(' 3) %5.2f      %i       # total magnitude     \n'%(mag,fitmag))
    outfile.write(' 4) %8.2f       %i       #     R_e              [Pixels] \n'%(rad,fitrad))
    outfile.write(' 5) %5.2f       %i       # Sersic exponent (deVauc=4, expdisk=1)   \n'%(nsersic,int(fitn)))
    outfile.write(' 9) %5.2f       %i       # axis ratio (b/a)    \n'%(BA,int(fitBA)))
    outfile.write('10) %5.2f       %i       # position angle (PA)  [Degrees: Up=0, Left=90] \n'%(PA,int(fitPA)))
    outfile.write(" Z) 0                  # Output option (0 = residual, 1 = Don't subtract)  \n")
    

    outfile.write(' \n')
    outfile.write('# Object number: 2 \n')
    outfile.write(' 0) sky             # Object type \n')
    outfile.write(' 1) %8.1f   1  # sky background at center of fitting region [ADUs] \n'%(sky))
    outfile.write(' 2) 0      0       # dsky/dx (sky gradient in x)    \n')
    outfile.write(' 3) 0      0       # dsky/dy (sky gradient in y) \n')
    outfile.write(" Z) 0                  # Output option (0 = residual, 1 = Don't subtract)  \n")
                  
    
    outfile.close()


if __name__ == '__main__':
    
    # run this from /mnt/astrophysics

    # TODO: move galfit output to a new destination
    # /mnt/astrophysics/rfinn/muchogalfit-output
    topdir = '/mnt/astrophysics/rfinn/muchogalfit-output/'
    os.chdir()
    
    # take as input the galaxy name
    galname = sys.argv[1]

    # move to muchogalfit-output directory
    output_dir = topdir+galname+'/'
    if not os.path.exists(output_dir):
        print('WARNING: {} does not exist\n Be sure to run setup_galfit.py first')
        os.chdir(topdir)
        sys.exit()
    
    os.chdir(output_dir)

    # read in input file to get galname, objname, ra, dec, and bandpass
    sourcefile = open(galname+'sourcelist','r')
    galaxies = sourcefile.readlines()
    if len(galaxies) > 1:
        # set the flag to have more than one galaxy in the galfit input file
        multiflag = True
    elif len(galaxies) == 1:
        multiflag = False
    else:
        print('Problem reading sourcelist for {}'.format(galname))
        print('Please check the setup directory')
        os.chdir(topdir)
        sys.exit()

    # parse information from file
    vfid, objname, ra, dec, bandpass = galaxies[0].rstrip().split()
    ra = float(ra)
    dec = float(dec)
    

    # set up path name for image directory
    # directory where galaxy images are
    data_dir = '/mnt/astrophysics/virgofilaments-data/{}/{}/'.format(int(ra),objname)



    # TODO: add code to generate galfit input for first run, no convolution, generic starting point
    write_galfit_input(data_dir, output_dir, ra, dec, bandpass)
    
    # code to run galfit
    print('running galfit')
    os.system(f"galfit galfit.input1")

    # TODO: read galfit output, and create new input to run with convolution
    write_galfit_input(data_dir, output_dir, bandpass, firstpass=False)

    print('running galfit second time')
    os.system(f"galfit galfit.input2")
    
    os.chdir(topdir)
