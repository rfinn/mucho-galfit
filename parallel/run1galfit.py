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
python ~/github/virgowise/run1galfit.py galdir bandpass

This assumes that the galfit input file starts with galname and ends with input

Each galaxy goes in its own folder, and the parallel program 
will move to each folder and run this script

'''
import os
import sys
import glob
### DICTIONARIES

pixel_scale = {'FUV':2.5,'NUV':2.5,'g':0.262,'r':0.262,'z':0.262,'W1':2.5,'W2':2.5,'W3':2.5,'W4':2.5}
psf_oversampling = {'FUV':1,'NUV':1,'g':1,'r':1,'g':1,'W1':1,'W2':1,'W3':1,'W4':1}
mag_zeropoint = {'FUV':22.5,'NUV':22.5,'g':22.5,'r':22.5,'g':22.5,'W1':22.5,'W2':22.5,'W3':22.5,'W4':22.5}

# TODO: set up a dictionary for the radius to use for the first guess of the sersic profile

### FUNCTIONS
def funpack_image(input,output):
    from astropy.io import fits
    hdu = fits.open(input)
    print('input file = ',input)
    print(hdu.info())
    try:
        fits.writeto(output,data=hdu[1].data, header=hdu[1].header, overwrite=True)
    except IndexError:
        fits.writeto(output,data=hdu.data, header=hdu.header, overwrite=True)        
    hdu.close()
    
def parse_galfit_1comp(galfit_outimage):
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
    hdu = fits.open(image)
    data_shape = hdu[1].data.shape
    hdu.close()
    return data_shape

def write_galfit_input(galdir, output_dir, bandpass, firstpass=True):
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

    

    # funpack the .fz images so galfit can read them
    # save them in output_dir

    funpack_image(os.path.join(galdir,image),os.path.join(output_dir,image.replace('.fz','')))
    funpack_image(os.path.join(galdir,psf_image),os.path.join(output_dir,psf_image.replace('.fz','')),nhdu=0)
                  
    # prepend output directory to all images

    image = image.replace('.fz','')
    psf_image = psf_image.replace('.fz','')
    
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
        fit_parameters = parse_galfit_output(output_image.replace('out1','out2'))
        # header_keywords=['1_XC','1_YC','1_MAG','1_RE','1_N','1_AR','1_PA','2_SKY','CHI2NU']
        xc, yc, mag, rad, nsersic, BA, PA, sky = fit_parameters[np.arange(0,15,2)]
        fitBA = 1
        fitPA = 1
        fitn = 1
        fitmag = 1
        fitrad = 1
        # get convolution size - set to cutout size?
    # make of values for xminfit, etc for now
    # get image size
    xmax,ymax = get_image_size(image) # am I mixing x and y dimensions here?
    print('image dimensions = ',xmax,ymax)
    xminfit = 1
    xmaxfit = xmax 
    yminfit = 1
    ymaxfit = ymax
    xobj = int(xmaxfit/2)
    yobj = int(ymaxfit/2)
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

    topdir = os.getcwd()
    
    # TODO: take as input the galaxy directory
    # directory of galaxy to work with, this will be 100/NGC2722A for example
    galdir = sys.argv[1]

    # wavelength of images, like W3
    bandpass = sys.argv[2]
    
    # fix this later when we figure how to extract two columns from data directory file
    galname = galdir.split('/')[1] # or pass VFID on second
    
    output_dir = '/mnt/astrophysics/rfinn/muchogalfit-output/'+galname+'/'
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    os.chdir(output_dir)
                      
    
    # set up path name for image directory
    # directory where galaxy images are
    data_dir = '/mnt/astrophysics/virgofilaments-data/'+galdir



    # TODO: add code to generate galfit input for first run, no convolution, generic starting point

    write_galfit_input(data_dir, output_dir, bandpass)
    # code to run galfit
    os.chdir(output_dir)
    print('running galfit')
    os.system(f"galfit galfit.input1")


    # TODO: read galfit output, and create new input to run with convolution


    os.chdir(topdir)
