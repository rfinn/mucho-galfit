#!/usr/bin/env python
'''
GOAL:
* create web page to inspect the galfit results
* adapting code from build_web_cutouts2.py
* create one page per galaxy, and show results for
  * g, r, z
  * W1, w/conv, fixedBA
  * W2
  * W3, w/conv, fixedBA
  * W4

USAGE:
* run from cutouts directory

NOTES:
* using John Moustakas's code as a reference (https://github.com/moustakas/legacyhalos/blob/main/py/legacyhalos/virgofilaments.py#L1131-L1202)
* https://docs.astropy.org/en/stable/visualization/normalization.html#:~:text=The%20astropy.visualization%20module%20provides%20a%20framework%20for%20transforming,x%20represents%20the%20values%20in%20the%20original%20image%3A

* on draco, the results are in /mnt/astrophysics/muchogalfit-output/



'''

import os
import sys
import numpy as np
import glob

from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse

from scipy.stats import scoreatpercentile
from astropy.io import fits
from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy.visualization import simple_norm
from astropy.visualization import SqrtStretch, PercentileInterval
from astropy.visualization import ImageNormalize
from astropy.visualization import LinearStretch,SinhStretch
from astropy import units as u
from astropy.nddata import Cutout2D
from astropy.stats import sigma_clip

import multiprocessing as mp

from PIL import Image

homedir = os.getenv("HOME")

os.sys.path.append(homedir+'/github/virgowise/')
import rungalfit as rg #This code has galfit defined functions 

#from build_web_coadds import get_galaxies_fov, plot_vf_gals
from build_web_common import *


###########################################################
####  GLOBAL VARIABLES
###########################################################


VFMAIN_PATH = homedir+'/research/Virgo/tables-north/v1/vf_north_v1_main.fits'
VFMAIN_PATH = homedir+'/research/Virgo/tables-north/v2/vf_v2_main.fits'
VFEPHOT_PATH = homedir+'/research/Virgo/tables-north/v2/vf_v2_legacy_ephot.fits'
haimaging_path = os.path.join(homedir,'github/HalphaImaging/python3/')
#sys.path.append(haimaging_path)

#vfmain = fits.getdata(VFMAIN_PATH)
residual_stretch = LinearStretch(slope=0.5, intercept=0.5) + SinhStretch() + \
    LinearStretch(slope=2, intercept=-1)
###########################################################
####  FUNCTIONS
###########################################################
def get_params_from_name(image_name):
    #print(t)
    tels = ['BOK','HDI','INT','MOS']
    for t in tels:
        if t in image_name:
            telescope = t
            break
    t = os.path.basename(image_name).split('-')
    for item in t:
        if item.startswith('20'):
            dateobs = item
            break
    pointing = t[-1]

    return telescope,dateobs,pointing

def get_galfit_results(galfit_output_image,ngal=1):

    """
    INPUT:
    * galfit_output_image : will have model fits in header 2
    * ngal : number of galaxies modeled in galfit image.  the default is 1

    RETURN:
    results : [list], output assuming two galaxies would be:
              first entry: 'XC','YC','MAG','RE','N','AR','PA','NUMERROR' for gal 1
              second entry: 'XC','YC','MAG','RE','N','AR','PA','NUMERROR' for gal 2
              third entry: SKY
              fourth entry: CHISQNU


    """
    # read in image header
    hdu = fits.open(galfit_output_image)
    # extension 2 has the model info
    imheader = hdu[2].header
    #print(imheader)
    hdu.close()
    

    
    header=['XC','YC','MAG','RE','N','AR','PA','SKY','CHI2NU']


    # create empty list for each header
    results = []
    # loop over galaxies
    for i in range(ngal):
        # build header keywords
        galparams = []
        for h in header[:-2]:
            hkey = f"{i+1}_{h}"
            # query galheader
            galparams.append(imheader[hkey])
        numerror = False
        for g in galparams:
            if '*' in g:
                numerror = True
        galparams.append(numerror)
        results.append(galparams)
    # get sky
    nsky = ngal + 1
    #print("HEYYYY: nsky = ",nsky)
    t = imheader[f"{nsky}_SKY"]
    results.append(t)
    # get chisqy
    results.append(imheader['CHI2NU'])

    return results
    
def buildone(subdir,outdir,flist):
    print(subdir)

    #telescope,dateobs,pointing = get_params_from_name(subdir)
    #run = dateobs+'-'+pointing
    #if os.path.isdir(subdir) & (subdir.startswith('pointing')) & (subdir.find('-') > -1):
    if os.path.isdir(subdir):
        print('##########################################')
        print('##########################################')        
        print('WORKING ON DIRECTORY: ',subdir)
        print('##########################################')
        print('##########################################')
        
        #try:
        # move to subdirectory
        # adding the telescope and run so that we don't write over
        # images if galaxy was observed more than once

        # moving to subdirectory
        if args.verbose:
            print("moving to directory ",subdir)
        os.chdir(subdir)

        gal_outdir = os.path.join(outdir,subdir+"")
        print('out directory for this galaxy = ',gal_outdir)
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        if not os.path.exists(gal_outdir):
            os.mkdir(gal_outdir)

        p = galfit_dir(cutoutdir=subdir,outdir=gal_outdir)
        p.runall()
        #p.get_galfit_model(band='r')
        i = flist.index(args.oneimage)            
        # define previous gal for html links
        if i > 0:
            previous = (flist[i-1])
            #print('previous = ',previous)
        else:
            previous = None
        if i < len(flist)-1:
            next = flist[i+1]
            #print('next = ',next)
        else:
            next = None
        h = build_html_cutout(p,gal_outdir,previous=previous,next=next,tel=None,run=None)
        h.build_html()
        #except:
        #    print('WARNING: problem building webpage for ',subdir)
    

def display_image(image,percentile1=.5,percentile2=99.5,stretch='asinh',mask=None,sigclip=True,zoom=None):
    if zoom is not None:
        print("who's zoomin' who?")
        # display central region of image
        
        # get image dimensions and center
        xmax,ymax = image.shape
        xcenter = int(xmax/2)
        ycenter = int(ymax/2)
        
        # calculate new size to display based on zoom factor
        new_xradius = int(xmax/2/(float(zoom)))
        new_yradius = int(ymax/2/(float(zoom)))
        
        # calculate pixels to keep based on zoom factor
        x1 = xcenter - new_xradius
        x2 = xcenter + new_xradius
        y1 = ycenter - new_yradius
        y2 = ycenter + new_yradius
        
        # check to make sure limits are not outsize image dimensions
        if (x1 < 1):
            x1 = 1
        if (y1 < 1):
            y1 = 1
        if (x2 > xmax):
            x2 = xmax
        if (y2 > ymax):
            y2 = ymax

        # cut images to new size
        image = image[x1:x2,y1:y2]
        if mask is not None:
            mask = mask[x1:x2,y1:y2]
    # use inner 80% of image
    xdim,ydim = image.shape
    xmin = int(.1*xdim)
    xmax = int(.9*xdim)    
    ymin = int(.1*ydim)
    ymax = int(.9*ydim)
    if mask is not None:
        imdata = np.ma.array(image,mask=mask)
        
    else:
        imdata = image
    v1 = scoreatpercentile(imdata,percentile1)    
    v2 = scoreatpercentile(imdata,percentile2)
    
    if mask is not None:
        statim = image[~mask]
    else:
        statim = image

    if sigclip:
        if mask is not None:
            clipped_data = sigma_clip(image[xmin:xmax,ymin:ymax][~mask[xmin:xmax,ymin:ymax]],sigma_lower=1.5,sigma_upper=1.5,grow=10,stdfunc='mad_std')
        else:
            clipped_data = sigma_clip(image[xmin:xmax,ymin:ymax],sigma_lower=1.5,sigma_upper=1.5,grow=10,stdfunc='mad_std')            
    else:
        clipped_data = image[xmin:xmax,ymin:ymax]

    norm = simple_norm(clipped_data, stretch=stretch,max_percent=percentile2,min_percent=percentile1)

    plt.imshow(image, norm=norm,cmap='gray_r',origin='lower')#,vmin=v1,vmax=v2)
    

def make_png(fitsimage,outname,mask=None,ellipseparams=None):
    imdata,imheader = fits.getdata(fitsimage,header=True)
    fig = plt.figure(figsize=(6,6))
    ax = plt.subplot(projection=wcs.WCS(imheader))
    plt.subplots_adjust(top=.95,right=.95,left=.2,bottom=.15)
    display_image(imdata,sigclip=True,mask=mask)
    plt.xlabel('RA (deg)',fontsize=16)
    plt.ylabel('DEC (deg)',fontsize=16)
    if ellipseparams is not None:
        ax = plt.gca()
        plot_ellipse(ax,ellipseparams)
    plt.savefig(outname)        
    plt.close(fig)

def plot_ellipse(ax,ellipseparams):

    xc,yc,r,BA,PA = ellipseparams

    #print("just checking - adding ellipse drawing ",self.ellipseparams)
    
    # need to reset b to be consistent with galfit ellipticity
    BA = float(BA)
    PA = float(PA)
    #print('THETA inside phot wrapper',THETA, BA)
    b = BA*r
    eps = 1 - BA
    #print(self.b,self.eps,self.sma,BA)
    t = PA
    if t < 0:
        theta = (180. + t)
    else:
        theta = (t) # orientation in radians
    # EllipticalAperture gives rotation angle in radians from +x axis, CCW
    # matplotlib uses total width and height, not semi-major /minor axes
    ellipse = Ellipse((xc,yc), 2*r, 2*b, angle=theta,facecolor='None',edgecolor='r',lw=2)
    ax.add_patch(ellipse)

    
def display_galfit_model(galfile,percentile1=.5,percentile2=99.5,p1residual=5,p2residual=99,cmap='viridis',zoom=None,outdir=None,mask=None,ellipseparams=None,suffix=''):
      '''
      ARGS:
      galfile = galfit output image (with image, model, residual)
      percentile1 = min percentile for stretch of image and model
      percentile2 = max percentile for stretch of image and model
      p1residual = min percentile for stretch of residual
      p2residual = max percentile for stretch of residual
      cmap = colormap, default is viridis
      mask = bad pixel mask, with bad values set to True
      '''
      # model name

      #print("inside display_galfit_model, mask = ",mask)
      image,h = fits.getdata(galfile,1,header=True)
      model = fits.getdata(galfile,2)
      residual = fits.getdata(galfile,3)

      if zoom is not None:
         print("who's zoomin' who?")
         # display central region of image

         # get image dimensions and center
         ymax,xmax = image.shape
         xcenter = int(xmax/2)
         ycenter = int(ymax/2)

         # calculate new size to display based on zoom factor
         new_xradius = int(xmax/2/(float(zoom)))
         new_yradius = int(ymax/2/(float(zoom)))

         # calculate pixels to keep based on zoom factor
         x1 = xcenter - new_xradius
         x2 = xcenter + new_xradius
         y1 = ycenter - new_yradius
         y2 = ycenter + new_yradius
         
         # check to make sure limits are not outsize image dimensions
         if (x1 < 1):
            x1 = 1
         if (y1 < 1):
            y1 = 1
         if (x2 > xmax):
            x2 = xmax
         if (y2 > ymax):
            y2 = ymax

         # cut images to new size
         # python is data[row,col]
         image = image[y1:y2,x1:x2]
         model = model[y1:y2,x1:x2]
         residual = residual[y1:y2,x1:x2]         
         pass
      imwcs = wcs.WCS(h)
      images = [image,model,residual]
      titles = ['image','model','residual']
      if mask is not None:
          print("\nshape of image = ",image.shape)
          print("shape of mask = ",mask.shape)
          print()
          try:
              im = image[~mask]
              res = residual[~mask]
          except IndexError:
              im = image
              res = residual
              
              print("no mask for galfit ",galfile)
          norms = [simple_norm(im,'asinh',max_percent=percentile2),
                   simple_norm(im,'asinh',max_percent=percentile2),
                   simple_norm(res,'asinh',max_percent=percentile2,min_percent=20)]

      else:
          norms = [simple_norm(image,'asinh',max_percent=percentile2),
                   simple_norm(image,'asinh',max_percent=percentile2),
                   simple_norm(residual,'asinh',max_percent=percentile2)]

      outim = [f'galfit_image{suffix}.png',\
               f'galfit_model{suffix}.png',\
               f'galfit_residual{suffix}.png']
      
      if outdir is not None:
          outim = [os.path.join(outdir,f) for f in outim]
      for i,im in enumerate(images):
          fig = plt.figure(figsize=(6,6))          
          plt.subplot(1,1,1,projection=imwcs)
          plt.subplots_adjust(top=.95,right=.95,left=.2,bottom=.15)
          plt.imshow(im,origin='lower',cmap=cmap,norm=norms[i])
          #plt.colorbar(fraction=.08)
          plt.xlabel('RA (deg)',fontsize=16)
          plt.ylabel('DEC (deg)',fontsize=16)
          #plt.title(titles[i],fontsize=16)
          
          # TODO add ellipse to the residual image
          if (i == 2) and (ellipseparams is not None):
              # plot the ellipse
              plot_ellipse(plt.gca(),ellipseparams)
          plt.savefig(outim[i])

      plt.close(fig)

###########################################################
####  CLASSES
###########################################################

    
    
class galfit_dir():

    def __init__(self,cutoutdir=None,outdir=None):
        '''
        INPUT:
        * directory containing galfit results
        * output directory for png images

        This creates the png images for different cutouts
        '''
        # this returns the VFID
        # the galname should be the NED name that JM uses
        self.vfid = cutoutdir

        if args.verbose:
            print('inside galfit_dir, vfid = ',self.vfid)
        #print('cutoutdir = ',cutoutdir)
        #print('outdir = ',outdir)
        
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        self.outdir = outdir
        if args.verbose:
            print("self.outdir = ",self.outdir)
        self.cutoutdir = cutoutdir
        if args.verbose:
            print("self.cutoutdir = ",self.cutoutdir)
        
    def runall(self):
        self.get_ngal()        
        self.get_ned_name()
        self.get_file_names()
        self.get_ellipse_params()
        self.make_png_mask()

        
    def get_ngal(self):
        """check to see how many galaxies in FOV  """
        infile = open('galsFOV.txt','r')
        mylines = infile.readlines()
        self.ngal = len(mylines)

        if args.verbose:
            print(f"number of galaxies = {self.ngal}")
            
    def get_ned_name(self):
        """ get galaxy NED name by grabbing an image """

        # running from within galfit output directory
        t = glob.glob('*-custom-image-r.fits')[0]
        self.gname = os.path.basename(t).split('-custom')[0]
        print("galaxy nedname = ",self.gname)
        
                
    def get_file_names(self):
        search_string = '*custom-image-r.fits'
        #print(search_string)
        t = glob.glob(search_string)
        #print(t)
        
        
        self.maskimage = self.gname+'-custom-image-r-mask.fits'
        self.wisemaskimage = self.maskimage.replace('r-mask.fits','wise-mask.fits')

        if not os.path.exists(self.maskimage):
            print(f"WARNING: can not find mask image {self.maskimage}!!!")
            
    def get_ellipse_params(self):
        """ get ellipse parameters from the header of the mask image  """
        try:
            header = fits.getheader(self.maskimage)
            ellipseparams = []
            keywords = ['ELLIP_XC','ELLIP_YC','ELLIP_A','ELLIP_BA','ELLIP_PA']
            # get ellipse params
            for k in keywords:
                ellipseparams.append(header[k])
            self.ellipseparams = ellipseparams
        except FileNotFoundError:
            print("\nWARNING: mask file was not found!!! ",self.maskimage)
            self.ellipseparams = None
        


    def make_png_mask(self):
        # fitsimages and pngimages should be dictionaries
        # so I am not relying on where they are in the list
        imnames = ['mask','wisemask']
        self.image_keys = imnames
        # build dictionaries to store fits and png images,
        # setting to None if image is not available
        self.fitsimages = {i:None for i in imnames}
        self.pngimages = {i:None for i in imnames}
            
        self.fitsimages['mask'] = self.maskimage
        self.fitsimages['wisemask'] = self.wisemaskimage        

        mask = fits.getdata(self.maskimage)
        mask = mask > 0
            
        for i,f in enumerate(self.fitsimages): # loop over keys

            try:
                pngfile = os.path.join(self.outdir,os.path.basename(self.fitsimages[f]).replace('.fits','.png'))
            except TypeError:
                continue
            try:
                if i < 4:
                    make_png(self.fitsimages[f],pngfile,mask=None)
                elif i == (len(self.fitsimages)-2): # add ellipse to mask image
                    if self.ellipseparams is not None:
                        make_png(self.fitsimages[f],pngfile,ellipseparams=self.ellipseparams)
                    else:
                        make_png(self.fitsimages[f],pngfile)
                else:
                    make_png(self.fitsimages[f],pngfile)                    
                self.pngimages[f] = pngfile
            except FileNotFoundError:
                print('WARNING: can not find ',self.fitsimages[f])

            except TypeError:
                print('WARNING: problem making png for ',self.fitsimages[f])
                

    def get_galfit_model(self,band='r'):
        ''' read in galfit model and make png '''
        self.galfit = f"{self.gname}-{band}-out2.fits"
        print("looking for galfit file {self.galfit}")
        if os.path.exists(self.galfit):
            # store fit results
            if args.verbose:
                print("making png and storing results from galfit model ",self.galfit)

            if 'W' in band:
                mask = fits.getdata(self.wisemaskimage)
            else:
                mask = fits.getdata(self.maskimage)
            mask = mask > 0


            display_galfit_model(self.galfit,outdir=self.outdir,mask=mask,ellipseparams=self.ellipseparams,suffix=f"_{band}")

            outim = [f'galfit_image_{band}.png',f'galfit_model_{band}.png',f'galfit_residual_{band}.png']
        
            self.galimage = os.path.join(self.outdir,outim[0])
            self.galmodel = os.path.join(self.outdir,outim[1])
            self.galresidual = os.path.join(self.outdir,outim[2])        

            # store fitted parameters

            self.results = get_galfit_results(self.galfit,ngal=self.ngal)
            # make png of mask

            #t = rg.parse_galfit_1comp(self.galfit)
        
            #header_keywords=['1_XC','1_YC','1_MAG','1_RE','1_N','1_AR','1_PA','2_SKY','ERROR','CHI2NU']
            #self.xc, self.xc_err = t[0]
            #self.yc, self.yc_err = t[1]
            #self.mag, self.mag_err = t[2]
            #self.re, self.re_err = t[3]
            #self.nsersic, self.nsersic_err = t[4]
            #self.BA, self.BA_err = t[5]
            #self.PA, self.PA_err = t[6]
            #self.sky, self.sky_err = t[7]
            #self.error = t[8]
            #self.chi2nu = t[9]
        else:
            self.galimage = None
        
        pass
    

    
class build_html_cutout():

    def __init__(self,cutoutdir,outdir,previous=None,next=None,tel=None,run=None):
        ''' pass in instance of cutout_dir class and output directory '''
        #print("in build_html_cutout!")
        
        # instance of class cutout_dir
        self.cutout = cutoutdir

        outfile = os.path.join(outdir,self.cutout.vfid+'.html')
        if args.verbose:
            print("outfile = ",outfile)
        vfindices = np.arange(len(vfmain))
        self.vfindex = vfindices[vfmain['VFID'] == self.cutout.vfid]
        #print('inside build html')
        #print('coutdir = ',coutdir)
        #print('outfile = ',outfile)        
        self.html = open(outfile,'w')
        self.htmlhome = 'index.html'
        self.next = next
        self.previous = previous

        
        # for reference, this is the order of the png images
        #self.fitsimages = [self.rimage,self.haimage,self.csimage,\
        #              self.legacy_g,self.legacy_r,self.legacy_z,\
        #              self.w1,self.w2,self.w3,self.w4]

        
        #self.build_html()
    def build_html(self):
        self.write_header()
        self.write_navigation_links()
        # adding this here so we can inspect the masks quickly
        # can remove once we are done with masks

        bands = ['r','g','W1','W2','W3','W4']
        bands = ['g','r','W1','W1-fixBA','W2','W3','W3-fixBA','W4']
        bands = ['r','W1','W1-fixBA','W3','W3-fixBA']                
        for b in bands:
            self.cutout.get_galfit_model(band=b)
        
            self.write_galfit_images(band=b)
            self.write_galfit_table(band=b)
        #self.write_image_stats()
        #if self.cutout.legacy_flag:
        #    self.write_legacy_images()

        #self.write_sfr_images()
        #if self.cutout.wise_flag:
        #    self.write_wise_images()
        #self.write_halpha_images()
        #if self.cutout.galimage is not None:
        #    self.write_galfit_images()
        #    self.write_galfit_table()
        #try:
        #    self.write_phot_profiles()
        #except AttributeError:
        #    pass
        #self.write_mag_table()
        #self.write_morph_table()        
        self.write_navigation_links()
        self.close_html()
    def write_header(self):
        # title
        # home link
        # previous
        # next
        self.html.write('<html><body>\n')
        self.html.write('<style type="text/css">\n')
        self.html.write('table, td, th {padding: 5px; text-align: center; border: 1px solid black}\n')
        self.html.write('</style>\n')

    def write_navigation_links(self):
        # Top navigation menu--
        self.html.write(f'<h1>{self.cutout.vfid}-{self.cutout.gname}</h1>\n')

        self.html.write('<a href="../{}">Home</a>\n'.format(self.htmlhome))
        self.html.write('<br />\n')
        if self.previous is not None:
            previoushtml = "{}.html".format(self.previous)
            self.html.write('<a href="../{}/{}">Previous ({})</a>\n'.format(self.previous,previoushtml, self.previous))
            self.html.write('<br />\n')        
        if self.next is not None:
            nexthtml = "{}.html".format(self.next)
            self.html.write('<a href="../{}/{}">Next ({})</a>\n'.format(self.next,nexthtml,self.next))
            self.html.write('<br />\n')

    

    def write_wise_images(self):
        ''' w1 - w4 images '''
        self.html.write('<h2>WISE Images</h2>\n')
        pngimages = [self.cutout.pngimages['w1'],self.cutout.pngimages['w2'],\
                     self.cutout.pngimages['w3'],self.cutout.pngimages['w4']]
        wlabels = ['W1','W2','W3','W4']
        images=[]
        labels=[]
        for i,im in enumerate(pngimages):
            if im is not None:
                images.append(os.path.basename(im))
                labels.append(wlabels[i])

        write_table(self.html,images=images,labels=labels)
    

    def write_galfit_images(self,band='r'):
        ''' display galfit model and fit parameters for r-band image '''
        if self.cutout.galimage is not None:
            self.html.write(f'<h2>GALFIT {band} Modeling </h2>\n')
            if 'W' in band:
                maskpng = self.cutout.gname+'-custom-image-wise-mask.png'
            else:
                maskpng = self.cutout.gname+'-custom-image-r-mask.png'
            images = [self.cutout.galimage,self.cutout.galmodel,self.cutout.galresidual,\
                      self.cutout.pngimages['mask']]
            images = [os.path.basename(i) for i in images]        
            labels = ['Image', 'Model', 'Residual','Mask']
            write_table(self.html,images=images,labels=labels)

    
    def write_galfit_table(self,band='r'):
        ''' display galfit model and fit parameters for r-band image '''
        self.html.write(f'<h4>GALFIT Sersic Parameters for {band}</h4>\n')                
        labels=['XC','YC','MAG','RE','N','AR','PA','ERROR','SKY','CHI2NU']


        self.html.write('<table width="90%"; table-layout: fixed>\n')
        self.html.write('<tr>')
        for l in labels:
            self.html.write('<th>{}</th>'.format(l))
        self.html.write('</tr></p>\n')        

        sky = self.cutout.results[self.cutout.ngal]
        chisqnu = self.cutout.results[self.cutout.ngal+1]
        for i in range(self.cutout.ngal):
            data = self.cutout.results[i]
            self.html.write('<tr>')
            for d in data:
                self.html.write('<td>{}</td>'.format(d))
            self.html.write('<td>{}</td>'.format(sky))
            self.html.write('<td>{}</td>'.format(chisqnu))            
            self.html.write('</tr>\n')            

        self.html.write('</table>\n')
        

    def close_html(self):
        self.html.close()
# wrap

if __name__ == '__main__':


    import argparse

    parser = argparse.ArgumentParser(description ='create psf image from image that contains stars')

    #parser.add_argument('--table-path', dest = 'tablepath', default = '/Users/rfinn/github/Virgo/tables/', help = 'path to github/Virgo/tables')
    parser.add_argument('--cutoutdir',dest = 'cutoutdir', default=None, help='set to cutout directory. default is the current directory, like you are running from the cutouts/ directory')
    parser.add_argument('--verbose',dest = 'verbose',default=False, action='store_true', help='set for additional print statements')    
    parser.add_argument('--oneimage',dest = 'oneimage',default=None, help='give directory for one image')
    parser.add_argument('--outdir',dest = 'outdir',default='/data-pool/Halpha/html_dev/galfit/', help='output directory.  default is /data-pool/Halpha/html_dev/galfit/')    
     
    args = parser.parse_args()

    # get tables, define as a global variable
    vfmain = fits.getdata(homedir+'/research/Virgo/tables-north/v2/vf_v2_main.fits')
    ephot = fits.getdata(homedir+'/research/Virgo/tables-north/v2/vf_v2_legacy_ephot.fits')

    topdir = os.getcwd()

    # get directory list to use with Previous and Next links
    rfiles = os.listdir()
    rfiles.sort()


    
    outdir = args.outdir
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    # always run with one image
    # make sure that the image exists
    if not os.path.exists(args.oneimage):
        print(f"Could not find {args.oneimage} - please check the cutout directory name you provided")
        sys.exit()
    # find index in rfiles that corresponds to image
    coadd_index = rfiles.index(args.oneimage)
    indices = [np.arange(len(rfiles))[coadd_index]]

    buildone(args.oneimage,outdir,rfiles)
    

    os.chdir(topdir)
