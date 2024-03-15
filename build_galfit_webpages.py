#!/usr/bin/env python
'''
GOAL:
* create web page to inspect the cutouts

USAGE:
* run from cutouts directory

NOTES:
* using John Moustakas's code as a reference (https://github.com/moustakas/legacyhalos/blob/main/py/legacyhalos/virgofilaments.py#L1131-L1202)
* https://docs.astropy.org/en/stable/visualization/normalization.html#:~:text=The%20astropy.visualization%20module%20provides%20a%20framework%20for%20transforming,x%20represents%20the%20values%20in%20the%20original%20image%3A

TO DO: (I think these are all fixed!)
* DONE need to figure out how to handle repeated observations
  - don't overwrite directory

* fix how I combine unwise images when multiple images are returned
* same for galex

2024-03-10 : 
* add panel to show color-based continuum subtract (added this to Halpha panel, instead of showing CS with two stretches)
* stellar mass, sfr and ssfr images

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

def buildone(subdir,outdir,flist):
    print(subdir)

    telescope,dateobs,pointing = get_params_from_name(subdir)
    run = dateobs+'-'+pointing
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
        gal_outdir = os.path.join(outdir,subdir+"")
        #print('out directory for this galaxy = ',gal_outdir)
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        if not os.path.exists(gal_outdir):
            os.mkdir(gal_outdir)

        p = cutout_dir(cutoutdir=subdir,outdir=gal_outdir)
        p.runall()

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
        h = build_html_cutout(p,gal_outdir,previous=previous,next=next,tel=telescope,run=run)
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
          im = image[~mask]
          res = residual[~mask]
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
          #plt.close(fig)

###########################################################
####  CLASSES
###########################################################

    
    
class cutout_dir():

    def __init__(self,cutoutdir=None,outdir=None):
        '''
        INPUT:
        * directory containing cutouts
        * output directory for png images

        This creates the png images for different cutouts
        '''
        self.gname = os.path.basename(os.path.dirname(cutoutdir))
        if len(self.gname) < 1:
            self.gname = cutoutdir
        self.vfid = self.gname.split('-')[0]
        #print('inside cutoutdir, gname = ',self.gname)
        #print('cutoutdir = ',cutoutdir)
        #print('outdir = ',outdir)        
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        self.outdir = outdir
        self.cutoutdir = cutoutdir
    def runall(self):
        self.get_file_names()
        self.get_ellipse_params()
        
        self.get_galfit_model()
    def get_file_names(self):
        search_string = os.path.join(self.cutoutdir,'*custom-image-r.fits')
        #print(search_string)
        t = glob.glob(search_string)
        #print(t)
        
        self.rimage = t[0]
        self.galname = self.rimage.split('-custom')[0]
        self.gimage = self.rimage.replace('-r.fits','-g.fits')

        self.w1 = self.rimage.replace('-r.fits','-W1.fits')
        self.w2 = self.rimage.replace('-r.fits','-W2.fits')
        self.w3 = self.rimage.replace('-r.fits','-W3.fits')
        self.w4 = self.rimage.replace('-r.fits','-W4.fits')
        
        self.maskimage = self.rimage.replace('.fits','-mask.fits')


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
        
    def get_legacy_names(self):
        ''' get names of legacy images  '''
        legdir = os.path.join(self.cutoutdir,'legacy')
        self.legacy_g = glob.glob(legdir+'/*-g.fits')[0]
        self.legacy_r = glob.glob(legdir+'/*-r.fits')[0]
        self.legacy_z = glob.glob(legdir+'/*-z.fits')[0]
    def get_legacy_jpg(self):
        ''' copy jpg to local directory '''
        legdir = os.path.join(self.cutoutdir,'legacy')        
        legacy_jpg = glob.glob(legdir+'/*.jpg')[0]
        jpeg_data = Image.open(legacy_jpg)
        fig = plt.figure(figsize=(6,6))


        #hdu = fits.open(filename)[0]
        header = fits.getheader(self.legacy_g)
        imwcs = wcs.WCS(header)
        plt.subplot(projection=imwcs)
        plt.subplots_adjust(left=.2,bottom=.15,top=.95,right=.95)        
        plt.imshow(jpeg_data, origin='lower')
        plt.xlabel('RA (deg)',fontsize=16)
        plt.ylabel('Dec (deg)',fontsize=16)
        self.legacy_jpg = os.path.join(self.outdir,os.path.basename(legacy_jpg))        
        plt.savefig(self.legacy_jpg)
        plt.close(fig)



    def make_png_plots(self):
        # fitsimages and pngimages should be dictionaries
        # so I am not relying on where they are in the list
        imnames = ['r','ha','cs','csgr','legacy_g','legacy_r','legacy_z',\
                   'w1','w2','w3','w4',\
                   'mask','nuv']
        self.image_keys = imnames
        # build dictionaries to store fits and png images,
        # setting to None if image is not available
        self.fitsimages = {i:None for i in imnames}
        self.pngimages = {i:None for i in imnames}
        keys = imnames[0:4]
        imlist = [self.rimage,self.haimage,self.csimage,self.csgrimage]
        for i,k in enumerate(keys):
            self.fitsimages[k] = imlist[i]
        
        if self.legacy_flag:
            keys = imnames[4:7]
            imlist = [self.legacy_g,self.legacy_r,self.legacy_z]
            for i,k in enumerate(keys):
                self.fitsimages[k] = imlist[i]
            
        if self.wise_flag:
            keys = imnames[7:11]            
            imlist = [self.w1,self.w2,self.w3,self.w4]
            for i,k in enumerate(keys):
                self.fitsimages[k] = imlist[i]
            
        self.fitsimages['mask'] = self.maskimage
        if self.nuv_flag:
            self.fitsimages['nuv'] = self.nuv

        mask = fits.getdata(self.maskimage)
        mask = mask > 0
            
        for i,f in enumerate(self.fitsimages): # loop over keys

            try:
                pngfile = os.path.join(self.outdir,os.path.basename(self.fitsimages[f]).replace('.fits','.png'))
            except TypeError:
                continue
            try:
                if i < 4:
                    make_png(self.fitsimages[f],pngfile,mask=mask)
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


    def make_cs_png(self,gr=False):
        if gr:
            csdata,csheader = fits.getdata(self.csgrimage,header=True)
            imx,imy,keepflag = get_galaxies_fov(self.csimage,vfmain['RA'],vfmain['DEC'])
        else:
            csdata,csheader = fits.getdata(self.csimage,header=True)
            imx,imy,keepflag = get_galaxies_fov(self.csimage,vfmain['RA'],vfmain['DEC'])
            
        mask = fits.getdata(self.maskimage)
        mask = mask > 0
        #galsize=60/(abs(csheader['CD1_1'])*3600)        
        p2 = [99.5,99.9]
        stretchs = ['asinh','linear']
        for i,s in enumerate(stretchs):
            fig = plt.figure(figsize=(6,6))
            plt.subplot(projection = wcs.WCS(csheader))
            plt.subplots_adjust(bottom=.15,left=.2,right=.95,top=.95)
            ax = plt.gca()
            #clipped_data = sigma_clip(image[xmin:xmax,ymin:ymax],sigma_lower=1.5,sigma_upper=1.5,grow=10,stdfunc='mad_std')            
            display_image(csdata,stretch=s,percentile1=.15,percentile2=p2[i],mask=mask)
            # mark VF galaxies
            #plot_vf_gals(imx,imy,keepflag,vfmain,ax,galsize=galsize)
            suffix = "-{}.png".format(p2[i])
            if gr:
                pngfile = os.path.join(self.outdir,os.path.basename(self.csgrimage).replace('.fits',suffix))
            else:
                pngfile = os.path.join(self.outdir,os.path.basename(self.csimage).replace('.fits',suffix))
            plt.xlabel('RA (deg)',fontsize=16)
            plt.ylabel('DEC (deg)',fontsize=16)        
            plt.savefig(pngfile)
            plt.close(fig)
            if i == 0:
                if gr:
                    self.csgr_png1 = pngfile
                else:
                    self.cs_png1 = pngfile
            elif i == 1:
                if gr:
                    self.csgr_png2 = pngfile
                else:
                    self.cs_png2 = pngfile 
    def get_galfit_model(self,band):
        ''' read in galfit model and make png '''
        self.galfit = self.rimage.replace('custom-image-r.fits','-1Comp-galfit-out.fits')
        if os.path.exists(self.galfit):
            # store fit results

            mask = fits.getdata(self.maskimage)
            mask = mask > 0


            display_galfit_model(self.galfit,outdir=self.outdir,mask=mask,ellipseparams=self.ellipseparams)

            outim = ['galfit_image.png','galfit_model.png','galfit_residual.png']
        
            self.galimage = os.path.join(self.outdir,outim[0])
            self.galmodel = os.path.join(self.outdir,outim[1])
            self.galresidual = os.path.join(self.outdir,outim[2])        

            # store fitted parameters

            t = rg.parse_galfit_1comp(self.galfit)
        
            #header_keywords=['1_XC','1_YC','1_MAG','1_RE','1_N','1_AR','1_PA','2_SKY','ERROR','CHI2NU']
            self.xc, self.xc_err = t[0]
            self.yc, self.yc_err = t[1]
            self.mag, self.mag_err = t[2]
            self.re, self.re_err = t[3]
            self.nsersic, self.nsersic_err = t[4]
            self.BA, self.BA_err = t[5]
            self.PA, self.PA_err = t[6]
            self.sky, self.sky_err = t[7]
            self.error = t[8]
            self.chi2nu = t[9]
        else:
            self.galimage = None
        
        pass
    

    def get_phot_tables(self):
        ''' read in phot tables and make plot of flux and sb vs sma '''


        # open data files
        cs_galfit_phot = self.csimage.replace('.fits','-GAL-phot.fits')
        cs_gphot = fits.getdata(cs_galfit_phot)
        #cs_phot = self.csimage.replace('.fits','-phot.fits')        
        r_galfit_phot = self.rimage.replace('.fits','-GAL-phot.fits')
        r_gphot = fits.getdata(r_galfit_phot)

        # photutils flux
        cs_galfit_phot = self.csimage.replace('.fits','-phot.fits')
        cs_phot = fits.getdata(cs_galfit_phot)
        #cs_phot = self.csimage.replace('.fits','-phot.fits')        
        r_galfit_phot = self.rimage.replace('.fits','-phot.fits')
        r_phot = fits.getdata(r_galfit_phot)
        

        # define colors - need this for plotting line and fill_between in the same color
        mycolors = plt.rcParams['axes.prop_cycle'].by_key()['color']

        # plot enclosed flux        
        fig = plt.figure(figsize=(6,6))
        plt.subplots_adjust(left=.15,bottom=.1,right=.95,top=.95)
        tabs = [r_gphot,cs_gphot,r_phot,cs_phot]
        labels = ['galfit r','galfit Halphax100','photutil r','photutil Halphax100']
        alphas = [1,.4,.6,.4]
        for i,t in enumerate(tabs):
            y0 = t['flux_erg']            
            y1 = t['flux_erg']+t['flux_erg_err']
            y2 = t['flux_erg']-t['flux_erg_err']

            if (i == 1) + (i == 3):
                y0=y0*100
                y1 = y1*100
                y2 = y2*100
            plt.fill_between(t['sma_arcsec'],y1,y2,label=labels[i],alpha=alphas[i],color=mycolors[i])
            # also plot line because you can't see the result when the error is small
            # this should fix issue #18 in Virgo github
            plt.plot(t['sma_arcsec'],y0,'-',lw=2,color=mycolors[i])

        plt.xlabel('SMA (arcsec)',fontsize=16)
        plt.ylabel('Flux (erg/s/cm^2/Hz)',fontsize=16)
        plt.gca().set_yscale('log')
        plt.gca().set_xscale('log')
        plt.legend(loc='lower right')        
        self.efluxsma_png = os.path.join(self.outdir,self.gname+'-enclosed-flux.png')
        plt.savefig(self.efluxsma_png)
        plt.close(fig)
        
        # plot mag sma
        fig = plt.figure(figsize=(6,6))
        plt.subplots_adjust(left=.15,bottom=.1,right=.95,top=.95)
        #tabs = [r_gphot,cs_gphot]
        #labels = ['r','Halpha']
        tabs = [r_gphot,cs_gphot,r_phot,cs_phot]
        labels = ['galfit r','galfit Halpha','photutil r','photutil Halpha']
        ncolor = 0
        for i,t in enumerate(tabs):
            y0 = t['mag']
            y1 = t['mag']+t['mag_err']
            y2 = t['mag']-t['mag_err']
            if (i == 1) + (i == 3):
                alpha=.4
            else:
                alpha=1

            plt.fill_between(t['sma_arcsec'],y1,y2,label=labels[i],alpha=alphas[i],color=mycolors[i])
            # also plot line because you can't see the result when the error is small
            # this should fix issue #18 in Virgo github
            plt.plot(t['sma_arcsec'],y0,'-',lw=2,color=mycolors[i])

            
        plt.xlabel('SMA (arcsec)',fontsize=16)
        plt.ylabel('magnitude (AB)',fontsize=16)
        plt.gca().set_xscale('log')
        plt.gca().invert_yaxis()        
        self.emagsma_png = os.path.join(self.outdir,self.gname+'-mag-sma.png')
        plt.legend(loc='lower right')        
        plt.savefig(self.emagsma_png)
        plt.close(fig)
        
        # plot sb erg vs sma
        fig = plt.figure(figsize=(6,6))
        plt.subplots_adjust(left=.15,bottom=.1,right=.95,top=.95)
        tabs = [r_gphot,cs_gphot,r_phot,cs_phot]
        labels = ['galfit r','galfit Halphax100','photutil r','photutil Halphax100']

        for i,t in enumerate(tabs):
            y0 = t['sb_erg_sqarcsec']
            y1 = t['sb_erg_sqarcsec']+t['sb_erg_sqarcsec_err']
            y2 = t['sb_erg_sqarcsec']-t['sb_erg_sqarcsec_err']
            if (i == 1) + (i == 3):
                y0 = y0*100
                y1 = y1*100
                y2 = y2*100
            plt.fill_between(t['sma_arcsec'],y1,y2,label=labels[i],alpha=alphas[i],color=mycolors[i])
            # also plot line because you can't see the result when the error is small
            # this should fix issue #18 in Virgo github
            plt.plot(t['sma_arcsec'],y0,'-',lw=2,color=mycolors[i])

                
        plt.xlabel('SMA (arcsec)',fontsize=16)
        plt.ylabel('SB (erg/s/cm^2/Hz/arcsec^2)',fontsize=16)
        plt.gca().set_yscale('log')
        plt.gca().set_xscale('log')
        plt.legend()
        self.sbfluxsma_png = os.path.join(self.outdir,self.gname+'-sb-sma.png')
        plt.savefig(self.sbfluxsma_png)
        plt.close(fig)
                             
        # plot mag sma
        fig = plt.figure(figsize=(6,6))
        plt.subplots_adjust(left=.15,bottom=.1,right=.95,top=.95)
        tabs = [r_gphot,cs_gphot,r_phot,cs_phot]
        labels = ['galfit r','galfit Halpha','photutil r','photutil Halpha']
        
        for i,t in enumerate(tabs):
            y0 = t['sb_mag_sqarcsec']
            y1 = t['sb_mag_sqarcsec']+t['sb_mag_sqarcsec_err']
            y2 = t['sb_mag_sqarcsec']-t['sb_mag_sqarcsec_err']

            plt.fill_between(t['sma_arcsec'],y1,y2,label=labels[i],alpha=alphas[i],color=mycolors[i])
            # also plot line because you can't see the result when the error is small
            # this should fix issue #18 in Virgo github
            plt.plot(t['sma_arcsec'],y0,'-',lw=2,color=mycolors[i])
            
        plt.xlabel('SMA (arcsec)',fontsize=16)
        plt.ylabel('Surface Brightness (mag/arcsec^2)',fontsize=16)
        plt.gca().set_xscale('log')
        plt.gca().invert_yaxis()
        plt.legend()        
        self.sbmagsma_png = os.path.join(self.outdir,self.gname+'-sbmag-sma.png')
        plt.savefig(self.sbmagsma_png)
        plt.close(fig)
    
class build_html_cutout():

    def __init__(self,cutoutdir,outdir,previous=None,next=None,tel=None,run=None):
        ''' pass in instance of cutout_dir class and output directory '''
        #print("in build_html_cutout!")
        self.cutout = cutoutdir


        outfile = os.path.join(outdir,self.cutout.gname+'.html')
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

        self.telescope = tel
        self.run=run
        
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
        self.write_galfit_images()
        self.write_image_stats()
        if self.cutout.legacy_flag:
            self.write_legacy_images()

        self.write_sfr_images()
        if self.cutout.wise_flag:
            self.write_wise_images()
        self.write_halpha_images()
        if self.cutout.galimage is not None:
            self.write_galfit_images()
            self.write_galfit_table()
        try:
            self.write_phot_profiles()
        except AttributeError:
            pass
        self.write_mag_table()
        self.write_morph_table()        
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
        self.html.write('<h1>{}</h1>\n'.format(self.cutout.gname))

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

    def write_image_stats(self):
        self.html.write('<h2>Image Statistics</h2>\n')        
        labels=['Telescope','Run','Pointing','R FWHM <br> (arcsec)','H&alpha; FWHM <br> (arcsec)','Filter Ratio','Filter Correction']
        if 'BOK' in self.cutout.cutoutdir:
            tel = 'BOK'
        elif 'HDI' in self.cutout.cutoutdir:
            tel = 'HDI'
        elif 'INT' in self.cutout.cutoutdir:
            tel = 'INT'
        elif 'MOS' in self.cutout.cutoutdir:
            tel = 'MOS'


        # add code to match to the correct line in the full halpha table
        telescope,dateobs,pointing = get_params_from_name(self.cutout.cutoutdir+".fits")

        # find matching rows
        vfid_match_index = np.arange(len(fullha))[fullha['VFID'] == self.cutout.vfid]
        if len(vfid_match_index) > 1:
            print("found duplicate",telescope,dateobs)
            print("directory name = ",self.cutout.cutoutdir)
            # duplicates - need to find the correct match

            # check if telescope matches directory
            for j in vfid_match_index:
                #print(fullha['TEL'][j],telescope, fullha['DATE-OBS'][j],dateobs)
                if fullha['TEL'][j] in telescope:
                    #print("Matched telescope")
                    if (int(fullha['DATE-OBS'][j]) == int(dateobs)):
                        #print("found a match!")
                        fullhaindex = j
                        break
        else:
            fullhaindex = vfid_match_index[0]
        self.fullhaindex = fullhaindex            
        myrow = fullha[fullhaindex]
        #print('HALPHA FILE')
        #print(self.vfindex)
        #print(myrow)
        colnames = ['POINTING','R_FWHM','H_FWHM','FILTER_RATIO','FILT_COR']


        
        # the pointing name from the vfha table might not be correct
        # b/c if the galaxy was observed multiple times, only one pointing will be included
        # I should instead construct the pointing from the image names
        pointing = myrow['POINTING']

        
        ## get telescope name to use to split on
        ##print('I think the cutout directory name is ',self.cutout.cutoutdir)
        #t = self.cutout.cutoutdir.split(tel)        
        #matchstring =  f"{tel}{t[1]}"
        ## check in coadd list for matching file name
        #coadd_list = open('../virgo-coadds-fullpath.txt')
        #
        #for line in coadd_list:
        #    if matchstring in line:
        #        pointing = line.rstrip()
        #        pointing = os.path.basename(pointing).replace("-r.fits","").replace('-r-shifted.fits','').replace('-R.fits','')
        #        #print("found matching coadd")
        #        break

        #data = [f'<a href="http://facultyweb.siena.edu/~rfinn/virgo/coadds/{pointing}/{pointing}.html">{pointing}</a>', \
        #        "{:.2f}".format(myrow['R_FWHM'][0]),\
        #        "{:.2f}".format(myrow['H_FWHM'][0]),\
        #        "{:.4f}".format(myrow['FILTER_RATIO'][0]),\
        #        "{:.2f}".format(myrow['FILT_COR'][0])]
        data = [f'<a href="../../coadds/{pointing}/{pointing}.html">{pointing}</a>', \
                "{:.2f}".format(myrow['R_FWHM']),\
                "{:.2f}".format(myrow['H_FWHM']),\
                "{:.4f}".format(myrow['FILTER_RATIO']),\
                "{:.2f}".format(myrow['FILT_COR'])]

        #print('self.run = ',self.run)
        data.insert(0,self.run)
        data.insert(0,self.telescope)
        write_text_table(self.html,labels,data)        
    def write_sfr_images(self):
        ''' legacy jpeg, galex nuv, halpha, w4 '''
        self.html.write('<h2>Star-Formation Indicators</h2>\n')
        #self.fitsimages = [self.rimage,self.haimage,self.csimage,\ # 0,1,2
        #                   self.legacy_g,self.legacy_r,self.legacy_z,\ # 3,4,5
        #                   self.w1,self.w2,self.w3,self.w4] # 6,7,8,9

        if self.cutout.nuv_flag:
            images = [self.cutout.pngimages['nuv'],\
                      self.cutout.cs_png1,\
                      self.cutout.pngimages['w3'],self.cutout.pngimages['w4']]
            labels = ['NUV','H&alpha;','W3','W4']
        else:
            try:
                images = [self.cutout.legacy_jpg,\
                          self.cutout.pngimages['cs'],\
                          self.cutout.pngimages['w3'],self.cutout.pngimages['w4']]                      

                labels = ['Legacy','H&alpha;','W3','W4']
            except IndexError:
                print("WARNING: problem plotting sfr images")
                return
            except AttributeError:
                print("WARNING: problem plotting sfr images, probably with legacy jpg")
                return
        newim = []
        for i in images:
            if i is not None:
                newim.append(os.path.basename(i))
            else:
                newim.append(i)
        images = newim
        #images = [os.path.basename(i) for i in images]            
        write_table(self.html,images=images,labels=labels)
        pass

    
    def write_halpha_images(self):
        '''  r, halpha, cs, and mask images '''
        self.html.write('<h2>Halpha Images</h2>\n')        
        #images = [self.cutout.pngimages['r'],self.cutout.pngimages['ha'],self.cutout.cs_png1,self.cutout.cs_png2]
        images = [self.cutout.pngimages['r'],self.cutout.pngimages['ha'],self.cutout.cs_png1,self.cutout.csgr_png1]
        # just changing order to see if halpha image is still the biggest in the table, re issue #15
        # the second was still the biggest
        # so what if we also change the label
        # seems to scale with label
        #images = [self.cutout.pngimages['ha'],self.cutout.pngimages['r'],self.cutout.cs_png1,self.cutout.cs_png2]        
        images = [os.path.basename(i) for i in images]

        #labels = ['R-band Image','H&alpha;+Cont','CS, stretch 1','CS, stretch 2']

        labels = ['R-band Image','H&alpha;+Cont','CS from ZP ratio','CS from ZP and g-r cor']        
        #labels = ['Halpha+Cont','R','CS, stretch 1','CS, stretch 2']        
        write_table(self.html,images=images,labels=labels)

    def write_mstar_sfr_images(self):
        '''  TODO : add panel for stellar mass, sfr and ssfr '''
        self.html.write('<h2>Halpha Images</h2>\n')        
        images = [self.cutout.pngimages['r'],self.cutout.pngimages['ha'],self.cutout.cs_png1,self.cutout.cs_png2]
        # just changing order to see if halpha image is still the biggest in the table, re issue #15
        # the second was still the biggest
        # so what if we also change the label
        # seems to scale with label
        #images = [self.cutout.pngimages['ha'],self.cutout.pngimages['r'],self.cutout.cs_png1,self.cutout.cs_png2]        
        images = [os.path.basename(i) for i in images]

        labels = ['R-band Image','H&alpha;+Cont','CS, stretch 1','CS, stretch 2']
        #labels = ['Halpha+Cont','R','CS, stretch 1','CS, stretch 2']        
        write_table(self.html,images=images,labels=labels)
        
    def write_legacy_images(self):
        ''' jpg, g,r,z legacy images '''
        self.html.write('<h2>Legacy Images</h2>\n')

        images = [self.cutout.pngimages['legacy_g'],self.cutout.pngimages['legacy_r'],self.cutout.pngimages['legacy_z']]
        images = [os.path.basename(i) for i in images]        
        images.insert(0,os.path.basename(self.cutout.legacy_jpg))        
        labels = ['grz','g','r','z']
        write_table(self.html,images=images,labels=labels)

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
    
    def write_galex_images(self):
        ''' right now just nuv '''
        self.html.write('<h2>GALEX Images</h2>\n')                
        pass

    def write_galfit_images(self):
        ''' display galfit model and fit parameters for r-band image '''
        if self.cutout.galimage is not None:
            self.html.write('<h2>GALFIT r-band Modeling </h2>\n')                
            images = [self.cutout.galimage,self.cutout.galmodel,self.cutout.galresidual,\
                      self.cutout.pngimages['mask']]
            images = [os.path.basename(i) for i in images]        
            labels = ['Image', 'Model', 'Residual','Mask']
            write_table(self.html,images=images,labels=labels)

    
    def write_galfit_table(self):
        ''' display galfit model and fit parameters for r-band image '''
        self.html.write('<h4>GALFIT Sersic Parameters</h4>\n')                
        labels=['1_XC','1_YC','1_MAG','1_RE','1_N','1_AR','1_PA','2_SKY','ERROR','CHI2NU']
        data = ['{:.1f}+/-{:.1f}'.format(self.cutout.xc, self.cutout.xc_err),
                '{:.1f}+/-{:.1f}'.format(self.cutout.yc, self.cutout.yc_err),
                '{:.2f}+/-{:.2f}'.format(self.cutout.mag, self.cutout.mag_err), 
                '{:.2f}+/-{:.2f}'.format(self.cutout.re, self.cutout.re_err),
                '{:.2f}+/-{:.2f}'.format(self.cutout.nsersic, self.cutout.nsersic_err),
                '{:.2f}+/-{:.2f}'.format(self.cutout.BA, self.cutout.BA_err),
                '{:.1f}+/-{:.1f}'.format(self.cutout.PA, self.cutout.PA_err),
                '{:.2f}+/-{:.2f}'.format(self.cutout.sky, self.cutout.sky_err),
                '{}'.format(self.cutout.error),
                '{:.2f}'.format(self.cutout.chi2nu)]
        
        #print(data)
        write_text_table(self.html,labels,data)
        pass

    def write_phot_profiles(self):
        ''' photometry table with galfit and photutil results '''
        self.html.write('<h2>Elliptical Photometry</h2>\n')
        self.html.write('<p>using galfit and photutil geometry</p>\n')                        
        images = [self.cutout.efluxsma_png,self.cutout.emagsma_png,self.cutout.sbfluxsma_png,self.cutout.sbmagsma_png]
        images = [os.path.basename(i) for i in images]
        labels = ['Enclosed Flux','Enclosed Magnitude','Surface Brightness','Surface Brightness']
        write_table(self.html,images=images,labels=labels)

    def write_mag_table(self):
        self.html.write('<h2>R-band AB Magnitudes</h2>\n')        
        labels=['mag24','mag25', 'mag26','mag_petro','mag_galfit','mag_photutil'] 
        #myrow = vfha[self.vfindex]
        myrow = fullha[self.fullhaindex]
        colnames = ['M24','M25','M26','PETRO_MAG','GAL_MAG','ELLIP_SUM_MAG']
        #data = ["{:.2f}".format(myrow[c][0]) for c in colnames]
        data = ["{:.2f}".format(myrow[c]) for c in colnames]
        write_text_table(self.html,labels,data)        

        self.html.write('<h2>Halpha AB Magnitudes/SFR</h2>\n')        
        labels=['mag16','mag17','mag_petro','SSFR_IN','SSFR_OUT'] 
        #myrow = vfha[self.vfindex]
        myrow = fullha[self.fullhaindex]
        colnames = ['HM16','HM17','HPETRO_MAG','SSFR_IN','SSFR_OUT']
        try:
            data = ["{:.2f}".format(myrow[c][0]) for c in colnames]
        except IndexError:
            data = ["{:.2f}".format(myrow[c]) for c in colnames]            
        write_text_table(self.html,labels,data)        
        
    def write_morph_table(self):
        self.html.write('<h2>Morphology Parameters</h2>\n')        
        labels=['Band','Gini','M20','Asym','C30','Petro Con']
        
        #myrow = vfha[self.vfindex]
        myrow = fullha[self.fullhaindex]
        colnames = ['ELLIP_GINI','ELLIP_M20','ELLIP_ASYM','C30','PETRO_CON']
        colnames2 = ['ELLIP_HGINI','ELLIP_HM20','ELLIP_HASYM','HC30','HPETRO_CON']
        try: # this is required to work with fullha instead of vfha
            data = ["{:.2f}".format(myrow[c][0]) for c in colnames]
        except IndexError:
            data = ["{:.2f}".format(myrow[c]) for c in colnames]            
        data.insert(0,'r')
        try:
            data2 = ["{:.2f}".format(myrow[c][0]) for c in colnames2]
        except IndexError:
            data2 = ["{:.2f}".format(myrow[c]) for c in colnames2]            
        data2.insert(0,'Halpha')
        write_text_table(self.html,labels,data,data2=data2)
    def close_html(self):
        self.html.close()
# wrap

if __name__ == '__main__':


    import argparse

    parser = argparse.ArgumentParser(description ='create psf image from image that contains stars')

    #parser.add_argument('--table-path', dest = 'tablepath', default = '/Users/rfinn/github/Virgo/tables/', help = 'path to github/Virgo/tables')
    parser.add_argument('--cutoutdir',dest = 'cutoutdir', default=None, help='set to cutout directory. default is the current directory, like you are running from the cutouts/ directory')
    parser.add_argument('--oneimage',dest = 'oneimage',default=None, help='give directory for one image')
    parser.add_argument('--outdir',dest = 'outdir',default='/data-pool/Halpha/html_dev/cutouts/', help='output directory.  default is /data-pool/Halpha/html_dev/cutouts/')    
     
    args = parser.parse_args()

    # get tables, define as a global variable
    vfmain = fits.getdata(homedir+'/research/Virgo/tables-north/v2/vf_v2_main.fits')
    vfha = fits.getdata(homedir+'/research/Virgo/tables-north/v2/vf_v2_halpha.fits')
    #fullha = fits.getdata('../halphagui-output-combined-2023-Aug-27.fits')
    # adding new filename
    fullha = fits.getdata('../halphagui-output-combined-2024-Jan-04.fits')        

    
    if args.cutoutdir is not None:
        os.chdir(args.cutoutdir)

    # get directory list to use with Previous and Next links
    rfiles = os.listdir()
    rfiles.sort()
    outdir = args.outdir
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    if args.oneimage is not None:
        # make sure that the image exists
        if not os.path.exists(args.oneimage):
            print(f"Could not find {args.oneimage} - please check the cutout directory name you provided")
            sys.exit()
        # find index in rfiles that corresponds to image
        #try:
        coadd_index = rfiles.index(args.oneimage)
        indices = [np.arange(len(rfiles))[coadd_index]]
        #print('when selecting one image, indices = ',indices,rfiles[indices[0]])
        buildone(args.oneimage,outdir,rfiles)
        #except ValueError:
        #    rfiles = [args.oneimage]
        #    indices = np.arange(len(rfiles))
    else:
        indices = np.arange(len(rfiles))

        image_pool = mp.Pool(mp.cpu_count())
        #image_pool = mp.pool.ThreadPool(mp.cpu_count())    
        myresults = [image_pool.apply_async(buildone,args=(subdir,outdir,rfiles,)) for subdir in rfiles]
    
        image_pool.close()
        image_pool.join()
        #image_results = [r.get() for r in myresults]

    
