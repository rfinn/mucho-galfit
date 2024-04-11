from astropy.io import fits
from matplotlib import pyplot as plt
import numpy as np
from astropy import wcs
import os

def write_coadd_prop_table(html,filter,zp,fwhm_arcsec):
    html.write('<h3>Image Characteristics</h3>\n')
    html.write('<table>\n')
    html.write('<tr>')
    html.write('<th ">Filter</th>\n')
    html.write('<th ">ZP<br />mag</th>\n')
    html.write('<th ">PSF FWHM <br />arcsec</th>\n')
    html.write('</tr>\n')
    html.write('<tr><td>{}</td><td>{:.2f}</td><td>{:.2f}</td>\n'.format(filter,zp,fwhm_arcsec))
    html.write('</tr>\n')
    html.write('</table>\n')

def write_table(html,images=None,labels=None,images2=None):    
    html.write('<table width="90%"; table-layout: fixed>\n')
    html.write('<tr>')
    for l in labels:
        html.write('<th>{}</th>'.format(l))
    html.write('</tr></p>\n')        
    html.write('<tr>')
    for i in images:
        html.write('<td><a href="{0}"><img src="{1}" alt="Missing file {0}" height="auto" width="100%"></a></td>'.format(i,i))
    html.write('</tr>\n')
    if images2 is not None:
        html.write('<tr>')
        for i in images2:
            html.write('<td><a href="{0}"><img src="{1}" alt="Missing file {0}" height="auto" width="100%"></a></td>'.format(i,i))            

        html.write('</tr>\n')            
    
    html.write('</table>\n')

def write_text_table(html,labels,data,data2=None):    
    html.write('<table width="90%"; table-layout: fixed>\n')
    html.write('<tr>')
    for l in labels:
        html.write('<th>{}</th>'.format(l))
    html.write('</tr></p>\n')        
    html.write('<tr>')
    for d in data:
        html.write('<td>{}</td>'.format(d))
    html.write('</tr>\n')            
    if data2 is not None:
        html.write('<tr>')
        for d in data2:
            html.write('<td>{}</td>'.format(d))
        html.write('</tr>\n')            

    html.write('</table>\n')

def get_galaxies_fov(imagename,RA,DEC):
    '''
    get the galaxies in the FOV  

    INPUT:
    * imagename
    * RA - array of RA values, like v.main['RA']
    * DEC - array of DEC values, like v.main['DEC']

    RETURNS:
    * imx, imy - RA and DEC, translated to image pixel coordinates
    * keepflag - boolean array, same length as RA and DEC
    * true for galaxies that fall within FOV of image

    '''
    # get image header
    header = fits.getheader(imagename)
    # define instance of wcs
    imwcs = wcs.WCS(header)

    # convert RA and DEC from catalog to pixels coordinates
    imx,imy = imwcs.wcs_world2pix(RA,DEC,1)

    # get image dimensions
    xmax = header['NAXIS1']
    ymax = header['NAXIS2']    

    # keep galaxies that fall within image boundary
    keepflag = (imx > 1) & (imx < xmax) & (imy > 1) & (imy < ymax)

    # should also check the weight image and remove galaxies with weight=0
    # this won't take care of images with partial exposures, but we can deal with that later...
    # TODO - how to handle images with partial exposures, meaning only part of galaxy is in FOV?
    if imagename.find('shifted.fits') > -1:
        weightimage = imagename.replace('-r-shifted.fits','-r.weight-shifted.fits')
    else:
        weightimage = imagename.replace('.fits','.weight.fits')
    if 'MOS' not in imagename:
        if os.path.exists(weightimage):
            print()
            print("cross checking object locations with weight image")
            print()
            whdu = fits.open(weightimage)
            # just check center position?
            int_imx = np.array(imx,'i')
            int_imy = np.array(imy,'i')        
            centerpixvals = whdu[0].data[int_imy[keepflag],int_imx[keepflag]]
            # weight image will have value > 0 if there is data there
            weightflag = centerpixvals > 0
            keepflag[keepflag] = keepflag[keepflag] & weightflag
    return imx, imy, keepflag

def plot_vf_gals(imx,imy,keepflag,cat,ax,galsize=120):
    ''' plot galaxies on the coadd png images

    INPUT:
    * imx, imy : pixel locations of the galaxies
    * keepflag : flag that tells which galaxies to plot
    * cat : vf main catalog, 
    * ax : plot axis

    OPTIONAL INPUT:
    * galsize : number or array, the size of the rectangle size; default is 120 pixels

    RETURNS:
    * Null
    '''
    gindex=np.arange(len(imx))[keepflag]
    galsizes = cat['radius']/.4*2

    for j in gindex:
        galsize = galsizes[j]
        rect= plt.Rectangle((imx[j]-galsize/2.,imy[j]-galsize/2.), galsize, galsize,fill=False, color='b',lw=2)
        ax.add_artist(rect)
        s='{}\n vr={:.0f}'.format(cat['VFID'][j],cat['vr'][j])
        plt.text(imx[j],imy[j]+galsize/2.,s,fontsize=10,clip_on=True,horizontalalignment='center',verticalalignment='bottom',bbox=dict(facecolor='c',alpha=.3))
        plt.text(imx[j],imy[j]-galsize/2.,cat['NEDname'][j],fontsize=10,clip_on=True,horizontalalignment='center',verticalalignment='top',bbox=dict(facecolor='c',alpha=.3))
        #if cat['COflag'][j]:
        #    size=galsize+10
        #    rect= plt.Rectangle((imx[j]-size/2.,imy[j]-size/2.), size, size,fill=False, color='g',lw=1.5)
        #    ax.add_artist(rect)
        if cat['COflag'][j]:
            size=galsize
            #rect= plt.Circle((ran-size/2.,decn-size/2.), size,fill=False, color='g')
            rect= plt.Circle((imx[j],imy[j]), size/2,fill=False, color='g',lw=1.5)
            ax.add_artist(rect)
        if cat['A100flag'][j]:
            size=galsize+20
            rect= plt.Rectangle((imx[j]-size/2.,imy[j]-size/2.), size, size,fill=False, color='c',lw=1)
            ax.add_artist(rect)
        #if cat['HAobsflag'][j]:
        #    size=galsize+2*.005
        #    rect= plt.Rectangle((imx[i]-size/2.,imy[j]-size/2.), size, size,fill=False, color='r')
        #    ax.add_artist(rect)
    pass
