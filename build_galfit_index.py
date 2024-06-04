#!/usr/bin/env python

'''
GOAL:
* create index web page that points to html pages for all cutouts 

USAGE:
* run from html-dev/cutouts directory

NOTES:
* using John Moustakas's code as a reference (https://github.com/moustakas/legacyhalos/blob/main/py/legacyhalos/virgofilaments.py#L1131-L1202)

'''

import os
import numpy as np
import glob

from astropy.io import fits
homedir = os.getenv("HOME")

###########################################################
####  FUNCTIONS
###########################################################

    
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

def write_table(html,images,labels):    
    html.write('<table width="90%">\n')
    html.write('<tr>')
    for l in labels:
        html.write('<th>{}</th>'.format(l))
    html.write('</tr></p>\n')        
    html.write('<tr>')
    for i in images:
        html.write('<td><a href="{0}"><img src="{1}" alt="Missing file {0}" height="auto" width="100%"></a></td>'.format(i,i))
    html.write('</tr>\n')            
    html.write('</table>\n')

def write_text_table(html,labels,data):    
    html.write('<table width="90%">\n')
    html.write('<tr>')
    for l in labels:
        html.write('<th>{}</th>'.format(l))
    html.write('</tr></p>\n')        
    html.write('<tr>')
    for d in data:
        html.write('<td>{}</td>'.format(d))
    html.write('</tr>\n')            
    html.write('</table>\n')


def legacy_link(ra,dec):
    return "https://www.legacysurvey.org/viewer?ra={:.4f}&dec={:.4f}&layer=ls-dr9&zoom=13".format(ra,dec)    
###########################################################
####  CLASSES
###########################################################
    
class build_html_coadd():

    def __init__(self,gallist,outdir):
        ''' pass in instance of cutout_dir class and output directory '''


        outfile = os.path.join(outdir,'index.html')
        self.outdir = outdir

        #print('inside build html')
        #print('coutdir = ',coutdir)
        #print('outfile = ',outfile)        
        self.html = open(outfile,'w')
        self.galnames = gallist
        self.build_html()
        
    def build_html(self):
        self.write_header()
        self.write_coadd_list()
        #self.write_navigation_links()
        self.close_html()
    def write_header(self):
        self.html.write('<html><body>\n')
        self.html.write('<style type="text/css">\n')
        self.html.write('table, td, th {padding: 5px; text-align: center; border: 1px solid black}\n')
        self.html.write('</style>\n')
        self.html.write('<h1>VF Coadds</h1>\n')        
    def write_coadd_list(self):
        self.html.write('<table width="50%">\n')
        self.html.write('<tr>')
        colnames = ['Index','COADD']
        for i,l in enumerate(colnames):
            if i == 1:
                colspan=2
            else:
                colspan=1
            self.html.write('<th colspan="{}">{}</th>'.format(colspan,l))
        self.html.write('</tr></p>\n')            

        #vfindices = np.arange(len(vfmain))
        # write one row for each galaxy
        #print('galnames = ',self.galnames)
        galindex = 1
        for i,g in enumerate(self.galnames):
            print(g)
            jpg_path = os.path.join(self.outdir,'r-coadd.png')
            
            self.html.write('<tr>')
            self.html.write('<td>{}</td>'.format(galindex))            
            htmlpage = "{}/{}.html".format(g,g)
            self.html.write('<td><a href="{}">{}</td>'.format(htmlpage,g.replace('-noback-coadd','')))

            self.html.write('</tr>\n')
            galindex += 1
        self.html.write('</table>\n')
    
    
    def close_html(self):
        self.html.close()
# wrap
if __name__ == '__main__':
    # work through coadd directory
    #global vmain

    
    #outdir = homedir+'/research/Virgo-dev/html-dev/coadds/'
    outdir = '/data-pool/Halpha/html_dev/galfit/'    
    #outdir = '/home/rfinn/Virgo-dev/html-dev/coadds/'    

    # this should contain a list of all the galaxy folders
    flist1 = os.listdir(outdir)
    flist1.sort()
    #print(flist1)
    galnames=[]
    for i,subdir in enumerate(flist1): # loop through list
        

        #if os.path.isdir(subdir) & (subdir.startswith('pointing')) & (subdir.find('-') > -1):
        if (os.path.isdir(subdir)):# & (subdir.startswith('VF')):
            #print('adding ',subdir)
            galnames.append(subdir)
    #print('galnames = ',galnames)
    h = build_html_coadd(galnames,outdir)
