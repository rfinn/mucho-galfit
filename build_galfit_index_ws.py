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
        self.html.write('<h1>WISESize GALFIT Analysis</h1>\n')        
    
    def write_coadd_list(self):
        self.html.write('<table width="50%">\n')
        self.html.write('<tr>')
        colnames = ['Index','Legacy Image','Galaxy']
        for i,l in enumerate(colnames):
            if i == 1:
                colspan=2
            else:
                colspan=1
            self.html.write('<th colspan="{}">{}</th>'.format(colspan,l))
        self.html.write('</tr></p>\n')            

        galindex = 1
        for i,g in enumerate(self.galnames):
            print(g)
            self.html.write('<tr>')
            self.html.write('<td>{}</td>'.format(galindex))  
            jpg_path = os.path.join(f'{g}/{g}-im-LS.jpg')
            
            self.html.write(f'<td><a href="{jpg_path}"><img src="{jpg_path}" alt="Missing file {jpg_path}" height="auto" width="50%"></a></td>')
            
            htmlpage = "{}/{}.html".format(g,g)
            
            #do not think the 'replace' is necessary...but I will keep it in case.
            self.html.write('<td><a href="{}">{}</td>'.format(htmlpage,g.replace('-noback-coadd','')))

            self.html.write('</tr>\n')
            galindex += 1
        self.html.write('</table>\n')
    
    
    def close_html(self):
        self.html.close()

# wrap
if __name__ == '__main__':
    # work through coadd directory
    
    param_file = homedir+'/github/mucho-galfit/paramfile.txt'
        
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
    
    outdir = param_dict['html_outdir']
    
    # this should contain a list of all the galaxy folders
    temp = os.listdir(outdir)
    temp.sort()
    flist1 = []
    for t in temp:
        if 'OBJID' in t:
            flist1.append(t)
    #print(flist1)
    galnames=[]
    for i,subdir in enumerate(flist1): # loop through list
        

        if (os.path.isdir(subdir)):
            #print('adding ',subdir)
            galnames.append(subdir)
    #print('galnames = ',galnames)
    h = build_html_coadd(galnames,outdir)
