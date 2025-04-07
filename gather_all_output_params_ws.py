#!/usr/bin/env  python
"""
GOAL: 

to gather all results from galfit fitting

USAGE: 

from /mnt/astrophysics/wisesize/mg_output_wisesize/

python3 ~/github/mucho-galfit/gather_all_output_params.py W3

where W3 is the bandpass and can  be [g,r,z,W1,W2,W3,W4]

OUTPUT:
creates a table in the current directory [/mnt/astrophysics/wisesize/mg_output_wisesize/] called

wisesize_galfit_{band}.fits

"""

import os
import sys
from astropy.table import Table
from astropy.io import fits
import glob
import numpy as np

HOME = os.getenv("HOME")
sys.path.append(HOME+'/github/halphagui/')
import rungalfit as rg


band = sys.argv[1]

try:
    flag = int(sys.argv[2])
    fixBA = '-fixBA'
    
except:
    fixBA = ''

print("fixBA = ",fixBA)

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

nobj = len(Table.read(param_dict['main_catalog']))

# set up table to store galfit output
header=['XC','YC','MAG','RE','N','AR','PA','SKY','CHI2NU']
header_err = [f"{h}_ERR" for h in header[:-1]]
# for convolution - second pass of galfit
cheader = ["C"+i for i in header]
cheader_err = [f"{h}_ERR" for h in cheader[:-1]]

hheader_with_err = []
cheader_with_err = []
for i in range(len(header)):
    hheader_with_err.append(header[i])
    if i < len(header)-1:
        hheader_with_err.append(header_err[i])    
for i in range(len(cheader)):
    cheader_with_err.append(cheader[i])
    if i < len(cheader) -1:
        cheader_with_err.append(cheader_err[i])
col1 = ['OBJID']


colnames = col1 + hheader_with_err + cheader_with_err
colnames.append('Numerical_Error')
colnames.append('CNumerical_Error')


dtype=['S8',\
           float,float,float,float,float,float,float,float,float,\
           float,float,float,float,float,float,float,float,\
           float,float,float,float,float,float,float,float,float,\
           float,float,float,float,float,float,float,float,\
           'bool','bool']
#print('length of dtype = ',len(dtype))

#print('length of colnames = ',len(colnames))
outtab = Table(np.zeros((nobj,len(dtype))),dtype=dtype,names=colnames)

for i in range(nobj):
    objid = f"OBJID{i:05d}"
    outtab['OBJID'][i] = objid

# create a dictionary to link the OBJID and row in table?
topdir = os.getcwd()
os.chdir(topdir)

dirlist = glob.glob('OBJID?????')
dirlist.sort()  #numerical order :-)

for d in dirlist:

    os.chdir(d)

    gfile = f"galsFOV-{band}.txt"
    if os.path.exists(gfile):
        objids = []
        xgal = []
        ygal = []

        gals = open(gfile,'r')
        for line in gals:
            t = line.strip().split(',')
            objids.append(t[0])
            xgal.append(t[1])
            ygal.append(t[2])

    #print("number of galaxies = ",len(xgal))

    #for no convolution and convolution cases...
    for k in range(2):
        if k == 0:
            #output from galfit is NAME-g-out1.fits = output from no convolution
            infile1 = glob.glob(f"*{band}{fixBA}-out1.fits")
            prefix = ''
        else:
            infile1 = glob.glob(f"*{band}{fixBA}-out2.fits")
            prefix = 'C'
        if len(infile1) < 1:
            print(f"no out{k+1} file for ",d)
            continue

        # read in image header
        hdu = fits.open(infile1[0])
        # extension 2 has the model info
        imheader = hdu[2].header
        #print(imheader)
        hdu.close()
        
        #print(infile1[0])
        for i in range(len(xgal)):
            table_index = int(objids[i].replace('OBJID',''))
            for h in header[:-2]:
                hkey = f"{i+1}_{h}"

                t = imheader[hkey]
                #print(hkey,t)
                if t.find('*') > -1:
                    #print("numerical error...")
                    #print(outtab[prefix+'Numerical_Error'][table_index])
                    outtab[prefix+'Numerical_Error'][table_index] = True
                    t=t.replace('*','')
                    #print("t w/out error = ",t)

                    
                if t.find('[') > -1:
                    t=t.replace('[','')
                    t=t.replace(']','')
                    t=t.split('+/-')
                    a,b=(float(t[0]),0.)# fit and error
                else:
                    a,b = t.split('+/-')
                outtab[prefix+h][table_index] = float(a)
                outtab[prefix+h+"_ERR"][table_index] = float(b)
            nsky = len(xgal) + 1
            t = imheader[f"{nsky}_SKY"]
            a,b = t.split(' +/- ')
            outtab[prefix+'SKY'][table_index] = float(a)
            outtab[prefix+'SKY_ERR'][table_index] = float(b)

    os.chdir(topdir)

# write output table
os.chdir(topdir)
outfilename = f"wisesize_galfit_{band}{fixBA}.fits"

outtab.write(outfilename,format='fits',overwrite=True)
