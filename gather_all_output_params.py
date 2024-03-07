#!/usr/bin/env  python
"""
GOAL: 

to gather all results from galfit fitting

USAGE: 

from /mnt/astrophysics/rfinn/muchogalfit-output/

python3 ~/github/mucho-galfit/gather_output_params_all.py W3

where W3 is the bandpass and can  be [g,r,z,W1,W2,W3,W4]

OUTPUT:
creates a table in the current directory called

vf_v2_galfit_{band}.fits

or should we combine all wavelengths into one table?

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


# read in VFID catalog
#tabledir = "/mnt/astrophysics/rfinn/catalogs/Virgo/v2/"

#vmain = Table.read(tabledir+"vf_v2_main.fits")

# don't actually need to read in the main table just to get the VFIDs
# they just run from 1 - 6780
#
# or from 0 - 6779

nvirgo = 6780
# set up table to store galfit output
header=['XC','YC','MAG','RE','N','AR','PA','SKY','CHI2NU']
header_err = [f"{h}_ERR" for h in header[:-1]]
# TODO - make an empty table
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
col1 = ['VFID']


colnames = col1 + hheader_with_err + cheader_with_err
colnames.append('Numerical_Error')
colnames.append('CNumerical_Error')


dtype=['S8',\
           float,float,float,float,float,float,float,float,float,\
           float,float,float,float,float,float,float,float,\
           float,float,float,float,float,float,float,float,float,\
           float,float,float,float,float,float,float,float,\
           'bool','bool']
print('length of dtype = ',len(dtype))

print('length of colnames = ',len(colnames))
outtab = Table(np.zeros((nvirgo,len(dtype))),dtype=dtype,names=colnames)




for i in range(nvirgo):
    vfid = f"VFID{i:04d}"
    outtab['VFID'][i] = vfid
    #outtab['VFID'][i] = vfid

# create a dictionary to link the VFID and row in table?
topdir = os.getcwd()
os.chdir(topdir)
# TODO get list of directories
dirlist = glob.glob('VFID????')
dirlist.sort()
#print(dirlist)
# go in each directory
for d in dirlist:
    if os.path.isdir(d):
        print(d)
    os.chdir(d)

    for k in range(2):
        if k == 0:
            #output from galfit is NAME-g-out1.fits = output from no convolution
            infile1 = glob.glob(f"*{band}-out1.fits")
            prefix = ''
        else:
            infile1 = glob.glob(f"*{band}-out2.fits")
            prefix = 'C'
        if len(infile1) < 1:
            print("no out1 file for ",vfid)
            continue

        # TODO : figure out how to extract results for groups, and how to match with the corresponding VFID!!!
        if 'GROUP' in infile1[0]:
            # read in galsFOV-{band}.txt to get VFIDs of sources in field
            gfile = f"galsFOV-{band}.txt"
            if os.path.exits(gfile):
                vfids = []
                xgal = []
                ygal = []

                vgals = open(gfile,'r')
                for line in vgals:
                    t = line.strip().split(',')
                    vfids.append(t[0])
                    xgal.append(t[1])
                    ygal.append(t[2])
        else:
            xgal = [0]
            vfids = [d]
        hdu = fits.open(infile1[0])
 
        imheader = hdu[2].header
        #print(imheader)
        hdu.close()
        print(infile1[0])
        for i in range(len(xgal)):
            table_index = int(vfids[i].replace('VFID',''))
            for h in header[:-2]:
                hkey = f"{i+1}_{h}"

                t = imheader[hkey]
                print(hkey,t)
                if '*' in t:
                    #print(outtab[prefix+'Numerical_Error'][table_index])
                    outtab[prefix+'Numerical_Error'][table_index] = True
                    t.replace('*','')
                else:
                    #print(table_index)
                    #print(prefix+'Numerical_Error')
                    #print(outtab[prefix+'Numerical_Error'][table_index])
                    outtab[prefix+'Numerical_Error'][table_index] = False
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
            print("HEYYYY: nsky = ",nsky)
            t = imheader[prefix+f"{nsky}_SKY"]
            a,b = t.split(' +/- ')
            outtab[prefix+'SKY'][table_index] = float(a)
            outtab[prefix+'SKY_ERR'][table_index] = float(b)

os.chdir(topdir)
# write output table
os.chdir(topdir)
outfilename = f"vf_v2_galfit_{band}.fits"

outtab.write(outfilename,format='fits',overwrite=True)
