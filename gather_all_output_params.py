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
header=['XC','YC','MAG','RE','N','BA','PA','SKY','ERROR','CHI2NU']
cheader = ["C"+i for i in header]
colnames = ['VFID']+header+cheader

dtype=['S8',\
           float,float,float,float,float,float,float,float,float,float,\
           float,float,float,float,float,float,float,float,float,float]
outtab = Table(dtype=dtype,names=colnames)

topdir = os.getcwd()
# for each VFID, look in directory
for i in range(nvirgo):

    vfid = f"VFID{i:04d}"
    outtab.add_row(('',0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
    outtab['VFID'][i] = vfid
    
    # move to directory for this galaxy
    os.chdir(vfid)

        

    # output from galfit is NAME-g-out1.fits = output from no convolution
    infile1 = glob.glob(f"*{band}-out1.fits")
    if len(infile1) < 1:
        print("no output file for ",vfid)
    else:
        t = rg.parse_galfit_results(infile1[0])
        for j in range(len(t)):
            if j < 8:
                outtab[cheader[j]][i] = t[j][0]
            else:
                outtab[cheader[j]][i] = t[j]
                
            
    # output from galfit is NAME-g-out2.fits = output with convolution        
    infile1 = glob.glob(f"*{band}-out2.fits")
    if len(infile1) < 1:
        print("no output file for ",vfid)
        continue
    else:
        t = rg.parse_galfit_results(infile1[0])
        for j in range(len(t)):
            if j < 8:
                outtab[cheader[j]][i] = t[j][0]
            else:
                outtab[cheader[j]][i] = t[j]
    os.chdir(topdir)

# write output table

outfilename = f"vf_v2_galfit_{band}.fits"

outtab.write(outfilename,format='fits',overwrite=True)
