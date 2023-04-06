#!/usr/bin/env python

"""
GOAL: 
* figure out which galaxies are in John's groups
* get list of groups

USAGE:
* run from muchogalfit-output directory



"""
import numpy as np
import os
from astropy.table import Table
import glob


topdir = '/mnt/astrophysics/rfinn/muchogalfit-output/'
try:
    os.chdir(topdir)
except FileNotFoundError: # assuming that we are running on virgo vms
    topdir = '/mnt/astrophysics/muchogalfit-output/'
    os.chdir(topdir)

# look through output directories
dirlist = []
t = open('Dirs.txt','r')
for line in t:
    dirlist.append(line.rstrip())
t.close()

# set up arrays to track group galaxies and central galaxies
groupflag = np.zeros(len(dirlist),'bool')
primaryflag = np.zeros(len(dirlist),'bool')

groupdirs_list = []

for i,dir in enumerate(dirlist):
    print(i,os.getcwd(),dir)
    # search for galaxies with no fits files
    os.chdir(dir)
    fitsfiles = glob.glob("*.fits")

    if len(fitsfiles) > 0:
        os.chdir(topdir)
        continue
    elif len(fitsfiles) == 0:
        groupflag[i] = True

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

        # check to see if this is the primary galaxy

        group_dir = f'/mnt/astrophysics/virgofilaments-data/{ra}/{objname}_GROUP'
        if os.path.exists(group_dir):
            primaryflag[i] = True
            groupdirs_list.append(group_dir)
        else:
            continue

        os.chdir(topdir)
# write out list of groupdirs
gtab = Table([groupdirs_list])
gtab.write('groupDirs.txt',format='ascii')

# save information about group
tab = Table([dirlist,groupflag,primaryflag],names=['VFID','groupflag','primary'])
tab.write('groupinfo.fits',format='fits')
