#!/usr/bin/env python
'''

call as
setup_galfit.py W3

where W3 can be any wavelength band [g,r,z,W1 through W4]

'''
import os
import sys

wavelength = sys.argv[1]
from astropy.io import fits

# read in the main table

# read in table that is line-matched to v2 main table that says galaxy is in sample


# cut main catalog based on sample flag
finalcat = v.main[sampleflag]

# for each galaxy, create a directory and write sourcelist
for i in range(len(finalcat)):
    galpath = '/mnt/astrophysics/muchogalfit-output/VFID'+finalcat['VFID'][i]
    os.mkdir(galpath)

    sourcefile = galpath+'/VFID{}sourcelist'.format(finalcat['VFID'][i])
    sourcelist = open(sourcefile,'w')
    # write out one line with VFID, objname, RA, DEC, wavelength

    sourcelist.close()
