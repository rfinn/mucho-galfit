#!/usr/bin/env python
'''
GOAL:
* set up directories for running galfit on the Virgo sample

PROCEDURE:
* the program will create one directory per galaxy
* the program creates an input file in each directory
* the input file contains the VFID, HL objname, RA, DEC, (and bandpass, which isn't used anymore
* the input file is then read by run1galfit.py

USAGE:

* On the virgo vms, move to /mnt/astrophysics/muchogalfit-output
 
* Call as:

python ~/github/muchogalfit/setup_galfit.py W3

where W3 can be any wavelength band [g,r,z,W1 through W4]

-will create repositories for sample galaxies, each containing a sourcelist with which to run galfit


NOTES: 
* 2023-06-20 : RF updating the program to run on the full virgo sample.  The program 
was originally created to run on Kim's sample. 


'''
import os
import sys

wavelength = sys.argv[1]

from astropy.table import Table

# filepath which contains VF tables 
filepath = '/mnt/astrophysics/muchogalfit-input-cats/'


# read in the main table
main_name = 'vf_v2_main.fits'
main_table = Table.read(filepath+main_name)


# RF - editing this so that we will run galfit on the entire VF sample
# basically resetting finalcat to the main table
finalcat = main_table


# read in table that is line-matched to v2 main table that says galaxy is in sample
#line_matched_name = 'v2_snrcoadd.fits'
#matched_table = Table.read(filepath+line_matched_name)


# cut main catalog based on sample flag
#sampleflag = matched_table['sgacut_flag']
#finalcat = matched_table[sampleflag]


# for each galaxy, create a directory and write sourcelist
for i in range(len(finalcat)):
    
    galpath = '/mnt/astrophysics/muchogalfit-output/'+finalcat['VFID'][i]

    # make directory if it doesn't already exist
    if not os.path.exists(galpath):
        os.mkdir(galpath)

    sourcefile = galpath+'/{}sourcelist'.format(finalcat['VFID'][i])
    sourcelist = open(sourcefile,'w')
    
    # write out one line with VFID, objname, RA, DEC, wavelength
    output_string = finalcat['VFID'][i] + ' ' + finalcat['objname'][i] + ' ' + str(finalcat['RA'][i]) + ' ' + str(finalcat['DEC'][i]) + ' ' + str(wavelength) + ' \n'.format()
    sourcelist.write(output_string)
    sourcelist.close()
