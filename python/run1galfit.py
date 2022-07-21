#!/usr/bin/env python

'''
GOAL: this program will run galfit on one galaxy

USAGE:
This should be run from a directory that contains
- galfit input file
- image
- noise image

To run:
python ~/github/virgowise/run1galfit.py galname

This assumes that the galfit input file starts with galname and ends with inpu

Each galaxy goes in its own folder, and the parallel program 
will move to each folder and run this script

'''
import os
import sys
import glob

# galaxy to work with
galid = sys.argv[1]

HOME = os.getenv("HOME")
# directory where galaxy folders are
data_dir = os.getcwd()+'/'

# move into directory for this galaxy
# again, this assumes that galid is the same as the directory name
os.chdir(data_dir+'/'+galid)


# code to run galfit

print('running galfit')
os.system(f"galfit {galfit_input}")
#os.system('rm *.lbr')
    
os.chdir(data_dir)
