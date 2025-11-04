#!/usr/bin/env python3
"""
GOAL:
* Set up directories for running GALFIT on WISESize sample.
* This script runs `setup_one_galaxy` from `setup_galfit_ws.py` in PARALLEL via GNU parallel.
"""

import os
import sys
from astropy.table import Table

# --- Import utility functions ---
sys.path.insert(0, '/mnt/astrophysics/wisesize/github/mucho-galfit/utils')
from setup_galfit_ws import setup_one_galaxy
from getprimarydirs import write_primarydirs, find_primaries

sys.path.insert(0, '/mnt/astrophysics/wisesize/github/mucho-galfit')
from merge_ns_catalogs import create_OBJIDs

# --- Read parameter file ---
param_file = '/mnt/astrophysics/wisesize/github/mucho-galfit/paramfile.txt'
param_dict = {}

print(os.getcwd())
print('#'*20)


with open(param_file) as f:
    for line in f:
        
        #allow comments and skip empty lines
        if not line.strip() or line.startswith('#'):
            continue
        
        try:
            key, val = line.split()[:2]
            param_dict[key] = val
        
        except ValueError:
            continue  #skip malformed lines

# --- Load main catalog ---
main_catalog_path = param_dict['main_catalog']
maintab = Table.read(main_catalog_path)

# --- Add OBJIDs if missing ---
if param_dict['objid_col'] not in maintab.colnames:
    maintab = create_OBJIDs(maintab)
    maintab.write(main_catalog_path, overwrite=True)

# --- Write primary directories file if missing ---
primary_dirs_file = os.path.join(param_dict['main_dir'], param_dict['path_to_images'], 'PrimaryDirs.txt')

if not os.path.exists(primary_dirs_file):
    primary_galaxies = find_primaries(maintab, param_dict['primary_group_col'])
    write_primarydirs(primary_galaxies, primary_dirs_file)

# --- Run GALFIT setup for a single galaxy ---
if len(sys.argv) < 2:
    sys.exit("Usage: run1directory.py <OBJID>")

objid = sys.argv[1]
setup_one_galaxy(objid, maintab, param_dict)
