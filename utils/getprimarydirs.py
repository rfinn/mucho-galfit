'''
AIM: generate primaryDirs.txt file for running scripts in parallel!
    
    * the primary galaxies are those which will have directories.
'''
import os
import numpy as np

def find_primaries(maintab, primary_group_col):
    
    #trim to only include primary galaxies; if no such flag exists, assume all galaxies are primary.
    try:
        primary_flag = maintab[primary_group_col]
    except:
        primary_flag = np.ones(len(maintab),dtype=bool)   #all true

    maintab = maintab[primary_flag]
    
    primary_objids = maintab[objid_col]
    
    return primary_objids

def write_primarydirs(primary_galaxies, destination):
    with open(destination, 'w') as f:
        for galaxy in primary_galaxies:
            f.write(str(galaxy) + '\n')