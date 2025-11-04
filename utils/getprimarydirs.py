'''
AIM: generate primaryDirs.txt file for running scripts in parallel!
    
    * the primary galaxies are those which will have directories.
'''
import os
from astropy.table import Table
import numpy as np

def find_primaries(maintab, primary_group_col, objid_col):
    
    #trim to only include primary galaxies; if no such flag exists, assume all galaxies are primary.
    try:
        primary_flag = maintab[primary_group_col]
    except:
        primary_flag = np.ones(len(maintab),dtype=bool)   #all true

    maintab = maintab[primary_flag]
    
    primary_objids = maintab[objid_col]
    
    return primary_objids

def write_primarydirs(primary_galaxies, destination):
    """Write one OBJID per line to the destination file."""
    os.makedirs(os.path.dirname(destination), exist_ok=True)
    with open(destination, 'w') as f:
        for galaxy in primary_galaxies:
            f.write(str(galaxy) + '\n')
            
if __name__ == '__main__':
    
    ######################
    ### Parameter File ###
    ######################
    
    param_file = '/mnt/astrophysics/wisesize/github/mucho-galfit/paramfile.txt'
        
    #create dictionary with keyword and values from param textfile
    param_dict = {}
    with open(param_file) as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            try:
                key, val = line.split()[:2]
                param_dict[key] = val
            except ValueError:
                continue
    
    maintab = Table.read(param_dict['main_catalog'])
    primary_group_col = param_dict['primary_group_col']
    objid_col = param_dict['objid_col']
    
    destination = os.path.join(param_dict['main_dir'], param_dict['path_to_images'], 'PrimaryDirs.txt')
    
    primary_galaxies = find_primaries(maintab, primary_group_col, objid_col)
    write_primarydirs(primary_galaxies, destination)
    print(f'Successfully written {destination}!')