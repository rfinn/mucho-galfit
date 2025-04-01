'''
Brief code to merge North and South NED-LVS sample catalogs from JM

PROCEDURE:
*read in catalogs, stack and remove duplicates, save output to location specified in paramfile.txt
'''

def stack_cats(north,south,data_root_dit):
    north=Table.read(data_root_dir+north)
    south=Table.read(data_root_dir+south)
    stacked_catalogs = vstack([north,south])
    return stacked_catalogs

def unique_cats(stacked_cat,create_objid=True):
    unique_cat = unique(stacked_cat,keys='OBJNAME')
    if create_objid:
        unique_cat = create_OBJIDs(unique_cat)
    return unique_cat
    
def create_OBJIDs(unique_cat):
    formatted_strings = np.array([f"OBJID{num:05}" for num in np.arange(0,len(unique_cat))])
    unique_cat['OBJID'] = formatted_strings  #add the new column
    return unique_cat
    
def save_table(unique_cat,main_catalog_name,overwrite_flag=True):
    unique_cat.write(main_catalog_name,overwrite=overwrite_flag)
    
#GALFIT requires an ephot tab. in cases where ephot tab does not exist, create a mock version
#where all group flags are true.
def make_ephot_mockcat(main_catalog,ephot_name,objid_col,primary_group_col,group_mult_col,overwrite_flag=True):
    
    #read main catalog; extract objid, RA, DEC columns
    maintab = Table.read(main_catalog)
    objid = maintab[objid_col]
    RA = maintab['RA']
    DEC = maintab['DEC']
    
    #define placeholder group columns
    primary_group = np.ones(len(maintab),dtype=bool)
    group_mult = np.ones(len(maintab),dtype=int)
    
    #make ephot mock catalog
    ephot_tab = Table([objid, RA, DEC, primary_group, group_mult, group_name],
                      names=[objid_col, 'RA', 'DEC', primary_group_col, group_mult_col])
    
    #save to ephot_name path
    ephot_tab.write(ephot_name,overwrite=overwrite_flag)
    

if __name__ == '__main__':

    from astropy.table import Table, vstack, unique
    import numpy as np
    import os
    
    param_file = '/mnt/astrophysics/kconger_wisesize/github/wisesize/mucho-galfit/paramfile.txt'
        
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
    
    data_root_dir = param_dict['data_root_dir']
    main_catalog = param_dict['main_catalog']    #where result will be saved
    north_cat = param_dict['north_cat_name']
    south_cat = param_dict['south_cat_name']
    
    ephot_tab = param_dict['phot_catalog']
    
    objid_col = param_dict['objid_col']
    primary_group_col = param_dict['primary_group_col']
    group_mult_col = param_dict['group_mult_col']
    
    stacked_catalog = stack_cats(north_cat, south_cat, data_root_dir)
    unique_catalog = unique_cats(stacked_catalog)
    save_table(unique_catalog, main_catalog, overwrite_flag=True)   
    
    #if no ephot tab, make one (or rather, a "placeholder" for when running GALFIT)
    if not os.path.exists(ephot_tab):
        print(f'ephot catalog not found. writing placeholder to {ephot_tab}')
        make_ephot_mockcat(main_catalog, ephot_tab, objid_col, primary_group_col, group_mult_col)