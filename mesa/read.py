import numpy as np
from astropy.table import Table, Column, unique

from qol.mesa.table.MesaTable import MesaTable

# Attributes to look for and store
model_attr_data_raw = \
    {'version_number': 'str',
             'M/Msun': 'float',
       'model_number': 'int',
           'star_age': 'float',
          'initial_z': 'float',
           'n_shells': 'int',
           'net_name': 'str',
            'species': 'int',
               'Teff': 'float',
     'power_nuc_burn': 'float',
       'power_h_burn': 'float',
      'power_he_burn': 'float',
       'power_z_burn': 'float',
        'power_photo': 'float',
       'total_energy': 'float',
        'cumulative_energy_error': 'float',
  'cumulative_error/total_energy': 'float',
        'num_retries': 'int',
  'previous n_shells': 'int',
          'previous mass (grams)': 'float',
 'timestep (seconds)': 'float',
  'dt_next (seconds)': 'float',
           'R_center': 'float',
             'xmstar': 'float',
               'time': 'float'}

model_attr_data_underscores = {key.replace(' ', '_'): value for key, value in zip(model_attr_data_raw.keys(), model_attr_data_raw.values())}


def read_data(fname, remove_duplicates=True):
    """
    Function which takes in .data files written by MESA/GYRE and return
    an astropy table.
    
    ARGUMENTS
    ---------
    fname: str
      filename (e.g., 'summary.txt')
    
    header_start: int (default: 4)
      header start line, default 4 in current MESA/GYRE

    remove_duplicates: bool (default: True)
      if True, removes duplicate rows. this can happen in MESA runs if runs are restarted from photo.
    
    RETURNS
    -------
    tab: astropy.table.Table
      summary table
    """
    f = open(fname, 'r')
    text = f.read().split('\n')
    f.close()

    # Get data table
    colname_row = text[5]
    data_rows = text[6:-1]

    # Remove all duplicate spaces
    colname_row = ' '.join(colname_row.split())
    data_rows = [' '.join(data_rows[ii].split()) for ii in range(len(data_rows))]

    colnames = colname_row.split()

    # Build each column, assemble table, save
    col_list = []
    existing_colnames = []

    for ii in range(len(colnames)):
        coldata = [float(data_rows[jj].split()[ii]) for jj in range(len(data_rows))]

        colname = colnames[ii]
        append_num = 0
        while colname in existing_colnames:
            colname = f'{colnames[ii]}{append_num}'
            append_num += 1

        existing_colnames.append(colname)
        
        col = Column(coldata, name=colname)
        col_list.append(col)

    tab = Table(col_list)
    if remove_duplicates:
        tab = unique(tab)

    # Get attributes and classify table
    attr_names = text[1].split()
    attr_vals = text[2].split()
    attr = {attr_names[ii]: attr_val.replace('"', '') if '"' in attr_val \
                       else float(attr_val.lower().replace('d', 'e'))    \
                       for ii, attr_val in enumerate(attr_vals)}

    tabtype = 'profile' if 'model_number' in attr_names else 'history'

    mesatab = MesaTable(tab, attr, tabtype=tabtype)

    return mesatab

def read_index(fname):
    """
    Reads profiles.index
    """
    header_start = 1
    
    f = open(fname, 'r')
    text = f.read().split('\n')
    f.close()
    
    data_rows = text[header_start:-1]
    
    # Remove all duplicate spaces
    data_rows = [' '.join(data_rows[ii].split()) for ii in range(len(data_rows))]
    
    colnames = ['model_number', 'priority', 'profile_number']
    
    # Build each column, assemble table, save
    col_list = []
    
    for ii in range(len(colnames)):
        coldata = [int(data_rows[jj].split()[ii]) for jj in range(len(data_rows))]
        
        col = Column(coldata, name=colnames[ii])
        col_list.append(col)
    
    tab = Table(col_list)

    mesatab = MesaTable(tab, attr={}, tabtype='index')
    
    return mesatab

def read_mod(fname):
    """
    Reads model files
    """
    with open(fname, 'r') as f:
        lines = np.array(f.readlines())

    # Loop through lines, extract needed info and then remove them
    # if not explicitly part of model data table
    attr = {}

    have_colnames = False
    table_data = []

    for line in lines:
        # remove comments
        line = line.split('!')[0]

        items = line.split()
        
        # remove empty lines
        if (len(items) == 0) or ('previous model' in line):
            continue

        # log setting for column data if not logged yet
        if ('model_setting' not in attr.keys()) and ('--' in items):
            model_setting = int(items[0])
            attr['model_setting'] = model_setting
            continue
        
        # log other attributes
        # to avoid issues, first replace all spaces in attr names with underscores
        # and then re-split the line to get items again
        for key in model_attr_data_raw.keys():
            line = line.replace(key, key.replace(' ', '_'))
        items = line.split()

        if np.in1d(items, list(model_attr_data_underscores.keys())).any():
            # loop through items
            # if one matches the name of an attr,
            # convert the next one into the desired datatype and
            # log it as an attribute
            for ii, item in enumerate(items):
                if item in model_attr_data_underscores.keys():
                    attr_name = item
                    attr_entry = items[ii+1]
                    attr_type = model_attr_data_underscores[attr_name]

                    match attr_type:
                        case 'str':
                            attr_entry = str(attr_entry.replace('"', '').replace("'", ''))
                            
                        case 'float':
                            attr_entry = attr_entry.lower()
                            attr_entry = attr_entry.replace('d', 'e')

                            # kludge: MESA sometimes overwrites the d if the exponential is too long
                            # just replace with 0, since this usually means a very large negative exponent
                            if (attr_entry.count('-') > 1):
                                attr_entry = 0
                            
                            attr_entry = float(attr_entry)
                            
                        case 'int':
                            attr_entry = int(attr_entry)
                    
                    attr[attr_name] = attr_entry

            continue

        # assume first unaltered row stores the column names
        if not have_colnames:
            colnames = items
            have_colnames = True
            continue

        # log data row
        # exclude first item, since it labels the row number
        entries = [float(item.lower().replace('d', 'e')) for item in items[1:]]
        table_data.append(entries)
    
    table_data = np.array(table_data)
    cols = [Column(table_data[:,ii], name=colname) for ii, colname in enumerate(colnames)]
    tab = Table(cols)
    
    mesatab = MesaTable(tab, attr, tabtype='model')

    return mesatab

