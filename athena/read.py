import numpy as np
from astropy.table import Table, Column, unique

from qol.athena.table.AthenaTable import AthenaTable

import pyvista as pv

def read_tab(fname, 
             colnames=None, # sometimes need to specify this manually, since header is not accurate...
            ):
    """
    Takes in .tab files written by Athena++ and returns an AthenaTable.
    
    ARGUMENTS
    ---------
    fname: str
      filename (e.g., 'summary.txt')

    RETURNS
    -------
    athenatab: AthenaTable
    """
    with open(fname, 'r') as f:
        text = f.read().split('\n')

    # Get data table
    attr_row = text[0].replace('# ', '')
    data_rows = text[2:-1]

    # Remove all duplicate and leading/trailing spaces
    attr_row = ' '.join(attr_row.split()).strip()
    data_rows = [' '.join(data_rows[ii].split()).strip() for ii in range(len(data_rows))]

    # Pull out attributes
    time = float(attr_row.split('time=')[1].split(' ')[0])
    cycle = int(attr_row.split('cycle=')[1].split(' ')[0])
    variables = attr_row.split('variables=')[1]

    # Get column names
    if colnames is None:
        colname_row = text[1].replace('# ', '')
        colname_row = ' '.join(colname_row.split()).strip()
        colnames = colname_row.split()
    assert len(colnames) == len(data_rows[0].split()) # ensure number of colnames matches columns

    # Build each column, assemble table, save
    col_list = []
    existing_colnames = []

    for ii in range(len(colnames)):
        # make sure to store 'i' column as int
        if ii == 0:
            coldata = [int(data_rows[jj].split()[ii]) for jj in range(len(data_rows))]
        else:
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

    athenatab = AthenaTable(tab, file_type='tab', time=time, cycle=cycle, variables=variables)

    return athenatab

def read_vtk(fname):
    """
    Read Athena++ output as VTK table
    """
    # Read file and access colnames
    mesh = pv.read(fname)
    colnames = mesh.cell_data.keys()

    # Convert to AthenaTable format
    col_list = []
    existing_colnames = []

    # add column for x1v, x2v, x3v
    cell_centers = mesh.cell_centers().points
    Ndim = cell_centers.shape[1]

    col_list.append(Column(cell_centers[:,0], name='x1v'))
    if Ndim >= 2:
        col_list.append(Column(cell_centers[:,1], name='x2v'))
    if Ndim >= 3:
        col_list.append(Column(cell_centers[:,2], name='x3v'))

    # read other columns
    for ii, colname in enumerate(colnames):
        coldata = mesh.cell_data[colname]

        colname = colnames[ii]
        append_num = 0
        while colname in existing_colnames:
            colname = f'{colnames[ii]}{append_num}'
            append_num += 1

        existing_colnames.append(colname)
        
        col = Column(coldata, name=colname)
        col_list.append(col)

    tab = Table(col_list)

    athenatab = AthenaTable(tab, file_type='vtk')

    return athenatab
