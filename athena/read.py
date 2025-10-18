import numpy as np
from astropy.table import Table, Column, unique

from qol.athena.table.AthenaTable import AthenaTable


def read_tab(fname):
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
    colname_row = text[1].replace('# ', '')
    data_rows = text[2:-1]

    # Remove all duplicate and leading/trailing spaces
    attr_row = ' '.join(attr_row.split()).strip()
    colname_row = ' '.join(colname_row.split()).strip()
    data_rows = [' '.join(data_rows[ii].split()).strip() for ii in range(len(data_rows))]

    # Pull out attributes
    time = float(attr_row.split('time=')[1].split(' ')[0])
    cycle = int(attr_row.split('cycle=')[1].split(' ')[0])
    variables = attr_row.split('variables=')[1]

    colnames = colname_row.split()

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

