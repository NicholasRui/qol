from astropy.table import Table, Column, unique

import warnings

class AthenaTable(Table):
    """
    Stores information about Athena++ output files. Inherits from astropy.table.Table.

    __init__ args not shared by astropy.table.Table:
    file_type: output data type, e.g., 'tab'
    time,cycle,variables: variables stored in the header of the 'tab' file

    Note: Right now, can ONLY accommodate 'tab' data_type for one-dimensional data
    """
    def __init__(self, data=None, file_type='tab', time=None, cycle=None, variables=None,
                 masked=False, names=None, dtype=None, meta=None, copy=True, rows=None, copy_indices=True, units=None, descriptions=None, **kwargs):
        super().__init__(data=data, masked=masked, names=names, dtype=dtype, meta=meta, copy=copy, rows=rows, copy_indices=copy_indices, units=units, descriptions=descriptions, **kwargs)

        self.file_type = file_type
        self.time = time
        self.cycle = cycle
        self.variables = variables # either 'prim' or 'cons'
