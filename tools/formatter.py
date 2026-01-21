# Handles conversion between formats
import warnings

def to_fortran(value):
    """
    Convert from Python to Fortran format
    """
    match value:
        case bool(): # important: need to check before int, since bool is a subclass of int
            return '.true.' if value else '.false.'
        case int():
            return f'{value:d}'
        case float():
            return format(value, '.15e').replace('e', 'D')
        case str():
            return f"'{value}'"
        case None:
            warnings.warn('to_fortran: tried to convert None, returning empty string')
            return ''
        case _:
            raise ValueError(f'to_fortran: invalid type {type(value)}')

def mesa_num_to_float(mesa_num):
    """
    Read in MESA number and return float
    """
    # replace D and d with e
    mesa_num = mesa_num.lower().replace('d', 'e')

    # there is a bug where exponents of -100 or lower overwrite the D
    # since this is vanishingly small, just set it to 0
    if ('-' in mesa_num[1:]) and ('e' not in mesa_num[1:]):
        return 0
    
    return float(mesa_num)
