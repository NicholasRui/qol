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
