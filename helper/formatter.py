# Handles conversion between formats

def to_fortran(value):
    """
    Convert from Python to Fortran format
    """
    match value:
        case int():
            return f'{value:d}'
        case float():
            return format(value, '.15e').replace('e', 'D')
        case bool():
            return '.true.' if value else '.false.'
        case str():
            return value
