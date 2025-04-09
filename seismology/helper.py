"""
Helper functions for Seismology class and associated attributes
"""
import numpy as np

# TODO: add methods to do additional conversions

def get_ω(ω=None, ν=None, ω_uHz=None, ν_uHz=None, P=None):
    """
    Method which allows other methods to take a frequency in multiple formats,
        and asserts that the user has only specified one of them
    """
    args = [ω, ν, ω_uHz, ν_uHz, P]
    provided_args = sum([1 if arg is not None else 0 for arg in args])
    assert provided_args == 1

    ω = 2 * np.pi * ν if ω is None and ν is not None else ω
    ω = 1e-6 * ω_uHz if ω is None and ω_uHz is not None else ω
    ω = 2 * np.pi * 1e-6 * ν_uHz if ω is None and ν_uHz is not None else ω
    ω = 2 * np.pi / P if ω is None and P is not None else ω

    return ω

def get_magnetic_acrit(l, m=None, acrit_ref='F+15'):
    """
    acrit_ref: sets uncertain dimensionless parameter for magnetic g-mode suppression,
              with 'a' defined in RF+23
       'F+15': Fuller+2015 normalization
    """
    match acrit_ref: # critical a_c parameter, which is fixed by certain papers but not known
        case 'F+15':
            acrit = 0.5 / np.sqrt(l * (l + 1))
        case 'L+17':
            raise NotImplementedError('acrit_ref for L+17 option not implemented yet') # TODO
        case 'RF+23':
            raise NotImplementedError('acrit_ref for RF+23 option not implemented yet') # TODO
        case _:
            raise ValueError(f'ac_ref is not valid: {acrit_ref}')

    return acrit


