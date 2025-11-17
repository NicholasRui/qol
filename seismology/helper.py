"""
Helper functions for Seismology class and associated attributes

Helpers are outside functions that don't depend on the state of the Seismology class.
"""

import numpy as np
from math import gcd

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

def get_δω(δω=None, δν=None, δω_uHz=None, δν_uHz=None, δP=None,
           P=None):
    """
    Method which allows other methods to take a frequency uncertainty in multiple formats,
        and asserts that the user has only specified one of them

        (frequency is also required to convert uncertainty, if period)

    Warning: This could fail in the strange case where you specify δP but not P in some function which calls this.
    """
    args = [δω, δν, δω_uHz, δν_uHz, δP]
    provided_args = sum([1 if arg is not None else 0 for arg in args])
    assert provided_args == 1

    if δP is None: # if uncertainty was in frequency, can just use get_ω
        δω = get_ω(ω=δω, ν=δν, ω_uHz=δω_uHz, ν_uHz=δν_uHz)
    else:
        assert P is not None

        δω = 2 * np.pi * δP / P**2
    
    return δω

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
            raise ValueError(f'acrit_ref is not valid: {acrit_ref}')

    return acrit

def get_magnetic_aasym(l, m1, m2, m3, asym_ref='RFH+25'):
    """
    Gets "good guess" for uncertain magnetic asymmetry, called 'a' in Li+2022 and other works for dipole modes.
    Definition for a_asym for arbitrary (l, m) triplets is given in RFH+25.

    asym_ref: sets "good guess" for uncertain dimensionless asymmetry parameter,
        'RFH+25': Rui, Fuller, Hermes+2025 definition, which assumes relative asymmetry parameters ~ oblique dipole
                  averaged over all possible angles in the RMS sense.
    """
    assert isinstance(l, int) and isinstance(m1, int) and isinstance(m2, int) and isinstance(m3, int)
    assert (np.abs(m1) <= l) and (np.abs(m2) <= l) and (np.abs(m3) <= l)
    assert len(set([m1, m2, m3])) == 3 # must be distinct

    # Flip all signs if majority are negative, since asymmetry is invariant under this
    if (m1 < 0) + (m2 < 0) + (m3 < 0) >= 2:
        m1, m2, m3 = -m1, -m2, -m3
    
    m_tuple = tuple(sorted([m1, m2, m3]))

    # Use set type to ignore order of (m1, m2, m3)
    match asym_ref:
        case 'RFH+25':
            P2cosβ_rms = np.sqrt(2/5) # RMS value of cosβ over all dipole obliquity angles

            match l:
                case 1: # DIPOLE
                    match m_tuple: # TODO: tihs does'nt work
                        case (-1, 0, 1):
                            aasym = P2cosβ_rms * 2/5
                        case _:
                            raise ValueError(f'(m1, m2, m3) is not valid: {(m1, m2, m3)}')
                
                case 2: # QUADRUPOLE
                    match m_tuple:
                        case (-1, 0, 1) | (0, 1, 2):
                            aasym = P2cosβ_rms * -2/35
                        case (-2, 0, 1) | (-1, 0, 2) | (-1, 1, 2):
                            aasym = P2cosβ_rms * -6/35
                        case (-2, 1, 2):
                            aasym = P2cosβ_rms * -12/35
                        case (-2, 0, 2):
                            aasym = P2cosβ_rms * -8/35
                        case _:
                            raise ValueError(f'(m1, m2, m3) is not valid: {(m1, m2, m3)}')
                
                case _:
                    raise ValueError(f'asym_ref is not valid: {asym_ref}')
        
        case _:
            raise ValueError(f'asym_ref is not valid: {asym_ref}')
    
    return aasym

def compute_Δasym(l, m1, m2, m3,
                  ω1=None, ν1=None, ω1_uHz=None, ν1_uHz=None, P1=None,
                  ω2=None, ν2=None, ω2_uHz=None, ν2_uHz=None, P2=None,
                  ω3=None, ν3=None, ω3_uHz=None, ν3_uHz=None, P3=None,
                  return_uncertainty=False,
                  δω1=None, δν1=None, δω1_uHz=None, δν1_uHz=None, δP1=None,
                  δω2=None, δν2=None, δω2_uHz=None, δν2_uHz=None, δP2=None,
                  δω3=None, δν3=None, δω3_uHz=None, δν3_uHz=None, δP3=None):
    """
    Compute the dimensionful asymmetry parameter as described by RFH+25, Appendix A.2.

    This takes linear combinations of three frequencies ν1, ν2, ν3 defined as
    Δasym = c1 ν1 + c2 ν2 + c3 ν3

    such that both
    c1 + c2 + c3 = 0
    c1 m1 + c2 m2 + c3 m3 = 0

    We choose c1, c2, c3 such that the integers are the lowest possible, and at most one is negative.
    """
    assert isinstance(l, int) and isinstance(m1, int) and isinstance(m2, int) and isinstance(m3, int)
    assert (np.abs(m1) <= l) and (np.abs(m2) <= l) and (np.abs(m3) <= l)
    assert len(set([m1, m2, m3])) == 3 # must be distinct

    ω1 = get_ω(ω=ω1, ν=ν1, ω_uHz=ω1_uHz, ν_uHz=ν1_uHz, P=P1)
    ω2 = get_ω(ω=ω2, ν=ν2, ω_uHz=ω2_uHz, ν_uHz=ν2_uHz, P=P2)
    ω3 = get_ω(ω=ω3, ν=ν3, ω_uHz=ω3_uHz, ν_uHz=ν3_uHz, P=P3)
    ν1, ν2, ν3 = ω1 / (2 * np.pi), ω2 / (2 * np.pi), ω3 / (2 * np.pi)

    # Sort frequencies by m
    m_sorted, ν_sorted = zip(*sorted(zip((m1, m2, m3), (ν1, ν2, ν3))))
    m1, m2, m3 = m_sorted
    ν1, ν2, ν3 = ν_sorted

    # Compute coefficients
    # This is naturally preserved by the choice of coefficients:
    #  c1 = m3 - m2, c2 = m1 - m3, c3 = m2 - m1
    # and we will simply divide by the GCD to get the smallest integers
    c1, c2, c3 = m3 - m2, m1 - m3, m2 - m1
    g = gcd(gcd(abs(c1), abs(c2)), abs(c3))
    c1, c2, c3 = c1 // g, c2 // g, c3 // g

    # Compute Δasym
    Δasym = c1 * ν1 + c2 * ν2 + c3 * ν3

    if return_uncertainty:
        δω1 = get_δω(δω=δω1, δν=δν1, δω_uHz=δω1_uHz, δν_uHz=δν1_uHz, δP=δP1,
           P=P1)
        δω2 = get_δω(δω=δω2, δν=δν2, δω_uHz=δω2_uHz, δν_uHz=δν2_uHz, δP=δP2,
           P=P2)
        δω3 = get_δω(δω=δω3, δν=δν3, δω_uHz=δω3_uHz, δν_uHz=δν3_uHz, δP=δP3,
           P=P3)
        
        δν1, δν2, δν3 = δω1 / (2 * np.pi), δω2 / (2 * np.pi), δω3 / (2 * np.pi)

        δΔasym = np.sqrt(c1**2 * δν1**2 + c2**2 * δν2**2 + c3**2 * δν3**2)

        return Δasym, δΔasym
    
    else:
        return Δasym

