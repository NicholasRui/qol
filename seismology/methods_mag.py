# Methods for Seismology class
# Calculate basic quantities related to seismic magnetometry
from qol.seismology.helper import get_ω, get_magnetic_acrit
from qol.tools.integrate import trapz_cond

import numpy as np
from scipy.integrate import cumulative_trapezoid



def get_magnetic_K(self, l, ω=None, ν=None, ω_uHz=None, ν_uHz=None, P=None, normed=True):
    """
    Magnetic weight function N^3 / Rho R^3, set to 0 outside of g-mode propagating region

    if normed, normalize K to unit-integral
    """
    assert self.Rho is not None

    ω = get_ω(ω=ω, ν=ν, ω_uHz=ω_uHz, ν_uHz=ν_uHz, P=P)
    is_g = self.get_is_prop(l=l, ω=ω, proptype='g')

    magnetic_K = self.N ** 3 / self.Rho / self.R ** 3
    magnetic_K[~is_g] = 0

    if normed:
        magnetic_K /= trapz_cond(x=self.R, y=magnetic_K)
    
    return magnetic_K

def get_magnetic_cumK(self, l, ω=None, ν=None, ω_uHz=None, ν_uHz=None, P=None, normed=True):
    """
    Cumulative integral of magnetic weight function, starting from center
    """
    ω = get_ω(ω=ω, ν=ν, ω_uHz=ω_uHz, ν_uHz=ν_uHz, P=P)
    magnetic_K = self.get_magnetic_K(l=l, ω=ω, normed=normed)

    if self.increasing_R:
        magnetic_cumK = cumulative_trapezoid(x=self.R, y=magnetic_K)
    if self.increasing_R:
        magnetic_cumK = np.flip(cumulative_trapezoid(x=np.flip(self.R), y=np.flip(magnetic_K)))
    
    return magnetic_cumK

def get_magnetic_scriptI(self, l, ω=None, ν=None, ω_uHz=None, ν_uHz=None, P=None):
    """
    Magnetic sensitivity function \mathscr{I} = int_N3_div_Rho_R3_dr / int_N_div_r_dr
    """
    ω = get_ω(ω=ω, ν=ν, ω_uHz=ω_uHz, ν_uHz=ν_uHz, P=P)

    magnetic_K_unnorm = self.get_magnetic_K(l=l, ω=ω, normed=False)
    int_N3_div_Rho_R3_dr = trapz_cond(x=self.R, y=magnetic_K_unnorm)
    int_N_div_r_dr = self.get_int_N_div_r_dr(l, ω=ω)

    magnetic_scriptI = int_N3_div_Rho_R3_dr / int_N_div_r_dr

    return magnetic_scriptI

def get_Bcrit(self, l, m=None, ω=None, ν=None, ω_uHz=None, ν_uHz=None, P=None, acrit_ref='F+15'):
    """
    Calculate critical (radial) magnetic field Bcrit needed to suppress g modes,
        according to Fuller et al. 2015 (F+15)
    """
    assert self.Rho is not None

    ω = get_ω(ω=ω, ν=ν, ω_uHz=ω_uHz, ν_uHz=ν_uHz, P=P)
    acrit = get_magnetic_acrit(l=l, m=m, acrit_ref=acrit_ref)
    
    Bcrit = acrit * np.sqrt(4 * np.pi * self.Rho) * ω ** 2 * self.R / self.N

    # If not propagating as g mode, just set Bcrit to infinity
    is_g = self.get_is_prop(l=l, ω=ω, proptype='g')
    Bcrit[~is_g] = np.inf

    return Bcrit

def get_ωB(self, l, m=None, ω=None, ν=None, ω_uHz=None, ν_uHz=None, P=None, acrit_ref='F+15', units='ωB'):
    """
    Get maximum suppressed frequency ωB (this is the RF+23 definition, NOT the Li et al. 2022 one)

    units allows the output to be returned in different formats
    """
    assert self.Rho is not None
    assert self.Br is not None

    acrit = get_magnetic_acrit(l=l, m=m, acrit_ref=acrit_ref)
    ωB = np.sqrt(self.N * self.Br / acrit / np.sqrt(4 * np.pi * self.Rho) / self.R)

    # If not propagating as a g mode, just set ωB to 0
    is_g = self.get_is_prop(l=l, ω=ω, proptype='g')
    ωB[~is_g] = np.inf

    match units:
        case 'ωB':
            return ωB
        
        case 'ωB_uHz':
            ωB_uHz = 1e6 * ωB
            return ωB_uHz
        
        case 'νB':
            νB = ωB / 2 / np.pi
            return νB
        
        case 'νB_uHz':
            νB_uHz = 1e6 * ωB / 2 / np.pi
            return νB_uHz
        
        case _:
            raise ValueError(f'invalid units: {units}')

def get_avg_Br2(self, l, ω=None, ν=None, ω_uHz=None, ν_uHz=None, P=None):
    """
    Get average <Br^2> wrt. magnetic K function
    """
    assert self.Br is not None

    ω = get_ω(ω=ω, ν=ν, ω_uHz=ω_uHz, ν_uHz=ν_uHz, P=P)
    is_g = self.get_is_prop(l=l, ω=ω, proptype='g')
    magnetic_K_norm = self.get_magnetic_K(l=l, ω=ω, normed=True)

    avg_Br2 = trapz_cond(x=self.R, y=magnetic_K_norm * self.Br ** 2)

    return avg_Br2

def get_δω_mag(self, l, ω=None, ν=None, ω_uHz=None, ν_uHz=None, P=None, units='δω_mag'):
    """
    Get δω_mag (called ωB in Li+2022)
    """
    ω = get_ω(ω=ω, ν=ν, ω_uHz=ω_uHz, ν_uHz=ν_uHz, P=P)

    magnetic_scriptI = self.get_magnetic_scriptI(l=l, ω=ω)
    avg_Br2 = self.get_avg_Br2(l=l, ω=ω)

    δω_mag = magnetic_scriptI * avg_Br2 / 4 / np.pi / ω ** 3

    match units:
        case 'δω_mag':
            return δω_mag
        
        case 'δν_mag':
            δν_mag = δω_mag / 2 / np.pi
            return δν_mag
        
        case 'δω_mag_uHz':
            δω_mag_uHz = 1e6 * δω_mag
            return δω_mag_uHz
        
        case 'δν_mag_uHz':
            δν_mag_uHz = 1e6 * δω_mag / 2 / np.pi
            return δν_mag_uHz
        
        case _:
            raise ValueError(f'invalid units: {units}')

