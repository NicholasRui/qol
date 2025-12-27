# Methods for Seismology class
# Calculate basic quantities related to seismic magnetometry
from qol.seismology.helper import get_ω, get_δω, get_magnetic_acrit, get_magnetic_aasym, compute_Δasym
from qol.tools.integrate import trapz_cond
import qol.mesa.const as const

import numpy as np
from scipy.integrate import cumulative_trapezoid



def get_magnetic_K(self, l, ω=None, ν=None, ω_uHz=None, ν_uHz=None, P=None, normed=True):
    """
    Magnetic weight function N^3 / Rho R^3, set to 0 outside of g-mode propagating region

    if normed, normalize K to unit-integral
    """
    assert self.Rho is not None

    ω = get_ω(ω=ω, ν=ν, ω_uHz=ω_uHz, ν_uHz=ν_uHz, P=P)
    is_g = self.get_is_g(l=l, ω=ω)

    magnetic_K = self.N ** 3 / self.Rho / self.R ** 3
    magnetic_K[~is_g] = 0

    if normed:
        magnetic_K /= trapz_cond(x=self.R, y=magnetic_K, c=is_g)
    
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

def get_magnetic_scriptI(self, l, ω=None, ν=None, ω_uHz=None, ν_uHz=None, P=None, use_f_corr_RFH25=False):
    """
    Magnetic sensitivity function scriptI = int_N3_div_Rho_r3_dr / int_N_div_r_dr
    """
    ω = get_ω(ω=ω, ν=ν, ω_uHz=ω_uHz, ν_uHz=ν_uHz, P=P)

    magnetic_K_unnorm = self.get_magnetic_K(l=l, ω=ω, normed=False)
    int_N3_div_Rho_r3_dr = trapz_cond(x=self.R, y=magnetic_K_unnorm)

    int_N_div_r_dr = self.get_int_N_div_r_dr(l, ω=ω)

    magnetic_scriptI = int_N3_div_Rho_r3_dr / int_N_div_r_dr

    if use_f_corr_RFH25:
        f_corr_RFH25 = self.get_f_corr_RFH25(l=l, ω=ω)
        magnetic_scriptI /= f_corr_RFH25 ** 2

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
    is_g = self.get_is_g(l=l, ω=ω)
    Bcrit[~is_g] = np.inf

    return Bcrit

def get_ωB(self, l, m=None, ω=None, ν=None, ω_uHz=None, ν_uHz=None, P=None, acrit_ref='F+15', units='ωB'):
    """
    Get maximum suppressed frequency ωB (this is the RF23 definition, NOT the Li et al. 2022 one)

    units allows the output to be returned in different formats
    """
    assert self.Rho is not None
    assert self.Br is not None

    acrit = get_magnetic_acrit(l=l, m=m, acrit_ref=acrit_ref)
    ωB = np.sqrt(self.N * self.Br / acrit / np.sqrt(4 * np.pi * self.Rho) / self.R)

    # If not propagating as a g mode, just set ωB to 0
    is_g = self.get_is_g(l=l, ω=ω)
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
    is_g = self.get_is_g(l=l, ω=ω)
    magnetic_K_norm = self.get_magnetic_K(l=l, ω=ω, normed=True)

    avg_Br2 = trapz_cond(x=self.R, y=magnetic_K_norm * self.Br ** 2, c=is_g)

    return avg_Br2

def get_δω_mag(self, l, ω=None, ν=None, ω_uHz=None, ν_uHz=None, P=None, units='δω_mag', use_f_corr_RFH25=False):
    """
    Get δω_mag (called ωB in Li+2022)

    Note the extra factor Al = l(l+1)/2, included in the definition in some works as well as here
    """
    ω = get_ω(ω=ω, ν=ν, ω_uHz=ω_uHz, ν_uHz=ν_uHz, P=P)

    magnetic_scriptI = self.get_magnetic_scriptI(l=l, ω=ω, use_f_corr_RFH25=use_f_corr_RFH25)
    avg_Br2 = self.get_avg_Br2(l=l, ω=ω)

    Al = l * (l + 1) / 2
    δω_mag = Al * magnetic_scriptI * avg_Br2 / 4 / np.pi / ω ** 3

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

def get_f_corr_RFH25(self, l, ω=None, ν=None, ω_uHz=None, ν_uHz=None, P=None):
    """
    non-asymptotic factor from Equation B3 of Rui, Fuller, Hermes 2025 (ApJ)
    """
    ω = get_ω(ω=ω, ν=ν, ω_uHz=ω_uHz, ν_uHz=ν_uHz, P=P)

    is_g = self.get_is_g(l=l, ω=ω)
    Sl_out = self.Sl(l)[is_g][0] # outermost value of Sl in g-mode cavity

    fN = 1 / np.sqrt(2)
    fS = 3
    f_corr_RFH25 = fN + (fS - fN) * (ω / Sl_out) ** 10.

    return f_corr_RFH25

def get_Brshift_from_νB_RFH25(self, νB_uHz, l, ω=None, ν=None, ω_uHz=None, ν_uHz=None, P=None, use_f_corr_RFH25=False):
    """
    For a given νB^l, calculate Brshift, as in RFH25 -- Equation 11, Figs. 3 and 7

    Note that Equation 11 has a typo and is missing a factor of (2l+1) in the denominator under the root.
    However, this doesn't matter here because we aren't involving the asymmetry.
    """
    # TODO: change νB to δν_mag
    ω = get_ω(ω=ω, ν=ν, ω_uHz=ω_uHz, ν_uHz=ν_uHz, P=P)
    P = 2 * const.pi / ω
    νB = 1e-6 * νB_uHz

    Al = l * (l + 1) / 2
    magnetic_scriptI = self.get_magnetic_scriptI(l=l, ω=ω, use_f_corr_RFH25=use_f_corr_RFH25)
    # Brshift = 8 * const.pi ** 2.5 * np.sqrt((2*l+1) * νB / (Al * magnetic_scriptI * P ** 3))
    Brshift = 8 * const.pi ** 2.5 * np.sqrt(νB / (Al * magnetic_scriptI * P ** 3))

    return Brshift

def get_νB_uHz_from_ν1_ν2_ν3_RFH25(self, l, m1, m2, m3,
                  ω1=None, ν1=None, ω1_uHz=None, ν1_uHz=None, P1=None,
                  ω2=None, ν2=None, ω2_uHz=None, ν2_uHz=None, P2=None,
                  ω3=None, ν3=None, ω3_uHz=None, ν3_uHz=None, P3=None,
                  return_uncertainty=False,
                  δω1=None, δν1=None, δω1_uHz=None, δν1_uHz=None, δP1=None,
                  δω2=None, δν2=None, δω2_uHz=None, δν2_uHz=None, δP2=None,
                  δω3=None, δν3=None, δω3_uHz=None, δν3_uHz=None, δP3=None,
                  asym_ref='RFH+25'):
    """
    Convert three frequencies within a multiplet and their quantum numbers to νB, as in RFH25 -- Equation 11 (note typo).
    """
    aasym = get_magnetic_aasym(l, m1, m2, m3, asym_ref=asym_ref)

    if return_uncertainty:
        Δasym, δΔasym = compute_Δasym(l, m1, m2, m3,
                    ω1=ω1, ν1=ν1, ω1_uHz=ω1_uHz, ν1_uHz=ν1_uHz, P1=P1,
                    ω2=ω2, ν2=ν2, ω2_uHz=ω2_uHz, ν2_uHz=ν2_uHz, P2=P2,
                    ω3=ω3, ν3=ν3, ω3_uHz=ω3_uHz, ν3_uHz=ν3_uHz, P3=P3,
                    return_uncertainty=True,
                    δω1=δω1, δν1=δν1, δω1_uHz=δω1_uHz, δν1_uHz=δν1_uHz, δP1=δP1,
                    δω2=δω2, δν2=δν2, δω2_uHz=δω2_uHz, δν2_uHz=δν2_uHz, δP2=δP2,
                    δω3=δω3, δν3=δν3, δω3_uHz=δω3_uHz, δν3_uHz=δν3_uHz, δP3=δP3)

        νB_uHz = 1e6 * Δasym / aasym / (2 * l + 1)
        δνB_uHz = 1e6 * δΔasym / aasym / (2 * l + 1)

        return νB_uHz, δνB_uHz

    else:
        Δasym = compute_Δasym(l, m1, m2, m3,
                    ω1=ω1, ν1=ν1, ω1_uHz=ω1_uHz, ν1_uHz=ν1_uHz, P1=P1,
                    ω2=ω2, ν2=ν2, ω2_uHz=ω2_uHz, ν2_uHz=ν2_uHz, P2=P2,
                    ω3=ω3, ν3=ν3, ω3_uHz=ω3_uHz, ν3_uHz=ν3_uHz, P3=P3,
                    return_uncertainty=False)
        
        νB_uHz = 1e6 * Δasym / aasym / (2 * l + 1)

        return νB_uHz

def get_Brshift_from_ν1_ν2_ν3_RFH25(self, l, m1, m2, m3,
                  ω1=None, ν1=None, ω1_uHz=None, ν1_uHz=None, P1=None,
                  ω2=None, ν2=None, ω2_uHz=None, ν2_uHz=None, P2=None,
                  ω3=None, ν3=None, ω3_uHz=None, ν3_uHz=None, P3=None,
                  include_uncertainty=True,
                  Nsigma=2,
                  δω1=None, δν1=None, δω1_uHz=None, δν1_uHz=None, δP1=None,
                  δω2=None, δν2=None, δω2_uHz=None, δν2_uHz=None, δP2=None,
                  δω3=None, δν3=None, δω3_uHz=None, δν3_uHz=None, δP3=None,
                  asym_ref='RFH+25', use_f_corr_RFH25=False):
    """
    For three frequencies within a multiplet and their quantum numbers, calculate Brshift, as in RFH25 -- Equation 11.
    Note that Equation 11 has a typo and is missing a factor of (2l+1) in the denominator under the root.

    Nsigma describes the number of sigma to use when computing the most conservative bound on Brshift.
    Default is 2 following RFH+25.
    """
    # Get νB in uHz by calculating dimensionful asymmetry with assumed dimensionless asymmetry
    if include_uncertainty:
        νB_uHz, δνB_uHz = self.get_νB_uHz_from_ν1_ν2_ν3_RFH25(l, m1, m2, m3,
                  ω1=ω1, ν1=ν1, ω1_uHz=ω1_uHz, ν1_uHz=ν1_uHz, P1=P1,
                  ω2=ω2, ν2=ν2, ω2_uHz=ω2_uHz, ν2_uHz=ν2_uHz, P2=P2,
                  ω3=ω3, ν3=ν3, ω3_uHz=ω3_uHz, ν3_uHz=ν3_uHz, P3=P3,
                  return_uncertainty=True,
                  δω1=δω1, δν1=δν1, δω1_uHz=δω1_uHz, δν1_uHz=δν1_uHz, δP1=δP1,
                  δω2=δω2, δν2=δν2, δω2_uHz=δω2_uHz, δν2_uHz=δν2_uHz, δP2=δP2,
                  δω3=δω3, δν3=δν3, δω3_uHz=δω3_uHz, δν3_uHz=δν3_uHz, δP3=δP3,
                  asym_ref='RFH+25')
    
        δω1 = get_δω(δω=δω1, δν=δν1, δω_uHz=δω1_uHz, δν_uHz=δν1_uHz, δP=δP1,
            P=P1)
        δω2 = get_δω(δω=δω2, δν=δν2, δω_uHz=δω2_uHz, δν_uHz=δν2_uHz, δP=δP2,
            P=P1)
        δω3 = get_δω(δω=δω3, δν=δν3, δω_uHz=δω3_uHz, δν_uHz=δν3_uHz, δP=δP3,
            P=P3)
        δν1, δν2, δν3 = δω1 / (2 * np.pi), δω2 / (2 * np.pi), δω3 / (2 * np.pi)
        δνs = np.array([δν1, δν2, δν3])
        ws = 1 / δνs**2 # for linear fit -- chi square weighting

    else:
        νB_uHz = self.get_νB_uHz_from_ν1_ν2_ν3_RFH25(l, m1, m2, m3,
                  ω1=ω1, ν1=ν1, ω1_uHz=ω1_uHz, ν1_uHz=ν1_uHz, P1=P1,
                  ω2=ω2, ν2=ν2, ω2_uHz=ω2_uHz, ν2_uHz=ν2_uHz, P2=P2,
                  ω3=ω3, ν3=ν3, ω3_uHz=ω3_uHz, ν3_uHz=ν3_uHz, P3=P3,
                  return_uncertainty=False,
                  asym_ref='RFH+25')
        δνB_uHz = 0.
        ws = np.ones(3) # for linear fit -- equal fit
    
    # Take most conservative absolute value of νB_uHz, according to specified Nsigma
    νB_uHz = np.abs(νB_uHz) + Nsigma * δνB_uHz

    # Find linear fit to get value of νi at m=0 for some approximation of central value
    ω1 = get_ω(ω=ω1, ν=ν1, ω_uHz=ω1_uHz, ν_uHz=ν1_uHz, P=P1)
    ω2 = get_ω(ω=ω2, ν=ν2, ω_uHz=ω2_uHz, ν_uHz=ν2_uHz, P=P2)
    ω3 = get_ω(ω=ω3, ν=ν3, ω_uHz=ω3_uHz, ν_uHz=ν3_uHz, P=P3)
    ν1, ν2, ν3 = ω1 / (2 * np.pi), ω2 / (2 * np.pi), ω3 / (2 * np.pi)

    ms, νs = np.array([m1, m2, m3]), np.array([ν1, ν2, ν3])
    mbar = np.sum(ws * ms) / np.sum(ws)
    νbar = np.sum(ws * νs) / np.sum(ws)
    dν_dm = np.sum(ws * (ms - mbar) * (νs - νbar)) / np.sum(ws * (ms - mbar)**2)
    ν0 = νbar - dν_dm * mbar

    Brshift = self.get_Brshift_from_νB_RFH25(νB_uHz, l, ν=ν0, use_f_corr_RFH25=use_f_corr_RFH25)

    return Brshift

    



