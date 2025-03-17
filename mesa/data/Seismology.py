"""
Purpose: Class keeps track of seismic quantities and performs common calculations
"""
from qol.mesa.data.MesaTable import MesaTable
import qol.mesa.const as const
from qol.helper.integrate import trapz_cond

import numpy as np
from scipy.integrate import cumulative_trapezoid

from typing import Union
import warnings

############################
##### HELPER FUNCTIONS #####
############################

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






class Seismology:
    """
    """
    def __init__(self,
            mesa_table: MesaTable = None,
            R: np.ndarray = None, R_in_Rsun: np.ndarray = None,
            N: np.ndarray = None, N_in_uHz: np.ndarray = None, N_floor: float = 1e-30,
            Sl1: np.ndarray = None, Sl1_in_uHz: np.ndarray = None,

            Rho: np.ndarray = None,
            Br: np.ndarray = None, Br_kG: np.ndarray = None, Br_MG: np.ndarray = None,

            Ωrot: Union[np.ndarray, list, float] = None, νrot: Union[np.ndarray, list, float] = None,
            Ωrot_uHz: Union[np.ndarray, list, float] = None, νrot_uHz: Union[np.ndarray, list, float] = None,
            Prot: Union[np.ndarray, list, float] = None, Prot_d: Union[np.ndarray, list, float] = None
                 ):
        """
        At minimum, require that we have an N and an Sl (without these, seismology cannot be done)

        If these are specified explicitly, then use the value. Else, take it from the mesa_table.
        Make sure everything is the same radius.

        N_floor clips N to the floor, to avoid division by 0 in some calculations

        TODO:
        - maybe ask for N2 instead?
        - ask for csound and calculate acoustic radius?
        - shove writing of N, R, etc., etc., into their own attributes to clean up __init__

        - compute averages of Omegarot, and delta_nu_rot
        """
        self.mesa_table = mesa_table

        # Read in inputs
        self.initialize_R(R=R, R_in_Rsun=R_in_Rsun)
        self.initialize_N(N=N, N_in_uHz=N_in_uHz, N_floor=N_floor)
        self.initialize_Sl(Sl1=Sl1, Sl1_in_uHz=Sl1_in_uHz)
        self.initialize_Br(Br=Br, Br_kG=Br_kG, Br_MG=Br_MG)
        self.initialize_Ωrot(Ωrot=Ωrot, νrot=νrot, Ωrot_uHz=Ωrot_uHz, νrot_uHz=νrot_uHz, Prot=Prot, Prot_d=Prot_d)

        self.Rho = np.array(Rho) if Rho is not None else None
        
        # Assert that radii are purely monotonic (allow for zero differences just in case there are precision issues)
        if (np.diff(self.R) >= 0).all():
            self.increasing_R = True
        elif (np.diff(self.R) <= 0).all():
            self.increasing_R = False
        else:
            raise ValueError('R is not purely increasing or decreasing')
        
        # Require that arrays that are specified are consistent lengths
        assert len(self.R) == len(self.N)
        assert len(self.R) == len(self.Sl1)

        if self.Rho is not None:
            assert len(self.R) == len(self.Rho)
        if self.Br is not None:
            assert len(self.R) == len(self.Br)
        if type(self.Ωrot) in [np.ndarray, list]:
            assert len(self.R) == len(self.Ωrot)

    ##############################
    ##### READ IN QUANTITIES #####
    ##############################
    def initialize_R(self, R, R_in_Rsun):
        self.R = np.array(R) if R is not None else None
    
        # Some checks to warn user if redundant keywords have been specified
        if (self.R is not None) and (R_in_Rsun is not None):
            warnings.warn('Both R and R_in_Rsun specified: defaulting to former')

        self.R = const.Rsun * np.array(R_in_Rsun) if self.R is None else self.R

        if hasattr(self.mesa_table, 'R'):
            if self.R is not None:
                warnings.warn('R already explicitly defined, so not using column in mesa_table')
            else:
                self.R = self.mesa_table.R

        # Required field
        assert self.R is not None

    def initialize_N(self, N, N_in_uHz, N_floor):
        self.N = np.array(N) if N is not None else None
        self.N_floor = N_floor

        # Some checks to warn user if redundant keywords have been specified
        if (self.N is not None) and (N_in_uHz is not None):
            warnings.warn('Both N and N_in_uHz specified: defaulting to former')
        
        # Alternative inputs
        self.N = 1e-6 * np.array(N_in_uHz) if self.N is None else self.N

        if hasattr(self.mesa_table, 'N'):
            if self.N is not None:
                warnings.warn('N already explicitly defined, so not using column in mesa_table')
            else:
                self.N = self.mesa_table.N

        # Clip zero N to a finite value, so zero-division doesn't happen
        self.N = np.clip(self.N, a_min=N_floor, a_max=None)

        # Calculate other helper quantities
        self.N_div_2pi = self.N / (2 * const.pi)
        self.N_in_uHz = 1e6 * self.N
        self.N_div_2pi_in_uHz = 1e6 * self.N_div_2pi

        # Required field
        assert self.N is not None

    def initialize_Sl(self, Sl1, Sl1_in_uHz):
        self.Sl1 = np.array(Sl1) if Sl1 is not None else None

        # Some checks to warn user if redundant keywords have been specified
        if (self.Sl1 is not None) and (Sl1_in_uHz is not None):
            warnings.warn('Both Sl1 and Sl1_in_uHz specified: defaulting to former')

        # Alternative inputs
        self.Sl1 = 1e-6 * np.array(Sl1_in_uHz) if self.Sl1 is None else self.Sl1

        if hasattr(self.mesa_table, 'Sl1'):
            if self.Sl1 is not None:
                warnings.warn('Sl1 already explicitly defined, so not using column in mesa_table')
            else:
                self.Sl1 = self.mesa_table.Sl1

        # Calculate other helper quantities
        self.Sl2 = np.sqrt(3) * self.Sl1
        self.Sl3 = np.sqrt(6) * self.Sl1
        self.Sl = lambda l: np.sqrt(l * (l + 1) / 2) * self.Sl1 # lambda for arbitrary l

        self.Sl1_div_2pi = self.Sl1 / (2 * const.pi)
        self.Sl2_div_2pi = self.Sl2 / (2 * const.pi)
        self.Sl3_div_2pi = self.Sl3 / (2 * const.pi)
        self.Sl_div_2pi = lambda l: self.Sl(l) / (2 * const.pi)

        self.Sl1_in_uHz = 1e6 * self.Sl1
        self.Sl2_in_uHz = 1e6 * self.Sl2
        self.Sl3_in_uHz = 1e6 * self.Sl3
        self.Sl_in_uHz = lambda l: 1e6 * self.Sl(l)

        self.Sl1_div_2pi_in_uHz = 1e6 * self.Sl1_div_2pi
        self.Sl2_div_2pi_in_uHz = 1e6 * self.Sl2_div_2pi
        self.Sl3_div_2pi_in_uHz = 1e6 * self.Sl3_div_2pi
        self.Sl_div_2pi_in_uHz = lambda l: 1e6 * self.Sl_div_2pi(l)

        # Required field
        assert self.Sl1 is not None

    def initialize_Br(self, Br, Br_kG, Br_MG):
        self.Br = np.array(Br) if Br is not None else None

        # Some checks to warn user if redundant keywords have been specified
        if (self.Br is not None) and (Br_kG is not None):
            warnings.warn('Both Br and Br_kG specified: defaulting to former')
        if (self.Br is not None) and (Br_MG is not None):
            warnings.warn('Both Br and Br_MG specified: defaulting to former')
        if (Br_kG is not None) and (Br_MG is not None):
            warnings.warn('Both Br_kG and Br_MG specified: defaulting to former') 

        # Alternative inputs
        self.Br = 1e3 * np.array(self.Br_kG) if self.Br is None else self.Br
        self.Br = 1e6 * np.array(self.Br_MG) if self.Br is None else self.Br

    def initialize_Ωrot(self, Ωrot, νrot, Ωrot_uHz, νrot_uHz, Prot, Prot_d):
        """
        Ωrot can be either array-like or float
        """
        self.Ωrot = np.array(Ωrot) if type(Ωrot) in [np.ndarray, list] else Ωrot

        # Some checks to warn user if redundant keywords have been specified
        if (self.Ωrot is not None) and (νrot is not None):
            warnings.warn('Both Ωrot and νrot specified: defaulting to former')
        if (self.Ωrot is not None) and (Ωrot_uHz is not None):
            warnings.warn('Both Ωrot and Ωrot_uHz specified: defaulting to former')
        if (self.Ωrot is not None) and (νrot_uHz is not None):
            warnings.warn('Both Ωrot and νrot_uHz specified: defaulting to former')
        if (self.Ωrot is not None) and (Prot is not None):
            warnings.warn('Both Ωrot and Prot specified: defaulting to former')
        if (self.Ωrot is not None) and (Prot_d is not None):
            warnings.warn('Both Ωrot and Prot_d specified: defaulting to former')

        # If array-like, make into array
        Ωrot = np.array(Ωrot) if type(Ωrot) in [np.ndarray, list] else Ωrot
        νrot = np.array(νrot) if type(νrot) in [np.ndarray, list] else νrot
        Ωrot_uHz = np.array(Ωrot_uHz) if type(Ωrot_uHz) in [np.ndarray, list] else Ωrot_uHz
        νrot_uHz = np.array(νrot_uHz) if type(νrot_uHz) in [np.ndarray, list] else νrot_uHz
        Prot = np.array(Prot) if type(Prot) in [np.ndarray, list] else Prot
        Prot_d = np.array(Prot_d) if type(Prot_d) in [np.ndarray, list] else Prot_d

        # Alternative inputs
        self.Ωrot = 2 * np.pi * νrot if self.Ωrot is None else self.Ωrot
        self.Ωrot = 1e-6 * Ωrot_uHz if self.Ωrot is None else self.Ωrot
        self.Ωrot = 2 * np.pi * 1e-6 * νrot_uHz if self.Ωrot is None else self.Ωrot
        self.Ωrot = 1 / Prot if self.Ωrot is None else self.Ωrot
        self.Ωrot = 1 / (const.secday * Prot_d) if self.Ωrot is None else self.Ωrot

    ############################
    ##### HELPER FUNCTIONS #####
    ############################
    def get_ω(self, ω=None, ν=None, ω_uHz=None, ν_uHz=None, P=None):
        """
        Helper method which allows other methods to take a frequency in multiple formats,
          and asserts that the user has only specified one of them
        """
        provided_args = [ω, ν, ω_uHz, ν_uHz, P]
        assert provided_args == 1

        ω = 2 * np.pi * ν if ω is None else ω
        ω = 1e-6 * ω_uHz if ω is None else ω
        ω = 2 * np.pi * 1e-6 * ν_uHz if ω is None else ω
        ω = 2 * np.pi / P if ω is None else ω

        return ω

    ##############################################
    ##### CALCULATE BASIC SEISMIC QUANTITIES #####
    ##############################################
    def get_is_prop(self, l, ω=None, ν=None, ω_uHz=None, ν_uHz=None, P=None, proptype=None):
        """
        Return array of Booleans which says whether a mode with a given frequency is propagating or not

        if proptype is specified and is either 'p' or 'g', return boolean
          showing only whether a mode would be propagating with that particular character
        """
        ω = self.get_ω(ω=ω, ν=ν, ω_uHz=ω_uHz, ν_uHz=ν_uHz, P=P)

        above_Sl1 = (ω >= self.Sl(l))
        above_N = (ω >= self.N)
        
        match proptype:
            case None:
                return (above_Sl1 & above_N) & (~above_Sl1 & ~above_N)
            case 'p':
                return (above_Sl1 & above_N)
            case 'g':
                return (~above_Sl1 & ~above_N)
            case _:
                raise ValueError(f'invalid proptype: {proptype}')

    def get_int_N_div_r_dr(self, l, ω=None, ν=None, ω_uHz=None, ν_uHz=None, P=None):
        """
        Calculate int_N_div_r_dr over g-propagating regions (no cavity contiguity enforced)
        """
        ω = self.get_ω(ω=ω, ν=ν, ω_uHz=ω_uHz, ν_uHz=ν_uHz, P=P)
        is_g = self.get_is_prop(l=l, ω=ω, proptype='g')

        int_N_div_r_dr = trapz_cond(x=self.R, y=self.N / self.R, c=is_g)

        return int_N_div_r_dr
    
    def get_delta_Pg(self, l, ω=None, ν=None, ω_uHz=None, ν_uHz=None, P=None):
        """
        Calculate g-mode period spacing
        """
        int_N_div_r_dr = self.get_int_N_div_r_dr(l=l, ω=ω)
        delta_Pg = 2 * np.pi ** 2 / int_N_div_r_dr / np.sqrt(l * (l + 1))

        return delta_Pg
    
    ########################################
    ##### SEISMIC MAGNETOMETRY METHODS #####
    ########################################

    def get_magnetic_K(self, l, ω=None, ν=None, ω_uHz=None, ν_uHz=None, P=None, normed=True):
        """
        Magnetic weight function N^3 / Rho R^3, set to 0 outside of g-mode propagating region

        if normed, normalize K to unit-integral
        """
        assert self.Rho is not None

        ω = self.get_ω(ω=ω, ν=ν, ω_uHz=ω_uHz, ν_uHz=ν_uHz, P=P)
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
        ω = self.get_ω(ω=ω, ν=ν, ω_uHz=ω_uHz, ν_uHz=ν_uHz, P=P)
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
        ω = self.get_ω(ω=ω, ν=ν, ω_uHz=ω_uHz, ν_uHz=ν_uHz, P=P)

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

        ω = self.get_ω(ω=ω, ν=ν, ω_uHz=ω_uHz, ν_uHz=ν_uHz, P=P)
        acrit = get_magnetic_acrit(l=l, m=m, acrit_ref=acrit_ref)
        
        Bcrit = acrit * np.sqrt(4 * np.pi * self.Rho) * ω ** 2 * self.R / self.N

        # If not propagating as g mode, just set Bcrit to infinity
        is_g = self.get_is_prop(l=l, ω=ω, proptype='g')
        Bcrit[~is_g] = np.inf

        return Bcrit
    
    def get_ωB(self, l, m=None, ω=None, ν=None, ω_uHz=None, ν_uHz=None, P=None, acrit_ref='F+15'):
        """
        Get maximum suppressed frequency ωB (this is the RF+23 definition, NOT the Li et al. 2022 one)

        See documentation for get_Bcrit for a description of ac_ref
        """
        assert self.Rho is not None
        assert self.Br is not None

        acrit = get_magnetic_acrit(l=l, m=m, acrit_ref=acrit_ref)
        ωB = np.sqrt(self.N * self.Br / acrit / np.sqrt(4 * np.pi * self.Rho) / self.R)

        # If not propagating as a g mode, just set ωB to 0
        is_g = self.get_is_prop(l=l, ω=ω, proptype='g')
        ωB[~is_g] = np.inf

        return ωB

    def get_ωB_uHz(self, l, m=None, ω=None, ν=None, ω_uHz=None, ν_uHz=None, P=None, acrit_ref='F+15'):
        """
        Calls get_ωB but converts output to different units
        """
        ωB = self.get_ωB(l=l, m=m, ω=ω, ν=ν, ω_uHz=ω_uHz, ν_uHz=ν_uHz, P=P, acrit_ref=acrit_ref)
        ωB_uHz = 1e6 * ωB

        return ωB_uHz
    
    def get_νB(self, l, m=None, ω=None, ν=None, ω_uHz=None, ν_uHz=None, P=None, acrit_ref='F+15'):
        """
        Calls get_ωB but converts output to different units
        """
        ωB = self.get_ωB(l=l, m=m, ω=ω, ν=ν, ω_uHz=ω_uHz, ν_uHz=ν_uHz, P=P, acrit_ref=acrit_ref)
        νB = ωB / 2 / np.pi

        return νB
    
    def get_νB_uHz(self, l, m=None, ω=None, ν=None, ω_uHz=None, ν_uHz=None, P=None, acrit_ref='F+15'):
        """
        Calls get_ωB but converts output to different units
        """
        ωB = self.get_ωB(l=l, m=m, ω=ω, ν=ν, ω_uHz=ω_uHz, ν_uHz=ν_uHz, P=P, acrit_ref=acrit_ref)
        νB_uHz = 1e6 * ωB / 2 / np.pi

        return νB_uHz

    def get_avg_Br2(self, l, ω=None, ν=None, ω_uHz=None, ν_uHz=None, P=None):
        """
        Get average <Br^2> wrt. magnetic K function
        """
        assert self.Br is not None

        ω = self.get_ω(ω=ω, ν=ν, ω_uHz=ω_uHz, ν_uHz=ν_uHz, P=P)
        is_g = self.get_is_prop(l=l, ω=ω, proptype='g')
        magnetic_K_norm = self.get_magnetic_K(l=l, ω=ω, normed=True)

        avg_Br2 = trapz_cond(x=self.R, y=magnetic_K_norm * self.Br ** 2)

        return avg_Br2

    def get_δω_mag(self, l, ω=None, ν=None, ω_uHz=None, ν_uHz=None, P=None):
        """
        Get δω_mag (called ωB in Li+2022)
        """
        ω = self.get_ω(ω=ω, ν=ν, ω_uHz=ω_uHz, ν_uHz=ν_uHz, P=P)

        magnetic_scriptI = self.get_magnetic_scriptI(l=l, ω=ω)
        avg_Br2 = self.get_avg_Br2(l=l, ω=ω)

        δω_mag = magnetic_scriptI * avg_Br2 / 4 / np.pi / ω ** 3

        return δω_mag

    def get_δν_mag(self, l, ω=None, ν=None, ω_uHz=None, ν_uHz=None, P=None):
        """
        Calls get_ωB but converts output to different units
        """
        δω_mag = self.get_δω_mag(l=l, ω=ω, ν=ν, ω_uHz=ω_uHz, ν_uHz=ν_uHz, P=P)
        δν_mag = δω_mag / 2 / np.pi

        return δν_mag
    
    def get_δω_mag_uHz(self, l, ω=None, ν=None, ω_uHz=None, ν_uHz=None, P=None):
        """
        Calls get_ωB but converts output to different units
        """
        δω_mag = self.get_δω_mag(l=l, ω=ω, ν=ν, ω_uHz=ω_uHz, ν_uHz=ν_uHz, P=P)
        δω_mag_uHz = 1e6 * δω_mag_uHz

        return δω_mag_uHz

    def get_δν_mag_uHz(self, l, ω=None, ν=None, ω_uHz=None, ν_uHz=None, P=None):
        """
        Calls get_ωB but converts output to different units
        """
        δω_mag = self.get_δω_mag(l=l, ω=ω, ν=ν, ω_uHz=ω_uHz, ν_uHz=ν_uHz, P=P)
        δν_mag_uHz = 1e6 * δω_mag / 2 / np.pi

        return δν_mag_uHz














































