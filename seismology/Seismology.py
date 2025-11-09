"""
Purpose: Class keeps track of seismic quantities and performs common calculations
"""
from qol.mesa.table.MesaTable import MesaTable
import qol.mesa.const as const
from qol.tools.integrate import trapz_cond

import qol.seismology.methods_prop as methods_prop
import qol.seismology.methods_mag as methods_mag
import qol.seismology.methods_rot as methods_rot

import numpy as np

from typing import Union
import warnings



class Seismology:
    """
    """
    # Import attributes from other files
    get_is_prop = methods_prop.get_is_prop
    get_is_g = methods_prop.get_is_g
    get_is_p = methods_prop.get_is_p
    get_int_N_div_r_dr = methods_prop.get_int_N_div_r_dr
    get_delta_Pg = methods_prop.get_delta_Pg
    
    get_magnetic_K = methods_mag.get_magnetic_K
    get_magnetic_cumK = methods_mag.get_magnetic_cumK
    get_magnetic_scriptI = methods_mag.get_magnetic_scriptI
    get_Bcrit = methods_mag.get_Bcrit
    get_ωB = methods_mag.get_ωB
    get_avg_Br2 = methods_mag.get_avg_Br2
    get_δω_mag = methods_mag.get_δω_mag
    get_f_corr_RFH25 = methods_mag.get_f_corr_RFH25
    get_Brshift_RFH25 = methods_mag.get_Brshift_from_νB_RFH25
    get_νB_uHz_from_ν1_ν2_ν3_RFH25 = methods_mag.get_νB_uHz_from_ν1_ν2_ν3_RFH25
    get_Brshift_from_ν1_ν2_ν3_RFH25 = methods_mag.get_Brshift_from_ν1_ν2_ν3_RFH25

    get_Ωrot_g = methods_rot.get_Ωrot_g

    def __init__(self,
            mesa_table: MesaTable = None,
            R: np.ndarray = None, R_in_Rsun: np.ndarray = None,
            N: np.ndarray = None, N_in_uHz: np.ndarray = None, N_floor: float = 1e-30,
            Sl1: np.ndarray = None, Sl1_in_uHz: np.ndarray = None,

            Rho: np.ndarray = None,
            Br: np.ndarray = None, Br_kG: np.ndarray = None, Br_MG: np.ndarray = None,

            Ωrot: Union[np.ndarray, list, float] = None, νrot: Union[np.ndarray, list, float] = None,
            Ωrot_uHz: Union[np.ndarray, list, float] = None, νrot_uHz: Union[np.ndarray, list, float] = None,
            Prot: Union[np.ndarray, list, float] = None, Prot_d: Union[np.ndarray, list, float] = None,

            Rin_in_R: float = None
                 ):
        """
        At minimum, require that we have an N and an Sl (without these, seismology cannot be done)

        If these are specified explicitly, then use the value. Else, take it from the mesa_table.
        Make sure everything is the same radius.

        N_floor clips N to the floor, to avoid division by 0 in some calculations

        Rin_in_R defines the inner g-mode cavity, in fractional radii (default: 0)

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
        self.calculate_min_max_N_Sl()
        self.initialize_Br(Br=Br, Br_kG=Br_kG, Br_MG=Br_MG)
        self.initialize_Ωrot(Ωrot=Ωrot, νrot=νrot, Ωrot_uHz=Ωrot_uHz, νrot_uHz=νrot_uHz, Prot=Prot, Prot_d=Prot_d)
        self.initialize_Rho(Rho=Rho)
        
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
        
        if Rin_in_R is not None:
            self.Rin_in_R = Rin_in_R
        else:
            self.Rin_in_R = 0

    ##############################
    ##### READ IN QUANTITIES #####
    ##############################

    def initialize_Rho(self, Rho):
        self.Rho = np.array(Rho) if Rho is not None else None
    
        if hasattr(self.mesa_table, 'Rho'):
            if self.Rho is not None:
                warnings.warn('Rho already explicitly defined, so not using column in mesa_table')
            else:
                self.Rho = self.mesa_table.Rho #10 ** self.mesa_table['logRho']

    def initialize_R(self, R, R_in_Rsun):
        self.R = np.array(R) if R is not None else None
    
        # Some checks to warn user if redundant keywords have been specified
        if (self.R is not None) and (R_in_Rsun is not None):
            warnings.warn('Both R and R_in_Rsun specified: defaulting to former')

        self.R = const.Rsun * np.array(R_in_Rsun) if self.R is None and R_in_Rsun is not None else self.R

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
        self.N = 1e-6 * np.array(N_in_uHz) if self.N is None and N_in_uHz is not None else self.N

        if hasattr(self.mesa_table, 'N'):
            if self.N is not None:
                warnings.warn('N already explicitly defined, so not using column in mesa_table')
            else:
                self.N = self.mesa_table.N

        # Required field
        assert self.N is not None

        # Clip zero N to a finite value, so zero-division doesn't happen
        self.N = np.clip(self.N, a_min=N_floor, a_max=None)

        # Calculate other helper quantities
        self.N_div_2pi = self.N / (2 * const.pi)
        self.N_in_uHz = 1e6 * self.N
        self.N_div_2pi_in_uHz = 1e6 * self.N_div_2pi

    def initialize_Sl(self, Sl1, Sl1_in_uHz):
        self.Sl1 = np.array(Sl1) if Sl1 is not None else None

        # Some checks to warn user if redundant keywords have been specified
        if (self.Sl1 is not None) and (Sl1_in_uHz is not None):
            warnings.warn('Both Sl1 and Sl1_in_uHz specified: defaulting to former')

        # Alternative inputs
        self.Sl1 = 1e-6 * np.array(Sl1_in_uHz) if self.Sl1 is None and Sl1_in_uHz is not None else self.Sl1

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

    def calculate_min_max_N_Sl(self):
        """
        Method which calculates the pointwise min and max between N and Sl
        """
        self.min_N_Sl1 = np.min([self.N, self.Sl1], axis=0)
        self.min_N_Sl2 = np.min([self.N, self.Sl2], axis=0)
        self.min_N_Sl = lambda l: np.min([self.N, self.Sl(l)], axis=0)
        self.max_N_Sl1 = np.max([self.N, self.Sl1], axis=0)
        self.max_N_Sl2 = np.max([self.N, self.Sl2], axis=0)
        self.max_N_Sl = lambda l: np.max([self.N, self.Sl(l)], axis=0)

        self.min_N_Sl1_div_2pi = np.min([self.N_div_2pi, self.Sl1_div_2pi], axis=0)
        self.min_N_Sl2_div_2pi = np.min([self.N_div_2pi, self.Sl2_div_2pi], axis=0)
        self.min_N_Sl_div_2pi = lambda l: np.min([self.N_div_2pi, self.Sl_div_2pi(l)], axis=0)
        self.max_N_Sl1_div_2pi = np.max([self.N_div_2pi, self.Sl1_div_2pi], axis=0)
        self.max_N_Sl2_div_2pi = np.max([self.N_div_2pi, self.Sl2_div_2pi], axis=0)
        self.max_N_Sl_div_2pi = lambda l: np.max([self.N_div_2pi, self.Sl_div_2pi(l)], axis=0)

        self.min_N_Sl1_in_uHz = np.min([self.N_in_uHz, self.Sl1_in_uHz], axis=0)
        self.min_N_Sl2_in_uHz = np.min([self.N_in_uHz, self.Sl2_in_uHz], axis=0)
        self.min_N_Sl_in_uHz = lambda l: np.min([self.N_in_uHz, self.Sl_in_uHz(l)], axis=0)
        self.max_N_Sl1_in_uHz = np.max([self.N_in_uHz, self.Sl1_in_uHz], axis=0)
        self.max_N_Sl2_in_uHz = np.max([self.N_in_uHz, self.Sl2_in_uHz], axis=0)
        self.max_N_Sl_in_uHz = lambda l: np.max([self.N_in_uHz, self.Sl_in_uHz(l)], axis=0)

        self.min_N_Sl1_div_2pi_in_uHz = np.min([self.N_div_2pi_in_uHz, self.Sl1_div_2pi_in_uHz], axis=0)
        self.min_N_Sl2_div_2pi_in_uHz = np.min([self.N_div_2pi_in_uHz, self.Sl2_div_2pi_in_uHz], axis=0)
        self.min_N_Sl_div_2pi_in_uHz = lambda l: np.min([self.N_div_2pi_in_uHz, self.Sl_div_2pi_in_uHz(l)], axis=0)
        self.max_N_Sl1_div_2pi_in_uHz = np.max([self.N_div_2pi_in_uHz, self.Sl1_div_2pi_in_uHz], axis=0)
        self.max_N_Sl2_div_2pi_in_uHz = np.max([self.N_div_2pi_in_uHz, self.Sl2_div_2pi_in_uHz], axis=0)
        self.max_N_Sl_div_2pi_in_uHz = lambda l: np.max([self.N_div_2pi_in_uHz, self.Sl_div_2pi_in_uHz(l)], axis=0)

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
        self.Br = 1e3 * np.array(Br_kG) if self.Br is None and Br_kG is not None else self.Br
        self.Br = 1e6 * np.array(Br_MG) if self.Br is None and Br_MG is not None else self.Br

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
        self.Ωrot = 2 * np.pi * νrot if self.Ωrot is None and νrot is not None else self.Ωrot
        self.Ωrot = 1e-6 * Ωrot_uHz if self.Ωrot is None and Ωrot_uHz is not None else self.Ωrot
        self.Ωrot = 2 * np.pi * 1e-6 * νrot_uHz if self.Ωrot is None and νrot_uHz is not None else self.Ωrot
        self.Ωrot = 2 * np.pi / Prot if self.Ωrot is None and Prot is not None else self.Ωrot
        self.Ωrot = 2 * np.pi / (const.secday * Prot_d) if self.Ωrot is None and Prot_d is not None else self.Ωrot


