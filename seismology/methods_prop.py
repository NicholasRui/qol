# Methods for Seismology class
# Calculate basic quantities related to wave propagation
from qol.seismology.helper import get_ω
from qol.tools.integrate import trapz_cond

import numpy as np



def get_is_prop(self, l, ω=None, ν=None, ω_uHz=None, ν_uHz=None, P=None, proptype=None):
    """
    Return array of Booleans which says whether a mode with a given frequency is propagating or not

    if proptype is specified and is either 'p' or 'g', return boolean
        showing only whether a mode would be propagating with that particular character
    """
    ω = get_ω(ω=ω, ν=ν, ω_uHz=ω_uHz, ν_uHz=ν_uHz, P=P)

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

def get_is_g(self, l, ω=None, ν=None, ω_uHz=None, ν_uHz=None, P=None):
    return self.get_is_prop(l, ω=ω, ν=ν, ω_uHz=ω_uHz, ν_uHz=ν_uHz, P=P, proptype='g')

def get_is_p(self, l, ω=None, ν=None, ω_uHz=None, ν_uHz=None, P=None):
    return self.get_is_prop(l, ω=ω, ν=ν, ω_uHz=ω_uHz, ν_uHz=ν_uHz, P=P, proptype='p')

def get_int_N_div_r_dr(self, l, ω=None, ν=None, ω_uHz=None, ν_uHz=None, P=None):
    """
    Calculate int_N_div_r_dr over g-propagating regions (no cavity contiguity enforced)
    """
    ω = get_ω(ω=ω, ν=ν, ω_uHz=ω_uHz, ν_uHz=ν_uHz, P=P)
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