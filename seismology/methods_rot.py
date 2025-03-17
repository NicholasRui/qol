# Methods for Seismology class
# Calculate basic quantities related to rotation
from qol.seismology.helper import get_ω
import qol.mesa.const as const
from qol.tools.integrate import trapz_cond

import numpy as np



def get_Ωrot_g(self, l, ω=None, ν=None, ω_uHz=None, ν_uHz=None, P=None, units='Ωrot_g'):
    """
    Average of Ωrot over the g-mode cavity, as sensed by a given mode
    
    If Ωrot was given as a float, just return it as-is
    """
    assert self.Ωrot is not None

    if type(self.Ωrot) == float:
        return self.Ωrot
    
    ω = get_ω(ω=ω, ν=ν, ω_uHz=ω_uHz, ν_uHz=ν_uHz, P=P)
    is_g = self.get_is_prop(l=l, ω=ω, proptype='g')

    Ωrot_g = trapz_cond(x=self.R, y=self.N * self.Ωrot / self.R, c=is_g)


    match units:
        case 'Ωrot_g':
            return Ωrot_g
        
        case 'νrot_g':
            νrot_g = Ωrot_g / 2 / np.pi
            return νrot_g
        
        case 'Ωrot_g_uHz':
            Ωrot_g_uHz = 1e6 * Ωrot_g
            return Ωrot_g_uHz
        
        case 'νrot_g_uHz':
            νrot_g_uHz = 1e6 * Ωrot_g / 2 / np.pi
            return νrot_g_uHz
        
        case 'Prot_g':
            Prot_g = 2 * np.pi / Ωrot_g
            return Prot_g

        case 'Prot_g_d':
            Prot_g_d = 2 * np.pi / Ωrot_g / const.secday
            return Prot_g_d
        
        case _:
            raise ValueError(f'invalid units: {units}')

