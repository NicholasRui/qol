from qol.mesa.read import read_mod
import qol.mesa.const as const
from qol.mesa.data import MesaTable

from astropy.table import vstack

import numpy as np
from scipy.integrate import cumulative_trapezoid

import warnings


def create_shell_burning_remnant(write_mod_fname, core_mod_fname, env_mod_fname,
                                 interface_setting=None,
                                 readjust_setting=None):
    """
    Create and write model file of a shell-burning merger remnant,
    by stitching together a core and envelope model (.mod file).

    interface_setting: what to do with the low-density mass in the outskirts of the core model?
    readjust_setting: what to recalculate about the profile (especially if interface_setting is used)
    """
    # Read models and grab tables (astropy.table.vstack is not supported for MesaTable)
    core_model = read_mod(core_mod_fname)
    env_model = read_mod(env_mod_fname)

    core_tab = core_model.tab
    env_tab = env_model.tab

    # Enforce that some important quantities match between model files
    assert core_model.attr['version_number'] == env_model.attr['version_number']
    assert core_model.attr['model_setting'] == env_model.attr['model_setting']
    assert core_model.attr['net_name'] == env_model.attr['net_name']
    assert core_model.attr['species'] == env_model.attr['species']

    version_number = core_model.attr['version_number']
    model_setting = core_model.attr['model_setting']
    net_name = core_model.attr['net_name']
    species = core_model.attr['species']
    initial_z = core_model.attr['initial_z']

    # Handle low-density interface material
    lnRho_env_base = np.max(env_tab['lnd'])
    lnT_env_base = np.max(env_tab['lnT'])
    Rho_too_low = (core_tab['lnd'] < lnRho_env_base)
    match interface_setting:
        case 'excise': # cut out lower-density material
            core_tab = core_tab[~Rho_too_low]

        case 'bulkup': # artificially increase density/temperature of lower-density material
            core_tab['lnd'][Rho_too_low] = lnRho_env_base
            core_tab['lnd'][Rho_too_low] = lnT_env_base
        
        case _:
            warnings.warn(f'No interface_setting = {interface_setting}', category=UserWarning)


    # Recalculate dq (weight by contribution from each object)
    # Note: dq means the mass in simulated matter, so the quotient doesn't include M_center
    core_M_in_Msun = core_model.attr['M/Msun']
    env_M_in_Msun = env_model.attr['xmstar'] / const.Msun # mass in simulated matter
    total_M_in_Msun = env_model.attr['M/Msun'] # outer mass coordinate, including inner boundary condition

    core_tab['dq'] *= core_M_in_Msun / total_M_in_Msun
    env_tab['dq'] *= env_M_in_Msun / total_M_in_Msun

    # Stack tables together
    full_tab = vstack([env_tab, core_tab])

    # if monotonic_Rho: # force density to be monotonic, if desired
    #     full_tab['lnd'] = np.maximum.accumulate(full_tab['lnd'])

    # d_Rho_R3 = 3 * const.Msun * total_M_in_Msun * dq / (4 * const.pi)
    # Rho_R3 = np.flip(np.cumsum(np.flip(d_Rho_R3)))
    # lnR = (np.log(Rho_R3) - full_tab['lnd']) / 3

    match readjust_setting:
        case 'change_R': # recalculate radii
            Rho = const.eulernum ** full_tab['lnd']
            dq = full_tab['dq']

            dR3 = const.Msun * total_M_in_Msun * dq / (const.four_thirds_pi * Rho)

            R3 = np.flip(np.cumsum(np.flip(dR3)))
            lnR = np.log(R3) / 3

            full_tab['lnR'] = lnR

        case 'change_m': # recalculate masses
            Rho = const.eulernum ** full_tab['lnd']
            R3 = const.eulernum ** (3 * full_tab['lnR'])

            dR3 = -np.diff(np.append(R3, 0))
            dm = const.four_thirds_pi * Rho * dR3

            total_M_in_Msun = np.sum(dm) / const.Msun
            dq = dm / np.sum(dm)
            full_tab['dq'] = dq

        case _:
            warnings.warn(f'No readjust_setting = {readjust_setting}', category=UserWarning)

            
    # full_tab['lnR'][:len(env_model)] = np.log(R[:len(env_model)])

    # Recast as MesaTable object
    attr =  {'M/Msun': total_M_in_Msun,
          'initial_z': initial_z,
      'model_setting': model_setting,
     'version_number': version_number,
           'net_name': net_name,
            'species': species}
    
    full_mod = MesaTable(full_tab, attr=attr, tabtype='model')
    full_mod.write_model(write_mod_fname)

    return full_mod
    
