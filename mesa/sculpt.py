from qol.mesa.read import read_mod
import qol.mesa.const as const
from qol.mesa.table.MesaTable import MesaTable

import qol.info as info
import qol.config as config
import qol.tools.formatter as formatter

from qol.mesa.launch.MesaInlist import MesaInlist

import numpy as np
from astropy.table import vstack

import warnings

def create_env_inlist_from_core(run_path, task_name, core_mod_fname, M_env_Msun, subdir='data', absdir=None):
    """
    Given a core model, save an inlist which implements boundary conditions for corresponding envelope
    This is meant to be called when creating an inlist as a prereq, i.e., when M_center and R_center are not known in advance of the run
      but depend on the output of another run. If this is not the case, just use MesaInlist.relax_to_inner_BC by itself.

    note: core_mod_fname should be an ABSOLUTE path
    """
    # Read models
    core_model = read_mod(core_mod_fname)

    # Calculate enclosed radius and mass needed
    M_center_Msun = core_model.M_in_Msun[0]
    R_center_Rsun = core_model.R_in_Rsun[0]
    dlgR_per_step = 0.1 # by default: increase dlgR

    # Create inlist
    inlist = MesaInlist(task_name)
    inlist.relax_to_inner_BC(M_new_Msun=M_env_Msun+M_center_Msun,
                             R_center_Rsun=R_center_Rsun,
                             dlgR_per_step=dlgR_per_step)

    # Save inlist
    inlist.save(run_path, subdir=subdir, absdir=absdir)

def create_shell_burning_remnant(write_mod_fname, core_mod_fname, env_mod_fname,
                                 interface_setting=None, readjust_setting=None):
    """
    Create and write model file of a shell-burning merger remnant,
    by stitching together a core and envelope model (.mod file).

    interface_setting: what to do with the low-density mass in the outskirts of the core model?
    readjust_setting: what to recalculate about the profile (especially if interface_setting is used)
    """
    # Read models and grab tables (astropy.table.vstack is not supported for MesaTable)
    core_model = read_mod(core_mod_fname)
    env_model = read_mod(env_mod_fname)

    # Enforce that some important quantities match between model files
    assert core_model.attr['version_number'] == env_model.attr['version_number']
    assert core_model.attr['model_setting'] == env_model.attr['model_setting']
    assert core_model.attr['net_name'] == env_model.attr['net_name']
    assert core_model.attr['species'] == env_model.attr['species']

    assert 'R_center' in env_model.attr.keys()

    version_number = core_model.attr['version_number']
    model_setting = core_model.attr['model_setting']
    net_name = core_model.attr['net_name']
    species = core_model.attr['species']
    initial_z = core_model.attr['initial_z']

    core_M_in_Msun = core_model.attr['M/Msun']
    env_M_in_Msun = env_model.attr['xmstar'] / const.Msun # mass in simulated matter
    total_M_in_Msun = env_model.attr['M/Msun'] # outer mass coordinate, including inner boundary condition

    # Handle low-density interface material
    lnRho_env_base = np.max(env_model['lnd'])
    lnT_env_base = np.max(env_model['lnT'])
    Rho_too_low = (core_model['lnd'] < lnRho_env_base)
    match interface_setting:
        case 'default':
            pass

        case 'excise': # cut out lower-density material
            core_model = core_model[~Rho_too_low]

        case 'bulkup': # artificially increase density/temperature of lower-density material
            core_model['lnd'][Rho_too_low] = lnRho_env_base
            core_model['lnd'][Rho_too_low] = lnT_env_base

        case _:
            warnings.warn(f'No interface_setting = {interface_setting}', category=UserWarning)

    # Recalculate dq (weight by contribution from each object)
    # Note: dq means the mass in simulated matter, so the quotient doesn't include M_center
    core_model['dq'] *= core_M_in_Msun / total_M_in_Msun
    env_model['dq'] *= env_M_in_Msun / total_M_in_Msun

    # Stack tables together
    full_tab = vstack([env_model, core_model])

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

    # Recast as MesaTable object
    attr =  {'M/Msun': total_M_in_Msun,
          'initial_z': initial_z,
      'model_setting': model_setting,
     'version_number': version_number,
           'net_name': net_name,
            'species': species}
    
    full_model = MesaTable(full_tab, attr=attr, tabtype='model')
    full_model.write_model(write_mod_fname)

    return full_model

def replace_elements(write_mod_fname, init_mod_fname, new_species, *old_species_arr):
    """
    takes args list old_species_arr, replaces each one with new_species, and writes it out
    """
    init_model = read_mod(init_mod_fname)

    # loop through species and set their abundances to zero
    old_species_X = np.zeros(len(init_model))
    for old_species in old_species_arr:
        old_species_X += init_model[old_species]
        init_model[old_species] = 0.
    
    # add mass fraction to new species
    init_model[new_species] += old_species_X

    init_model.write_model(write_mod_fname)
