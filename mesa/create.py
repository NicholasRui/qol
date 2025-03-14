from qol.mesa.read import read_mod
import qol.mesa.const as const
from qol.mesa.data import MesaTable

import qol.paths as paths
import qol.config as config
import qol.helper.formatter as formatter

from astropy.table import vstack

import numpy as np
from scipy.integrate import cumulative_trapezoid

import os; import shutil

import warnings

class MesaInlistControl:
    """
    Simple class which stores a single MESA control
    """
    def __init__(self, section, control, value, category=None, comment=None):
        """
        section: e.g., star_job, controls, etc.
        control: name of inlist option
        value: value of option
        category: optional -- allows grouping of options together under comment string given by this input
        comment: optional -- add trailing comment
        """
        self.section = section
        self.control = control
        self.value = value
        self.category = category
        self.comment = comment

    def inlist_string(self):
        """
        Return the string which should appear in inlist
        """
        fortran_value = formatter.to_fortran(self.value)
        
        if self.comment is not None:
            comment_str = f' ! {self.comment}'
        else:
            comment_str = None

        return f'{self.control} = {fortran_value}{comment_str}'




class MesaInlistProjectFile:
    """
    Stores inlist information for MESA 
    """
    def __init__(self, inlist_name='inlist_project', mesa_version=config.mesa_version):
        self.inlist_name = inlist_name
        self.mesa_version = mesa_version

        self.inlist_controls = []

        # defaults
        self.make_pre_ms = False


        # TODO
        # - disable specific nuclear reactions (or groups, like pp chain)
        # - save gyre, model with profile, etc.
        # - mass_change? max or min star etc..
        # - add easy function for making a wind
        # - define metal fractions

        # - disable burning, or set max_abar_for_burning
        # - saving an inlist should add comments at the top about how this was generated, and what version
        # - terminate at max age
        # - limit_for_rel_error_in_energy_conservation


        # - writeout interval

        # - pgstar enable

        # - net!!
        # - opacities
        # - overshoot
        # - disable mixing
        # - 


    def add_control(self, section, control, value, category=None, skip_if_None=False):
        """
        skip_if_None: if True and value is None, do nothing (for optional keywords)
        """
        if skip_if_None and value is None:
            return

        inlist_control = MesaInlistControl(section=section, 
                control=control,
                value=value,
                category=value)
        
        self.inlist_controls.append(inlist_control)



    #############################################
    ###### Preset combinations of controls ######
    #############################################
    def load_model(self, fname):
        section = 'star_job'
        category = 'load initial model'

        self.add_control(section=section, category=category,
                control='load_saved_model', value=True)
        self.add_control(section=section, category=category,
                control='load_model_filename', value=fname)
        
    def save_final_model(self, fname):
        section = 'star_job'
        category = 'save final model'

        self.add_control(section=section, category=category,
                control='save_model_when_terminate', value=True)
        self.add_control(section=section, category=category,
                control='save_model_filename', value=fname)

    def enable_hydrodynamics(self):
        section = 'star_job'
        category = 'enable hydrodynamics'

        self.add_control(section=section, category=category,
                control='change_v_flag', value=True)
        self.add_control(section=section, category=category,
                control='change_initial_v_flag', value=True)
        self.add_control(section=section, category=category,
                control='new_v_flag', value=True)

    def set_Zbase(self, Zbase=0.02):
        section = 'kap'
        
        self.add_control(section=section,
                control='use_Type2_opacities', value=True)
        self.add_control(section=section,
                control='Zbase', value=Zbase)

    def change_net(self, net_name):
        section = 'star_job'
        category = 'reaction network'

        self.add_control(section=section, category=category,
                control='change_initial_net', value=True)
        self.add_control(section=section, category=category,
                control='new_net_name', value=net_name)

    def relax_to_inner_BC(self, M_new_Msun=None, R_center_Rsun=None, L_center_Lsun=None,
            dlgm_per_step=None, dlgR_per_step=None, dlgL_per_step=None,
            relax_M_center_dt=None, relax_R_center_dt=None, relax_L_center_dt=None):
        """
        Relax to inner boundary condition
        """
        section = 'star_job'

        if M_new_Msun is not None:
            category = 'relax to inner BC: M_center'
            comment = 'total mass in Msun'

            self.add_control(section=section, category=category,
                    control='relax_M_center', value=True)
            self.add_control(section=section, category=category,
                    control='new_mass', value=M_new_Msun, comment=comment)
            
            self.add_control(section=section, category=category, skip_if_None=True,
                    control='dlgm_per_step', value=dlgm_per_step)
            self.add_control(section=section, category=category, comment='sec', skip_if_None=True,
                    control='relax_M_center_dt', value=relax_M_center_dt)

        if R_center_Rsun is not None:
            category = 'relax to inner BC: R_center'
            R_center_cgs = R_center_Rsun * const.Rsun
            comment = f'cm = {formatter.to_fortran(R_center_Rsun)} ! Rsun'

            self.add_control(section=section, category=category,
                    control='relax_R_center', value=True)
            self.add_control(section=section, category=category,
                    control='new_R_center', value=R_center_cgs, comment=comment)
            
            self.add_control(section=section, category=category, skip_if_None=True,
                    control='dlgR_per_step', value=dlgR_per_step)
            self.add_control(section=section, category=category, comment='sec', skip_if_None=True,
                    control='relax_R_center_dt', value=relax_R_center_dt)    
        
        if L_center_Lsun is not None:
            category = 'relax to inner BC: L_center'
            L_center_cgs = L_center_Lsun * const.Lsun
            comment = f'erg s-1 = {formatter.to_fortran(L_center_Lsun)} ! Lsun'

            self.add_control(section=section, category=category,
                    control='relax_L_center', value=True)
            self.add_control(section=section, category=category,
                    control='new_L_center', value=L_center_cgs, comment=comment)
            
            self.add_control(section=section, category=category, skip_if_None=True,
                    control='dlgL_per_step', value=dlgL_per_step)
            self.add_control(section=section, category=category, comment='sec', skip_if_None=True,
                    control='relax_L_center_dt', value=relax_L_center_dt)

    def drag_for_HSE(self, drag_coefficient, use_drag_energy=False):
        """
        For hydrodynamical mode, add extra drag coefficient
        in order to relax model into HSE
        """
        section = 'controls'
        category = 'artificial drag to relax into HSE'

        self.add_control(section=section, category=category,
                control='use_drag_energy', value=use_drag_energy)
        self.add_control(section=section, category=category,
                control='drag_coefficient', value=drag_coefficient)

    def set_min_D_mix(self, min_D_mix,
            mass_lower_limit_for_min_D_mix=None, mass_upper_limit_for_min_D_mix=None):
        section = 'controls'
        category = 'mixing: minimum mixing coefficient'

        self.add_control(section=section, category=category,
                control='set_min_D_mix', value=True)
        self.add_control(section=section, category=category,
                control='min_D_mix', value=min_D_mix)
        
        self.add_control(section=section, category=category, skip_if_None=True,
                control='mass_lower_limit_for_min_D_mix', value=mass_lower_limit_for_min_D_mix)
        self.add_control(section=section, category=category, skip_if_None=True,
                control='mass_upper_limit_for_min_D_mix', value=mass_upper_limit_for_min_D_mix)

    def set_max_num_retries(self, value):
        self.add_control(section='controls', category='retry limit',
            control='max_number_retries', value=value)
        self.add_control(section='controls', category='retry limit',
            control='relax_max_number_retries', value=value)



    ###########################################
    ###### Shortcuts for common controls ######
    ###########################################
    def create_pre_main_sequence_model(self, value=True):
        self.add_control(section='star_job', category='create initial model',
                control='create_pre_main_sequence_model', value=value)
        
        self.make_pre_ms = True

    def use_gold_tolerances(self, value=True):
        self.add_control(section='controls', category='solver',
                control='use_gold_tolerances', value=value)
    
    def energy_eqn_option(self, value='dedt'):
        """
        either dedt or eps_grav
        """
        self.add_control(section='controls', category='solver',
                control='energy_eqn_option', value=value)

    def convergence_ignore_equL_residuals(self, value=True):
        self.add_control(section='controls', category='solver',
                control='convergence_ignore_equL_residuals', value=value)

    def mesh_delta_coeff(self, value):
        self.add_control(section='controls', category='resolution',
                control='mesh_delta_coeff', value=value)
    
    def min_timestep_limit(self, value):
        self.add_control(section='controls', category='timestepping',
                control='min_timestep_limit', value=value)

    def initial_y(self, value):
        self.add_control(section='controls', category='initial composition',
                control='initial_y', value=value)
        
        if not self.make_pre_ms:
            warnings.warn('not used yet, need to make model from scratch (e.g., pre-MS)')
    
    def initial_z(self, value):
        self.add_control(section='controls', category='initial composition',
                control='initial_z', value=value)


class MesaWorkDirectory:
    """
    Stores information for creating a custom MESA work directory
    """
    def __init__(self):
        ...
        # - inlist pgstar

        # - save shell for running, according to template



    def save_directory(path):
        """
        Create MESA directory
        """
        # copy working directory


        ...















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
    
