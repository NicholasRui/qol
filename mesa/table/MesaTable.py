## Stores class definition for generalized table with MESA properties
import numpy as np
from astropy.table import Table, Column, unique

import qol.mesa.const as const
import qol.info as info
import qol.tools.formatter as formatter

import os

class MesaTable(Table):
    """
    Generalized astropy.table.Table object which adds additional functionality

    __init__ args not shared by astropy.table.Table:
    attr: dict of other information
    tabtype: type of MESA data, either 'history', 'profile', 'index', 'model'
    """
    def __init__(self, data=None, attr={}, tabtype=None,
                 masked=False, names=None, dtype=None, meta=None, copy=True, rows=None, copy_indices=True, units=None, descriptions=None, **kwargs):
        super().__init__(data=data, masked=masked, names=names, dtype=dtype, meta=meta, copy=copy, rows=rows, copy_indices=copy_indices, units=units, descriptions=descriptions, **kwargs)

        self.attr = attr
        self.tabtype = tabtype
        self.update_attributes()
        # rewrite n_shells just in case
        if self.tabtype == 'model':
            self.attr['n_shells'] = len(self)

    def update_attributes(self):
        # save some additional attributes, depending on tabtype
        match self.tabtype:
            case 'history':
                pass
            case 'profile':
                self.update_attributes_for_profile()
            case 'index':
                pass
            case 'model':
                self.update_attributes_for_model()

    def update_attributes_for_history(self):
        # TODO: add other useful stuff

        if 'star_age' in self.colnames:
            self.star_age = self['star_age']
            self.star_age_in_kyr = 1e-3 * self['star_age']
            self.star_age_in_Myr = 1e-6 * self['star_age']
            self.star_age_in_Gyr = 1e-9 * self['star_age']

        if 'age' in self.colnames:
            self.age = self['age']
            self.age_in_kyr = 1e-3 * self['age']
            self.age_in_Myr = 1e-6 * self['age']
            self.age_in_Gyr = 1e-9 * self['age']

        if 'delta_nu' in self.colnames:
            self.delta_nu_in_uHz = self['delta_nu']
            self.delta_nu = 1e-6 * self.delta_nu_in_uHz

        if 'delta_Pg' in self.colnames:
            self.delta_Pg = self['delta_Pg']

        if 'nu_max' in self.colnames:
            self.nu_max_in_uHz = self['nu_max']
            self.nu_max = 1e-6 * self.nu_max_in_uHz

    def update_attributes_for_profile(self):
        if 'mass' in self.colnames:
            self.mass = self.M_in_Msun = self['mass']
            self.dM_in_Msun = -np.diff(self.M_in_Msun, append=0)

            self.M = const.Msun * self.M_in_Msun
            self.dM = const.Msun * self.dM_in_Msun
            self.dq = self.dM / self.M[0]

        if 'logR' in self.colnames:
            self.logR_in_Rsun = self['logR']
            self.lnR = const.ln10 * self.logR_in_Rsun + np.log(const.Rsun)

            self.R_in_Rsun = 10 ** self.logR_in_Rsun
            self.dR_in_Rsun = -np.diff(self.R_in_Rsun, append=0)

            self.R = const.Rsun * self.R_in_Rsun
            self.dR = const.Rsun * self.dR_in_Rsun

        if 'logT' in self.colnames:
            self.logT = self['logT']
            self.lnT = const.ln10 * self.logT
            self.T = 10 ** self.logT
        
        if 'logRho' in self.colnames:
            self.logRho = self['logRho']
            self.lnRho = self.lnd = const.ln10 * self.logRho
            self.Rho = 10 ** self.logRho
        
        if 'logP' in self.colnames:
            self.logP = self['logP']
            self.lnP = const.ln10 * self.logP
            self.P = 10 ** self.logP

        if 'brunt_N2' in self.colnames:
            self.N2 = self['brunt_N2']
            self.N = np.sqrt(np.clip(self.N2, a_min=0, a_max=None)) # Brunt zeroed out at imaginary values

            self.N_div_2pi = self.N / (2 * const.pi)
            self.N_in_uHz = 1e6 * self.N
            self.N_div_2pi_in_uHz = 1e6 * self.N_div_2pi

        if 'brunt_N2_structure_term' in self.colnames:
            self.N2_structure_term = self['brunt_N2_structure_term']
            self.N_structure_term = np.sqrt(np.clip(self.N2_structure_term, a_min=0, a_max=None)) # Brunt zeroed out at imaginary values

            self.N_structure_term_div_2pi = self.N_structure_term / (2 * const.pi)
            self.N_structure_term_in_uHz = 1e6 * self.N_structure_term
            self.N_structure_term_div_2pi_in_uHz = 1e6 * self.N_structure_term_div_2pi

        if 'brunt_N2_composition_term' in self.colnames:
            self.N2_composition_term = self['brunt_N2_composition_term']
            self.N_composition_term = np.sqrt(np.clip(self.N2_composition_term, a_min=0, a_max=None)) # Brunt zeroed out at imaginary values

            self.N_composition_term_div_2pi = self.N_composition_term / (2 * const.pi)
            self.N_composition_term_in_uHz = 1e6 * self.N_composition_term
            self.N_composition_term_div_2pi_in_uHz = 1e6 * self.N_composition_term_div_2pi
        
        if 'csound' in self.colnames:
            self.csound = self['csound']

        have_S = False
        if 'lamb_S2' in self.colnames:
            self.Sl1_2 = self['lamb_S2']
            self.Sl1 = np.sqrt(self.Sl1_2)

            have_S = True
        elif 'lamb_S' in self.colnames: # MESA writes this in Hz by default
            self.Sl1 = self['lamb_S']
            self.Sl1_2 = self.Sl1 ** 2

            have_S = True
        elif 'lamb_Sl1' in self.colnames: # MESA writes this in uHz by default, so need to convert
            self.Sl1 = 1e-6 * self['lamb_Sl1']
            self.Sl1_2 = self.Sl1 ** 2

            have_S = True
        elif ('csound' in self.colnames) and ('logR' in self.colnames):
            self.Sl1 = const.sqrt2 * self.csound / self.R
            self.Sl1_2 = self.Sl1 ** 2

            have_S = True

        if have_S: # if have Lamb frequency somehow, calculate some more stuff
            self.Sl2 = np.sqrt(3) * self.Sl1
            self.Sl3 = np.sqrt(6) * self.Sl1

            self.Sl1_div_2pi = self.Sl1 / (2 * const.pi)
            self.Sl2_div_2pi = self.Sl2 / (2 * const.pi)
            self.Sl3_div_2pi = self.Sl3 / (2 * const.pi)

            self.Sl1_in_uHz = 1e6 * self.Sl1
            self.Sl2_in_uHz = 1e6 * self.Sl2
            self.Sl3_in_uHz = 1e6 * self.Sl3

            self.Sl1_div_2pi_in_uHz = 1e6 * self.Sl1_div_2pi
            self.Sl2_div_2pi_in_uHz = 1e6 * self.Sl2_div_2pi
            self.Sl3_div_2pi_in_uHz = 1e6 * self.Sl3_div_2pi

    def update_attributes_for_model(self):
        if 'lnd' in self.colnames:
            self.lnd = self.lnRho = self['lnd']
            self.logRho = self.lnd / const.ln10
            self.Rho = const.eulernum ** self.lnd
        
        if 'lnT' in self.colnames:
            self.lnT = self['lnT']
            self.logT = self.lnT / const.ln10
            self.T = const.eulernum ** self.lnT
        
        if 'lnR' in self.colnames:
            self.lnR = self['lnR']
            self.logR_in_Rsun = (self.lnR - np.log(const.Rsun)) / const.ln10

            self.R_in_Rsun = 10 ** self.logR_in_Rsun
            self.dR_in_Rsun = -np.diff(self.R_in_Rsun, append=0)

            self.R = const.Rsun * self.R_in_Rsun
            self.dR = const.Rsun * self.dR_in_Rsun
        
        if 'L' in self.colnames:
            self.L = self['L']
        
        if 'dq' in self.colnames:
            self.dq = dq = self['dq']
            if 'M/Msun' in self.attr.keys():
                M = self.attr['M/Msun']

                # need to account correctly if xmstar is nonzero
                if 'xmstar' in self.attr.keys():
                    self.M_center = M_center = M - self.attr['xmstar'] / const.Msun
                else:
                    self.M_center = M_center = 0

                self.dM_in_Msun = dM_in_Msun = (M - M_center) * dq
                self.M_in_Msun = np.flip(np.cumsum(np.flip(dM_in_Msun))) + M_center

                self.dM = self.dM_in_Msun / const.Msun
                self.M = self.M_in_Msun / const.Msun
        
        if 'conv_vel' in self.colnames:
            self.conv_vel = self['conv_vel']
        if 'mlt_vc' in self.colnames:
            self.mlt_vc = self['mlt_vc']

    def write_model(self, fname, model_setting=36, model_number=1, star_age=0., time=0.):
        """
        model_setting is the integer option at the top of the model file
        """
        assert self.tabtype == 'model' # for now, only writing something which was already a model is supported

        # Update attributes, in case tab attribute has been modified by user
        self.update_attributes()

        match model_setting:
            case 36:
                template_path = os.path.join(info.qol_path, 'mesa/templates/model/model36.mod')

                version_number = self.attr['version_number']
                M_in_Msun = self.attr['M/Msun']
                initial_z = self.attr['initial_z']
                n_shells = previous_n_shells = len(self)
                net_name = self.attr['net_name']
                species = self.attr['species']
                previous_mass_grams = const.Msun * M_in_Msun
                timestep_seconds = dt_next_seconds = 100. * const.secyer # arbitrary

                format_attribute = lambda name, val: f'{name:>32}   {formatter.to_fortran(val)}' '\n'

                attributes_text = ''
                attributes_text += format_attribute('version_number', version_number)
                attributes_text += format_attribute('M/Msun', M_in_Msun)
                attributes_text += format_attribute('initial_z', initial_z)
                attributes_text += format_attribute('model_number', model_number)
                attributes_text += format_attribute('star_age', star_age)
                attributes_text += format_attribute('n_shells', n_shells)
                attributes_text += format_attribute('net_name', net_name)
                attributes_text += format_attribute('species', species)
                attributes_text += format_attribute('time', time)

                # Add some other important attributes if they exist
                if 'xmstar' in self.attr.keys():
                    xmstar = self.attr['xmstar']
                    attributes_text += format_attribute('xmstar', xmstar)
                if 'R_center' in self.attr.keys():
                    R_center = self.attr['R_center']
                    attributes_text += format_attribute('R_center', R_center)

                with open(template_path, 'r') as f:
                    text = f.read()
                
                text = text.replace('<<ATTRIBUTES_TEXT>>', attributes_text)
                text = text.replace('<<PREVIOUS_N_SHELLS>>', formatter.to_fortran(previous_n_shells))
                text = text.replace('<<PREVIOUS_MASS_GRAMS>>', formatter.to_fortran(previous_mass_grams))
                text = text.replace('<<TIMESTEP_SECONDS>>', formatter.to_fortran(timestep_seconds))
                text = text.replace('<<DT_NEXT_SECONDS>>', formatter.to_fortran(dt_next_seconds))

                # Format and write table
                model_table = ' ' + ' '.join(self.colnames) + ' \n'
                model_table += ''.join(f' {ii+1} ' + ' '.join([formatter.to_fortran(self[colname][ii]) \
                                                for colname in self.colnames]) + ' \n' \
                                                for ii in range(len(self)))
                text = text.replace('<<MODEL_TABLE>>', model_table)

            case _:
                raise ValueError(f'model_setting={model_setting} is not supported yet')

        # Write model file
        with open(fname, 'w') as f:
            f.write(text)

    def get_he_core_mass(self, he_core_boundary_h1_fraction=0.1, min_boundary_fraction=0.1):
        """
        MESA documentation says:
        If fraction >= 0, boundary is outermost location where h1 mass fraction is <= this value,
        and he4 mass fraction >= min_boundary_fraction (see below). If fraction < 0, boundary is
        outermost location where he4 is the most abundant species.

        (just handle the first case for now)
        """
        assert self.tabtype in ['profile', 'model']
        assert hasattr(self, 'M_in_Msun')

        bad_frac_thresh = he_core_boundary_h1_fraction
        assert np.isin(['h1', 'he4'], self.colnames).all()
        bad_frac, good_frac = self['h1'], self['he4']

        # print(type(bad_frac), type(bad_frac_thresh), type(good_frac), type(min_boundary_fraction))

        cond = (bad_frac <= bad_frac_thresh) & (good_frac >= min_boundary_fraction)
        if cond.any():
            k = np.where(cond)[0][0] # get qualifying index
            core_mass = float(self.M_in_Msun[k])
            return core_mass
        else:
            return 0.

    def get_co_core_mass(self, co_core_boundary_he4_fraction=0.1, min_boundary_fraction=0.1):
        """
        MESA documentation says:
        If fraction >= 0, boundary is outermost location where he4 mass fraction is <= this value,
        and c12+o16 mass fraction >= min_boundary_fraction (see below). If fraction < 0, boundary is
        outermost location where c12+o16 is more abundant than any other species.
        """
        assert self.tabtype in ['profile', 'model']
        assert hasattr(self, 'M_in_Msun')

        bad_frac_thresh = co_core_boundary_he4_fraction
        assert np.isin(['he4', 'c12', 'o16'], self.colnames).all()
        bad_frac, good_frac = self['he4'], self['c12']+self['o16']

        cond = (bad_frac <= bad_frac_thresh) & (good_frac >= min_boundary_fraction)
        if cond.any():
            k = np.where(cond)[0][0] # get qualifying index
            core_mass = float(self.M_in_Msun[k])
            return core_mass
        else:
            return 0.

    def get_one_core_mass(self, one_core_boundary_he4_c12_fraction=0.1, min_boundary_fraction=0.1):
        """
        MESA documentation says:
        If fraction >= 0, boundary is outermost location where combine he4+c12 mass fraction is <= this value,
        and o16+ne20 mass fraction >= min_boundary_fraction (see below). If fraction < 0, boundary is
        outermost location where o16+ne20 is more abundant than any other species.
        """
        assert self.tabtype in ['profile', 'model']
        assert hasattr(self, 'M_in_Msun')

        bad_frac_thresh = one_core_boundary_he4_c12_fraction
        assert np.isin(['he4', 'c12', 'o16', 'ne20'], self.colnames).all()
        bad_frac, good_frac = self['he4']+self['c12'], self['o16']+self['ne20']

        cond = (bad_frac <= bad_frac_thresh) & (good_frac >= min_boundary_fraction)
        if cond.any():
            k = np.where(cond)[0][0] # get qualifying index
            core_mass = float(self.M_in_Msun[k])
            return core_mass
        else:
            return 0.

    def __getitem__(self, item):
        """
        For model tables, redefine getitem (i.e., table[:5]) to keep attributes
        """
        result = super().__getitem__(item)

        # only if the output is a MesaTable (vs. Row, Column)
        if isinstance(result, MesaTable) and self.tabtype == 'model':
            result.attr = self.attr
            result.tabtype = 'model'
            result.update_attributes()
            result.attr['n_shells'] = len(result)

        return result






