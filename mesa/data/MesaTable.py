## Stores class definition for generalized table with MESA properties
import numpy as np
from astropy.table import Table, Column, unique

import qol.mesa.const as const
import qol.paths as paths
import qol.helper.formatter as formatter

class MesaTable(Table):
    """
    Generalized astropy.table.Table object which adds additional functionality
    """
    def __init__(self, tab, attr={}, tabtype=None):
        """
        tab: base astropy.table.Table object containing data of table
        attr: dict of other information
        tabtype: type of MESA data, either 'history', 'profile', 'index', 'model'
        """
        self.tab = tab
        self.attr = attr
        self.tabtype = tabtype

        super().__init__(tab)
        self.update_attributes()

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


        if 'delta_nu' in self.keys():
            self.delta_nu_in_uHz = self['delta_nu']
            self.delta_nu = 1e-6 * self.delta_nu_in_uHz

        if 'delta_Pg' in self.keys():
            self.delta_Pg = self['delta_Pg']

        if 'nu_max' in self.keys():
            self.nu_max_in_uHz = self['nu_max']
            self.nu_max = 1e-6 * self.nu_max_in_uHz

    def update_attributes_for_profile(self):
        if 'mass' in self.keys():
            self.mass = self.M_in_Msun = self['mass']
            self.dM_in_Msun = -np.diff(self.M_in_Msun, append=0)

            self.M = const.Msun * self.M_in_Msun
            self.dM = const.Msun * self.dM_in_Msun
            self.dq = self.dM / self.M[0]

        if 'logR' in self.keys():
            self.logR_in_Rsun = self['logR']
            self.lnR = const.ln10 * self.logR_in_Rsun + np.log(const.Rsun)

            self.R_in_Rsun = 10 ** self.logR_in_Rsun
            self.dR_in_Rsun = -np.diff(self.R_in_Rsun, append=0)

            self.R = const.Rsun * self.R_in_Rsun
            self.dR = const.Rsun * self.dR_in_Rsun

        if 'logT' in self.keys():
            self.logT = self['logT']
            self.lnT = const.ln10 * self.logT
            self.T = 10 ** self.logT
        
        if 'logRho' in self.keys():
            self.logRho = self['logRho']
            self.lnRho = self.lnd = const.ln10 * self.logRho
            self.Rho = 10 ** self.logRho
        
        if 'logP' in self.keys():
            self.logP = self['logRho']
            self.lnP = const.ln10 * self.logP
            self.P = 10 ** self.logP

        if 'brunt_N2' in self.keys():
            self.N2 = self['brunt_N2']
            self.N = np.clip(self.brunt_N2, a_min=0, a_max=None) # Brunt zeroed out at imaginary values

            self.N_div_2pi = self.N / (2 * const.pi)
            self.N_in_uHz = 1e6 * self.N
            self.N_div_2pi_in_uHz = 1e6 * self.N_div_2pi

        if 'brunt_N2_structure_term' in self.keys():
            self.N2_structure_term = self['brunt_N2_structure_term']
            self.N_structure_term = np.clip(self.brunt_N2_structure_term, a_min=0, a_max=None) # Brunt zeroed out at imaginary values

            self.N_structure_term_div_2pi = self.N_structure_term / (2 * const.pi)
            self.N_structure_term_in_uHz = 1e6 * self.N_structure_term
            self.N_structure_term_div_2pi_in_uHz = 1e6 * self.N_structure_term_div_2pi

        if 'brunt_N2_composition_term' in self.keys():
            self.N2_composition_term = self['brunt_N2_composition_term']
            self.N_composition_term = np.clip(self.brunt_N2_composition_term, a_min=0, a_max=None) # Brunt zeroed out at imaginary values

            self.N_composition_term_div_2pi = self.N_composition_term / (2 * const.pi)
            self.N_composition_term_in_uHz = 1e6 * self.N_composition_term
            self.N_composition_term_div_2pi_in_uHz = 1e6 * self.N_composition_term_div_2pi
        
        if 'csound' in self.keys():
            self.csound = self['csound']

        have_S = False
        if 'lamb_S2' in self.keys():
            self.Sl1_2 = self['lamb_S2']
            self.Sl1 = np.sqrt(self.Sl1_2)

            have_S = True
        elif 'lamb_S': # MESA writes this in Hz by default
            self.Sl1 = self['lamb_S']
            self.Sl1_2 = self.Sl1 ** 2

            have_S = True
        elif 'lamb_Sl1': # MESA writes this in uHz by default, so need to convert
            self.Sl1 = 1e-6 * self['lamb_Sl1']
            self.Sl1_2 = self.Sl1 ** 2

            have_S = True
        elif ('csound' in self.keys()) and ('logR' in self.keys()):
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
        if 'lnd' in self.keys():
            self.lnd = self.lnRho = self['lnd']
            self.logRho = self.lnd / const.ln10
            self.Rho = const.eulernum ** self.lnd
        
        if 'lnT' in self.keys():
            self.lnT = self['lnT']
            self.logT = self.lnT / const.ln10
            self.T = const.eulernum ** self.lnT
        
        if 'lnR' in self.keys():
            self.lnR = self['lnR']
            self.logR_in_Rsun = (self.lnR - np.log(const.Rsun)) / const.ln10

            self.R_in_Rsun = 10 ** self.logR_in_Rsun
            self.dR_in_Rsun = -np.diff(self.R_in_Rsun, append=0)

            self.R = const.Rsun * self.R_in_Rsun
            self.dR = const.Rsun * self.dR_in_Rsun
        
        if 'L' in self.keys():
            self.L = self['L']
        
        if 'dq' in self.keys():
            self.dq = dq = self['dq']
            
            if 'M/Msun' in self.attr.keys():
                M = self.attr['M/Msun']

                self.dM_in_Msun = dM_in_Msun = M * dq
                self.M_in_Msun = np.flip(np.cumsum(np.flip(dM_in_Msun)))

                self.dM = self.dM_in_Msun / const.Msun
                self.M = self.M_in_Msun / const.Msun
        
        if 'conv_vel' in self.keys():
            self.conv_vel = self['conv_vel']
        if 'mlt_vc' in self.keys():
            self.mlt_vc = self['mlt_vc']

    def write_model(self, fname, model_setting=36):
        """
        model_setting is the integer option at the top of the model file
        """
        assert self.tabtype == 'model' # for now, only writing something which was already a model is supported

        # Update attributes, in case tab attribute has been modified by user
        self.update_attributes()

        match model_setting:
            case 36:
                template_path = f'{paths.qol_path}/mesa/templates/model/model36.mod'

                version_number = self.attr['version_number']
                M_in_Msun = self.attr['M/Msun']
                initial_z = self.attr['initial_z']
                model_number = 1
                star_age = 0
                N_shells = previous_N_shells = len(self)
                net_name = self.attr['net_name']
                species = self.attr['species']
                time = 0
                previous_mass_grams = const.Msun * M_in_Msun
                timestep_seconds = dt_next_seconds = 100 * const.secyer # arbitrary

                with open(template_path, 'r') as f:
                    text = f.read()
                
                text = text.replace('<<VERSION_NUMBER>>', formatter.to_fortran(version_number))
                text = text.replace('<<M_IN_MSUN>>', formatter.to_fortran(M_in_Msun))
                text = text.replace('<<INITIAL_Z>>', formatter.to_fortran(initial_z))
                text = text.replace('<<MODEL_NUMBER>>', formatter.to_fortran(model_number))
                text = text.replace('<<STAR_AGE>>', formatter.to_fortran(star_age))
                text = text.replace('<<N_SHELLS>>', formatter.to_fortran(N_shells))
                text = text.replace('<<NET_NAME>>', formatter.to_fortran(net_name))
                text = text.replace('<<SPECIES>>', formatter.to_fortran(species))
                text = text.replace('<<TIME>>', formatter.to_fortran(time))
                text = text.replace('<<PREVIOUS_N_SHELLS>>', formatter.to_fortran(previous_N_shells))
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
