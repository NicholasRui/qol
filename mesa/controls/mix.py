# Methods handling controls about mixing

import re

def set_min_D_mix(self, min_D_mix,
        mass_lower_limit_for_min_D_mix=None, mass_upper_limit_for_min_D_mix=None):
    namelist = 'controls'
    category = 'mixing: minimum mixing coefficient'

    self.add_control(namelist=namelist, category=category,
            control='set_min_D_mix', value=True)
    self.add_control(namelist=namelist, category=category,
            control='min_D_mix', value=min_D_mix)
    
    self.add_control(namelist=namelist, category=category, optional=True,
            control='mass_lower_limit_for_min_D_mix', value=mass_lower_limit_for_min_D_mix)
    self.add_control(namelist=namelist, category=category, optional=True,
            control='mass_upper_limit_for_min_D_mix', value=mass_upper_limit_for_min_D_mix)
    
def disable_mixing(self):
    self.add_control(namelist='controls', category='mixing: disabled',
            control='mix_factor', value=0.)

def use_Ledoux_criterion(self):
    self.add_control(namelist='controls', category='mixing: composition',
            control='use_Ledoux_criterion', value=True)

def alpha_semiconvection(self, value):
    self.add_control(namelist='controls', category='mixing: composition',
            control='alpha_semiconvection', value=value)

def thermohaline_coeff(self, value):
    self.add_control(namelist='controls', category='mixing: composition',
            control='thermohaline_coeff', value=value)

def gravitational_settling(self, diffusion_class_representatives,
        diffusion_use_cgs_solver=True,
        show_diffusion_info=True,
        diffusion_steps_hard_limit=2000,
        diffusion_maxsteps_for_isolve=2000):
    """
    diffusion_class_representatives: list of species names

    assume that diffusion_class_A_max() is the atomic mass of the species, except for the last which is set to 10000
    """
    category = 'mixing: gravitational settling'

    self.add_control(namelist='controls', category=category,
            control='do_element_diffusion', value=True)
    
    self.add_control(namelist='controls', category=category,
            control='diffusion_use_cgs_solver', value=diffusion_use_cgs_solver)
    self.add_control(namelist='controls', category=category,
            control='show_diffusion_info', value=show_diffusion_info)
    self.add_control(namelist='controls', category=category,
            control='diffusion_steps_hard_limit', value=diffusion_steps_hard_limit)
    self.add_control(namelist='controls', category=category,
            control='diffusion_maxsteps_for_isolve', value=diffusion_maxsteps_for_isolve)

    self.add_control(namelist='controls', category=category,
            control='diffusion_num_classes', value=len(diffusion_class_representatives))
    
    for ii, rep in enumerate(diffusion_class_representatives):
        A = int(re.sub(r'\D', '', rep))
        if ii+1 == len(diffusion_class_representatives):
            A = 10000

        self.add_control(namelist='controls', category=category,
                control=f'diffusion_class_representative({ii+1})', value=rep)
        self.add_control(namelist='controls', category=category,
                control=f'diffusion_class_A_max({ii+1})', value=A)

