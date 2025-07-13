# Methods handling controls about mixing

import re
import warnings

def add_overshoot_zone(self, overshoot_scheme, overshoot_zone_type,
                       overshoot_zone_loc, overshoot_bdy_loc,
                       overshoot_f, overshoot_f0,
                       overshoot_D0=None, overshoot_Delta0=None,
                       overshoot_mass_full_on=None, overshoot_mass_full_off=None):
    """
    Add an overshooting zone

    Note: The MESA documentation says:
    "The overshooting parameter values corresponding to the first set of matching criteria will be used.
    Therefore, narrower criteria should precede more general ones (i.e have lower array indices)."

    qol is not going to automatically enforced that overshoot is added in the correct order,
    but will follow what the user has specified
    """
    assert overshoot_scheme in ['exponential', 'step', 'other']
    assert overshoot_zone_type in ['burn_H', 'burn_He', 'burn_Z', 'nonburn', 'any']
    assert overshoot_zone_loc in ['core', 'shell', 'any']
    assert overshoot_bdy_loc in ['bottom', 'top', 'any']

    self.num_overshoot_zones += 1
    category = f'mixing: overshoot zone {self.num_overshoot_zones}'

    self.add_to_controls(category=category,
            control=f'overshoot_scheme({self.num_overshoot_zones})', value=overshoot_scheme)
    self.add_to_controls(category=category,
            control=f'overshoot_zone_type({self.num_overshoot_zones})', value=overshoot_zone_type)
    self.add_to_controls(category=category,
            control=f'overshoot_zone_loc({self.num_overshoot_zones})', value=overshoot_zone_loc)
    self.add_to_controls(category=category,
            control=f'overshoot_bdy_loc({self.num_overshoot_zones})', value=overshoot_bdy_loc)
    
    self.add_to_controls(category=category,
            control=f'overshoot_f({self.num_overshoot_zones})', value=overshoot_f)
    self.add_to_controls(category=category,
            control=f'overshoot_f0({self.num_overshoot_zones})', value=overshoot_f0)
    
    self.add_to_controls(category=category, optional=True,
            control=f'overshoot_D0({self.num_overshoot_zones})', value=overshoot_D0)
    self.add_to_controls(category=category, optional=True,
            control=f'overshoot_Delta0({self.num_overshoot_zones})', value=overshoot_Delta0)
    self.add_to_controls(category=category, optional=True,
            control=f'overshoot_mass_full_on({self.num_overshoot_zones})', value=overshoot_mass_full_on)
    self.add_to_controls(category=category, optional=True,
            control=f'overshoot_mass_full_off({self.num_overshoot_zones})', value=overshoot_mass_full_off)

def add_predictive_mix_zone(self, predictive_zone_type, predictive_zone_loc, predictive_bdy_loc,
                            predictive_superad_thresh=None,
                            predictive_bdy_q_min=None, predictive_bdy_q_max=None,
                            predictive_avoid_reversal=None, predictive_limit_ingestion=None, predictive_ingestion_factor=None):
    """
    For predictive_superad_thresh, see MESA doc:
    Threshold for minimum superadiabaticity in the predictive mixing scheme;
    boundary expansion stops when gradr/grada-1 drops below this threshold.
    Default value is usually good for main-sequence evolution;
    for core He-burning, set to 0.005, 0.01 or larger to prevent splitting
    of the core convection zone and/or core breathing pulses.
    """
    assert predictive_zone_type in ['burn_H', 'burn_He', 'burn_Z', 'nonburn', 'any']
    assert predictive_zone_loc in ['core', 'shell', 'any']
    assert predictive_bdy_loc in ['bottom', 'top', 'any']

    self.num_predictive_mix_zones += 1
    category = f'mixing: predictive mixing zone {self.num_overshoot_zones}'

    self.add_to_controls(category=category,
            control=f'predictive_mix({self.num_predictive_mix_zones})', value=True)
    self.add_to_controls(category=category,
            control=f'predictive_zone_type({self.num_predictive_mix_zones})', value=predictive_zone_type)
    self.add_to_controls(category=category,
            control=f'predictive_zone_loc({self.num_predictive_mix_zones})', value=predictive_zone_loc)
    self.add_to_controls(category=category,
            control=f'predictive_bdy_loc({self.num_predictive_mix_zones})', value=predictive_bdy_loc)
    
    self.add_to_controls(category=category, optional=True,
            control=f'predictive_superad_thresh({self.num_predictive_mix_zones})', value=predictive_superad_thresh)
    self.add_to_controls(category=category, optional=True,
            control=f'predictive_bdy_q_min({self.num_predictive_mix_zones})', value=predictive_bdy_q_min)
    self.add_to_controls(category=category, optional=True,
            control=f'predictive_bdy_q_max({self.num_predictive_mix_zones})', value=predictive_bdy_q_max)
    self.add_to_controls(category=category, optional=True,
            control=f'predictive_avoid_reversal({self.num_predictive_mix_zones})', value=predictive_avoid_reversal)
    self.add_to_controls(category=category, optional=True,
            control=f'predictive_limit_ingestion({self.num_predictive_mix_zones})', value=predictive_limit_ingestion)
    self.add_to_controls(category=category, optional=True,
            control=f'predictive_ingestion_factor({self.num_predictive_mix_zones})', value=predictive_ingestion_factor)

def set_min_D_mix(self, min_D_mix,
        mass_lower_limit_for_min_D_mix=None, mass_upper_limit_for_min_D_mix=None):
    category = 'mixing: minimum mixing coefficient'

    self.add_to_controls(category=category,
            control='set_min_D_mix', value=True)
    self.add_to_controls(category=category,
            control='min_D_mix', value=min_D_mix)
    
    self.add_to_controls(category=category, optional=True,
            control='mass_lower_limit_for_min_D_mix', value=mass_lower_limit_for_min_D_mix)
    self.add_to_controls(category=category, optional=True,
            control='mass_upper_limit_for_min_D_mix', value=mass_upper_limit_for_min_D_mix)
    
def disable_mixing(self):
    self.add_to_controls(category='mixing: disabled',
            control='mix_factor', value=0.)

def use_Ledoux_criterion(self):
    self.add_to_controls(category='mixing: composition',
            control='use_Ledoux_criterion', value=True)

def alpha_semiconvection(self, value):
    self.add_to_controls(category='mixing: composition',
            control='alpha_semiconvection', value=value)

def thermohaline_coeff(self, value):
    self.add_to_controls(category='mixing: composition',
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

    self.add_to_controls(category=category,
            control='do_element_diffusion', value=True)
    
    self.add_to_controls(category=category,
            control='diffusion_use_cgs_solver', value=diffusion_use_cgs_solver)
    self.add_to_controls(category=category,
            control='show_diffusion_info', value=show_diffusion_info)
    self.add_to_controls(category=category,
            control='diffusion_steps_hard_limit', value=diffusion_steps_hard_limit)
    self.add_to_controls(category=category,
            control='diffusion_maxsteps_for_isolve', value=diffusion_maxsteps_for_isolve)

    self.add_to_controls(category=category,
            control='diffusion_num_classes', value=len(diffusion_class_representatives))
    
    for ii, rep in enumerate(diffusion_class_representatives):
        A = int(re.sub(r'\D', '', rep))
        if ii+1 == len(diffusion_class_representatives):
            A = 10000

        self.add_to_controls(category=category,
                control=f'diffusion_class_representative({ii+1})', value=rep)
        self.add_to_controls(category=category,
                control=f'diffusion_class_A_max({ii+1})', value=A)

def phase_separation(self, phase_separation_option,
                     do_phase_separation_heating=True, phase_separation_mixing_use_brunt=True,
                     do_WD_sedimentation_heating=True, min_xa_for_WD_sedimentation_heating=None):
    """
    Do phase separation for white dwarf

    phase_separation_option -> 'CO' or 'ONe'
    """
    assert phase_separation_option in ['CO', 'ONe']

    category = 'phase separation'

    self.add_to_controls(category=category,
            control='do_phase_separation', value=True)
    self.add_to_controls(category=category,
            control='phase_separation_option', value=phase_separation_option)
    self.add_to_controls(category=category,
            control='do_phase_separation_heating', value=do_phase_separation_heating)
    self.add_to_controls(category=category,
            control='phase_separation_mixing_use_brunt', value=phase_separation_mixing_use_brunt)
    
    if do_WD_sedimentation_heating:
        self.add_to_controls(category=category, comment='disabled to avoid double-counting with do_WD_sedimentation_heating',
                control='do_diffusion_heating', value=False)
        self.add_to_controls(category=category, optional=True,
                control='min_xa_for_WD_sedimentation_heating', value=min_xa_for_WD_sedimentation_heating)
    elif min_xa_for_WD_sedimentation_heating is not None:
        warnings.warn('since do_WD_sedimentation_heating==False, min_xa_for_WD_sedimentation_heating is unused.')

def okay_to_reduce_gradT_excess(self):
    """
    MLT++ (Paxton et al. 2013)
    """
    self.add_to_controls(category='mixing',
            control='okay_to_reduce_gradT_excess', value=True)
