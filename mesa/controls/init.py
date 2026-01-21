# Methods handling controls about model initialization
import warnings

def create_pre_main_sequence_model(self, value=True):
    self.add_to_star_job(category='create initial model',
            control='create_pre_main_sequence_model', value=value)
        
    self.make_new_model = True

def create_initial_model(self, M_in_Msun, R_in_Rsun):
    self.add_to_star_job(category='create initial model',
            control='create_initial_model', value=True)

    self.add_to_star_job(category='create initial model',
            control='mass_in_gm_for_create_initial_model', value=M_in_Msun * const.Msun)
    self.add_to_star_job(category='create initial model',
            control='radius_in_cm_for_create_initial_model', value=R_in_Rsun * const.Rsun)
    
    self.make_new_model = True

def initial_mass(self, value):
    self.add_to_controls(category='initial mass and composition',
            control='initial_mass', value=value)

def initial_y(self, value):
    self.add_to_controls(category='initial mass and composition',
            control='initial_y', value=value)
        
    if not self.make_new_model:
        warnings.warn('not used yet, need to make model from scratch (e.g., pre-MS)')

def initial_he3(self, value):
    self.add_to_controls(category='initial mass and composition',
            control='initial_he3', value=value)

def initial_z(self, value):
    self.add_to_controls(category='initial mass and composition',
            control='initial_z', value=value)

def relax_initial_mass(self, M_new_Msun, lg_max_abs_mdot=None):
    self.add_to_star_job(category='relax mass',
        control='relax_initial_mass', value=True)
    self.add_to_star_job(category='relax mass',
        control='new_mass', value=M_new_Msun)
        
    self.add_to_star_job(category='relax mass', optional=True,
        control='lg_max_abs_mdot', value=lg_max_abs_mdot)

def relax_initial_mass_to_remove_H_env(self, extra_mass_retained_by_remove_H_env=None):
    self.add_to_star_job(category='relax mass',
        control='relax_initial_mass_to_remove_H_env', value=True)    
    self.add_to_star_job(category='relax mass', optional=True,
        control='extra_mass_retained_by_remove_H_env', value=extra_mass_retained_by_remove_H_env)

def reset_age(self):
    self.add_to_star_job(category='reset age',
        control='set_initial_age', value=True)
    self.add_to_star_job(category='reset age',
        control='initial_age', value=0.)

def reset_model_number(self):
    self.add_to_star_job(category='reset age',
        control='set_initial_model_number', value=True)
    self.add_to_star_job(category='reset age',
        control='initial_model_number', value=0)

def relax_tau_factor(self, relax_to_this_tau_factor,
                     dlogtau_factor=None, relax_tau_factor_after_core_He_burn=None, relax_tau_factor_after_core_C_burn=None,
                     initial_only=True):
    """
    initial_only: if True, only do this at the beginning, and not when doing ./re
    """
    category = 'relax tau factor'

    if initial_only:
        self.add_to_star_job(category=category,
            control='relax_initial_tau_factor', value=True)
    else:
        self.add_to_star_job(category=category,
            control='relax_tau_factor', value=True)

    self.add_to_star_job(category=category,
        control='relax_to_this_tau_factor', value=True)
    
    # optional
    self.add_to_star_job(category=category, optional=True,
        control='dlogtau_factor', value=dlogtau_factor)
    self.add_to_star_job(category=category, optional=True,
        control='relax_tau_factor_after_core_He_burn', value=relax_tau_factor_after_core_He_burn)
    self.add_to_star_job(category=category, optional=True,
        control='relax_tau_factor_after_core_C_burn', value=relax_tau_factor_after_core_C_burn)

def remove_surface_at_cell_k(self, value):
    self.add_to_star_job(category='remove surface',
        control='remove_surface_at_cell_k', value=value)

def remove_surface_at_he_core_boundary(self, value):
    self.add_to_star_job(category='remove surface',
        control='remove_surface_at_he_core_boundary', value=value)
    
def remove_surface_by_optical_depth(self, value):
    self.add_to_star_job(category='remove surface',
        control='remove_surface_by_optical_depth', value=value)

def remove_surface_by_density(self, value):
    self.add_to_star_job(category='remove surface',
        control='remove_surface_by_density', value=value)

def remove_surface_by_pressure(self, value):
    self.add_to_star_job(category='remove surface',
        control='remove_surface_by_pressure', value=value)

def remove_surface_by_mass_fraction_q(self, value):
    self.add_to_star_job(category='remove surface',
        control='remove_surface_by_mass_fraction_q', value=value)

def remove_surface_by_mass_gm(self, value):
    self.add_to_star_job(category='remove surface',
        control='remove_surface_by_mass_gm', value=value)

def remove_surface_by_mass_Msun(self, value):
    self.add_to_star_job(category='remove surface',
        control='remove_surface_by_mass_Msun', value=value)

def remove_surface_by_radius_cm(self, value):
    self.add_to_star_job(category='remove surface',
        control='remove_surface_by_radius_cm', value=value)

def remove_surface_by_radius_Rsun(self, value):
    self.add_to_star_job(category='remove surface',
        control='remove_surface_by_radius_Rsun', value=value)

def remove_surface_by_v_surf_km_s(self, value):
    self.add_to_star_job(category='remove surface',
        control='remove_surface_by_v_surf_km_s', value=value)

def remove_surface_by_v_surf_div_cs(self, value):
    self.add_to_star_job(category='remove surface',
        control='remove_surface_by_v_surf_div_cs', value=value)

def remove_surface_by_v_surf_div_v_escape(self, value,
                    min_q_for_remove_surface_by_v_surf_div_v_escape=None,
                    max_q_for_remove_surface_by_v_surf_div_v_escape=None):
    self.add_to_star_job(category='remove surface',
        control='remove_surface_by_v_surf_div_v_escape', value=value)
    
    self.add_to_star_job(category='remove surface', optional=True,
        control='min_q_for_remove_surface_by_v_surf_div_v_escape', value=min_q_for_remove_surface_by_v_surf_div_v_escape)
    self.add_to_star_job(category='remove surface', optional=True,
        control='max_q_for_remove_surface_by_v_surf_div_v_escape', value=max_q_for_remove_surface_by_v_surf_div_v_escape)
