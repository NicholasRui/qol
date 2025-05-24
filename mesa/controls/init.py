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

