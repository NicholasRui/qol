# Methods handling controls about the solver and convergence

def use_gold_tolerances(self, value=True):
    self.add_to_controls(category='solver',
            control='use_gold_tolerances', value=value)

def energy_eqn_option(self, value='dedt'):
    """
    either dedt or eps_grav
    """
    self.add_to_controls(category='solver',
            control='energy_eqn_option', value=value)

def convergence_ignore_equL_residuals(self, value=True):
    self.add_to_controls(category='solver',
            control='convergence_ignore_equL_residuals', value=value)

def limit_for_rel_error_in_energy_conservation(self, value):
    self.add_to_controls(category='solver',
            control='limit_for_rel_error_in_energy_conservation', value=value)

def set_max_num_retries(self, value):
    self.add_to_controls(category='retry limit',
        control='max_number_retries', value=value)
    self.add_to_controls(category='retry limit',
        control='relax_max_number_retries', value=value)
