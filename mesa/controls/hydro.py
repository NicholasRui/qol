# Methods handling controls about hydrodynamics

def toggle_hydrodynamics(self):
    category = 'enable hydrodynamics'

    self.add_to_star_job(category=category,
            control='change_v_flag', value=True)
    self.add_to_star_job(category=category,
            control='change_initial_v_flag', value=True)
    self.add_to_star_job(category=category,
            control='new_v_flag', value=True)
        
def add_hydrodynamical_drag(self, drag_coefficient, use_drag_energy=False, min_q_for_drag=None):
    """
    For hydrodynamical mode, add extra drag coefficient
    in order to relax model into HSE
    """
    category = 'artificial drag'

    self.add_to_controls(category=category,
            control='use_drag_energy', value=use_drag_energy)
    self.add_to_controls(category=category,
            control='drag_coefficient', value=drag_coefficient)
    self.add_to_controls(category=category, optional=True,
            control='min_q_for_drag', value=min_q_for_drag)

