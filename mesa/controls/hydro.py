# Methods handling controls about hydrodynamics

def enable_hydrodynamics(self):
    namelist = 'star_job'
    category = 'enable hydrodynamics'

    self.add_control(namelist=namelist, category=category,
            control='change_v_flag', value=True)
    self.add_control(namelist=namelist, category=category,
            control='change_initial_v_flag', value=True)
    self.add_control(namelist=namelist, category=category,
            control='new_v_flag', value=True)
        
def add_drag_for_HSE(self, drag_coefficient, use_drag_energy=False):
    """
    For hydrodynamical mode, add extra drag coefficient
    in order to relax model into HSE
    """
    namelist = 'controls'
    category = 'artificial drag to relax into HSE'

    self.add_control(namelist=namelist, category=category,
            control='use_drag_energy', value=use_drag_energy)
    self.add_control(namelist=namelist, category=category,
            control='drag_coefficient', value=drag_coefficient)

