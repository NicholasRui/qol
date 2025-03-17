# Methods handling controls about mixing

def set_min_D_mix(self, min_D_mix,
        mass_lower_limit_for_min_D_mix=None, mass_upper_limit_for_min_D_mix=None):
    namelist = 'controls'
    category = 'mixing: minimum mixing coefficient'

    self.add_control(namelist=namelist, category=category,
            control='set_min_D_mix', value=True)
    self.add_control(namelist=namelist, category=category,
            control='min_D_mix', value=min_D_mix)
    
    self.add_control(namelist=namelist, category=category, skip_if_None=True,
            control='mass_lower_limit_for_min_D_mix', value=mass_lower_limit_for_min_D_mix)
    self.add_control(namelist=namelist, category=category, skip_if_None=True,
            control='mass_upper_limit_for_min_D_mix', value=mass_upper_limit_for_min_D_mix)
    
def disable_mixing(self):
    self.add_control(namelist='controls', category='mixing: disabled',
            control='mix_factor', value=0.)

        