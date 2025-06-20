import warnings

# TODO add hot wind

def get_update_wind_scaling_factor(self, scheme, scaling_factor):
    """
    Helper function which sorts wind options, makes sure they are valid, and that a scaling factor is not being overwritten
    """
    match scheme:
        case 'Reimers':
            current_scaling_factor = self.Reimers_scaling_factor
        case 'Blocker':
            current_scaling_factor = self.Blocker_scaling_factor
        case 'de Jager':
            current_scaling_factor = self.de_Jager_scaling_factor
        case 'van Loon':
            current_scaling_factor = self.van_Loon_scaling_factor
        case 'Nieuwenhuijzen':
            current_scaling_factor = self.Nieuwenhuijzen_scaling_factor
        case 'Dutch_scaling_factor':
            current_scaling_Factor = self.Dutch_scaling_factor
        case _:
            raise ValueError(f'wind schemne not found or supported yet: {scheme}')

    if (current_scaling_factor is not None) and (scaling_factor != current_scaling_factor):
        warnings.warn(f'Attempted to update scaling_factor for scheme {scheme} to {scaling_factor}, but it is already set to {current_scaling_factor}, so default to previous value')
        return False
    else:
        return True

def update_wind_scaling_factor(self, scheme, scaling_factor):
    """
    Helper function which updates the wind scaling factor
    """
    match scheme:
        case 'Reimers':
            self.Reimers_scaling_factor = scaling_factor
        case 'Blocker':
            self.Blocker_scaling_factor = scaling_factor
        case 'de Jager':
            self.de_Jager_scaling_factor = scaling_factor
        case 'van Loon':
            self.van_Loon_scaling_factor = scaling_factor
        case 'Nieuwenhuijzen':
            self.Nieuwenhuijzen_scaling_factor = scaling_factor
        case 'Dutch_scaling_factor':
            self.Dutch_scaling_factor = scaling_factor
        case _:
            raise ValueError(f'wind schemne not found or supported yet: {scheme}')

def get_wind_scaling_factor_control(self, scheme):
    """
    Helper function which retrieves the name of the MESA control for the desired wind scheme
    """
    match scheme:
        case 'Reimers': # suggested: 0.5
            scaling_factor_control = 'Reimers_scaling_factor'
        case 'Blocker': # suggested: 0.1
            scaling_factor_control = 'Blocker_scaling_factor'
        case 'de Jager':
            scaling_factor_control = 'de_Jager_scaling_factor'
        case 'van Loon':
            scaling_factor_control = 'van_Loon_scaling_factor'
        case 'Nieuwenhuijzen':
            scaling_factor_control = 'Nieuwenhuijzen_scaling_factor'
        case 'Dutch_scaling_factor':
            scaling_factor_control = 'Dutch_scaling_factor'
        case _:
            raise ValueError(f'wind schemne not found or supported yet: {scheme}')

    return scaling_factor_control

# Methods handling controls about winds

def cool_wind_RGB(self, scheme, scaling_factor):
    update_wind_scaling_factor = self.get_update_wind_scaling_factor(scheme, scaling_factor)
    scaling_factor_control = self.get_wind_scaling_factor_control(scheme)

    self.add_to_controls(category='wind',
            control='cool_wind_RGB_scheme', value=scheme)
    
    if update_wind_scaling_factor:
        scaling_factor_control = self.get_wind_scaling_factor_control(scheme)
        self.add_to_controls(category='wind',
                control=scaling_factor_control, value=scaling_factor)
        self.update_wind_scaling_factor(scheme, scaling_factor)

def cool_wind_AGB(self, scheme, scaling_factor):
    update_wind_scaling_factor = self.get_update_wind_scaling_factor(scheme, scaling_factor)
    scaling_factor_control = self.get_wind_scaling_factor_control(scheme)

    self.add_to_controls(category='wind',
            control='cool_wind_AGB_scheme', value=scheme)
    
    if update_wind_scaling_factor:
        scaling_factor_control = self.get_wind_scaling_factor_control(scheme)
        self.add_to_controls(category='wind',
                control=scaling_factor_control, value=scaling_factor)
        self.update_wind_scaling_factor(scheme, scaling_factor)

def RGB_to_AGB_wind_switch(self, value):
    self.add_to_controls(category='wind',
            control='RGB_to_AGB_wind_switch', value=value)

def cool_wind_full_on_T(self, value):
    self.add_to_controls(category='wind',
            control='cool_wind_full_on_T', value=value)

def hot_wind_full_on_T(self, value):
    self.add_to_controls(category='wind',
            control='hot_wind_full_on_T', value=value)

# Mass gain and loss is here
# TODO: add ability to use accrete_given_mass_fractions
def gain_mass(self, max_star_mass_for_gain, mass_change,
              accrete_same_as_surface=None, accretion_h1=None, accretion_h2=None, accretion_he3=None, accretion_he4=None, accretion_zfracs=None):
    # Note: force mass_change to be positive for this case
    assert mass_change > 0

    self.add_to_controls(category='mass change',
            control='max_star_mass_for_gain', value=max_star_mass_for_gain)
    self.add_to_controls(category='mass change',
            control='mass_change', value=mass_change)
    
    self.add_to_controls(category='mass change', optional=True,
            control='accrete_same_as_surface', value=accrete_same_as_surface)
    self.add_to_controls(category='mass change', optional=True,
            control='accretion_h1', value=accretion_h1)
    self.add_to_controls(category='mass change', optional=True,
            control='accretion_h2', value=accretion_h2)
    self.add_to_controls(category='mass change', optional=True,
            control='accretion_he3', value=accretion_he3)
    self.add_to_controls(category='mass change', optional=True,
            control='accretion_he4', value=accretion_he4)
    self.add_to_controls(category='mass change', optional=True,
            control='accretion_zfracs', value=accretion_zfracs)

def lose_mass(self, min_star_mass_for_loss, mass_change,
              accrete_same_as_surface=None, accretion_h1=None, accretion_h2=None, accretion_he3=None, accretion_he4=None, accretion_zfracs=None):
    # Note: force mass_change to be negative for this case
    if mass_change > 0:
        mass_change *= -1
        warnings.warn('For mass loss, mass_change should be negative. Flipping sign of mass_change automatically in inlist.')

    self.add_to_controls(category='mass change',
            control='min_star_mass_for_loss', value=min_star_mass_for_loss)
    self.add_to_controls(category='mass change',
            control='mass_change', value=mass_change)
    
    self.add_to_controls(category='mass change', optional=True,
            control='accrete_same_as_surface', value=accrete_same_as_surface)
    self.add_to_controls(category='mass change', optional=True,
            control='accretion_h1', value=accretion_h1)
    self.add_to_controls(category='mass change', optional=True,
            control='accretion_h2', value=accretion_h2)
    self.add_to_controls(category='mass change', optional=True,
            control='accretion_he3', value=accretion_he3)
    self.add_to_controls(category='mass change', optional=True,
            control='accretion_he4', value=accretion_he4)
    self.add_to_controls(category='mass change', optional=True,
            control='accretion_zfracs', value=accretion_zfracs)






