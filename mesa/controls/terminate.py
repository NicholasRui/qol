# Methods handling controls about termination conditions

def max_age(self, value):
    self.add_control(namelist='controls', category='termination conditions', comment='yr',
            control='max_age', value=value)

def Teff_upper_limit(self, value):
    self.add_control(namelist='controls', category='termination conditions',
            control='Teff_upper_limit', value=value)

def Teff_lower_limit(self, value):
    self.add_control(namelist='controls', category='termination conditions',
            control='Teff_lower_limit', value=value)

def max_model_number(self, value):
    self.add_control(namelist='controls', category='termination conditions',
            control='max_model_number', value=value)

def he_core_mass_limit(self, value):
    self.add_control(namelist='controls', category='termination conditions',
            control='he_core_mass_limit', value=value)

def co_core_mass_limit(self, value):
    self.add_control(namelist='controls', category='termination conditions',
            control='co_core_mass_limit', value=value)

def one_core_mass_limit(self, value):
    self.add_control(namelist='controls', category='termination conditions',
            control='one_core_mass_limit', value=value)

