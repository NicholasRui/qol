# Methods handling controls about timestepping

def min_timestep_limit(self, value):
    self.add_control(namelist='controls', category='timestepping',
            control='min_timestep_limit', value=value)

