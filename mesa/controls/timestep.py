# Methods handling controls about timestepping

def min_timestep_limit(self, value):
    self.add_to_controls(category='timestepping',
            control='min_timestep_limit', value=value)

