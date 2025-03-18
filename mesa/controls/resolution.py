# Methods handling controls about resolution

def mesh_delta_coeff(self, value):
    self.add_control(namelist='controls', category='resolution',
            control='mesh_delta_coeff', value=value)

def min_dq(self, value):
    self.add_control(namelist='controls', category='resolution',
            control='min_dq', value=value)

def max_surface_cell_dq(self, value):
    self.add_control(namelist='controls', category='resolution',
            control='max_surface_cell_dq', value=value)

