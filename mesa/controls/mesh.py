# Methods handling controls about resolution

def mesh_delta_coeff(self, value):
    self.add_to_controls(category='spatial resolution',
            control='mesh_delta_coeff', value=value)

def min_dq(self, value):
    self.add_to_controls(category='spatial resolution',
            control='min_dq', value=value)

def max_surface_cell_dq(self, value):
    self.add_to_controls(category='spatial resolution',
            control='max_surface_cell_dq', value=value)

def okay_to_remesh(self, value):
    self.add_to_controls(category='spatial resolution',
            control='okay_to_remesh', value=value)

