# Methods handling controls about nuclear reactions

def change_net(self, net_name):
    namelist = 'star_job'
    category = 'reaction network'

    self.add_control(namelist=namelist, category=category,
            control='change_initial_net', value=True)
    self.add_control(namelist=namelist, category=category,
            control='new_net_name', value=net_name)
    
def disable_nuclear_burning(self):
    self.add_control(namelist='controls', category='burning: disabled',
            control='max_abar_for_burning', value=-1)

def disable_dxdt_from_nuclear_burning(self):
    """
    disable composition changes from burning, but not energy release
    (i.e., artificially setting nuclear timescale to infinity)
    """
    self.add_control(namelist='controls', category='burning: dxdt disabled',
            control='dxdt_nuc_factor', value=0.)

