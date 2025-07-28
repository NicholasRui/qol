# Methods handling controls about nuclear reactions

def change_net(self, net_name):
    namelist = 'star_job'
    category = 'reaction network'

    self.add_to_star_job(category=category,
            control='change_initial_net', value=True)
    self.add_to_star_job(category=category,
            control='new_net_name', value=net_name)

def max_abar_for_burning(self, value):
    self.add_to_controls(category='burning: disabled',
            control='max_abar_for_burning', value=value)

def disable_nuclear_burning(self):
    self.max_abar_for_burning(-1)

def disable_dxdt_from_nuclear_burning(self):
    """
    disable composition changes from burning, but not energy release
    (i.e., artificially setting nuclear timescale to infinity)
    """
    self.add_to_controls(category='burning: dxdt disabled',
            control='dxdt_nuc_factor', value=0.)

def replace_one_element_with_another(self, chem_name1, chem_name2, replace_element_nzlo=None, replace_element_nzhi=None):
    self.add_to_star_job(category='replace species',
            control='replace_element', value=True)
    self.add_to_star_job(category='replace species',
            control='replace_initial_element', value=True)
    
    self.add_to_star_job(category='replace species',
            control='chem_name1', value=chem_name1)
    self.add_to_star_job(category='replace species',
            control='chem_name2', value=chem_name2)

    self.add_to_star_job(category='replace species', optional=True,
            control='replace_element_nzlo', value=replace_element_nzlo)
    self.add_to_star_job(category='replace species', optional=True,
            control='replace_element_nzhi', value=replace_element_nzhi)

