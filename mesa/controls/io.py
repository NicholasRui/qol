# Methods handling controls about writing or reading out data

def load_model(self, rel_path, absdir='data/'):
    """
    by default, stored at work/data/rel_path
    """
    category = 'load initial model'

    if absdir[-1] != '/':
        absdir += '/'

    self.add_to_star_job(category=category,
            control='load_saved_model', value=True)
    self.add_to_star_job(category=category,
            control='load_model_filename', value=f'{absdir}{rel_path}')
    
    self.data_prereqs += [rel_path]

def save_final_model(self, rel_path, absdir='data/'):
    category = 'save final model'

    self.add_to_star_job(category=category,
            control='save_model_when_terminate', value=True)
    self.add_to_star_job(category=category,
            control='save_model_filename', value=f'{absdir}{rel_path}')
    
    self.data_products += [rel_path]

def read_extra_inlist(self, namelist, rel_path, absdir='data/', category=None, comment=None):
    """
    read extra inlist (used specifically in the case that it is a prereq)
    by default, looks in work/data/rel_path

    if absdir specified, subdir will get ignored
    """
    match namelist:
        case 'star_job':
            self.num_extra_star_job_inlists += 1
            control_bool = f'read_extra_star_job_inlist({self.num_extra_star_job_inlists})'
            control_path = f'extra_star_job_inlist_name({self.num_extra_star_job_inlists})'

        case 'controls':
            self.num_extra_controls_inlists += 1
            control_bool = f'read_extra_controls_inlist({self.num_extra_controls_inlists})'
            control_path = f'extra_controls_inlist_name({self.num_extra_controls_inlists})'
    
        case 'pgstar':
            self.num_extra_pgstar_inlists += 1
            control_bool = f'read_extra_pgstar_inlist({self.num_extra_pgstar_inlists})'
            control_path = f'extra_pgstar_inlist_name({self.num_extra_pgstar_inlists})'

        case _:
            raise ValueError(f'invalid namelist: {namelist}')

    self.add_control(namelist=namelist, category=category, comment=comment,
            control=control_bool, value=True)
    
    self.add_control(namelist=namelist, category=category, comment=comment,
            control=control_path, value=f'{absdir}/{rel_path}')

    if rel_path not in self.data_prereqs: # only add Sprereq if not already there
        self.data_prereqs.append(rel_path)

def write_model_with_profile(self):
    self.add_control(namelist='controls', category='write-out',
            control='write_model_with_profile', value=True)

def write_gyre_data_with_profile(self):
    """
    for each profile, write input file to GYRE
    """
    self.add_to_controls(category='write-out',
            control='write_pulse_data_with_profile', value=True)
    self.add_to_controls(category='write-out',
            control='pulse_data_format', value='GYRE')

def surface_avg_abundance_dq(self, value):
    self.add_to_controls(category='surface abundance write-out',
            control='surface_avg_abundance_dq', value=value)



