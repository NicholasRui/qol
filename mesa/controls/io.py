# Methods handling controls about writing or reading out data

def load_model(self, rel_path):
    namelist = 'star_job'
    category = 'load initial model'

    self.add_control(namelist=namelist, category=category,
            control='load_saved_model', value=True)
    self.add_control(namelist=namelist, category=category,
            control='load_model_filename', value=rel_path)
    
    self.prereqs += [rel_path]

def save_final_model(self, rel_path):
    namelist = 'star_job'
    category = 'save final model'

    self.add_control(namelist=namelist, category=category,
            control='save_model_when_terminate', value=True)
    self.add_control(namelist=namelist, category=category,
            control='save_model_filename', value=rel_path)
    
    self.products += [rel_path]

def read_extra_inlist(self, namelist, rel_path, category=None, comment=None):
    """
    read extra inlist (used specifically in the case that it is a prereq)

    namelist is either 'star_job' or 'controls' for now
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
            raise ValueError('namelist not yet supported: pgstar')

        case _:
            raise ValueError(f'invalid namelist: {namelist}')

    self.add_control(namelist=namelist, category=category, comment=comment,
            control=control_bool, value=True)
    self.add_control(namelist=namelist, category=category, comment=comment,
            control=control_path, value=rel_path)

    self.prereqs.append(rel_path)

def write_model_with_profile(self):
    self.add_control(namelist='controls', category='write-out',
            control='write_model_with_profile', value=True)

