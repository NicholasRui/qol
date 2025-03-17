# Methods handling controls about pgstar

def show_pgstar(self, abs_path=None):
    self.add_control(namelist='star_job', category='enable pgstar',
            control='pgstar_flag', value=True)
    
    self.enable_pgstar(abs_path)

def enable_pgstar(self, abs_path=None):
    self.use_pgstar = True

    # Store information about pgstar path, if specified
    self.inlist_pgstar_path = abs_path

