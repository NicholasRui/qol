
from qol.mesa.launch.MesaControl import MesaControl
import qol.config as config
import qol.mesa.const as const
import qol.tools.formatter as formatter

import mesa.controls.io as io
import mesa.controls.init as init
import mesa.controls.hydro as hydro
import mesa.controls.kap as kap
import mesa.controls.mix as mix
import mesa.controls.net as net
import qol.mesa.controls.terminate as terminate
import mesa.controls.bc as bc
import mesa.controls.pgstar as pgstar
import mesa.controls.resolution as resolution
import mesa.controls.solver as solver
import mesa.controls.timestep as timestep
import mesa.controls.coredef as coredef

import numpy as np

import warnings


class MesaInlist:
    """
    Stores inlist information for MESA 
    """
    def __init__(self, rel_path, mesa_version=config.mesa_version, LOGS_dir='LOGS/'):
        """
        rel_path : relative path in work directory to save this inlist file
        """
        self.rel_path = rel_path
        self.mesa_version = mesa_version
        self.LOGS_dir = LOGS_dir

        self.inlist_controls = [] # controls
        self.prereqs = [] # model files and other things which are required for this to work
        self.products = [] # model files and other things which are saved by this inlist

        if LOGS_dir is not None: # unless LOGS_dir is explicitly set to None, automatically define it
            if LOGS_dir[-1] != '/':
                LOGS_dir += '/'
            self.add_control(namelist='controls', control='log_directory', value=LOGS_dir)

        # defaults
        self.make_new_model = False
        self.use_pgstar = False # if True, also copy in qol preset inlist_pgstar

        self.num_extra_star_job_inlists = 0
        self.num_extra_controls_inlists = 0
        # self.num_extra_pgstar_inlists = 0


        # TODO
        # - disable specific nuclear reactions (or groups, like pp chain)
        # - save gyre, model with profile, etc.
        # - mass_change? max or min star etc..
        # - add easy function for making a wind
        # - define metal fractions

        # automate settling, overshoot

        # - saving an inlist should add comments at the top about how this was generated, and what version
        # - limit_for_rel_error_in_energy_conservation

        # - writeout interval

        # - overshoot

    def rn_string(self):
        """
        rel_path = rel_path of file, with respect to the work directory (i.e., relative path)
        """
        return f'./do_one {self.rel_path} {self.LOGS_dir}'

    def add_control(self, namelist, control, value, category=None, comment=None, skip_if_None=False):
        """
        skip_if_None: if True and value is None, do nothing (for optional keywords)
        """
        if skip_if_None and value is None:
            return

        inlist_control = MesaControl(namelist=namelist, 
                control=control,
                value=value,
                category=category,
                comment=comment)
        
        self.inlist_controls.append(inlist_control)

    def save(self, run_path):
        """
        for saving inlist file and return text of inlist file
        run_path refers to a DIRECTORY, not the filename
        what is saved is run_path/rel_path

        if abs_path is None, still return the text but don't save anything
        """
        namelist_names = ['star_job', 'eos', 'kap', 'controls']
        
        if not self.enable_pgstar: # kludge to make a blank pgstar section
            namelist_names += ['pgstar']

        namelist_texts = []

        # grab attributes from existing inlist controls
        item_namelists = np.array([inlist_control.namelist for inlist_control in self.inlist_controls]).astype(str)
        item_categories = np.array([inlist_control.category for inlist_control in self.inlist_controls]).astype(str)

        for namelist_name in namelist_names:
            namelist_text = ''
            namelist_text += f'&{namelist_name}\n'
            namelist_text += f'  ! {namelist_name} options\n\n'

            # get all controls matching namelist
            in_namelist = (item_namelists == namelist_name)
            item_categories_in_namelist = np.unique(item_categories[in_namelist])

            # if None is in item_categories_in_namelist, sort it to the front
            # this is to take uncategorized controls and write them first
            if None in item_categories_in_namelist:
                item_categories_in_namelist = np.append(None, item_categories_in_namelist[item_categories_in_namelist != None])

            for category in item_categories_in_namelist:
                # get all controls matching both namelist and category
                in_category = (item_categories == category)
                indices_in_namelist_and_category, = np.where(in_namelist & in_category)

                # write brief comment corresponding to category and then add all controls
                if category != '':
                    namelist_text += f'  ! {category}\n'

                for ii in indices_in_namelist_and_category:
                    namelist_text += f'    {self.inlist_controls[ii].inlist_string()}\n'
                
                namelist_text += '\n'

            namelist_text += f'/ ! end of {namelist_name} namelist\n\n\n'
            namelist_texts.append(namelist_text)

        preface = f'! inlist created using qol package for MESA version {self.mesa_version}\n\n'
        inlist_text = preface + ''.join(namelist_texts)

        if self.use_pgstar and self.inlist_pgstar_path is not None:
            with open(self.inlist_pgstar_path, 'r') as f:
                pgstar_text = f.read()
            
            inlist_text += pgstar_text

        abs_path = f'{run_path}/{self.rel_path}'
        if abs_path is not None:
            with open(abs_path, 'w') as f:
                f.write(inlist_text)
        
        return inlist_text

    ###############################################
    ###### Methods for implementing controls ######
    ###############################################
    load_model = io.load_model
    save_final_model = io.save_final_model
    read_extra_inlist = io.read_extra_inlist
    write_model_with_profile = io.write_model_with_profile

    create_initial_model = init.create_initial_model
    create_pre_main_sequence_model = init.create_pre_main_sequence_model
    relax_initial_mass = init.relax_initial_mass
    initial_mass = init.initial_mass
    initial_y = init.initial_y
    initial_z = init.initial_z
    reset_age = init.reset_age
    reset_model_number = init.reset_model_number

    enable_hydrodynamics = hydro.enable_hydrodynamics
    add_drag_for_HSE = hydro.add_drag_for_HSE

    set_Zbase = kap.set_Zbase

    change_net = net.change_net
    disable_nuclear_burning = net.disable_nuclear_burning
    disable_dxdt_from_nuclear_burning = net.disable_dxdt_from_nuclear_burning

    set_min_D_mix = mix.set_min_D_mix
    disable_mixing = mix.disable_mixing

    relax_to_inner_BC = bc.relax_to_inner_BC

    set_max_num_retries = solver.set_max_num_retries
    use_gold_tolerances = solver.use_gold_tolerances
    energy_eqn_option = solver.energy_eqn_option
    convergence_ignore_equL_residuals = solver.convergence_ignore_equL_residuals
    limit_for_rel_error_in_energy_conservation = solver.limit_for_rel_error_in_energy_conservation

    mesh_delta_coeff = resolution.mesh_delta_coeff
    min_dq = resolution.min_dq

    min_timestep_limit = timestep.min_timestep_limit

    max_age = terminate.max_age
    Teff_upper_limit = terminate.Teff_upper_limit
    Teff_lower_limit = terminate.Teff_lower_limit
    max_model_number = terminate.max_model_number
    he_core_mass_limit = terminate.he_core_mass_limit
    co_core_mass_limit = terminate.co_core_mass_limit
    one_core_mass_limit = terminate.one_core_mass_limit

    show_pgstar = pgstar.show_pgstar
    enable_pgstar = pgstar.enable_pgstar

    he_core_boundary_h1_fraction = coredef.he_core_boundary_h1_fraction
    co_core_boundary_he4_fraction = coredef.co_core_boundary_he4_fraction
    one_core_boundary_he4_c12_fraction = coredef.one_core_boundary_he4_c12_fraction

