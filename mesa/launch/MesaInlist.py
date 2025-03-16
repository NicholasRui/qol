
from qol.mesa.launch.MesaControl import MesaControl
import qol.config as config
import qol.mesa.const as const
import qol.helper.formatter as formatter

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

        # - core definitions

        # automate settling, overshoot

        # - saving an inlist should add comments at the top about how this was generated, and what version
        # - limit_for_rel_error_in_energy_conservation

        # - writeout interval

        # - pgstar enable

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

    #############################################
    ###### Preset combinations of controls ######
    #############################################
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

    def enable_hydrodynamics(self):
        namelist = 'star_job'
        category = 'enable hydrodynamics'

        self.add_control(namelist=namelist, category=category,
                control='change_v_flag', value=True)
        self.add_control(namelist=namelist, category=category,
                control='change_initial_v_flag', value=True)
        self.add_control(namelist=namelist, category=category,
                control='new_v_flag', value=True)

    def set_Zbase(self, Zbase=0.02):
        namelist = 'kap'
        
        self.add_control(namelist=namelist,
                control='use_Type2_opacities', value=True)
        self.add_control(namelist=namelist,
                control='Zbase', value=Zbase)

    def change_net(self, net_name):
        namelist = 'star_job'
        category = 'reaction network'

        self.add_control(namelist=namelist, category=category,
                control='change_initial_net', value=True)
        self.add_control(namelist=namelist, category=category,
                control='new_net_name', value=net_name)

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

    def create_initial_model(self, M_in_Msun, R_in_Rsun):
        self.add_control(namelist='star_job', category='create initial model',
                control='create_initial_model', value=True)
        
        self.add_control(namelist='star_job', category='create initial model',
                control='mass_in_gm_for_create_initial_model', value=M_in_Msun * const.Msun)
        self.add_control(namelist='star_job', category='create initial model',
                control='radius_in_cm_for_create_initial_model', value=R_in_Rsun * const.Rsun)
        
        self.make_new_model = True

    def relax_to_inner_BC(self, M_new_Msun=None, R_center_Rsun=None, L_center_Lsun=None,
            dlgm_per_step=None, dlgR_per_step=None, dlgL_per_step=None,
            relax_M_center_dt=None, relax_R_center_dt=None, relax_L_center_dt=None):
        """
        Relax to inner boundary condition
        """
        namelist = 'star_job'

        if M_new_Msun is not None:
            category = 'relax to inner BC: mass'
            comment = 'total mass in Msun'

            self.add_control(namelist=namelist, category=category,
                    control='relax_M_center', value=True)
            self.add_control(namelist=namelist, category=category,
                    control='new_mass', value=M_new_Msun, comment=comment)
            
            self.add_control(namelist=namelist, category=category, skip_if_None=True,
                    control='dlgm_per_step', value=dlgm_per_step)
            self.add_control(namelist=namelist, category=category, comment='sec', skip_if_None=True,
                    control='relax_M_center_dt', value=relax_M_center_dt)

        if R_center_Rsun is not None:
            category = 'relax to inner BC: radius'
            R_center_cgs = R_center_Rsun * const.Rsun
            comment = f'cm = {formatter.to_fortran(R_center_Rsun)} Rsun'

            self.add_control(namelist=namelist, category=category,
                    control='relax_R_center', value=True)
            self.add_control(namelist=namelist, category=category,
                    control='new_R_center', value=R_center_cgs, comment=comment)
            
            self.add_control(namelist=namelist, category=category, skip_if_None=True,
                    control='dlgR_per_step', value=dlgR_per_step)
            self.add_control(namelist=namelist, category=category, comment='sec', skip_if_None=True,
                    control='relax_R_center_dt', value=relax_R_center_dt)    
        
        if L_center_Lsun is not None:
            category = 'relax to inner BC: luminosity'
            L_center_cgs = L_center_Lsun * const.Lsun
            comment = f'erg s-1 = {formatter.to_fortran(L_center_Lsun)} Lsun'

            self.add_control(namelist=namelist, category=category,
                    control='relax_L_center', value=True)
            self.add_control(namelist=namelist, category=category,
                    control='new_L_center', value=L_center_cgs, comment=comment)
            
            self.add_control(namelist=namelist, category=category, skip_if_None=True,
                    control='dlgL_per_step', value=dlgL_per_step)
            self.add_control(namelist=namelist, category=category, comment='sec', skip_if_None=True,
                    control='relax_L_center_dt', value=relax_L_center_dt)

    def add_drag_for_HSE(self, drag_coefficient, use_drag_energy=False):
        """
        For hydrodynamical mode, add extra drag coefficient
        in order to relax model into HSE
        """
        namelist = 'controls'
        category = 'artificial drag to relax into HSE'

        self.add_control(namelist=namelist, category=category,
                control='use_drag_energy', value=use_drag_energy)
        self.add_control(namelist=namelist, category=category,
                control='drag_coefficient', value=drag_coefficient)

    def set_min_D_mix(self, min_D_mix,
            mass_lower_limit_for_min_D_mix=None, mass_upper_limit_for_min_D_mix=None):
        namelist = 'controls'
        category = 'mixing: minimum mixing coefficient'

        self.add_control(namelist=namelist, category=category,
                control='set_min_D_mix', value=True)
        self.add_control(namelist=namelist, category=category,
                control='min_D_mix', value=min_D_mix)
        
        self.add_control(namelist=namelist, category=category, skip_if_None=True,
                control='mass_lower_limit_for_min_D_mix', value=mass_lower_limit_for_min_D_mix)
        self.add_control(namelist=namelist, category=category, skip_if_None=True,
                control='mass_upper_limit_for_min_D_mix', value=mass_upper_limit_for_min_D_mix)

    def set_max_num_retries(self, value):
        self.add_control(namelist='controls', category='retry limit',
            control='max_number_retries', value=value)
        self.add_control(namelist='controls', category='retry limit',
            control='relax_max_number_retries', value=value)

    def reset_age(self):
        self.add_control(namelist='star_job', category='reset age',
            control='set_initial_age', value=True)
        self.add_control(namelist='star_job', category='reset age',
            control='initial_age', value=0.)

    def reset_model_number(self):
        self.add_control(namelist='star_job', category='reset age',
            control='set_initial_model_number', value=True)
        self.add_control(namelist='star_job', category='reset age',
            control='initial_model_number', value=0)

    def relax_initial_mass(self, M_new_Msun, lg_max_abs_mdot=None):
        self.add_control(namelist='star_job', category='relax mass',
            control='relax_initial_mass', value=True)
        self.add_control(namelist='star_job', category='relax mass',
            control='new_mass', value=M_new_Msun)
        
        self.add_control(namelist='star_job', category='relax mass', skip_if_None=True,
            control='lg_max_abs_mdot', value=lg_max_abs_mdot)

    ###########################################
    ###### Shortcuts for common controls ######
    ###########################################
    def create_pre_main_sequence_model(self, value=True):
        self.add_control(namelist='star_job', category='create initial model',
                control='create_pre_main_sequence_model', value=value)
        
        self.make_new_model = True

    def use_gold_tolerances(self, value=True):
        self.add_control(namelist='controls', category='solver',
                control='use_gold_tolerances', value=value)
    
    def energy_eqn_option(self, value='dedt'):
        """
        either dedt or eps_grav
        """
        self.add_control(namelist='controls', category='solver',
                control='energy_eqn_option', value=value)

    def convergence_ignore_equL_residuals(self, value=True):
        self.add_control(namelist='controls', category='solver',
                control='convergence_ignore_equL_residuals', value=value)

    def limit_for_rel_error_in_energy_conservation(self, value):
        self.add_control(namelist='controls', category='solver',
                control='limit_for_rel_error_in_energy_conservation', value=value)

    def mesh_delta_coeff(self, value):
        self.add_control(namelist='controls', category='resolution',
                control='mesh_delta_coeff', value=value)
    
    def min_dq(self, value):
        self.add_control(namelist='controls', category='resolution',
                control='min_dq', value=value)

    def min_timestep_limit(self, value):
        self.add_control(namelist='controls', category='timestepping',
                control='min_timestep_limit', value=value)

    def initial_mass(self, value):
        self.add_control(namelist='controls', category='initial mass and composition',
                control='initial_mass', value=value)

    def initial_y(self, value):
        self.add_control(namelist='controls', category='initial mass and composition',
                control='initial_y', value=value)
        
        if not self.make_new_model:
            warnings.warn('not used yet, need to make model from scratch (e.g., pre-MS)')
    
    def initial_z(self, value):
        self.add_control(namelist='controls', category='initial mass and composition',
                control='initial_z', value=value)
    
    def disable_mixing(self):
        self.add_control(namelist='controls', category='mixing: disabled',
                control='mix_factor', value=0.)

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

    def max_age(self, value):
        self.add_control(namelist='controls', category='termination conditions', comment='yr',
                control='max_age', value=value)
    
    def Teff_upper_limit(self, value):
        self.add_control(namelist='controls', category='termination conditions',
                control='Teff_upper_limit', value=value)

    def Teff_lower_limit(self, value):
        self.add_control(namelist='controls', category='termination conditions',
                control='Teff_lower_limit', value=value)
    
    def max_model_number(self, value):
        self.add_control(namelist='controls', category='termination conditions',
                control='max_model_number', value=value)
    
    def he_core_mass_limit(self, value):
        self.add_control(namelist='controls', category='termination conditions',
                control='he_core_mass_limit', value=value)
    
    def co_core_mass_limit(self, value):
        self.add_control(namelist='controls', category='termination conditions',
                control='co_core_mass_limit', value=value)
    
    def one_core_mass_limit(self, value):
        self.add_control(namelist='controls', category='termination conditions',
                control='one_core_mass_limit', value=value)

    def he_core_boundary_h1_fraction(self, value):
        self.add_control(namelist='controls', category='core definition',
                control='he_core_boundary_h1_fraction', value=value)
    
    def co_core_boundary_he4_fraction(self, value):
        self.add_control(namelist='controls', category='core definition',
                control='co_core_boundary_he4_fraction', value=value)
    
    def one_core_boundary_he4_c12_fraction(self, value):
        self.add_control(namelist='controls', category='core definition',
                control='one_core_boundary_he4_c12_fraction', value=value)

    def show_pgstar(self, abs_path=None):
        self.add_control(namelist='star_job', category='enable pgstar',
                control='pgstar_flag', value=True)
        
        self.enable_pgstar(abs_path)
    
    def enable_pgstar(self, abs_path=None):
        self.use_pgstar = True

        # Store information about pgstar path, if specified
        self.inlist_pgstar_path = abs_path

    def write_model_with_profile(self):
        self.add_control(namelist='controls', category='write-out',
                control='write_model_with_profile', value=True)



