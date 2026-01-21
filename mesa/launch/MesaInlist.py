
from qol.mesa.launch.MesaControl import MesaControl
import qol.config as config
import qol.info as info

import qol.mesa.controls.io as io
import qol.mesa.controls.init as init
import qol.mesa.controls.hydro as hydro
import qol.mesa.controls.kap as kap
import qol.mesa.controls.mix as mix
import qol.mesa.controls.net as net
import qol.mesa.controls.stop as stop
import qol.mesa.controls.bc as bc
import qol.mesa.controls.pgstar as pgstar
import qol.mesa.controls.mesh as mesh
import qol.mesa.controls.solver as solver
import qol.mesa.controls.timestep as timestep
import qol.mesa.controls.coredef as coredef
import qol.mesa.controls.wind as wind

import os

import numpy as np
from collections.abc import Iterable

class MesaInlist:
    """
    Stores inlist information for MESA

    NOTE: could generalize this and MesaPythonScrip into a single Task object
    """
    def __init__(self, name, mesa_version=config.mesa_version, data_path='data/'):
        """
        name: an informative name

        relative path (rel_path) will be "inlist_{name}"
        """
        self.name = name
        self.rel_path = f'inlist_{name}'
        self.mesa_version = mesa_version
        self.data_path = data_path

        self.inlist_controls = [] # controls
        self.data_prereqs = [] # model files and other things which are required for this to work
        self.data_products = [] # model files and other things which are saved by this inlist

        self.LOGS_dir = LOGS_dir = f'LOGS/{name}/'
        self.add_to_controls(control='log_directory', value=LOGS_dir, optional=True)
        
        self.photos_dir = photos_dir = f'photos/{name}/'
        self.add_to_controls(control='photo_directory', value=photos_dir, optional=True)

        # defaults
        self.make_new_model = False
        self.use_pgstar = False # if True, also copy in qol preset inlist_pgstar

        self.num_extra_star_job_inlists = 0
        self.num_extra_controls_inlists = 0
        self.num_extra_pgstar_inlists = 0

        self.num_overshoot_zones = 0
        self.num_predictive_mix_zones = 0

        self.Reimers_scaling_factor = None
        self.Blocker_scaling_factor = None
        self.de_Jager_scaling_factor = None
        self.van_Loon_scaling_factor = None
        self.Nieuwenhuijzen_scaling_factor = None
        self.Dutch_scaling_factor = None

        self.Grid1_plot_coords = []


        # TODO
        # - disable specific nuclear reactions (or groups, like pp chain)
        # - define metal fractions
        # - writeout interval

    def rn_string(self):
        return f'./do_one tasks/{self.rel_path} {self.LOGS_dir}'

    def re_string(self):
        return f'./re_one tasks/{self.rel_path} {self.photos_dir} {self.LOGS_dir}'

    def add_control(self, namelist, control, value, category=None, comment=None, optional=False):
        """
        optional: if True and value is None, do nothing (for optional keywords)
        """
        if optional and value is None:
            return

        inlist_control = MesaControl(namelist=namelist, 
                control=control,
                value=value,
                category=category,
                comment=comment)
        
        self.inlist_controls.append(inlist_control)

    def add_list_control(self, namelist, control, values, category=None, comment=None, optional=False):
        """
        adds controls based off of list of arguments, rather than single argument

        e.g., control(1), control(2), etc.

        values MUST be iterable (or None, if optional=True)
        category and comment can either be iterable or not
        """
        if optional and values is None:
            return
        
        if isinstance(category, Iterable) and not isinstance(category, str):
            assert len(category) == len(values)
        else:
            category = len(values) * [category]
        
        if isinstance(comment, Iterable) and not isinstance(comment, str):
            assert len(comment) == len(values)
        else:
            comment = len(values) * [comment]
            
        for ii in range(len(values)):
            inlist_control = MesaControl(namelist=namelist, 
                    control=f'{control}({ii+1})',
                    value=values[ii],
                    category=category[ii],
                    comment=comment[ii])
        
            self.inlist_controls.append(inlist_control)

    # Aliases for add_control for specific namelists    
    def add_to_star_job(self, control, value, category=None, comment=None, optional=False):
        self.add_control(namelist='star_job', control=control, value=value, category=category, comment=comment, optional=optional)
    def add_to_eos(self, control, value, category=None, comment=None, optional=False):
        self.add_control(namelist='eos', control=control, value=value, category=category, comment=comment, optional=optional)
    def add_to_kap(self, control, value, category=None, comment=None, optional=False):
        self.add_control(namelist='kap', control=control, value=value, category=category, comment=comment, optional=optional)
    def add_to_controls(self, control, value, category=None, comment=None, optional=False):
        self.add_control(namelist='controls', control=control, value=value, category=category, comment=comment, optional=optional)
    def add_to_pgstar(self, control, value, category=None, comment=None, optional=False):
        self.add_control(namelist='pgstar', control=control, value=value, category=category, comment=comment, optional=optional)

    def add_list_to_star_job(self, control, values, category=None, comment=None, optional=False):
        self.add_read_extra_inlistlist_control(namelist='star_job', control=control, values=values, category=category, comment=comment, optional=optional)
    def add_list_to_eos(self, control, values, category=None, comment=None, optional=False):
        self.add_list_control(namelist='eos', control=control, values=values, category=category, comment=comment, optional=optional)
    def add_list_to_kap(self, control, values, category=None, comment=None, optional=False):
        self.add_list_control(namelist='kap', control=control, values=values, category=category, comment=comment, optional=optional)
    def add_list_to_controls(self, control, values, category=None, comment=None, optional=False):
        self.add_list_control(namelist='controls', control=control, values=values, category=category, comment=comment, optional=optional)
    def add_list_to_pgstar(self, control, values, category=None, comment=None, optional=False):
        self.add_list_control(namelist='pgstar', control=control, values=values, category=category, comment=comment, optional=optional)

    def save(self, run_path, subdir='tasks', absdir=None):
        """
        for saving inlist file and return text of inlist file
        run_path refers to a DIRECTORY, not the filename
        what is saved is run_path/tasks/rel_path

        ## if abs_path is None, still return the text but don't save anything

        if absdir is given, save file to directory at this absolute path instead
        """
        namelist_names = ['star_job', 'eos', 'kap', 'controls', 'pgstar']
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
            item_categories_in_namelist = np.sort(item_categories_in_namelist) # alphabetize for consistency

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

        inlist_text = f'! inlist created using the QoL package\n'
        inlist_text += f'! MESA version: {self.mesa_version}\n'
        inlist_text += f'! QoL commit hash: {info.qol_commit_hash}\n\n'

        inlist_text += ''.join(namelist_texts)

        if absdir is None:
            abs_path = os.path.join(run_path, subdir, self.rel_path)
        else:
            abs_path = os.path.join(absdir, self.rel_path)

        #if abs_path is not None:
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
    write_gyre_data_with_profile = io.write_gyre_data_with_profile
    surface_avg_abundance_dq = io.surface_avg_abundance_dq
    history_interval = io.history_interval
    profile_interval = io.profile_interval
    max_num_profile_models = io.max_num_profile_models

    create_initial_model = init.create_initial_model
    create_pre_main_sequence_model = init.create_pre_main_sequence_model
    relax_initial_mass = init.relax_initial_mass
    relax_initial_mass_to_remove_H_env = init.relax_initial_mass_to_remove_H_env
    initial_mass = init.initial_mass
    initial_y = init.initial_y
    initial_he3 = init.initial_he3
    initial_z = init.initial_z
    reset_age = init.reset_age
    reset_model_number = init.reset_model_number
    relax_tau_factor = init.relax_tau_factor
    remove_surface_at_cell_k = init.remove_surface_at_cell_k
    remove_surface_at_he_core_boundary = init.remove_surface_at_he_core_boundary
    remove_surface_by_optical_depth = init.remove_surface_by_optical_depth
    remove_surface_by_density = init.remove_surface_by_density
    remove_surface_by_pressure = init.remove_surface_by_pressure
    remove_surface_by_mass_fraction_q = init.remove_surface_by_mass_fraction_q
    remove_surface_by_mass_gm = init.remove_surface_by_mass_gm
    remove_surface_by_mass_Msun = init.remove_surface_by_mass_Msun
    remove_surface_by_radius_cm = init.remove_surface_by_radius_cm
    remove_surface_by_radius_Rsun = init.remove_surface_by_radius_Rsun
    remove_surface_by_v_surf_km_s = init.remove_surface_by_v_surf_km_s
    remove_surface_by_v_surf_div_cs = init.remove_surface_by_v_surf_div_cs
    remove_surface_by_v_surf_div_v_escape = init.remove_surface_by_v_surf_div_v_escape

    enable_hydrodynamics = hydro.enable_hydrodynamics
    disable_hydrodynamics = hydro.disable_hydrodynamics
    add_hydrodynamical_drag = hydro.add_hydrodynamical_drag

    set_Zbase = kap.set_Zbase

    change_net = net.change_net
    max_abar_for_burning = net.max_abar_for_burning
    disable_nuclear_burning = net.disable_nuclear_burning
    disable_dxdt_from_nuclear_burning = net.disable_dxdt_from_nuclear_burning
    replace_one_element_with_another = net.replace_one_element_with_another

    set_min_D_mix = mix.set_min_D_mix
    disable_mixing = mix.disable_mixing
    use_Ledoux_criterion = mix.use_Ledoux_criterion
    alpha_semiconvection = mix.alpha_semiconvection
    thermohaline_coeff = mix.thermohaline_coeff
    thermohaline_option = mix.thermohaline_option
    add_thermohaline_mixing = mix.add_thermohaline_mixing
    gravitational_settling = mix.gravitational_settling
    add_overshoot_zone = mix.add_overshoot_zone
    add_predictive_mix_zone = mix.add_predictive_mix_zone
    phase_separation = mix.phase_separation
    okay_to_reduce_gradT_excess = mix.okay_to_reduce_gradT_excess

    relax_to_inner_BC = bc.relax_to_inner_BC
    use_table_atmosphere = bc.use_table_atmosphere

    set_max_num_retries = solver.set_max_num_retries
    use_gold_tolerances = solver.use_gold_tolerances
    energy_eqn_option = solver.energy_eqn_option
    convergence_ignore_equL_residuals = solver.convergence_ignore_equL_residuals
    limit_for_rel_error_in_energy_conservation = solver.limit_for_rel_error_in_energy_conservation

    mesh_delta_coeff = mesh.mesh_delta_coeff
    min_dq = mesh.min_dq
    max_surface_cell_dq = mesh.max_surface_cell_dq
    okay_to_remesh = mesh.okay_to_remesh

    set_dX_limits = timestep.set_dX_limits
    set_delta_lg_star_mass_limits = timestep.set_delta_lg_star_mass_limits
    set_initial_dt = timestep.set_initial_dt

    max_age = stop.max_age
    Teff_upper_limit = stop.Teff_upper_limit
    Teff_lower_limit = stop.Teff_lower_limit
    max_model_number = stop.max_model_number
    he_core_mass_limit = stop.he_core_mass_limit
    co_core_mass_limit = stop.co_core_mass_limit
    one_core_mass_limit = stop.one_core_mass_limit
    fe_core_mass_limit = stop.fe_core_mass_limit
    neutron_rich_core_mass_limit = stop.neutron_rich_core_mass_limit
    stop_near_zams = stop.stop_near_zams
    stop_at_phase_PreMS = stop.stop_at_phase_PreMS
    stop_at_phase_ZAMS = stop.stop_at_phase_ZAMS
    stop_at_phase_IAMS = stop.stop_at_phase_IAMS
    stop_at_phase_TAMS = stop.stop_at_phase_TAMS
    stop_at_phase_He_Burn = stop.stop_at_phase_He_Burn
    stop_at_phase_ZACHeB = stop.stop_at_phase_ZACHeB
    stop_at_phase_TACHeB = stop.stop_at_phase_TACHeB
    stop_at_phase_TP_AGB = stop.stop_at_phase_TP_AGB
    stop_at_phase_C_Burn = stop.stop_at_phase_C_Burn
    stop_at_phase_Ne_Burn = stop.stop_at_phase_Ne_Burn
    stop_at_phase_O_Burn = stop.stop_at_phase_O_Burn
    stop_at_phase_Si_Burn = stop.stop_at_phase_Si_Burn
    stop_at_phase_WDCS = stop.stop_at_phase_WDCS
    add_xa_central_lower_limit = stop.add_xa_central_lower_limit
    gamma_center_limit = stop.gamma_center_limit
    eta_center_limit = stop.eta_center_limit
    log_center_density_upper_limit = stop.log_center_density_upper_limit
    log_center_density_lower_limit = stop.log_center_density_lower_limit
    log_center_temp_upper_limit = stop.log_center_temp_upper_limit
    log_center_temp_lower_limit = stop.log_center_temp_lower_limit
    surface_accel_div_grav_limit = stop.surface_accel_div_grav_limit
    log_max_temp_upper_limit = stop.log_max_temp_upper_limit
    log_max_temp_lower_limit = stop.log_max_temp_lower_limit
    center_entropy_upper_limit = stop.center_entropy_upper_limit
    center_entropy_lower_limit = stop.center_entropy_lower_limit
    max_entropy_upper_limit = stop.max_entropy_upper_limit
    max_entropy_lower_limit = stop.max_entropy_lower_limit
    HB_limit = stop.HB_limit
    star_mass_min_limit = stop.star_mass_min_limit
    remnant_mass_min_limit = stop.remnant_mass_min_limit
    ejecta_mass_max_limit = stop.ejecta_mass_max_limit
    envelope_mass_limit = stop.envelope_mass_limit
    envelope_fraction_left_limit = stop.envelope_fraction_left_limit
    xmstar_min_limit = stop.xmstar_min_limit
    xmstar_max_limit = stop.xmstar_max_limit
    he_layer_mass_lower_limit = stop.he_layer_mass_lower_limit
    abs_diff_lg_LH_lg_Ls_limit = stop.abs_diff_lg_LH_lg_Ls_limit
    photosphere_r_upper_limit = stop.photosphere_r_upper_limit
    photosphere_r_lower_limit = stop.photosphere_r_lower_limit
    photosphere_m_upper_limit = stop.photosphere_m_upper_limit
    photosphere_m_lower_limit = stop.photosphere_m_lower_limit
    photosphere_m_sub_M_center_limit = stop.photosphere_m_sub_M_center_limit
    log_Teff_upper_limit = stop.log_Teff_upper_limit
    log_Teff_lower_limit = stop.log_Teff_lower_limit
    log_Tsurf_upper_limit = stop.log_Tsurf_upper_limit
    log_Tsurf_lower_limit = stop.log_Tsurf_lower_limit
    log_Rsurf_upper_limit = stop.log_Rsurf_upper_limit
    log_Rsurf_lower_limit = stop.log_Rsurf_lower_limit
    log_L_upper_limit = stop.log_L_upper_limit
    log_L_lower_limit = stop.log_L_lower_limit
    log_g_upper_limit = stop.log_g_upper_limit
    log_g_lower_limit = stop.log_g_lower_limit
    log_Psurf_upper_limit = stop.log_Psurf_upper_limit
    log_Psurf_lower_limit = stop.log_Psurf_lower_limit
    log_Dsurf_upper_limit = stop.log_Dsurf_upper_limit
    log_Dsurf_lower_limit = stop.log_Dsurf_lower_limit
    power_nuc_burn_upper_limit = stop.power_nuc_burn_upper_limit
    power_h_burn_upper_limit = stop.power_h_burn_upper_limit
    power_he_burn_upper_limit = stop.power_he_burn_upper_limit
    power_z_burn_upper_limit = stop.power_z_burn_upper_limit
    power_nuc_burn_lower_limit = stop.power_nuc_burn_lower_limit
    power_h_burn_lower_limit = stop.power_h_burn_lower_limit
    power_he_burn_lower_limit = stop.power_he_burn_lower_limit
    power_z_burn_lower_limit = stop.power_z_burn_lower_limit
    max_abs_rel_run_E_err = stop.max_abs_rel_run_E_err
    max_number_retries = stop.max_number_retries
    min_timestep_limit = stop.min_timestep_limit
    center_Ye_lower_limit = stop.center_Ye_lower_limit
    center_R_lower_limit = stop.center_R_lower_limit
    fe_core_infall_limit = stop.fe_core_infall_limit
    non_fe_core_rebound_limit = stop.non_fe_core_rebound_limit
    v_div_csound_max_limit = stop.v_div_csound_max_limit
    v_div_csound_surf_limit = stop.v_div_csound_surf_limit
    v_surf_div_v_kh_upper_limit = stop.v_surf_div_v_kh_upper_limit
    v_surf_div_v_kh_lower_limit = stop.v_surf_div_v_kh_lower_limit
    v_surf_div_v_esc_limit = stop.v_surf_div_v_esc_limit
    v_surf_kms_limit = stop.v_surf_kms_limit
    Lnuc_div_L_zams_limit = stop.Lnuc_div_L_zams_limit
    Lnuc_div_L_upper_limit = stop.Lnuc_div_L_upper_limit
    Lnuc_div_L_lower_limit = stop.Lnuc_div_L_lower_limit
    set_gamma1_limit = stop.set_gamma1_limit
    gamma1_limit_max_v_div_vesc = stop.gamma1_limit_max_v_div_vesc
    set_Pgas_div_P_limit = stop.set_Pgas_div_P_limit
    peak_burn_vconv_div_cs_limit = stop.peak_burn_vconv_div_cs_limit
    omega_div_omega_crit_limit = stop.omega_div_omega_crit_limit
    delta_nu_lower_limit = stop.delta_nu_lower_limit
    delta_nu_upper_limit = stop.delta_nu_upper_limit
    shock_mass_upper_limit = stop.delta_nu_upper_limit
    mach1_mass_upper_limit = stop.mach1_mass_upper_limit
    delta_Pg_lower_limit = stop.delta_Pg_lower_limit
    delta_Pg_upper_limit = stop.delta_Pg_upper_limit
    stop_when_reach_this_cumulative_extra_heating = stop.stop_when_reach_this_cumulative_extra_heating

    enable_pgstar = pgstar.enable_pgstar
    use_qol_pgstar = pgstar.use_qol_pgstar
    pgstar_Grid1_enable = pgstar.pgstar_Grid1_enable
    pgstar_Grid1_add_plot = pgstar.pgstar_Grid1_add_plot
    pgstar_Grid1_qol_default = pgstar.pgstar_Grid1_qol_default
    save_pgstar = pgstar.save_pgstar

    he_core_boundary_h1_fraction = coredef.he_core_boundary_h1_fraction
    co_core_boundary_he4_fraction = coredef.co_core_boundary_he4_fraction
    one_core_boundary_he4_c12_fraction = coredef.one_core_boundary_he4_c12_fraction

    get_update_wind_scaling_factor = wind.get_update_wind_scaling_factor # helper
    update_wind_scaling_factor = wind.update_wind_scaling_factor # helper
    get_wind_scaling_factor_control = wind.get_wind_scaling_factor_control # helper
    cool_wind_RGB = wind.cool_wind_RGB
    cool_wind_AGB = wind.cool_wind_AGB
    RGB_to_AGB_wind_switch = wind.RGB_to_AGB_wind_switch
    cool_wind_full_on_T = wind.cool_wind_full_on_T
    hot_wind_full_on_T = wind.hot_wind_full_on_T
    gain_mass = wind.gain_mass
    lose_mass = wind.lose_mass
