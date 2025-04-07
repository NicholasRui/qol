from qol.mesa.launcher import *
import qol.paths as paths

def make_merger_MS_HeWD(
        run_path, # absolute path to save run
        MWD_in_Msun, # WD mass
        MMS_in_Msun, # MS mass
        T_WD,        # temperature of WD
        net_name='pp_cno_extras_o18_ne22.net', #'cno_extras_o18_to_mg26.net', # net
        ringdown_time_yr=1e5, # ringdown timescale to HSE
        enable_pgstar=False,
        rgb_wind=True,
        ):
    """
    Make merger between HeWD and MS
    This is reproducing what was done in m_plus_wd (Rui & Fuller 2024, OJAp) but going further
    """
    # generate tasks
    task_evolve_rg = helper_merger_MS_HeWD_evolve_rg(enable_pgstar=enable_pgstar, net_name=net_name, MWD_in_Msun=MWD_in_Msun)
    task_strip_rg = helper_merger_MS_HeWD_strip_rg(enable_pgstar=enable_pgstar, MWD_in_Msun=MWD_in_Msun)
    task_cool_he_wd = helper_merger_MS_HeWD_cool_he_wd(enable_pgstar=enable_pgstar, T_WD=T_WD)
    task_inner_bc = helper_merger_MS_HeWD_inner_bc(MMS_in_Msun=MMS_in_Msun)
    task_env_to_th_eq = helper_merger_MS_HeWD_env_to_th_eq(enable_pgstar=enable_pgstar, net_name=net_name, MMS_in_Msun=MMS_in_Msun)
    task_merge = helper_merger_MS_HeWD_merge()
    task_remnant_ringdown = helper_merger_MS_HeWD_remnant_ringdown(enable_pgstar=enable_pgstar, ringdown_time_yr=ringdown_time_yr)
    task_remnant_to_trgb = helper_merger_MS_HeWD_remnant_to_trgb(enable_pgstar=enable_pgstar, rgb_wind=rgb_wind)
    task_trgb_to_zacheb = helper_merger_MS_HeWD_trgb_to_zacheb(enable_pgstar=enable_pgstar, rgb_wind=rgb_wind)
    task_zacheb_to_co_wd = helper_merger_MS_HeWD_zacheb_to_co_wd(enable_pgstar=enable_pgstar)
    task_cool_co_wd = helper_merger_MS_HeWD_cool_co_wd(enable_pgstar=enable_pgstar)

    # create and save work directory
    work = MesaWorkingDirectory(run_path=run_path)
    work.copy_history_columns_list(f'{paths.qol_path}mesa/resources/r24.08.1/history_columns.list')
    work.copy_profile_columns_list(f'{paths.qol_path}mesa/resources/r24.08.1/profile_columns.list')
    work.load_qol_pgstar()

    work.add_task(task_evolve_rg)
    work.add_task(task_strip_rg)
    work.add_task(task_cool_he_wd)
    work.add_task(task_inner_bc)
    work.add_task(task_env_to_th_eq)
    work.add_task(task_merge)
    work.add_task(task_remnant_ringdown)
    work.add_task(task_remnant_to_trgb)
    work.add_task(task_trgb_to_zacheb)
    work.add_task(task_zacheb_to_co_wd)
    work.add_task(task_cool_co_wd)

    work.save_directory(grant_perms=True, slurm_job_name=f'M{MMS_in_Msun:.2f}+HeWD{MWD_in_Msun:.2f}TWD{T_WD/1000.:.1f}')


##########################################################
###### HELPER FUNCTIONS FOR CREATING REQUIRED TASKS ######
##########################################################

def helper_merger_MS_HeWD_evolve_rg(enable_pgstar, net_name, MWD_in_Msun):
    """
    make RG, evolve to desired core mass
    """
    inlist = MesaInlist('inlist_evolve_rg', LOGS_dir='LOGS/evolve_rg/', photos_dir='photos/evolve_rg/')
    if enable_pgstar:
        inlist.enable_pgstar()
    inlist.save_pgstar(write_path='Grid1/evolve_rg/')
    inlist.use_qol_pgstar()

    # 1.2 Msun RG
    inlist.initial_mass(1.2)
    inlist.set_Zbase(0.02)
    inlist.change_net(net_name)

    inlist.energy_eqn_option('eps_grav')

    # terminate when core reaches desired He WD mass
    inlist.he_core_boundary_h1_fraction(1e-3)
    inlist.he_core_mass_limit(MWD_in_Msun)
    inlist.save_final_model('rg.mod')

    return inlist

def helper_merger_MS_HeWD_strip_rg(enable_pgstar, MWD_in_Msun):
    """
    remove mass from RG
    """
    inlist = MesaInlist('inlist_strip_rg', LOGS_dir='LOGS/strip_rg/', photos_dir='photos/strip_rg/')
    if enable_pgstar:
        inlist.enable_pgstar()
    inlist.save_pgstar(write_path='Grid1/strip_rg/')

    # remove envelope by relaxation
    inlist.load_model('rg.mod')
    inlist.relax_initial_mass(M_new_Msun=MWD_in_Msun, lg_max_abs_mdot=-3)
    inlist.he_core_boundary_h1_fraction(1e-3)

    inlist.set_Zbase(0.02)

    inlist.energy_eqn_option('eps_grav')

    inlist.max_model_number(1)
    inlist.save_final_model('hot_he_wd.mod')

    return inlist

def helper_merger_MS_HeWD_cool_he_wd(enable_pgstar, T_WD):
    """
    cool He WD to desired temperature
    """
    inlist = MesaInlist('inlist_cool_he_wd', LOGS_dir='LOGS/cool_he_wd/', photos_dir='photos/cool_he_wd/')
    if enable_pgstar:
        inlist.enable_pgstar()
    inlist.save_pgstar(write_path='Grid1/cool_he_wd/')
    inlist.use_qol_pgstar()

    inlist.load_model('hot_he_wd.mod')
    inlist.set_Zbase(0.02)

    # resolution
    inlist.min_dq(1e-25)
    inlist.max_surface_cell_dq(1e-18)

    inlist.energy_eqn_option('eps_grav')

    # terminate when He WD reaches desired temperature
    inlist.Teff_lower_limit(T_WD)
    inlist.save_final_model('cool_he_wd.mod')

    return inlist

def helper_merger_MS_HeWD_inner_bc(MMS_in_Msun):
    """
    create envelope matching core model
    """
    script = MesaPythonScript('inner_bc.py',
            template=f'{paths.qol_path}mesa/templates/scripts/call_create_env_inlist_from_core.py',
              const_args=[MMS_in_Msun], prereqs=['cool_he_wd.mod'], products=['inlist_env_inner_bc'])
    
    return script

def helper_merger_MS_HeWD_env_to_th_eq(enable_pgstar, net_name, MMS_in_Msun):
    """
    run envelope model to thermal equilibrium (no dxdt_nuc)
    """
    inlist = MesaInlist('inlist_env_to_th_eq', LOGS_dir='LOGS/env_to_th_eq/', photos_dir='photos/env_to_th_eq/')
    if enable_pgstar:
        inlist.enable_pgstar()
    inlist.save_pgstar(write_path='Grid1/env_to_th_eq/')
    inlist.use_qol_pgstar()

    # set inner boundary condition
    inlist.read_extra_inlist(namelist='star_job', rel_path='inlist_env_inner_bc', category='relax inner BC to accommodate core model')

    # initialize as pre-MS
    inlist.create_pre_main_sequence_model(True)
    inlist.initial_mass(MMS_in_Msun)
    inlist.initial_y(0.28)
    inlist.initial_z(0.02)
    inlist.set_Zbase(0.02)
    inlist.change_net(net_name)

    # relax some convergence conditions, disable dx/dt from burning
    inlist.min_timestep_limit(1e-12)
    inlist.energy_eqn_option('dedt')
    inlist.use_gold_tolerances(False)
    inlist.convergence_ignore_equL_residuals(True)
    inlist.set_max_num_retries(3000)
    inlist.limit_for_rel_error_in_energy_conservation(-1.)
    inlist.disable_dxdt_from_nuclear_burning()

    inlist.max_age(1e10)
    inlist.save_final_model('env_th_eq.mod')

    return inlist


def helper_merger_MS_HeWD_merge():
    """
    stitch core and envelope together
    """
    script = MesaPythonScript(rel_path='merge.py',
            template=f'{paths.qol_path}mesa/templates/scripts/call_create_shell_burning_remnant.py',
            prereqs=['cool_he_wd.mod', 'env_th_eq.mod'],
            products=['remnant_init.mod'])
    
    return script

def helper_merger_MS_HeWD_remnant_ringdown(enable_pgstar, ringdown_time_yr):
    """
    run remnant into HSE
    """
    inlist = MesaInlist('inlist_remnant_ringdown', LOGS_dir='LOGS/remnant_ringdown/', photos_dir='photos/remnant_ringdown/')
    if enable_pgstar:
        inlist.enable_pgstar()
    inlist.save_pgstar(write_path='Grid1/remnant_ringdown/')
    inlist.use_qol_pgstar()

    inlist.load_model('remnant_init.mod')
    inlist.set_Zbase(0.02)

    # enable hydro with drag
    # still disable dx/dt from burning
    inlist.enable_hydrodynamics()
    inlist.add_drag_for_HSE(drag_coefficient=1.)
    inlist.energy_eqn_option('eps_grav')
    inlist.min_timestep_limit(1e-12)
    inlist.use_gold_tolerances(False)
    inlist.convergence_ignore_equL_residuals(True)
    inlist.disable_dxdt_from_nuclear_burning()
    inlist.disable_mixing()

    # set arbitrary ringdown time
    inlist.max_age(ringdown_time_yr)
    inlist.save_final_model('remnant_hse.mod')

    return inlist

def helper_merger_MS_HeWD_remnant_to_trgb(enable_pgstar, rgb_wind):
    """
    run remnant to tRGB
    """
    inlist = MesaInlist('inlist_remnant_to_trgb', LOGS_dir='LOGS/evolve_remnant/', photos_dir='photos/evolve_remnant/')
    if enable_pgstar:
        inlist.enable_pgstar()
    inlist.save_pgstar(write_path='Grid1/remnant_to_trgb/')
    inlist.use_qol_pgstar()

    inlist.load_model('remnant_hse.mod')
    inlist.set_Zbase(0.02)

    inlist.energy_eqn_option('eps_grav')

    # wind
    if rgb_wind:
        inlist.cool_wind_RGB(scheme='Reimers', scaling_factor=0.5)
        inlist.cool_wind_AGB(scheme='Blocker', scaling_factor=0.1)
        inlist.RGB_to_AGB_wind_switch(1e-4)

    # evolve remnant up to tRGB
    inlist.stop_at_phase_He_Burn()
    inlist.save_final_model('remnant_trgb.mod')

    return inlist

def helper_merger_MS_HeWD_trgb_to_zacheb(enable_pgstar, rgb_wind):
    """
    run remnant through He flash to ZACHeB
    """
    inlist = MesaInlist('inlist_trgb_to_zacheb', LOGS_dir='LOGS/trgb_to_zacheb/', photos_dir='photos/trgb_to_zacheb/')
    if enable_pgstar:
        inlist.enable_pgstar()
    inlist.save_pgstar(write_path='Grid1/trgb_to_zacheb/')
    inlist.use_qol_pgstar()

    inlist.load_model('remnant_trgb.mod')
    inlist.he_core_boundary_h1_fraction(1e-3)

    inlist.set_Zbase(0.02)

    # relax some convergence conditions
    inlist.min_timestep_limit(1e-12)
    inlist.energy_eqn_option('eps_grav')
    inlist.use_gold_tolerances(False)
    inlist.convergence_ignore_equL_residuals(True)
    inlist.set_max_num_retries(3000)

    # wind
    if rgb_wind:
        inlist.cool_wind_RGB(scheme='Reimers', scaling_factor=0.5)
        inlist.cool_wind_AGB(scheme='Blocker', scaling_factor=0.1)
        inlist.RGB_to_AGB_wind_switch(1e-4)

    # stop at ZACHEB
    inlist.stop_at_phase_ZACHeB()
    inlist.save_final_model('remnant_zacheb.mod')

    return inlist

def helper_merger_MS_HeWD_zacheb_to_co_wd(enable_pgstar):
    """
    run remnant from ZACHEB to CO WD
    """
    inlist = MesaInlist('inlist_zacheb_to_co_wd', LOGS_dir='LOGS/zacheb_to_co_wd/', photos_dir='photos/zacheb_to_co_wd/')
    if enable_pgstar:
        inlist.enable_pgstar()
    inlist.save_pgstar(write_path='Grid1/zacheb_to_co_wd/')
    inlist.use_qol_pgstar()

    inlist.load_model('remnant_zacheb.mod')
    inlist.he_core_boundary_h1_fraction(1e-3)
    inlist.co_core_boundary_he4_fraction(1e-3)

    inlist.set_Zbase(0.02)

    inlist.energy_eqn_option('eps_grav')

    # wind
    inlist.cool_wind_RGB(scheme='Reimers', scaling_factor=0.5)
    inlist.cool_wind_AGB(scheme='Blocker', scaling_factor=0.1)
    inlist.RGB_to_AGB_wind_switch(1e-4)
    inlist.cool_wind_full_on_T(9.99e9) # force usage of cool wind only
    inlist.hot_wind_full_on_T(1.e10)

    # stop at ZACHEB
    inlist.stop_at_phase_WDCS()
    inlist.save_final_model('hot_co_wd.mod')

    return inlist

def helper_merger_MS_HeWD_cool_co_wd(enable_pgstar):
    """
    cool leftover CO WD for a long time
    """
    inlist = MesaInlist('inlist_cool_co_wd', LOGS_dir='LOGS/cool_co_wd/', photos_dir='photos/cool_co_wd/')
    if enable_pgstar:
        inlist.enable_pgstar()
    inlist.save_pgstar(write_path='Grid1/cool_co_wd/')
    inlist.use_qol_pgstar()

    inlist.load_model('hot_co_wd.mod')
    inlist.he_core_boundary_h1_fraction(1e-3)
    inlist.co_core_boundary_he4_fraction(1e-3)

    inlist.set_Zbase(0.02)

    inlist.energy_eqn_option('eps_grav')

    # resolution
    inlist.min_dq(1e-25)
    inlist.max_surface_cell_dq(1e-18)

    # wind
    inlist.cool_wind_RGB(scheme='Reimers', scaling_factor=0.5)
    inlist.cool_wind_AGB(scheme='Blocker', scaling_factor=0.1)
    inlist.RGB_to_AGB_wind_switch(1e-4)
    inlist.cool_wind_full_on_T(9.99e9) # force usage of cool wind only
    inlist.hot_wind_full_on_T(1.e10)

    # mixing
    inlist.use_Ledoux_criterion()
    inlist.alpha_semiconvection(4e-2)
    inlist.thermohaline_coeff(2.)

    inlist.gravitational_settling(diffusion_class_representatives=['h1', 'he3', 'he4', 'c12', 'o16', 'ne20', 'ne22', 'mg26'],
            diffusion_use_cgs_solver=True,
            show_diffusion_info=True,
            diffusion_steps_hard_limit=2000,
            diffusion_maxsteps_for_isolve=2000)

    # phase separation
    inlist.phase_separation(phase_separation_option='CO', do_phase_separation_heating=True, phase_separation_mixing_use_brunt=True)

    # stop after a long time, if needed... okay to fail here, if sufficiently cooled
    inlist.max_age(1e10)
    inlist.save_final_model('cool_co_wd.mod')

    return inlist

