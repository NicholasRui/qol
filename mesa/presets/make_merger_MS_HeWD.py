from qol.mesa.launcher import *
import qol.info as info

def make_merger_MS_HeWD(
        root_path, # absolute path in which to write directory
        MWD_in_Msun, # WD mass
        MMS_in_Msun, # MS mass
        T_WD,        # temperature of WD
        net_name='pp_cno_extras_o18_ne22.net', #'cno_extras_o18_to_mg26.net', # net
        ringdown_time_yr=1e5, # ringdown timescale to HSE
        disable_hydro_after_ringdown=False, # if True, disable hydro mode on after ringdown, else keep it enabled
        enable_pgstar=False,
        rgb_wind=True,
        alpha_semiconvection=0., #4e-2, # semiconvection
        thermohaline_coeff=2., # thermohaline -- probably more important
        source_sdk=True, # manually activate sdk, since Caltech HPC doesn't seem to like it
        mesh_delta_coeff=1.,
        ):
    """
    Make merger between HeWD and MS
    This is reproducing what was done in m_plus_wd (Rui & Fuller 2024, OJAp) but going further
    """
    run_name = f'MS{MMS_in_Msun:.3f}+HeWD{MWD_in_Msun:.3f}TWD{T_WD/1000.:.1f}_sc{alpha_semiconvection:.4f}_th{thermohaline_coeff:.4f}_mdc{mesh_delta_coeff:.2f}_hydro{int(not disable_hydro_after_ringdown)}'
    run_path = f'{root_path}/{run_name}'

    # generate tasks
    task_evolve_rg = helper_merger_MS_HeWD_evolve_rg(enable_pgstar=enable_pgstar, net_name=net_name, MWD_in_Msun=MWD_in_Msun, mesh_delta_coeff=mesh_delta_coeff)
    task_strip_rg = helper_merger_MS_HeWD_strip_rg(enable_pgstar=enable_pgstar, MWD_in_Msun=MWD_in_Msun, mesh_delta_coeff=mesh_delta_coeff)
    task_cool_he_wd = helper_merger_MS_HeWD_cool_he_wd(enable_pgstar=enable_pgstar, T_WD=T_WD, mesh_delta_coeff=mesh_delta_coeff)
    task_inner_bc = helper_merger_MS_HeWD_inner_bc(MMS_in_Msun=MMS_in_Msun)
    task_env_to_th_eq = helper_merger_MS_HeWD_env_to_th_eq(enable_pgstar=enable_pgstar, net_name=net_name, MMS_in_Msun=MMS_in_Msun, mesh_delta_coeff=mesh_delta_coeff)
    task_merge = helper_merger_MS_HeWD_merge()
    task_remnant_ringdown = helper_merger_MS_HeWD_remnant_ringdown(enable_pgstar=enable_pgstar, ringdown_time_yr=ringdown_time_yr, mesh_delta_coeff=mesh_delta_coeff)
    task_remnant_to_trgb = helper_merger_MS_HeWD_remnant_to_trgb(enable_pgstar=enable_pgstar, rgb_wind=rgb_wind, mesh_delta_coeff=mesh_delta_coeff, disable_hydro_after_ringdown=disable_hydro_after_ringdown)
    task_trgb_to_zacheb = helper_merger_MS_HeWD_trgb_to_zacheb(enable_pgstar=enable_pgstar, rgb_wind=rgb_wind, mesh_delta_coeff=mesh_delta_coeff)
    task_zacheb_to_co_wd = helper_merger_MS_HeWD_zacheb_to_co_wd(enable_pgstar=enable_pgstar, mesh_delta_coeff=mesh_delta_coeff)
    task_cool_co_wd_early = helper_merger_MS_HeWD_cool_co_wd_early(enable_pgstar=enable_pgstar, alpha_semiconvection=alpha_semiconvection, thermohaline_coeff=thermohaline_coeff, mesh_delta_coeff=mesh_delta_coeff)
    task_cool_co_wd_late = helper_merger_MS_HeWD_cool_co_wd_late(enable_pgstar=enable_pgstar, alpha_semiconvection=alpha_semiconvection, thermohaline_coeff=thermohaline_coeff, mesh_delta_coeff=mesh_delta_coeff)
    
    # create and save work directory
    work = MesaWorkingDirectory(run_path=run_path)
    work.copy_history_columns_list(f'{info.qol_path}mesa/resources/r24.08.1/history_columns_hewd_ms.list')
    work.copy_profile_columns_list(f'{info.qol_path}mesa/resources/r24.08.1/profile_columns_qol.list')
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
    work.add_task(task_cool_co_wd_early)
    work.add_task(task_cool_co_wd_late)

    work.save_directory(slurm_job_name=run_name, grant_perms=True, source_sdk=source_sdk)

    return work


##########################################################
###### HELPER FUNCTIONS FOR CREATING REQUIRED TASKS ######
##########################################################

def helper_merger_MS_HeWD_evolve_rg(enable_pgstar, net_name, MWD_in_Msun, mesh_delta_coeff):
    """
    make RG, evolve to desired core mass
    """
    inlist = MesaInlist(name='evolve_rg')
    if enable_pgstar:
        inlist.enable_pgstar()
    inlist.save_pgstar(write_path='Grid1/evolve_rg/')
    inlist.use_qol_pgstar()

    # 1.2 Msun RG
    inlist.initial_mass(1.2)
    inlist.set_Zbase(0.02)
    inlist.change_net(net_name)

    inlist.energy_eqn_option('eps_grav')
    inlist.mesh_delta_coeff(mesh_delta_coeff)

    # average composition of outer layers for write-out
    inlist.surface_avg_abundance_dq(1e-2)

    # terminate when core reaches desired He WD mass
    inlist.he_core_boundary_h1_fraction(1e-3)
    inlist.he_core_mass_limit(MWD_in_Msun)
    inlist.save_final_model('rg.mod')

    return inlist

def helper_merger_MS_HeWD_strip_rg(enable_pgstar, MWD_in_Msun, mesh_delta_coeff):
    """
    remove mass from RG
    """
    inlist = MesaInlist(name='strip_rg')
    if enable_pgstar:
        inlist.enable_pgstar()
    inlist.save_pgstar(write_path='Grid1/strip_rg/')

    # average composition of outer layers for write-out
    inlist.surface_avg_abundance_dq(1e-2)

    # remove envelope by relaxation
    inlist.load_model('rg.mod')
    inlist.relax_initial_mass(M_new_Msun=MWD_in_Msun, lg_max_abs_mdot=-3)
    inlist.he_core_boundary_h1_fraction(1e-3)

    inlist.set_Zbase(0.02)

    inlist.energy_eqn_option('eps_grav')
    inlist.mesh_delta_coeff(mesh_delta_coeff)

    inlist.max_model_number(1)
    inlist.save_final_model('hot_he_wd.mod')

    return inlist

def helper_merger_MS_HeWD_cool_he_wd(enable_pgstar, T_WD, mesh_delta_coeff):
    """
    cool He WD to desired temperature
    """
    inlist = MesaInlist(name='cool_he_wd')
    if enable_pgstar:
        inlist.enable_pgstar()
    inlist.save_pgstar(write_path='Grid1/cool_he_wd/')
    inlist.use_qol_pgstar()

    inlist.load_model('hot_he_wd.mod')
    inlist.set_Zbase(0.02)

    # average composition of outer layers for write-out
    inlist.surface_avg_abundance_dq(1e-2)

    # resolution
    inlist.min_dq(1e-25)
    inlist.max_surface_cell_dq(1e-18)

    inlist.energy_eqn_option('eps_grav')
    inlist.mesh_delta_coeff(mesh_delta_coeff)

    # terminate when He WD reaches desired temperature
    inlist.Teff_lower_limit(T_WD)
    inlist.save_final_model('cool_he_wd.mod')

    return inlist

def helper_merger_MS_HeWD_inner_bc(MMS_in_Msun):
    """
    create envelope matching core model
    """
    script = MesaPythonScript(name='inner_bc',
            template=f'{info.qol_path}mesa/templates/scripts/call_create_env_inlist_from_core.py',
              const_args=[MMS_in_Msun], prereqs=['cool_he_wd.mod'], products=['inlist_env_inner_bc'])
    
    return script

def helper_merger_MS_HeWD_env_to_th_eq(enable_pgstar, net_name, MMS_in_Msun, mesh_delta_coeff):
    """
    run envelope model to thermal equilibrium (no dxdt_nuc)
    """
    inlist = MesaInlist(name='env_to_th_eq')
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

    # average composition of outer layers for write-out
    inlist.surface_avg_abundance_dq(1e-2)

    # relax some convergence conditions, disable dx/dt from burning
    inlist.min_timestep_limit(1e-12)
    inlist.energy_eqn_option('dedt')
    inlist.use_gold_tolerances(False)
    inlist.convergence_ignore_equL_residuals(True)
    inlist.set_max_num_retries(3000)
    inlist.limit_for_rel_error_in_energy_conservation(-1.)
    inlist.disable_dxdt_from_nuclear_burning()

    inlist.mesh_delta_coeff(mesh_delta_coeff)

    inlist.max_age(1e10)
    inlist.save_final_model('env_th_eq.mod')

    return inlist

def helper_merger_MS_HeWD_merge():
    """
    stitch core and envelope together
    """
    script = MesaPythonScript(name='merge',
            template=f'{info.qol_path}mesa/templates/scripts/call_create_shell_burning_remnant.py',
            prereqs=['cool_he_wd.mod', 'env_th_eq.mod'],
            products=['remnant_init.mod'])
    
    return script

def helper_merger_MS_HeWD_remnant_ringdown(enable_pgstar, ringdown_time_yr, mesh_delta_coeff):
    """
    run remnant into HSE
    """
    inlist = MesaInlist(name='remnant_ringdown')
    if enable_pgstar:
        inlist.enable_pgstar()
    inlist.save_pgstar(write_path='Grid1/remnant_ringdown/')
    inlist.use_qol_pgstar()

    # Write GYRE model files
    inlist.write_gyre_data_with_profile()

    inlist.load_model('remnant_init.mod')
    inlist.set_Zbase(0.02)

    # average composition of outer layers for write-out
    inlist.surface_avg_abundance_dq(1e-2)

    # enable hydro with drag
    # still disable dx/dt from burning
    inlist.enable_hydrodynamics()
    inlist.add_hydrodynamical_drag(drag_coefficient=1.)
    inlist.energy_eqn_option('eps_grav')
    inlist.min_timestep_limit(1e-12)
    inlist.use_gold_tolerances(False)
    inlist.convergence_ignore_equL_residuals(True)
    inlist.disable_dxdt_from_nuclear_burning()
    inlist.disable_mixing()

    inlist.mesh_delta_coeff(mesh_delta_coeff)

    # set arbitrary ringdown time
    inlist.max_age(ringdown_time_yr)
    inlist.save_final_model('remnant_hse.mod')

    return inlist

def helper_merger_MS_HeWD_remnant_to_trgb(enable_pgstar, rgb_wind, mesh_delta_coeff, disable_hydro_after_ringdown):
    """
    run remnant to tRGB
    """
    inlist = MesaInlist(name='remnant_to_trgb')
    if enable_pgstar:
        inlist.enable_pgstar()
    inlist.save_pgstar(write_path='Grid1/remnant_to_trgb/')
    inlist.use_qol_pgstar()

    # Write GYRE model files
    inlist.write_gyre_data_with_profile()

    inlist.load_model('remnant_hse.mod')
    inlist.set_Zbase(0.02)

    inlist.energy_eqn_option('eps_grav')
    inlist.mesh_delta_coeff(mesh_delta_coeff)

    # average composition of outer layers for write-out
    inlist.surface_avg_abundance_dq(1e-2)

    if disable_hydro_after_ringdown:
        inlist.disable_hydrodynamics()

    # add artificial damping to outermost layers
    # inlist.add_hydrodynamical_drag(drag_coefficient=1., min_q_for_drag=0.95)

    # wind
    if rgb_wind:
        inlist.cool_wind_RGB(scheme='Reimers', scaling_factor=0.5)
        inlist.cool_wind_AGB(scheme='Blocker', scaling_factor=0.1)
        inlist.RGB_to_AGB_wind_switch(1e-4)

    # evolve remnant up to tRGB
    inlist.stop_at_phase_He_Burn()
    inlist.save_final_model('remnant_trgb.mod')

    return inlist

def helper_merger_MS_HeWD_trgb_to_zacheb(enable_pgstar, rgb_wind, mesh_delta_coeff):
    """
    run remnant through He flash to ZACHeB
    """
    inlist = MesaInlist(name='trgb_to_zacheb')
    if enable_pgstar:
        inlist.enable_pgstar()
    inlist.save_pgstar(write_path='Grid1/trgb_to_zacheb/')
    inlist.use_qol_pgstar()

    # Write GYRE model files
    inlist.write_gyre_data_with_profile()

    inlist.load_model('remnant_trgb.mod')
    inlist.he_core_boundary_h1_fraction(1e-3)

    inlist.set_Zbase(0.02)

    # relax some convergence conditions
    inlist.min_timestep_limit(1e-12)
    inlist.energy_eqn_option('eps_grav')
    inlist.use_gold_tolerances(False)
    inlist.convergence_ignore_equL_residuals(True)
    inlist.set_max_num_retries(3000)

    inlist.mesh_delta_coeff(mesh_delta_coeff)

    # predictive mixing
    inlist.add_predictive_mix_zone(predictive_zone_type='any', predictive_zone_loc='core', predictive_bdy_loc='top',
                            predictive_superad_thresh=0.01, predictive_avoid_reversal='he4')

    # average composition of outer layers for write-out
    inlist.surface_avg_abundance_dq(1e-2)

    # add artificial damping to outermost layers
    # inlist.add_hydrodynamical_drag(drag_coefficient=1., min_q_for_drag=0.95)

    # wind
    if rgb_wind:
        inlist.cool_wind_RGB(scheme='Reimers', scaling_factor=0.5)
        inlist.cool_wind_AGB(scheme='Blocker', scaling_factor=0.1)
        inlist.RGB_to_AGB_wind_switch(1e-4)

    # stop at ZACHEB
    inlist.stop_at_phase_ZACHeB()
    inlist.save_final_model('remnant_zacheb.mod')

    return inlist

def helper_merger_MS_HeWD_zacheb_to_co_wd(enable_pgstar, mesh_delta_coeff):
    """
    run remnant from ZACHEB to CO WD
    """
    inlist = MesaInlist(name='zacheb_to_co_wd')
    if enable_pgstar:
        inlist.enable_pgstar()
    inlist.save_pgstar(write_path='Grid1/zacheb_to_co_wd/')
    inlist.use_qol_pgstar()

    # Write GYRE model files
    inlist.write_gyre_data_with_profile()

    inlist.load_model('remnant_zacheb.mod')
    inlist.he_core_boundary_h1_fraction(1e-3)
    inlist.co_core_boundary_he4_fraction(1e-3)

    inlist.set_Zbase(0.02)

    inlist.energy_eqn_option('eps_grav')
    inlist.mesh_delta_coeff(mesh_delta_coeff)

    # average composition of outer layers for write-out
    inlist.surface_avg_abundance_dq(1e-2)

    # add artificial damping to outermost layers
    # inlist.add_hydrodynamical_drag(drag_coefficient=1., min_q_for_drag=0.95)

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

def helper_merger_MS_HeWD_cool_co_wd_early(enable_pgstar, alpha_semiconvection, thermohaline_coeff, mesh_delta_coeff):
    """
    cool leftover CO WD through "early" stages -- include elemental diffusion but not phase separation
    """
    inlist = MesaInlist(name='cool_co_wd_early')
    if enable_pgstar:
        inlist.enable_pgstar()
    inlist.save_pgstar(write_path='Grid1/cool_co_wd_early/')
    inlist.use_qol_pgstar()

    # Write GYRE model files
    inlist.write_gyre_data_with_profile()

    inlist.load_model('hot_co_wd.mod')
    inlist.he_core_boundary_h1_fraction(1e-3)
    inlist.co_core_boundary_he4_fraction(1e-3)

    inlist.set_Zbase(0.02)

    inlist.energy_eqn_option('eps_grav')
    inlist.mesh_delta_coeff(mesh_delta_coeff)

    # resolution
    inlist.min_dq(1e-25)
    inlist.max_surface_cell_dq(1e-18)

    # MLT++ to get through late thermal pulses
    inlist.okay_to_reduce_gradT_excess()

    # average composition of outer layers for write-out
    inlist.surface_avg_abundance_dq(1e-2)

    # add artificial damping to outermost layers
    # inlist.add_hydrodynamical_drag(drag_coefficient=1., min_q_for_drag=0.95)

    # wind
    inlist.cool_wind_RGB(scheme='Reimers', scaling_factor=0.5)
    inlist.cool_wind_AGB(scheme='Blocker', scaling_factor=0.1)
    inlist.RGB_to_AGB_wind_switch(1e-4)
    inlist.cool_wind_full_on_T(9.99e9) # force usage of cool wind only
    inlist.hot_wind_full_on_T(1.e10)

    # mixing
    inlist.use_Ledoux_criterion()

    inlist.alpha_semiconvection(alpha_semiconvection)
    inlist.thermohaline_coeff(thermohaline_coeff)

    inlist.gravitational_settling(diffusion_class_representatives=['h1', 'he3', 'he4', 'c12', 'o16', 'ne20', 'ne22', 'mg26'],
            diffusion_use_cgs_solver=True,
            show_diffusion_info=True,
            diffusion_steps_hard_limit=2000,
            diffusion_maxsteps_for_isolve=2000)

    # stop after cool down enough
    inlist.log_L_lower_limit(0.)
    inlist.save_final_model('cool_co_wd_early.mod')

    return inlist

def helper_merger_MS_HeWD_cool_co_wd_late(enable_pgstar, alpha_semiconvection, thermohaline_coeff, mesh_delta_coeff):
    """
    cool leftover CO WD through "late" stages -- include phase separation but not elemental diffusion
    """
    inlist = MesaInlist(name='cool_co_wd_late')
    if enable_pgstar:
        inlist.enable_pgstar()
    inlist.save_pgstar(write_path='Grid1/cool_co_wd_late/')
    inlist.use_qol_pgstar()

    # Write GYRE model files
    inlist.write_gyre_data_with_profile()

    inlist.load_model('cool_co_wd_early.mod')
    inlist.he_core_boundary_h1_fraction(1e-3)
    inlist.co_core_boundary_he4_fraction(1e-3)

    inlist.set_Zbase(0.02)

    inlist.energy_eqn_option('eps_grav')
    inlist.mesh_delta_coeff(mesh_delta_coeff)

    # resolution
    inlist.min_dq(1e-25)
    inlist.max_surface_cell_dq(1e-18)

    # MLT++ to get through late thermal pulses
    inlist.okay_to_reduce_gradT_excess()

    # average composition of outer layers for write-out
    inlist.surface_avg_abundance_dq(1e-2)

    # add artificial damping to outermost layers
    # inlist.add_hydrodynamical_drag(drag_coefficient=1., min_q_for_drag=0.95)

    # wind
    inlist.cool_wind_RGB(scheme='Reimers', scaling_factor=0.5)
    inlist.cool_wind_AGB(scheme='Blocker', scaling_factor=0.1)
    inlist.RGB_to_AGB_wind_switch(1e-4)
    inlist.cool_wind_full_on_T(9.99e9) # force usage of cool wind only
    inlist.hot_wind_full_on_T(1.e10)

    # mixing
    inlist.use_Ledoux_criterion()

    inlist.alpha_semiconvection(alpha_semiconvection)
    inlist.thermohaline_coeff(thermohaline_coeff)

    # phase separation
    inlist.phase_separation(phase_separation_option='CO', do_phase_separation_heating=True, phase_separation_mixing_use_brunt=True)

    # stop after a long time, if needed... okay to fail here, if sufficiently cooled
    inlist.max_age(1e10)
    inlist.save_final_model('cool_co_wd_late.mod')

    return inlist
