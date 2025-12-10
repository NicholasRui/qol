from qol.mesa.launcher import *
import qol.info as info

import os

def make_merger_MS_HeWD(
        root_path, # absolute path in which to write directory
        MWD_in_Msun, # WD mass
        MMS_in_Msun, # MS mass
        T_WD,        # temperature of WD
        net_name='cno_extras_o18_to_mg26.net', # net
        ringdown_time_yr=1e4, # ringdown timescale to HSE
        disable_hydro_after_ringdown=True, # if True, disable hydro mode on after ringdown, else keep it enabled
        enable_pgstar=True,
        rgb_wind=True,
        alpha_semiconvection=0., #4e-2, # semiconvection
        thermohaline_coeff=1., # thermohaline -- probably more important
        thermohaline_option='Brown_Garaud_Stellmach_13', # can be 'Kippenhahn', 'Traxler_Garaud_Stellmach_11', or 'Brown_Garaud_Stellmach_13'
        source_sdk=True, # manually activate sdk, since Caltech HPC doesn't seem to like it
        mesh_delta_coeff=1.,
        include_late=False, # include "late" WD cooling phase with crystallization but no settling
        include_overshoot=True, # include overshoot during COWD phase
        save_directory=True, # if False, don't save directory
        data_path='data/',
        ):
    """
    Make merger between HeWD and MS
    This is reproducing what was done in m_plus_wd (Rui & Fuller 2024, OJAp) but going further
    """
    assert thermohaline_option in ['Kippenhahn', 'Traxler_Garaud_Stellmach_11', 'Brown_Garaud_Stellmach_13'], "thermohaline_option must be one of 'Kippenhahn', 'Traxler_Garaud_Stellmach_11', or 'Brown_Garaud_Stellmach_13'"
    tho_string_dict = {'Kippenhahn': 'K80', 'Traxler_Garaud_Stellmach_11': 'TGS11', 'Brown_Garaud_Stellmach_13': 'BGS13'}

    run_name = f'MS{MMS_in_Msun:.3f}+HeWD{MWD_in_Msun:.3f}TWD{T_WD/1000.:.1f}_sc{alpha_semiconvection:.4f}_th{thermohaline_coeff:.4f}_tho{tho_string_dict[thermohaline_option]}_rgbw{int(rgb_wind)}_mdc{mesh_delta_coeff:.2f}_hydro{int(not disable_hydro_after_ringdown)}_il{int(include_late)}_ov{int(include_overshoot)}'
    run_path = os.path.join(root_path, run_name)

    argdict = {'root_path': root_path,
               'MWD_in_Msun': MWD_in_Msun,
               'MMS_in_Msun': MMS_in_Msun,
               'T_WD': T_WD,
               'net_name': net_name,
               'ringdown_time_yr': ringdown_time_yr,
               'disable_hydro_after_ringdown': disable_hydro_after_ringdown,
               'enable_pgstar': enable_pgstar,
               'rgb_wind': rgb_wind,
               'alpha_semiconvection': alpha_semiconvection,
               'thermohaline_coeff': thermohaline_coeff,
               'thermohaline_option': thermohaline_option,
               'source_sdk': source_sdk,
               'mesh_delta_coeff': mesh_delta_coeff,
               'include_late': include_late,
               'include_overshoot': include_overshoot,
               'data_path': data_path,
               'id_str': '-', # # appended to by each function to keep track of used inputs
               }
    # create and save work directory
    work = MesaWorkingDirectory(run_path=run_path)
    work.copy_history_columns_list(os.path.join(info.qol_path,
         'mesa/resources/r24.08.1/history_columns_sbr.list'))
    work.copy_profile_columns_list(os.path.join(info.qol_path,
         'mesa/resources/r24.08.1/profile_columns_sbr.list'))
    work.load_qol_pgstar()

    # create HeWD
    work.add_task(helper_merger_MS_HeWD_evolve_rg(argdict))
    work.add_task(helper_merger_MS_HeWD_strip_rg(argdict))
    work.add_task(helper_merger_MS_HeWD_cool_he_wd(argdict))

    # create envelope
    work.add_task(helper_merger_MS_HeWD_inner_bc(argdict))
    work.add_task(helper_merger_MS_HeWD_env_to_th_eq(argdict))

    # create and evolve merger remnant
    work.add_task(helper_merger_MS_HeWD_merge(argdict))
    work.add_task(helper_merger_MS_HeWD_remnant_ringdown(argdict))
    work.add_task(helper_merger_MS_HeWD_remnant_to_trgb(argdict))
    work.add_task(helper_merger_MS_HeWD_trgb_to_zacheb(argdict))
    work.add_task(helper_merger_MS_HeWD_zacheb_to_co_wd(argdict))
    work.add_task(helper_merger_MS_HeWD_cool_co_wd_early(argdict))
    if include_late:
        work.add_task(helper_merger_MS_HeWD_cool_co_wd_late(argdict))

    if save_directory:
        work.save_directory(slurm_job_name=run_name, grant_perms=True, source_sdk=source_sdk, data_path=data_path)

    return work


##########################################################
###### HELPER FUNCTIONS FOR CREATING REQUIRED TASKS ######
##########################################################

def helper_merger_MS_HeWD_evolve_rg(argdict):
    """
    make RG, evolve to desired core mass
    """
    enable_pgstar = argdict['enable_pgstar']
    net_name = argdict['net_name']
    MWD_in_Msun = argdict['MWD_in_Msun']
    mesh_delta_coeff = argdict['mesh_delta_coeff']
    data_path = argdict['data_path']

    argdict['id_str'] += f'{net_name}_MWD{MWD_in_Msun:.3f}_mdc{mesh_delta_coeff:.2f}'
    id_str = argdict['id_str']

    inlist = MesaInlist(name='evolve_rg', data_path=data_path)
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
    inlist.save_final_model(f'rg.mod{id_str}')

    return inlist

def helper_merger_MS_HeWD_strip_rg(argdict):
    """
    remove mass from RG
    """
    enable_pgstar = argdict['enable_pgstar']
    MWD_in_Msun = argdict['MWD_in_Msun']
    mesh_delta_coeff = argdict['mesh_delta_coeff']
    data_path = argdict['data_path']
    
    id_str = argdict['id_str']

    inlist = MesaInlist(name='strip_rg', data_path=data_path)
    if enable_pgstar:
        inlist.enable_pgstar()
    inlist.save_pgstar(write_path='Grid1/strip_rg/')

    # average composition of outer layers for write-out
    inlist.surface_avg_abundance_dq(1e-2)

    # remove envelope by relaxation
    inlist.load_model(f'rg.mod{id_str}')
    inlist.relax_initial_mass(M_new_Msun=MWD_in_Msun, lg_max_abs_mdot=-3)
    inlist.he_core_boundary_h1_fraction(1e-3)

    inlist.set_Zbase(0.02)

    inlist.energy_eqn_option('eps_grav')
    inlist.mesh_delta_coeff(mesh_delta_coeff)

    inlist.max_model_number(1)
    inlist.save_final_model(f'hot_he_wd.mod{id_str}')

    return inlist

def helper_merger_MS_HeWD_cool_he_wd(argdict):
    """
    cool He WD to desired temperature
    """
    enable_pgstar = argdict['enable_pgstar']
    T_WD = argdict['T_WD']
    mesh_delta_coeff = argdict['mesh_delta_coeff']
    data_path = argdict['data_path']

    id_str = argdict['id_str']

    inlist = MesaInlist(name='cool_he_wd', data_path=data_path)
    if enable_pgstar:
        inlist.enable_pgstar()
    inlist.save_pgstar(write_path='Grid1/cool_he_wd/')
    inlist.use_qol_pgstar()

    inlist.load_model(f'hot_he_wd.mod{id_str}')
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

    argdict['id_str'] += f'TWD{T_WD/1000.:.1f}'
    id_str = argdict['id_str0'] = argdict['id_str']

    inlist.save_final_model(f'cool_he_wd.mod{id_str}')

    return inlist

def helper_merger_MS_HeWD_inner_bc(argdict):
    """
    create envelope matching core model
    """
    MMS_in_Msun = argdict['MMS_in_Msun']
    data_path = argdict['data_path']
    
    argdict['id_str'] += f'MMS{MMS_in_Msun:.3f}'
    id_str = argdict['id_str']
    id_str0 = argdict['id_str0']

    script = MesaPythonScript(name='inner_bc',
            template=os.path.join(info.qol_path, 'mesa/templates/scripts/call_create_env_inlist_from_core.py'),
            const_args=[MMS_in_Msun, data_path],
            prereqs=[f'cool_he_wd.mod{id_str0}'],
            products=[f'inlist_env_inner_bc{id_str}'],
            data_path=data_path)

    return script

def helper_merger_MS_HeWD_env_to_th_eq(argdict):
    """
    run envelope model to thermal equilibrium (no dxdt_nuc)
    """
    enable_pgstar = argdict['enable_pgstar']
    net_name = argdict['net_name']
    MMS_in_Msun = argdict['MMS_in_Msun']
    mesh_delta_coeff = argdict['mesh_delta_coeff']
    data_path = argdict['data_path']

    id_str = argdict['id_str']

    inlist = MesaInlist(name='env_to_th_eq', data_path=data_path)
    if enable_pgstar:
        inlist.enable_pgstar()
    inlist.save_pgstar(write_path='Grid1/env_to_th_eq/')
    inlist.use_qol_pgstar()

    # set inner boundary condition
    inlist.read_extra_inlist(namelist='star_job', rel_path=f'inlist_env_inner_bc{id_str}',
                             category='relax inner BC to accommodate core model')

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
    inlist.save_final_model(f'env_th_eq.mod{id_str}')

    return inlist

def helper_merger_MS_HeWD_merge(argdict):
    """
    stitch core and envelope together
    """
    data_path = argdict['data_path']

    id_str = argdict['id_str']
    id_str0 = argdict['id_str0']

    script = MesaPythonScript(name='merge',
            template=os.path.join(info.qol_path, 'mesa/templates/scripts/call_create_shell_burning_remnant.py'),
            const_args=['excise', 'change_m', data_path],
            prereqs=[f'cool_he_wd.mod{id_str0}', f'env_th_eq.mod{id_str}'],
            products=[f'remnant_init.mod{id_str}'],
            data_path=data_path)

    return script

def helper_merger_MS_HeWD_remnant_ringdown(argdict):
    """
    run remnant into HSE
    """
    enable_pgstar = argdict['enable_pgstar']
    ringdown_time_yr = argdict['ringdown_time_yr']
    mesh_delta_coeff = argdict['mesh_delta_coeff']
    data_path = argdict['data_path']

    id_str = argdict['id_str']

    inlist = MesaInlist(name='remnant_ringdown', data_path=data_path)
    if enable_pgstar:
        inlist.enable_pgstar()
    inlist.save_pgstar(write_path='Grid1/remnant_ringdown/')
    inlist.use_qol_pgstar()

    # Write GYRE model files
    inlist.write_gyre_data_with_profile()

    inlist.load_model(f'remnant_init.mod{id_str}')
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
    argdict['id_str'] += f'_trd{ringdown_time_yr:.0f}'
    id_str = argdict['id_str']

    inlist.max_age(ringdown_time_yr)
    inlist.save_final_model(f'remnant_hse.mod{id_str}')

    return inlist

def helper_merger_MS_HeWD_remnant_to_trgb(argdict):
    """
    run remnant to tRGB
    """
    enable_pgstar = argdict['enable_pgstar']
    rgb_wind = argdict['rgb_wind']
    mesh_delta_coeff = argdict['mesh_delta_coeff']
    disable_hydro_after_ringdown = argdict['disable_hydro_after_ringdown']
    data_path = argdict['data_path']

    id_str = argdict['id_str']    
    
    inlist = MesaInlist(name='remnant_to_trgb', data_path=data_path)
    if enable_pgstar:
        inlist.enable_pgstar()
    inlist.save_pgstar(write_path='Grid1/remnant_to_trgb/')
    inlist.use_qol_pgstar()

    # Write GYRE model files
    inlist.write_gyre_data_with_profile()

    inlist.load_model(f'remnant_hse.mod{id_str}')
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

    argdict['id_str'] += f'_rgbw{int(rgb_wind)}_dhar{int(disable_hydro_after_ringdown)}'
    id_str = argdict['id_str'] 

    inlist.save_final_model(f'remnant_trgb.mod{id_str}')

    return inlist

def helper_merger_MS_HeWD_trgb_to_zacheb(argdict):
    """
    run remnant through He flash to ZACHeB
    """
    enable_pgstar = argdict['enable_pgstar']
    rgb_wind = argdict['rgb_wind']
    mesh_delta_coeff = argdict['mesh_delta_coeff']
    data_path = argdict['data_path']

    id_str = argdict['id_str']    
    
    inlist = MesaInlist(name='trgb_to_zacheb', data_path=data_path)
    if enable_pgstar:
        inlist.enable_pgstar()
    inlist.save_pgstar(write_path='Grid1/trgb_to_zacheb/')
    inlist.use_qol_pgstar()

    # Write GYRE model files
    inlist.write_gyre_data_with_profile()

    inlist.load_model(f'remnant_trgb.mod{id_str}')
    inlist.he_core_boundary_h1_fraction(1e-3)

    inlist.set_Zbase(0.02)

    # relax some convergence conditions
    inlist.min_timestep_limit(1e-12)
    inlist.energy_eqn_option('eps_grav')
    inlist.use_gold_tolerances(False)
    inlist.convergence_ignore_equL_residuals(True)
    inlist.set_max_num_retries(3000)

    inlist.mesh_delta_coeff(mesh_delta_coeff)

    # average composition of outer layers for write-out
    inlist.surface_avg_abundance_dq(1e-2)

    # wind
    if rgb_wind:
        inlist.cool_wind_RGB(scheme='Reimers', scaling_factor=0.5)
        inlist.cool_wind_AGB(scheme='Blocker', scaling_factor=0.1)
        inlist.RGB_to_AGB_wind_switch(1e-4)

    # stop at ZACHeB
    inlist.stop_at_phase_ZACHeB()
    inlist.save_final_model(f'remnant_zacheb.mod{id_str}')

    return inlist

def helper_merger_MS_HeWD_zacheb_to_co_wd(argdict):
    """
    run remnant from ZACHeB to CO WD
    """
    enable_pgstar = argdict['enable_pgstar']
    mesh_delta_coeff = argdict['mesh_delta_coeff']
    data_path = argdict['data_path']
    
    id_str = argdict['id_str']

    inlist = MesaInlist(name='zacheb_to_co_wd', data_path=data_path)
    if enable_pgstar:
        inlist.enable_pgstar()
    inlist.save_pgstar(write_path='Grid1/zacheb_to_co_wd/')
    inlist.use_qol_pgstar()

    # Write GYRE model files
    inlist.write_gyre_data_with_profile()

    inlist.load_model(f'remnant_zacheb.mod{id_str}')
    inlist.he_core_boundary_h1_fraction(1e-3)
    inlist.co_core_boundary_he4_fraction(1e-3)

    inlist.set_Zbase(0.02)

    inlist.energy_eqn_option('eps_grav')
    inlist.mesh_delta_coeff(mesh_delta_coeff)

    # average composition of outer layers for write-out
    inlist.surface_avg_abundance_dq(1e-2)

    # predictive mixing
    inlist.add_predictive_mix_zone(predictive_zone_type='any', predictive_zone_loc='core', predictive_bdy_loc='top',
                            predictive_superad_thresh=0.01, predictive_avoid_reversal='he4')

    # wind
    inlist.cool_wind_RGB(scheme='Reimers', scaling_factor=0.5)
    inlist.cool_wind_AGB(scheme='Blocker', scaling_factor=0.1)
    inlist.RGB_to_AGB_wind_switch(1e-4)
    inlist.cool_wind_full_on_T(9.99e9) # force usage of cool wind only
    inlist.hot_wind_full_on_T(1.e10)

    # stop at WD
    inlist.stop_at_phase_WDCS()
    inlist.save_final_model(f'hot_co_wd.mod{id_str}')

    return inlist

def helper_merger_MS_HeWD_cool_co_wd_early(argdict):
    """
    cool leftover CO WD through "early" stages
    """
    enable_pgstar = argdict['enable_pgstar']
    alpha_semiconvection = argdict['alpha_semiconvection']
    thermohaline_coeff = argdict['thermohaline_coeff']
    thermohaline_option = argdict['thermohaline_option']
    mesh_delta_coeff = argdict['mesh_delta_coeff']
    include_late = argdict['include_late']
    include_overshoot = argdict['include_overshoot']
    data_path = argdict['data_path']
    
    id_str = argdict['id_str']
    
    inlist = MesaInlist(name='cool_co_wd_early', data_path=data_path)
    if enable_pgstar:
        inlist.enable_pgstar()
    inlist.save_pgstar(write_path='Grid1/cool_co_wd_early/')
    inlist.use_qol_pgstar()

    # Write GYRE model files
    inlist.write_gyre_data_with_profile()

    inlist.load_model(f'hot_co_wd.mod{id_str}')
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
    inlist.add_thermohaline_mixing(thermohaline_coeff=thermohaline_coeff, thermohaline_option=thermohaline_option)

    inlist.gravitational_settling(diffusion_class_representatives=['h1', 'he3', 'he4', 'c12', 'o16', 'ne20', 'ne22', 'mg26'],
            diffusion_use_cgs_solver=True,
            show_diffusion_info=True,
            diffusion_steps_hard_limit=2000,
            diffusion_maxsteps_for_isolve=2000)

    # overshoot for very late thermal pulses
    if include_overshoot:
        inlist.add_overshoot_zone(overshoot_scheme='exponential', overshoot_zone_type='any',
                        overshoot_zone_loc='any', overshoot_bdy_loc='any',
                        overshoot_f=0.015, overshoot_f0=0.005)

    # stop after cool down enough, if include_late
    if include_late:
        inlist.log_Teff_lower_limit(4.4)
    else:
        inlist.max_age(1e10)

    # disable hydro, since this seems to cause issues during the thermal pulse
    inlist.disable_hydrodynamics()

    # convergence?
    inlist.use_gold_tolerances(False)

    tho_string_dict = {'Kippenhahn': 'K80', 'Traxler_Garaud_Stellmach_11': 'TGS11', 'Brown_Garaud_Stellmach_13': 'BGS13'}
    argdict['id_str'] += f'_sc{alpha_semiconvection:.2f}_th{thermohaline_coeff:.1f}_tho{tho_string_dict[thermohaline_option]}_il{int(include_late)}'
    id_str = argdict['id_str']

    inlist.save_final_model(f'cool_co_wd_early.mod{id_str}')

    # increase write-out interval to get the moment of crystallization
    inlist.profile_interval(1)

    # also stop at "crystallization"
    # inlist.gamma_center_limit(175)

    return inlist

def helper_merger_MS_HeWD_cool_co_wd_late(argdict):
    """
    cool leftover CO WD through "late" stages -- use table atmosphere for cool WD phase
    """
    enable_pgstar = argdict['enable_pgstar']
    alpha_semiconvection = argdict['alpha_semiconvection']
    thermohaline_coeff = argdict['thermohaline_coeff']
    thermohaline_option = argdict['thermohaline_option']
    mesh_delta_coeff = argdict['mesh_delta_coeff']
    include_overshoot = argdict['include_overshoot']
    data_path = argdict['data_path']
    
    id_str = argdict['id_str']

    inlist = MesaInlist(name='cool_co_wd_late', data_path=data_path)
    if enable_pgstar:
        inlist.enable_pgstar()
    inlist.save_pgstar(write_path='Grid1/cool_co_wd_late/')
    inlist.use_qol_pgstar()

    # Use tabulated WD atmosphere (Rohrmann+2011)
    inlist.use_table_atmosphere('WD_tau_25')

    # Write GYRE model files
    inlist.write_gyre_data_with_profile()

    inlist.load_model(f'cool_co_wd_early.mod{id_str}')
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

    # include settling
    # inlist.gravitational_settling(diffusion_class_representatives=['h1', 'he3', 'he4', 'c12', 'o16', 'ne20', 'ne22', 'mg26'],
    #         diffusion_use_cgs_solver=True,
    #         show_diffusion_info=True,
    #         diffusion_steps_hard_limit=2000,
    #         diffusion_maxsteps_for_isolve=2000)

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
    inlist.add_thermohaline_mixing(thermohaline_coeff=thermohaline_coeff, thermohaline_option=thermohaline_option)

    # phase separation
    # inlist.phase_separation(phase_separation_option='CO', do_phase_separation_heating=True, phase_separation_mixing_use_brunt=True)

    # overshoot for very late thermal pulses
    if include_overshoot:
        inlist.add_overshoot_zone(overshoot_scheme='exponential', overshoot_zone_type='any',
                        overshoot_zone_loc='any', overshoot_bdy_loc='any',
                        overshoot_f=0.015, overshoot_f0=0.005)

    # convergence?
    inlist.use_gold_tolerances(False)

    # increase write-out interval to get the moment of crystallization
    inlist.history_interval(1)
    inlist.profile_interval(1)

    # stop after a long time, if needed... okay to fail here, if sufficiently cooled
    inlist.max_age(1e10)
    inlist.save_final_model(f'cool_co_wd_late.mod{id_str}')

    return inlist
