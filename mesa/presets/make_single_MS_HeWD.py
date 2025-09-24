# Purpose: make a single-star model with similar parameters to those used in the MS+HeWD project, for comparison purposes
# 
#

from qol.mesa.launcher import *
import qol.info as info

import os

def make_single_MS_HeWD(
        root_path, # absolute path in which to write directory
        MMS_in_Msun, # ZAMS mass
        net_name='cno_extras_o18_to_mg26.net', # net
        enable_pgstar=True,
        rgb_wind=True,
        alpha_semiconvection=0., #4e-2, # semiconvection
        thermohaline_coeff=1., # thermohaline -- probably more important
        thermohaline_option='Brown_Garaud_Stellmach_13', # can be 'Kippenhahn', 'Traxler_Garaud_Stellmach_11', or 'Brown_Garaud_Stellmach_13'
        source_sdk=True, # manually activate sdk, since Caltech HPC doesn't seem to like it
        mesh_delta_coeff=1.,
        include_late=False, # include "late" WD cooling phase with crystallization but no settling
        save_directory=True, # if False, don't save directory
        data_path='data/',
        ):
    assert thermohaline_option in ['Kippenhahn', 'Traxler_Garaud_Stellmach_11', 'Brown_Garaud_Stellmach_13'], "thermohaline_option must be one of 'Kippenhahn', 'Traxler_Garaud_Stellmach_11', or 'Brown_Garaud_Stellmach_13'"
    tho_string_dict = {'Kippenhahn': 'K80', 'Traxler_Garaud_Stellmach_11': 'TGS11', 'Brown_Garaud_Stellmach_13': 'BGS13'}

    if data_path != 'data/': # TODO
        raise NotImplementedError()

    run_name = f'MS{MMS_in_Msun:.3f}_sc{alpha_semiconvection:.4f}_th{thermohaline_coeff:.4f}_tho{tho_string_dict[thermohaline_option]}_rgbw{int(rgb_wind)}_mdc{mesh_delta_coeff:.2f}'
    run_path = os.path.join(root_path, run_name)

    argdict = {'root_path': root_path,
               'MMS_in_Msun': MMS_in_Msun,
               'net_name': net_name,
               'enable_pgstar': enable_pgstar,
               'rgb_wind': rgb_wind,
               'alpha_semiconvection': alpha_semiconvection,
               'thermohaline_coeff': thermohaline_coeff,
               'thermohaline_option': thermohaline_option,
               'source_sdk': source_sdk,
               'mesh_delta_coeff': mesh_delta_coeff,
               'include_late': include_late,
               'data_path': data_path,}
    
    # create and save work directory
    work = MesaWorkingDirectory(run_path=run_path)
    work.copy_history_columns_list(os.path.join(info.qol_path,
            'mesa/resources/r24.08.1/history_columns_hewd_ms.list'))
    work.copy_profile_columns_list(os.path.join(info.qol_path,
            'mesa/resources/r24.08.1/profile_columns_qol.list'))
    work.load_qol_pgstar()

    work.add_task(helper_single_MS_HeWD_zams_to_tams(argdict))
    work.add_task(helper_single_MS_HeWD_tams_to_trgb(argdict))
    work.add_task(helper_single_MS_HeWD_trgb_to_zacheb(argdict))
    work.add_task(helper_single_MS_HeWD_zacheb_to_co_wd(argdict))
    work.add_task(helper_single_MS_HeWD_cool_co_wd_early(argdict))
    if include_late:
        work.add_task(helper_single_MS_HeWD_cool_co_wd_late(argdict))

    if save_directory:
        work.save_directory(slurm_job_name=run_name, grant_perms=True, source_sdk=source_sdk, data_path=data_path)

    return work

def helper_single_MS_HeWD_zams_to_tams(argdict):
    enable_pgstar = argdict['enable_pgstar']
    MMS_in_Msun = argdict['MMS_in_Msun']
    net_name = argdict['net_name']
    mesh_delta_coeff = argdict['mesh_delta_coeff']
    
    inlist = MesaInlist(name='zams_to_tams')
    if enable_pgstar:
        inlist.enable_pgstar()
    inlist.save_pgstar(write_path='Grid1/zams_to_tams/')
    inlist.use_qol_pgstar()

    # 1.2 Msun RG
    inlist.initial_mass(MMS_in_Msun)
    inlist.set_Zbase(0.02)
    inlist.change_net(net_name)

    inlist.energy_eqn_option('eps_grav')
    inlist.mesh_delta_coeff(mesh_delta_coeff)

    # average composition of outer layers for write-out
    inlist.surface_avg_abundance_dq(1e-2)

    # terminate when core reaches desired He WD mass
    inlist.he_core_boundary_h1_fraction(1e-3)
    inlist.stop_at_phase_TAMS()
    inlist.save_final_model('tams.mod')

    return inlist

def helper_single_MS_HeWD_tams_to_trgb(argdict):
    enable_pgstar = argdict['enable_pgstar']
    rgb_wind = argdict['rgb_wind']
    mesh_delta_coeff = argdict['mesh_delta_coeff']
    data_path = argdict['data_path']

    inlist = MesaInlist(name='tams_to_trgb')
    if enable_pgstar:
        inlist.enable_pgstar()
    inlist.save_pgstar(write_path='Grid1/tams_to_trgb/')
    inlist.use_qol_pgstar()

    # Write GYRE model files
    inlist.write_gyre_data_with_profile()

    inlist.load_model('tams.mod')
    inlist.set_Zbase(0.02)

    inlist.energy_eqn_option('eps_grav')
    inlist.mesh_delta_coeff(mesh_delta_coeff)

    # average composition of outer layers for write-out
    inlist.surface_avg_abundance_dq(1e-2)

    # wind
    if rgb_wind:
        inlist.cool_wind_RGB(scheme='Reimers', scaling_factor=0.5)
        inlist.cool_wind_AGB(scheme='Blocker', scaling_factor=0.1)
        inlist.RGB_to_AGB_wind_switch(1e-4)

    # evolve remnant up to tRGB
    inlist.stop_at_phase_He_Burn()
    inlist.save_final_model('trgb.mod')

    return inlist

def helper_single_MS_HeWD_trgb_to_zacheb(argdict):
    """
    run remnant through He flash to ZACHeB
    """
    enable_pgstar = argdict['enable_pgstar']
    rgb_wind = argdict['rgb_wind']
    mesh_delta_coeff = argdict['mesh_delta_coeff']
    data_path = argdict['data_path']

    inlist = MesaInlist(name='trgb_to_zacheb')
    if enable_pgstar:
        inlist.enable_pgstar()
    inlist.save_pgstar(write_path='Grid1/trgb_to_zacheb/')
    inlist.use_qol_pgstar()

    # Write GYRE model files
    inlist.write_gyre_data_with_profile()

    inlist.load_model('trgb.mod')
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

    # wind
    if rgb_wind:
        inlist.cool_wind_RGB(scheme='Reimers', scaling_factor=0.5)
        inlist.cool_wind_AGB(scheme='Blocker', scaling_factor=0.1)
        inlist.RGB_to_AGB_wind_switch(1e-4)

    # stop at ZACHEB
    inlist.stop_at_phase_ZACHeB()
    inlist.save_final_model('zacheb.mod')

    return inlist

def helper_single_MS_HeWD_zacheb_to_co_wd(argdict):
    """
    run remnant from ZACHEB to CO WD
    """
    enable_pgstar = argdict['enable_pgstar']
    mesh_delta_coeff = argdict['mesh_delta_coeff']
    data_path = argdict['data_path']

    inlist = MesaInlist(name='zacheb_to_co_wd')
    if enable_pgstar:
        inlist.enable_pgstar()
    inlist.save_pgstar(write_path='Grid1/zacheb_to_co_wd/')
    inlist.use_qol_pgstar()

    # Write GYRE model files
    inlist.write_gyre_data_with_profile()

    inlist.load_model('zacheb.mod')
    inlist.he_core_boundary_h1_fraction(1e-3)
    inlist.co_core_boundary_he4_fraction(1e-3)

    inlist.set_Zbase(0.02)

    inlist.energy_eqn_option('eps_grav')
    inlist.mesh_delta_coeff(mesh_delta_coeff)

    # average composition of outer layers for write-out
    inlist.surface_avg_abundance_dq(1e-2)

    # wind
    inlist.cool_wind_RGB(scheme='Reimers', scaling_factor=0.5)
    inlist.cool_wind_AGB(scheme='Blocker', scaling_factor=0.1)
    inlist.RGB_to_AGB_wind_switch(1e-4)
    inlist.cool_wind_full_on_T(9.99e9) # force usage of cool wind only
    inlist.hot_wind_full_on_T(1.e10)

    # stop at WD
    inlist.stop_at_phase_WDCS()
    inlist.save_final_model('hot_co_wd.mod')

    return inlist

def helper_single_MS_HeWD_cool_co_wd_early(argdict):
    """
    cool leftover CO WD through "early" stages -- include elemental diffusion but not phase separation
    """
    enable_pgstar = argdict['enable_pgstar']
    alpha_semiconvection = argdict['alpha_semiconvection']
    thermohaline_coeff = argdict['thermohaline_coeff']
    thermohaline_option = argdict['thermohaline_option']
    mesh_delta_coeff = argdict['mesh_delta_coeff']
    include_late = argdict['include_late']
    data_path = argdict['data_path']
    
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

    # stop after cool down enough
    if include_late:
        inlist.log_L_lower_limit(0.)
    else:
        inlist.max_age(1e10)
    inlist.save_final_model('cool_co_wd_early.mod')

    return inlist

def helper_single_MS_HeWD_cool_co_wd_late(argdict):
    """
    cool leftover CO WD through "late" stages -- include phase separation but not elemental diffusion
    """
    enable_pgstar = argdict['enable_pgstar']
    alpha_semiconvection = argdict['alpha_semiconvection']
    thermohaline_coeff = argdict['thermohaline_coeff']
    mesh_delta_coeff = argdict['mesh_delta_coeff']
    data_path = argdict['data_path']
    
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

    # wind
    inlist.cool_wind_RGB(scheme='Reimers', scaling_factor=0.5)
    inlist.cool_wind_AGB(scheme='Blocker', scaling_factor=0.1)
    inlist.RGB_to_AGB_wind_switch(1e-4)
    inlist.cool_wind_full_on_T(9.99e9) # force usage of cool wind only
    inlist.hot_wind_full_on_T(1.e10)

    # mixing
    inlist.use_Ledoux_criterion()

    inlist.alpha_semiconvection(alpha_semiconvection)
    inlist.add_thermohaline_mixing(thermohaline_coeff=thermohaline_coeff, thermohaline_option='Brown_Garaud_Stellmach_13')

    # phase separation
    inlist.phase_separation(phase_separation_option='CO', do_phase_separation_heating=True, phase_separation_mixing_use_brunt=True)

    # stop after a long time, if needed... okay to fail here, if sufficiently cooled
    inlist.max_age(1e10)
    inlist.save_final_model('cool_co_wd_late.mod')

    return inlist
