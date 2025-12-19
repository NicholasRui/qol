# Purpose: Make a RG+HeWD merger remnant
from qol.mesa.launcher import *
import qol.info as info

import os
import warnings

# NOTE
# re: id_str: there HAS to be a more elegant way to do this......
# e.g., allow inlist to store dependent variables and autogenerate the id_str
# also consider .inlist as a suffix to inlists

def make_merger_RG_HeWD(
        root_path, # absolute path in which to write directory
        MWD_in_Msun, # WD mass
        Mcore_in_Msun, # RG core mass
        Menv_in_Msun, # RG envelope mass
        T_WD,        # temperature of WD
        net_name='cno_extras_o18_to_mg26.net', # 'pp_cno_extras_o18_ne22.net' # net
        ringdown_time_yr=1e3, # ringdown timescale to HSE
        disable_hydro_after_ringdown=True, # if True, disable hydro mode on after ringdown, else keep it enabled
        enable_pgstar=True,
        rgb_wind=True,
        alpha_semiconvection=0., #4e-2, # semiconvection
        thermohaline_coeff=2., # thermohaline -- probably more important
        thermohaline_option='Kippenhahn', # can be 'Kippenhahn', 'Traxler_Garaud_Stellmach_11', or 'Brown_Garaud_Stellmach_13'
        source_sdk=True, # manually activate sdk, since Caltech HPC doesn't seem to like it
        mesh_delta_coeff=1.,
        include_late=False, # include "late" WD cooling phase with crystallization but no settling
        include_overshoot=True, # include overshoot during COWD phase
        save_directory=True, # if False, don't save directory
        data_path='data/',
        ):
    """
    Make merger between HeWD and RG
    This is reproducing what was done in Zhang & Jeffery 2013, MNRAS
    """
    assert thermohaline_option in ['Kippenhahn', 'Traxler_Garaud_Stellmach_11', 'Brown_Garaud_Stellmach_13'], "thermohaline_option must be one of 'Kippenhahn', 'Traxler_Garaud_Stellmach_11', or 'Brown_Garaud_Stellmach_13'"
    tho_string_dict = {'Kippenhahn': 'KRT80', 'Traxler_Garaud_Stellmach_11': 'TGS11', 'Brown_Garaud_Stellmach_13': 'BGS13'}
    if thermohaline_option == 'Brown_Garaud_Stellmach_13':
        warnings.warn('BGS13 prescription does not take into account electron conduction and will give wrong results for COWD phase.')

    run_name = f'RGc{Mcore_in_Msun:.3f}e{Menv_in_Msun:.3f}+HeWD{MWD_in_Msun:.3f}TWD{T_WD/1000.:.1f}_sc{alpha_semiconvection:.4f}_th{thermohaline_coeff:.4f}_tho{tho_string_dict[thermohaline_option]}_rgbw{int(rgb_wind)}_mdc{mesh_delta_coeff:.2f}_hydro{int(not disable_hydro_after_ringdown)}_il{int(include_late)}_ov{int(include_overshoot)}'
    run_path = os.path.join(root_path, run_name)

    argdict = {'root_path': root_path,
               'MWD_in_Msun': MWD_in_Msun,
               'Mcore_in_Msun': Mcore_in_Msun,
               'Menv_in_Msun': Menv_in_Msun,
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
               'id_str': '-', # appended to by each function to keep track of used inputs
               }    
    # create and save work directory
    work = MesaWorkingDirectory(run_path=run_path)
    work.copy_history_columns_list(os.path.join(info.qol_path,
         'mesa/resources/r24.08.1/history_columns_sbr.list'))
    work.copy_profile_columns_list(os.path.join(info.qol_path,
         'mesa/resources/r24.08.1/profile_columns_sbr.list'))
    work.load_qol_pgstar()

    # forming and cooling the HeWD
    work.add_task(helper_merger_RG_HeWD_evolve_rg(argdict))
    work.add_task(helper_merger_RG_HeWD_strip_rg(argdict))
    work.add_task(helper_merger_RG_HeWD_clean_he_wd(argdict))
    work.add_task(helper_merger_RG_HeWD_cool_he_wd(argdict))

    # artificially increasing the density of the outer layers of the HeWD
    work.add_task(helper_merger_MS_HeWD_inner_bc_for_compress(argdict))
    work.add_task(helper_merger_MS_HeWD_env_for_compress(argdict))
    work.add_task(helper_merger_MS_HeWD_merge_for_compress(argdict))
    work.add_task(helper_merger_MS_HeWD_ringdown_for_compress(argdict))
    work.add_task(helper_merger_RG_HeWD_disable_hydro_for_he_wd_with_env(argdict))
    work.add_task(helper_merger_MS_HeWD_remove_compress_env(argdict))
    work.add_task(helper_merger_RG_HeWD_clean_compressed_he_wd(argdict))

    # create hot outer core of merger remnant (progenitor: RG core)
    work.add_task(helper_merger_RG_HeWD_outer_core_inner_bc(argdict))
    work.add_task(helper_merger_RG_HeWD_make_outer_core(argdict))
    work.add_task(helper_merger_RG_HeWD_clean_outer_core(argdict))

    # add hydrogen to outer core and merge with HeWD
    work.add_task(helper_merger_RG_HeWD_accrete_env(argdict))
    work.add_task(helper_merger_RG_HeWD_merge(argdict))
    work.add_task(helper_merger_RG_HeWD_remnant_ringdown(argdict))

    # evolve merger remnant
    work.add_task(helper_merger_RG_HeWD_ignite_he(argdict))
    work.add_task(helper_merger_RG_HeWD_zacheb_to_co_wd(argdict))
    work.add_task(helper_merger_RG_HeWD_cool_co_wd_early(argdict))
    if include_late:
        work.add_task(helper_merger_RG_HeWD_cool_co_wd_late(argdict))

    if save_directory:
       work.save_directory(slurm_job_name=run_name, grant_perms=True, source_sdk=source_sdk, data_path=data_path)
    
    return work

def helper_merger_RG_HeWD_evolve_rg(argdict):
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

def helper_merger_RG_HeWD_strip_rg(argdict):
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

def helper_merger_RG_HeWD_clean_he_wd(argdict):
    """
    remove "contaminant" species
    """
    data_path = argdict['data_path']

    id_str = argdict['id_str']

    script = MesaPythonScript(name='clean_he_wd',
        template=os.path.join(info.qol_path, 'mesa/templates/scripts/call_replace_elements.py'),
        const_args=['he4', data_path,
                    'h1', 'he3', 'n13', 'o14', 'o15', 'f17', 'f18', 'ne18'],
        prereqs=[f'hot_he_wd.mod{id_str}'],
        products=[f'hot_he_wd_clean.mod{id_str}'],
        data_path=data_path)

    return script

def helper_merger_RG_HeWD_cool_he_wd(argdict):
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

    inlist.load_model(f'hot_he_wd_clean.mod{id_str}')
    inlist.set_Zbase(0.02)

    # average composition of outer layers for write-out
    inlist.surface_avg_abundance_dq(1e-2)

    # resolution
    inlist.min_dq(1e-25)
    inlist.max_surface_cell_dq(1e-18)

    inlist.energy_eqn_option('eps_grav')
    inlist.mesh_delta_coeff(mesh_delta_coeff)

    argdict['id_str'] += f'_TWD{T_WD/1000.:.1f}'
    id_str = argdict['id_str']

    # terminate when He WD reaches desired temperature
    inlist.Teff_lower_limit(T_WD)
    inlist.save_final_model(f'cool_he_wd.mod{id_str}')

    return inlist

def helper_merger_MS_HeWD_inner_bc_for_compress(argdict):
    """
    boundary condition for cool HeWD in order to compress the outer layers of HeWD
    """
    data_path = argdict['data_path']

    id_str = argdict['id_str']

    script = MesaPythonScript(name='inner_bc_for_compress',
              template=os.path.join(info.qol_path, 'mesa/templates/scripts/call_create_env_inlist_from_core.py'),
              const_args=[1.0, data_path],
              prereqs=[f'cool_he_wd.mod{id_str}'],
              products=[f'inlist_inner_bc_for_compress{id_str}'],
              data_path=data_path)

    return script

def helper_merger_MS_HeWD_env_for_compress(argdict):
    """
    create artificial H envelope to compress the outer layers of the HeWD
    """
    enable_pgstar = argdict['enable_pgstar']
    net_name = argdict['net_name']
    mesh_delta_coeff = argdict['mesh_delta_coeff']
    data_path = argdict['data_path']

    id_str = argdict['id_str']

    inlist = MesaInlist(name='env_for_compress', data_path=data_path)
    if enable_pgstar:
        inlist.enable_pgstar()
    inlist.save_pgstar(write_path='Grid1/env_for_compress/')
    inlist.use_qol_pgstar()

    # set inner boundary condition
    inlist.read_extra_inlist(namelist='star_job', rel_path=f'inlist_inner_bc_for_compress{id_str}',
                             category='relax inner BC to accommodate core model')

    # initialize as pre-MS
    inlist.create_pre_main_sequence_model(True)
    inlist.initial_mass(1.0) # arbitrarily chosen envelope mass
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
    inlist.save_final_model(f'env_for_compress.mod{id_str}')

    return inlist

def helper_merger_MS_HeWD_merge_for_compress(argdict):
    """
    merge envelope with cool HeWD in order to compress the HeWD's outer layers
    """
    data_path = argdict['data_path']

    id_str = argdict['id_str']

    script = MesaPythonScript(name='merge_for_compress',
            template=os.path.join(info.qol_path, 'mesa/templates/scripts/call_create_shell_burning_remnant.py'),
            const_args=['excise', 'change_m', data_path],
            prereqs=[f'cool_he_wd.mod{id_str}', f'env_for_compress.mod{id_str}'],
            products=[f'cool_he_wd_with_env.mod{id_str}'],
            data_path=data_path)

    return script

def helper_merger_MS_HeWD_ringdown_for_compress(argdict):
    """
    run remnant into HSE
    """
    enable_pgstar = argdict['enable_pgstar']
    mesh_delta_coeff = argdict['mesh_delta_coeff']
    data_path = argdict['data_path']

    id_str = argdict['id_str']
    
    inlist = MesaInlist(name='ringdown_for_compress', data_path=data_path)
    if enable_pgstar:
        inlist.enable_pgstar()
    inlist.save_pgstar(write_path='Grid1/ringdown_for_compress/')
    inlist.use_qol_pgstar()

    # Write GYRE model files
    inlist.write_gyre_data_with_profile()

    inlist.load_model(f'cool_he_wd_with_env.mod{id_str}')
    inlist.set_Zbase(0.02)

    # average composition of outer layers for write-out
    inlist.surface_avg_abundance_dq(1e-2)

    # enable hydro with drag
    # disable burning
    inlist.enable_hydrodynamics()
    inlist.add_hydrodynamical_drag(drag_coefficient=1.)
    inlist.energy_eqn_option('eps_grav')
    inlist.min_timestep_limit(1e-12)
    inlist.use_gold_tolerances(False)
    inlist.convergence_ignore_equL_residuals(True)
    inlist.disable_nuclear_burning()
    inlist.disable_mixing()

    inlist.mesh_delta_coeff(mesh_delta_coeff)

    inlist.max_age(1e4) # chosen somewhat arbitrarily
    inlist.save_final_model(f'cool_he_wd_with_env_rungdown.mod{id_str}')

    return inlist

def helper_merger_RG_HeWD_disable_hydro_for_he_wd_with_env(argdict):
    """
    simply read in compressed HeWD model with envelope and save it again but with hydro turned off
    """
    enable_pgstar = argdict['enable_pgstar']
    data_path = argdict['data_path']

    id_str = argdict['id_str']

    inlist = MesaInlist(name='disable_hydro_for_he_wd_with_env', data_path=data_path)
    if enable_pgstar:
        inlist.enable_pgstar()
    inlist.save_pgstar(write_path='Grid1/disable_hydro_for_he_wd_with_env/')
    inlist.use_qol_pgstar()
    
    inlist.disable_hydrodynamics()
    inlist.load_model(f'cool_he_wd_with_env_rungdown.mod{id_str}')
    inlist.set_Zbase(0.02)

    inlist.max_model_number(1)
    inlist.save_final_model(f'cool_he_wd_with_env_hse.mod{id_str}')

    return inlist

def helper_merger_MS_HeWD_remove_compress_env(argdict):
    """
    excise artificial envelope used to compress outer layers of cool HeWD
    """
    data_path = argdict['data_path']

    id_str = argdict['id_str']

    bad_frac_thresh = 1e-3
    min_boundary_fraction = 0.1

    script = MesaPythonScript(name='remove_compress_env',
            template=os.path.join(info.qol_path, 'mesa/templates/scripts/excise_envelope.py'),
            const_args=['he', bad_frac_thresh, min_boundary_fraction, data_path],
            prereqs=[f'cool_he_wd_with_env_hse.mod{id_str}'],
            products=[f'cool_he_wd_compress_dirty.mod{id_str}'],
            data_path=data_path)

    return script

def helper_merger_RG_HeWD_clean_compressed_he_wd(argdict):
    """
    remove "contaminant" species
    """
    data_path = argdict['data_path']

    id_str = argdict['id_str0'] = argdict['id_str']

    script = MesaPythonScript(name='clean_compressed_he_wd',
        template=os.path.join(info.qol_path, 'mesa/templates/scripts/call_replace_elements.py'),
        const_args=['he4', data_path,
                    'h1', 'he3', 'n13', 'o14', 'o15', 'f17', 'f18', 'ne18'],
        prereqs=[f'cool_he_wd_compress_dirty.mod{id_str}'],
        products=[f'cool_he_wd_compress.mod{id_str}'],
        data_path=data_path)

    return script

def helper_merger_RG_HeWD_outer_core_inner_bc(argdict):
    """
    we will try to construct the core model by creating a helium pre-MS model to "warm" degeneracy
    and putting it on the HeWD
    """
    Mcore_in_Msun = argdict['Mcore_in_Msun']
    data_path = argdict['data_path']

    id_str = argdict['id_str']

    # Generate envelope boundary conditions
    script = MesaPythonScript(name='outer_core_inner_bc',
                template=os.path.join(info.qol_path, 'mesa/templates/scripts/call_create_env_inlist_from_core.py'),
                const_args=[Mcore_in_Msun, data_path], prereqs=[f'cool_he_wd_compress.mod{id_str}'], products=[f'inlist_outer_core_inner_bc{id_str}'],
                data_path=data_path)

    return script

def helper_merger_RG_HeWD_make_outer_core(argdict):
    """
    create outer core model (shell-burning He object)
    """
    enable_pgstar = argdict['enable_pgstar']
    net_name = argdict['net_name']
    Mcore_in_Msun = argdict['Mcore_in_Msun']
    mesh_delta_coeff = argdict['mesh_delta_coeff']
    data_path = argdict['data_path']

    id_str = argdict['id_str']

    inlist = MesaInlist(name='make_outer_core', data_path=data_path)
    if enable_pgstar:
        inlist.enable_pgstar()
    inlist.save_pgstar(write_path='Grid1/make_outer_core/')
    inlist.use_qol_pgstar()

    # set inner boundary condition
    inlist.read_extra_inlist(namelist='star_job', rel_path=f'inlist_outer_core_inner_bc{id_str}',
                             category='relax inner BC to accommodate core model')

    # initialize as pre-MS
    inlist.create_pre_main_sequence_model(True)
    inlist.initial_mass(Mcore_in_Msun)
    inlist.initial_y(0.98)
    inlist.initial_z(0.02)
    inlist.set_Zbase(0.02)
    inlist.change_net(net_name)

    # convergence
    inlist.energy_eqn_option('dedt')
    inlist.use_gold_tolerances(False)
    inlist.convergence_ignore_equL_residuals(True)
    inlist.set_max_num_retries(3000)
    inlist.limit_for_rel_error_in_energy_conservation(-1.)
    inlist.mesh_delta_coeff(mesh_delta_coeff)

    # stop after some preset age
    # inlist.max_model_number(1)
    inlist.reset_age()
    inlist.max_age(1)

    argdict['id_str'] += f'_Mcore{Mcore_in_Msun:.3f}'
    id_str = argdict['id_str']

    inlist.save_final_model(f'hot_outer_core_dirty.mod{id_str}')

    return inlist

def helper_merger_RG_HeWD_clean_outer_core(argdict):
    """
    remove "contaminant" species
    """
    data_path = argdict['data_path']
    
    id_str = argdict['id_str']

    script = MesaPythonScript(name='clean_outer_core',
        template=os.path.join(info.qol_path, 'mesa/templates/scripts/call_replace_elements.py'),
        const_args=['he4', data_path,
                    'h1', 'he3', 'n13', 'o14', 'o15', 'f17', 'f18', 'ne18'],
        prereqs=[f'hot_outer_core_dirty.mod{id_str}'],
        products=[f'hot_outer_core.mod{id_str}'],
        data_path=data_path)

    return script

def helper_merger_RG_HeWD_accrete_env(argdict):
    """
    accrete envelope onto clean outer core before doing any "merging"
    """
    enable_pgstar = argdict['enable_pgstar']
    Menv_in_Msun = argdict['Menv_in_Msun']
    Mcore_in_Msun = argdict['Mcore_in_Msun']
    MWD_in_Msun = argdict['MWD_in_Msun']
    mesh_delta_coeff = argdict['mesh_delta_coeff']
    data_path = argdict['data_path']

    id_str = argdict['id_str']

    inlist = MesaInlist(name='accrete_env', data_path=data_path)
    if enable_pgstar:
        inlist.enable_pgstar()
    inlist.save_pgstar(write_path='Grid1/accrete_env/')
    inlist.use_qol_pgstar()

    inlist.load_model(f'hot_outer_core.mod{id_str}')

    inlist.set_Zbase(0.02)

    # disable composition changes
    inlist.disable_mixing()
    inlist.disable_dxdt_from_nuclear_burning()
    inlist.max_abar_for_burning(3) # allow pp chain to happen but disable helium burning (so the core has to become marginally degenerate)

    # add material with solar composition
    # note that MESA sets He3/Y = 1e-4
    He3_div_Y = 1e-4
    Z = 0.02
    Y = 0.24 + 2 * Z

    Mdot = 1e-4 # set to allow convergence somehow

    inlist.gain_mass(max_star_mass_for_gain=Menv_in_Msun+Mcore_in_Msun+MWD_in_Msun, mass_change=Mdot,
                 accrete_same_as_surface=False, accretion_h1=1. - Y - Z, accretion_h2=0.,
                 accretion_he3=He3_div_Y * Y, accretion_he4=(1-He3_div_Y) * Y, accretion_zfracs=3)

    # convergence
    # inlist.energy_eqn_option('dedt')
    # inlist.use_gold_tolerances(False)
    # inlist.convergence_ignore_equL_residuals(True)
    # inlist.set_max_num_retries(3000)
    # inlist.limit_for_rel_error_in_energy_conservation(-1.)
    inlist.mesh_delta_coeff(mesh_delta_coeff)

    # cool until core is at the threshold of degeneracy
    inlist.eta_center_limit(1.)

    argdict['id_str'] += f'_Menv{Menv_in_Msun:.3f}'
    id_str = argdict['id_str']

    inlist.save_final_model(f'hot_outer_core_and_env.mod{id_str}')

    return inlist

def helper_merger_RG_HeWD_merge(argdict):
    data_path = argdict['data_path']

    id_str = argdict['id_str']
    id_str0 = argdict['id_str0']

    script = MesaPythonScript(name='merge',
        template=os.path.join(info.qol_path, 'mesa/templates/scripts/call_create_shell_burning_remnant.py'),
        const_args=['default', 'change_m', data_path],
        prereqs=[f'cool_he_wd_compress.mod{id_str0}', f'hot_outer_core_and_env.mod{id_str}'],
        products=[f'remnant_init.mod{id_str}'],
        data_path=data_path)

    return script

def helper_merger_RG_HeWD_remnant_ringdown(argdict):
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
    inlist.enable_hydrodynamics()
    inlist.add_hydrodynamical_drag(drag_coefficient=1.)
    inlist.energy_eqn_option('eps_grav')
    inlist.min_timestep_limit(1e-12)
    inlist.use_gold_tolerances(False)
    inlist.convergence_ignore_equL_residuals(True)
    inlist.disable_mixing()
    # inlist.disable_nuclear_burning()
    inlist.max_abar_for_burning(3) # allow pp chain to happen but disable helium burning
    inlist.disable_dxdt_from_nuclear_burning()

    inlist.limit_for_rel_error_in_energy_conservation(-1.)

    inlist.mesh_delta_coeff(mesh_delta_coeff)

    # set arbitrary ringdown time
    argdict['id_str'] += f'_trd{ringdown_time_yr:.0f}'
    id_str = argdict['id_str']

    inlist.max_age(ringdown_time_yr)
    inlist.save_final_model(f'remnant_hse.mod{id_str}')

    # initial timestep
    inlist.set_initial_dt(years_for_initial_dt=1e-8)

    return inlist

def helper_merger_RG_HeWD_ignite_he(argdict):
    """
    helium flash: make this a separate inlist to disable remeshing during this stage
    """
    enable_pgstar = argdict['enable_pgstar']
    rgb_wind = argdict['rgb_wind']
    disable_hydro_after_ringdown = argdict['disable_hydro_after_ringdown']
    mesh_delta_coeff = argdict['mesh_delta_coeff']
    data_path = argdict['data_path']

    id_str = argdict['id_str']

    inlist = MesaInlist(name='ignite_he', data_path=data_path)
    if enable_pgstar:
        inlist.enable_pgstar()
    inlist.save_pgstar(write_path='Grid1/ignite_he/')
    inlist.use_qol_pgstar()

    if disable_hydro_after_ringdown:
        inlist.disable_hydrodynamics()

    # Write GYRE model files
    inlist.write_gyre_data_with_profile()

    inlist.load_model(f'remnant_hse.mod{id_str}')
    inlist.he_core_boundary_h1_fraction(1e-3)
    inlist.co_core_boundary_he4_fraction(1e-3)

    inlist.set_Zbase(0.02)

    inlist.mesh_delta_coeff(mesh_delta_coeff)

    # important! disable remeshing and initialize timestep to something low
    # new note: seems to not be necessary if you allow some thermalization across the entropy jump in the core during ringdown
    # inlist.okay_to_remesh(False)
    # inlist.set_initial_dt(1e-8)

    # relax some convergence conditions to get through degenerate helium ignition
    inlist.min_timestep_limit(1e-12)
    inlist.energy_eqn_option('eps_grav')
    inlist.use_gold_tolerances(False)

    inlist.convergence_ignore_equL_residuals(True)
    inlist.set_max_num_retries(3000)
    inlist.limit_for_rel_error_in_energy_conservation(-1.)

    # average composition of outer layers for write-out
    inlist.surface_avg_abundance_dq(1e-2)

    # wind
    if rgb_wind:
        inlist.cool_wind_RGB(scheme='Reimers', scaling_factor=0.5)
        inlist.cool_wind_AGB(scheme='Blocker', scaling_factor=0.1)
        inlist.RGB_to_AGB_wind_switch(1e-4)

    # stop at ZACHeB
    inlist.stop_at_phase_ZACHeB()

    argdict['id_str'] += f'_rgbw{int(rgb_wind)}_dhar{int(disable_hydro_after_ringdown)}'
    id_str = argdict['id_str']

    inlist.save_final_model(f'zacheb.mod{id_str}')

    return inlist

def helper_merger_RG_HeWD_zacheb_to_co_wd(argdict):
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

    inlist.load_model(f'zacheb.mod{id_str}')
    inlist.he_core_boundary_h1_fraction(1e-3)
    inlist.co_core_boundary_he4_fraction(1e-3)

    inlist.set_Zbase(0.02)

    inlist.energy_eqn_option('eps_grav')
    inlist.mesh_delta_coeff(mesh_delta_coeff)

    # predictive mixing
    inlist.add_predictive_mix_zone(predictive_zone_type='any', predictive_zone_loc='core', predictive_bdy_loc='top',
                            predictive_superad_thresh=0.01, predictive_avoid_reversal='he4')

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
    inlist.save_final_model(f'hot_co_wd.mod{id_str}')

    return inlist

def helper_merger_RG_HeWD_cool_co_wd_early(argdict):
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

    tho_string_dict = {'Kippenhahn': 'KRT80', 'Traxler_Garaud_Stellmach_11': 'TGS11', 'Brown_Garaud_Stellmach_13': 'BGS13'}
    argdict['id_str'] += f'_sc{alpha_semiconvection:.2f}_th{thermohaline_coeff:.1f}_tho{tho_string_dict[thermohaline_option]}_il{int(include_late)}_ov{int(include_overshoot)}'
    id_str = argdict['id_str']

    # convergence?
    inlist.use_gold_tolerances(False)

    # stop after cool down enough, if include_late
    if include_late:
        inlist.log_Teff_lower_limit(4.2)
    else:
        inlist.max_age(1e10)
    inlist.save_final_model(f'cool_co_wd_early.mod{id_str}')

    # increase write-out interval to get the moment of crystallization
    inlist.profile_interval(1)

    # also stop at "crystallization"
    # inlist.gamma_center_limit(175)

    return inlist

def helper_merger_RG_HeWD_cool_co_wd_late(argdict):
    """
    cool leftover CO WD through "late" stages

    disable:
    - Ledoux criterion -- this fixes an unphysical helium mixing event at cooler WD temperatures
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
    # inlist.use_table_atmosphere('WD_tau_25')

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

    # average composition of outer layers for write-out
    inlist.surface_avg_abundance_dq(1e-2)

    # wind
    inlist.cool_wind_RGB(scheme='Reimers', scaling_factor=0.5)
    inlist.cool_wind_AGB(scheme='Blocker', scaling_factor=0.1)
    inlist.RGB_to_AGB_wind_switch(1e-4)
    inlist.cool_wind_full_on_T(9.99e9) # force usage of cool wind only
    inlist.hot_wind_full_on_T(1.e10)

    # mixing
    inlist.gravitational_settling(diffusion_class_representatives=['h1', 'he3', 'he4', 'c12', 'o16', 'ne20', 'ne22', 'mg26'],
            diffusion_use_cgs_solver=True,
            show_diffusion_info=True,
            diffusion_steps_hard_limit=2000,
            diffusion_maxsteps_for_isolve=2000)

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
    inlist.max_num_profile_models(1000)

    # stop after a long time, if needed... okay to fail here, if sufficiently cooled
    inlist.max_age(1e10)
    inlist.log_Teff_lower_limit(3.7)
    inlist.save_final_model(f'cool_co_wd_late.mod{id_str}')

    return inlist
