from qol.mesa.launcher import *
import qol.info as info

# controlling dX limits during main sequence
dX_div_X_limit_min_X = 1e-4
dX_div_X_limit = 1e-2
dX_nuc_drop_min_X_limit = 1e-4
dX_nuc_drop_limit = 1e-2
# controlling change in log mstar
delta_lg_star_mass_limit = delta_lg_star_mass_hard_limit = 1e-5

def make_MS_single(root_path, # absolute path in which to write directory
                   M_initial_in_Msun,
                   initial_z=0.02, # metallicity
                   overshoot_f=0.015, # overshoot parameter f
                   overshoot_f0=0.005, # overshoot parameter f0
                   enable_pgstar=True,
                   source_sdk=True, # manually activate sdk, since Caltech HPC doesn't seem to like it
                   to_lower_rgb=False, # add task evolving to lower RGB)
                   alpha_semiconvection=0., # semiconvection parameter
                   save_directory=True, # if False, don't save directory
                   data_path='data/',
                   ):
    """
    single MS star model for comparison purposes
    """
    if data_path != 'data/': # TODO
        raise NotImplementedError()

    run_name = f'M{M_initial_in_Msun:.2f}_z{initial_z:.4f}_osf{overshoot_f:.4f}_osf0{overshoot_f0:.4f}'
    run_path = f'{root_path}/{run_name}'

    argdict = {'root_path': root_path,
               'M_initial_in_Msun': M_initial_in_Msun,
               'Xcen_accrete': None,
               'initial_z': initial_z,
               'overshoot_f': overshoot_f,
               'overshoot_f0': overshoot_f0,
               'enable_pgstar': enable_pgstar,
               'source_sdk': source_sdk,
               'to_lower_rgb': to_lower_rgb,
               'alpha_semiconvection': alpha_semiconvection,
               'data_path': data_path,}

    # Put it together
    work = MesaWorkingDirectory(run_path=run_path)
    work.copy_history_columns_list(f'{info.qol_path}mesa/resources/r24.08.1/history_columns_qol.list')
    work.copy_profile_columns_list(f'{info.qol_path}mesa/resources/r24.08.1/profile_columns_qol.list')
    work.load_qol_pgstar()

    work.add_task(helper_MS_mass_changer_zams_to_mt(argdict))
    if to_lower_rgb:
        work.add_task(helper_MS_mass_changer_tams_to_lower_rgb(argdict, last_task='zams_to_mt'))

    if save_directory:
        work.save_directory(slurm_job_name=run_name, grant_perms=True, source_sdk=source_sdk, data_path=data_path)

    return work


def make_MS_mass_changer(root_path, # absolute path in which to write directory
                         M_initial_in_Msun, # initial mass
                         M_final_in_Msun, # final mass
                         Xcen_accrete, # core hydrogen fraction at which time mass change starts
                         log_abs_Mdot_accrete=-5, # mass change rate in Msun per yr
                         initial_z=0.02, # metallicity
                         overshoot_f=0.015, # overshoot parameter f
                         overshoot_f0=0.005, # overshoot parameter f0
                         enable_pgstar=True,
                         source_sdk=True, # manually activate sdk, since Caltech HPC doesn't seem to like it
                         to_lower_rgb=False, # add task evolving to lower RGB
                         alpha_semiconvection=0., # semiconvection parameter
                         save_directory=True, # if False, don't save directory
                         data_path='data/',
                         ):
    """
    Evolve an MS star for a bit before adding / removing mass from it
    """
    if data_path != 'data/': # TODO
        raise NotImplementedError()

    run_name = f'M{M_initial_in_Msun:.2f}to{M_final_in_Msun:.2f}atX{Xcen_accrete:.2f}_dM{log_abs_Mdot_accrete:.1f}_z{initial_z:.4f}_osf{overshoot_f:.4f}_osf0{overshoot_f0:.4f}_sc{alpha_semiconvection:.4f}'
    run_path = f'{root_path}/{run_name}'

    argdict = {'root_path': root_path,
               'M_initial_in_Msun': M_initial_in_Msun,
               'M_final_in_Msun': M_final_in_Msun,
               'Xcen_accrete': Xcen_accrete,
               'log_abs_Mdot_accrete': log_abs_Mdot_accrete,
               'initial_z': initial_z,
               'overshoot_f': overshoot_f,
               'overshoot_f0': overshoot_f0,
               'enable_pgstar': enable_pgstar,
               'source_sdk': source_sdk,
               'to_lower_rgb': to_lower_rgb,
               'alpha_semiconvection': alpha_semiconvection,
               'data_path': data_path,}

    # Put it together
    work = MesaWorkingDirectory(run_path=run_path)
    work.copy_history_columns_list(f'{info.qol_path}mesa/resources/r24.08.1/history_columns_qol.list')
    work.copy_profile_columns_list(f'{info.qol_path}mesa/resources/r24.08.1/profile_columns_qol.list')
    work.load_qol_pgstar()

    work.add_task(helper_MS_mass_changer_zams_to_mt(argdict))
    work.add_task(helper_MS_mass_changer_mt_to_tams(argdict))
    if to_lower_rgb:
        work.add_task(helper_MS_mass_changer_tams_to_lower_rgb(argdict, last_task='mt_to_tams'))

    if save_directory:
        work.save_directory(slurm_job_name=run_name, grant_perms=True, source_sdk=source_sdk, data_path=data_path)

    return work

def helper_MS_mass_changer_zams_to_mt(argdict):
    """
    Evolve up to the mass transfer time
    """
    enable_pgstar = argdict['enable_pgstar']
    M_initial_in_Msun = argdict['M_initial_in_Msun']
    Xcen_accrete = argdict['Xcen_accrete']
    initial_z = argdict['initial_z']
    overshoot_f = argdict['overshoot_f']
    overshoot_f0 = argdict['overshoot_f0']
    alpha_semiconvection = argdict['alpha_semiconvection']
    data_path = argdict['data_path']

    inlist = MesaInlist('zams_to_mt')
    if enable_pgstar:
        inlist.enable_pgstar()
    inlist.save_pgstar(write_path='Grid1/zams_to_mt/')
    inlist.use_qol_pgstar()

    # Write GYRE
    inlist.write_gyre_data_with_profile()

    # Initial parameters
    inlist.initial_mass(M_initial_in_Msun)
    inlist.initial_z(initial_z)
    inlist.set_Zbase(initial_z)

    inlist.alpha_semiconvection(alpha_semiconvection)

    # Timestepping
    inlist.set_dX_limits(dX_div_X_limit_min_X=dX_div_X_limit_min_X,
                     dX_div_X_limit=dX_div_X_limit, # adjust composition only very slowly
                     dX_nuc_drop_min_X_limit=dX_nuc_drop_min_X_limit,
                     dX_nuc_drop_limit=dX_nuc_drop_limit,
                     )
    inlist.set_delta_lg_star_mass_limits(delta_lg_star_mass_limit=delta_lg_star_mass_limit,
                                  delta_lg_star_mass_hard_limit=delta_lg_star_mass_hard_limit)

    # Termination conditions
    inlist.stop_at_phase_TAMS()
    if Xcen_accrete is not None:
        inlist.add_xa_central_lower_limit(xa_central_lower_limit_species='h1',
                                        xa_central_lower_limit=Xcen_accrete) # stop this evolution when core drops to desired hydrogen fraction
    inlist.save_final_model('zams_to_mt.mod')

    # Overshoot
    inlist.add_overshoot_zone(overshoot_scheme='exponential',
                            overshoot_zone_type='any', overshoot_zone_loc='core', overshoot_bdy_loc='top',
                            overshoot_f=overshoot_f, overshoot_f0=overshoot_f0)

    return inlist

def helper_MS_mass_changer_mt_to_tams(argdict):
    """
    Perform mass transfer and evolve to TAMS
    """
    enable_pgstar = argdict['enable_pgstar']
    M_initial_in_Msun = argdict['M_initial_in_Msun']
    M_final_in_Msun = argdict['M_final_in_Msun']
    log_abs_Mdot_accrete = argdict['log_abs_Mdot_accrete']
    initial_z = argdict['initial_z']
    overshoot_f = argdict['overshoot_f']
    overshoot_f0 = argdict['overshoot_f0']
    alpha_semiconvection = argdict['alpha_semiconvection']
    data_path = argdict['data_path']
    
    inlist = MesaInlist('mt_to_tams')
    if enable_pgstar:
        inlist.enable_pgstar()
    inlist.save_pgstar(write_path='Grid1/mt_to_tams/')
    inlist.use_qol_pgstar()

    # Initial conditions
    inlist.load_model('zams_to_mt.mod')
    inlist.set_Zbase(initial_z)

    inlist.alpha_semiconvection(alpha_semiconvection)

    # Timestepping
    inlist.set_dX_limits(dX_div_X_limit_min_X=dX_div_X_limit_min_X,
                     dX_div_X_limit=dX_div_X_limit, # adjust composition only very slowly
                     dX_nuc_drop_min_X_limit=dX_nuc_drop_min_X_limit,
                     dX_nuc_drop_limit=dX_nuc_drop_limit,
                     )
    inlist.set_delta_lg_star_mass_limits(delta_lg_star_mass_limit=delta_lg_star_mass_limit,
                                  delta_lg_star_mass_hard_limit=delta_lg_star_mass_hard_limit)

    # Write GYRE
    inlist.write_gyre_data_with_profile()

    # Mass transfer
    if M_final_in_Msun > M_initial_in_Msun:
        inlist.gain_mass(max_star_mass_for_gain=M_final_in_Msun, mass_change=10 ** log_abs_Mdot_accrete)
    elif M_final_in_Msun < M_initial_in_Msun:
        inlist.lose_mass(min_star_mass_for_loss=M_final_in_Msun, mass_change=-10 ** log_abs_Mdot_accrete)

    # Termination conditions
    inlist.stop_at_phase_TAMS()
    inlist.save_final_model('mt_to_tams.mod')

    # Overshoot
    inlist.add_overshoot_zone(overshoot_scheme='exponential',
                            overshoot_zone_type='any', overshoot_zone_loc='core', overshoot_bdy_loc='top',
                            overshoot_f=overshoot_f, overshoot_f0=overshoot_f0)

    return inlist

def helper_MS_mass_changer_tams_to_lower_rgb(argdict, last_task):
    enable_pgstar = argdict['enable_pgstar']
    initial_z = argdict['initial_z']
    alpha_semiconvection = argdict['alpha_semiconvection']
    data_path = argdict['data_path']

    inlist = MesaInlist('tams_to_lower_rgb')
    if enable_pgstar:
        inlist.enable_pgstar()
    inlist.save_pgstar(write_path='Grid1/tams_to_lower_rgb/')
    inlist.use_qol_pgstar()

    # Initial conditions
    inlist.set_Zbase(initial_z)
    match last_task:
        case 'zams_to_mt':
            inlist.load_model('zams_to_mt.mod')
        case 'mt_to_tams':
            inlist.load_model('mt_to_tams.mod')
        case _:
            raise ValueError(f"invalid last_task {last_task}: must be either 'zams_to_mt' or 'mt_to_tams'")

    inlist.alpha_semiconvection(alpha_semiconvection)

    # Write GYRE
    inlist.write_gyre_data_with_profile()

    # Termination conditions
    inlist.he_core_mass_limit(0.20)
    inlist.he_core_boundary_h1_fraction(1e-3)
    inlist.save_final_model('tams_to_lower_rgb.mod')

    return inlist
