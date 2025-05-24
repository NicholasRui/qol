from qol.mesa.launcher import *
import qol.info as info

def make_MS_mass_changer(root_path, # absolute path in which to write directory
                         M_initial_in_Msun, # initial mass
                         M_final_in_Msun, # final mass
                         Xcen_accrete, # core hydrogen fraction at which time mass change starts
                         log_abs_Mdot_accrete=-5, # mass change rate in Msun per yr
                         initial_z=0.02, # metallicity
                         overshoot_f=0.015, # overshoot parameter f
                         overshoot_f0=0.005, # overshoot parameter f0
                         enable_pgstar=False,
                         source_sdk=True, # manually activate sdk, since Caltech HPC doesn't seem to like it
                         ):
    """
    Evolve an MS star for a bit before adding / removing mass from it
    """
    run_name = f'MS{M_initial_in_Msun:.2f}to{M_final_in_Msun:.2f}atX{Xcen_accrete:.2f}_dM{log_abs_Mdot_accrete:.1f}_z{initial_z:.4f}_osf{overshoot_f:.4f}_osf0{overshoot_f0:.4f}'
    run_path = f'{root_path}/{run_name}'

    task_zams_to_mt = helper_MS_mass_changer_zams_to_mt(enable_pgstar=enable_pgstar, M_initial_in_Msun=M_initial_in_Msun, Xcen_accrete=Xcen_accrete, initial_z=initial_z, overshoot_f=overshoot_f, overshoot_f0=overshoot_f0)
    task_mt_to_tams = helper_MS_mass_changer_mt_to_tams(enable_pgstar=enable_pgstar, M_initial_in_Msun=M_initial_in_Msun, M_final_in_Msun=M_final_in_Msun, log_abs_Mdot_accrete=log_abs_Mdot_accrete, initial_z=initial_z, overshoot_f=overshoot_f, overshoot_f0=overshoot_f0)

    # Put it together
    work = MesaWorkingDirectory(run_path=run_path)
    work.copy_history_columns_list(f'{info.qol_path}mesa/resources/r24.08.1/history_columns.list')
    work.copy_profile_columns_list(f'{info.qol_path}mesa/resources/r24.08.1/profile_columns.list')
    work.load_qol_pgstar()

    work.add_task(task_zams_to_mt)
    work.add_task(task_mt_to_tams)

    work.save_directory(slurm_job_name=run_name, grant_perms=True, source_sdk=source_sdk)

    return work

def helper_MS_mass_changer_zams_to_mt(enable_pgstar, M_initial_in_Msun, Xcen_accrete, initial_z, overshoot_f, overshoot_f0):
    """
    Evolve up to the mass transfer time
    """    
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

    # Termination conditions
    inlist.stop_at_phase_TAMS()
    inlist.add_xa_central_lower_limit(xa_central_lower_limit_species='h1',
                                    xa_central_lower_limit=Xcen_accrete) # stop this evolution when core drops to desired hydrogen fraction
    inlist.save_final_model('zams_to_mt.mod')

    # Overshoot
    inlist.add_overshoot_zone(overshoot_scheme='exponential',
                            overshoot_zone_type='any', overshoot_zone_loc='core', overshoot_bdy_loc='top',
                            overshoot_f=overshoot_f, overshoot_f0=overshoot_f0)

    return inlist

def helper_MS_mass_changer_mt_to_tams(enable_pgstar, M_initial_in_Msun, M_final_in_Msun, log_abs_Mdot_accrete, initial_z, overshoot_f, overshoot_f0):
    """
    Perform mass transfer and evolve to TAMS
    """
    inlist = MesaInlist('mt_to_tams')
    if enable_pgstar:
        inlist.enable_pgstar()
    inlist.save_pgstar(write_path='Grid1/mt_to_tams/')
    inlist.use_qol_pgstar()

    # Initial conditions
    inlist.load_model('zams_to_mt.mod')
    inlist.set_Zbase(initial_z)

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
