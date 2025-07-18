# This is the only file that should be modified by the user, except in development
import os
import warnings

####################
###  MESA PATHS  ###
####################
# TODO: come up with a better way to store this info and handle multiple MESA setups
# e.g., could change MESA_DIR in the ./mk and ./rn files each time we run MESA
mesa_paths = {}

if 'MESA_DIR' in os.environ:
    MESA_DIR = os.environ['MESA_DIR']

    if os.path.exists(f'{MESA_DIR}/data/version_number'):
        with open(f'{MESA_DIR}/data/version_number', 'r') as f:
            mesa_version = f.read().replace('\n', '')
        mesa_paths[mesa_version] = MESA_DIR
else:
    warnings.warn('No $MESA_DIR found!') # TODO: handle exception better later

# slurm job defaults
slurm_job_mail_user = 'nrui.mailing.list@gmail.com'
slurm_job_time_default = '7-00:00:00'
slurm_job_mem_per_cpu_default = '10G'
slurm_job_nodes_default = 10
slurm_job_ntasks_default = 2 * slurm_job_nodes_default
mesa_OMP_NUM_THREADS = slurm_job_ntasks_default
# slurm_job_nodes_default = 1
# slurm_job_ntasks_default = 1
# mesa_OMP_NUM_THREADS = 1
