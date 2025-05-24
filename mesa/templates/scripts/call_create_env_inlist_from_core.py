# Wrapper script to call create_env_inlist_from_core within a MESA directory
import sys
from qol.mesa.create import create_env_inlist_from_core

# skip argv[0], which is the script path itself
M_env_Msun = sys.argv[1]
core_mod_fname = sys.argv[2]
task_name = sys.argv[3]

run_path = '.'

try:
    M_env_Msun = float(M_env_Msun.lower().replace('d', 'e'))
except:
    raise ValueError(f'M_env_Msun -- cannot interpret as float: {M_env_Msun}')

create_env_inlist_from_core(run_path=run_path,
    name=task_name,
    core_mod_fname=core_mod_fname,
    M_env_Msun=M_env_Msun)
