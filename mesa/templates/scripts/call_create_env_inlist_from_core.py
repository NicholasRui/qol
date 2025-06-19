# Wrapper script to call create_env_inlist_from_core within a MESA directory
import sys
from qol.mesa.create import create_env_inlist_from_core

# skip argv[0], which is the script path itself
M_env_Msun = sys.argv[1]
core_mod_fname = sys.argv[2]
inlist_fname = sys.argv[3]

run_path = '.'

# Trim off "inlist_" from beginning of inlist_fname, since it will be readded by the MesaInlist __init__()
assert len(inlist_fname) > 7
assert inlist_fname[:7] == 'inlist_'
task_name = inlist_fname[7:]

try:
    M_env_Msun = float(M_env_Msun.lower().replace('d', 'e'))
except:
    raise ValueError(f'M_env_Msun -- cannot interpret as float: {M_env_Msun}')

create_env_inlist_from_core(run_path=run_path,
    task_name=task_name,
    core_mod_fname=core_mod_fname,
    M_env_Msun=M_env_Msun)
