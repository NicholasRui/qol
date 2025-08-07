# Wrapper script to call create_shell_burning_remant within a MESA directory
import sys
from qol.mesa.sculpt import create_shell_burning_remnant

# skip argv[0], which is the script path itself
assert len(sys.argv) == 6

core_mod_fname = sys.argv[1]
env_mod_fname = sys.argv[2]
write_mod_fname = sys.argv[3]
interface_setting = sys.argv[4]
readjust_setting = sys.argv[5]

# Append 'data/' to the front of all of these
core_mod_fname = f'data/{core_mod_fname}'
env_mod_fname = f'data/{env_mod_fname}'
write_mod_fname = f'data/{write_mod_fname}'

create_shell_burning_remnant(write_mod_fname, core_mod_fname, env_mod_fname,
            interface_setting=interface_setting, readjust_setting=readjust_setting)
