# Wrapper script to call create_shell_burning_remant within a MESA directory
import sys
from qol.mesa.create import create_shell_burning_remnant

# skip argv[0], which is the script path itself
core_mod_fname = sys.argv[1]
env_mod_fname = sys.argv[2]
write_mod_fname = sys.argv[3]

interface_setting = 'excise'
readjust_setting = 'change_m'

create_shell_burning_remnant(write_mod_fname, core_mod_fname, env_mod_fname,
            interface_setting=interface_setting, readjust_setting=readjust_setting)
