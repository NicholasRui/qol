# Wrapper script to call create_env_inlist_from_core within a MESA directory
import sys
from qol.mesa.sculpt import replace_elements

# skip argv[0], which is the script path itself
init_mod_fname = sys.argv[1]
write_mod_fname = sys.argv[2]
new_species = sys.argv[3]
old_species_arr = sys.argv[4:]

run_path = '.'

# Append 'data/' to the front of all of these
write_mod_fname = f'data/{write_mod_fname}'
init_mod_fname = f'data/{init_mod_fname}'

replace_elements(write_mod_fname,
                 init_mod_fname,
                 new_species, *old_species_arr)
