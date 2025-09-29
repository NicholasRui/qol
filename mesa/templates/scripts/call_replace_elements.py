# Wrapper script to call replace_elements within a MESA directory
import sys
from qol.mesa.sculpt import replace_elements

import os

# skip argv[0], which is the script path itself
init_mod_fname = sys.argv[1]
write_mod_fname = sys.argv[2]
new_species = sys.argv[3]
absdir = sys.argv[4]
old_species_arr = sys.argv[5:]

run_path = '.'

# Append data path to the front of all of these
write_mod_fname = os.path.join(absdir, write_mod_fname)
init_mod_fname = os.path.join(absdir, init_mod_fname)

replace_elements(write_mod_fname,
                 init_mod_fname,
                 new_species, *old_species_arr)
