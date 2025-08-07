"""
Purpose: remove the envelope shells off of a model.
This is different than something like relax_mass because we are not doing anything
to keep the model a valid model, and are simply just slicing off the shells.
"""
import sys
from qol.mesa.read import read_mod

# skip argv[0], which is the script path itself
orig_mod_fname = sys.argv[1]
write_mod_fname = sys.argv[2]
coretype = sys.argv[3]

# bad_frac_thresh is he_core_boundary_h1_fraction, co_core_boundary_he4_fraction, or one_core_boundary_he4_c12_fraction
bad_frac_thresh = sys.argv[4] if len(sys.argv) >= 5 else 0.1
min_boundary_fraction = sys.argv[5] if len(sys.argv) >= 6 else 0.1

assert coretype in ['he', 'co', 'one']
orig_mod_fname = f'data/{orig_mod_fname}'
write_mod_fname = f'data/{write_mod_fname}'

orig_mod = read_mod(orig_mod_fname)

match coretype:
    case 'he':
        core_mass = orig_mod.get_he_core_mass(he_core_boundary_h1_fraction=bad_frac_thresh, min_boundary_fraction=min_boundary_fraction)
    case 'co':
        core_mass = orig_mod.get_co_core_mass(co_core_boundary_he4_fraction=bad_frac_thresh, min_boundary_fraction=min_boundary_fraction)
    case 'one':
        core_mass = orig_mod.get_one_core_mass(one_core_boundary_he4_c12_fraction=bad_frac_thresh, min_boundary_fraction=min_boundary_fraction)
    case _:
        raise ValueError("Only coretype 'he', 'co', and 'one' supported for now.")

write_mod = orig_mod[orig_mod.M_in_Msun <= core_mass]
write_mod.write_model(write_mod_fname)
