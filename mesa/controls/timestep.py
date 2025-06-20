# Methods handling controls about timestepping
from collections.abc import Iterable

def set_dX_limits(self,
                  # list arguments
                  dX_limit_species=['h1', 'he4'],
                  dX_limit_min_X=None,
                  dX_limit=None,
                  dX_hard_limit=None,
                  dX_div_X_limit_min_X=None,
                  dX_div_X_limit=None,
                  dX_div_X_hard_limit=None,
                  dX_div_X_at_high_T_limit=None,
                  dX_div_X_at_high_T_hard_limit=None,
                  dX_div_X_at_high_T_limit_lgT_min=None,
                  dX_decreases_only=None, # bool
                  # float args
                  dX_nuc_drop_min_X_limit=None,
                  dX_nuc_drop_limit_at_high_T=None,
                  dX_nuc_drop_limit=None,
                  dX_nuc_drop_hard_limit=None,
                  dX_nuc_drop_min_yrs_for_dt=None,
                  # int args
                  dX_nuc_drop_max_A_limit=None,
                  ):
    """
    controls all dX related timestepping controls

    assert that all array arguments MUST be the same length
    if a single value is given, substitute for an array of the same length
    """
    list_args = {'dX_limit_min_X': dX_limit_min_X,
                 'dX_limit': dX_limit,
                 'dX_hard_limit': dX_hard_limit,
                 'dX_div_X_limit_min_X': dX_div_X_limit_min_X,
                 'dX_div_X_limit': dX_div_X_limit,
                 'dX_div_X_hard_limit': dX_div_X_hard_limit,
                 'dX_div_X_at_high_T_limit': dX_div_X_at_high_T_limit,
                 'dX_div_X_at_high_T_hard_limit': dX_div_X_at_high_T_hard_limit,
                 'dX_div_X_at_high_T_limit_lgT_min': dX_decreases_only,
                 } # compare these against dX_limit_species
    category = 'dX limits'
    
    num_limit_species = len(dX_limit_species)

    for control, list_arg in list_args.items():
        if list_arg is not None:
            if isinstance(list_arg, Iterable) and not isinstance(list_arg, str):
                assert len(list_arg) == num_limit_species
            elif list_arg is not None:
                list_args[control] = num_limit_species * [list_arg]

    # Now add all list arguments, if no exceptions triggered
    self.add_list_to_controls(category=category,
            control='dX_limit_species', values=dX_limit_species)

    for control, list_arg in list_args.items():
        if list_arg is not None:
            self.add_list_to_controls(category=category, comment=dX_limit_species, # label species for convenience
                    control=control, values=list_arg)

    # single value arguments
    self.add_to_controls(category=category, optional=True,
            control='dX_nuc_drop_min_X_limit', value=dX_nuc_drop_min_X_limit)
    self.add_to_controls(category=category, optional=True,
            control='dX_nuc_drop_limit_at_high_T', value=dX_nuc_drop_limit_at_high_T)
    self.add_to_controls(category=category, optional=True,
            control='dX_nuc_drop_limit', value=dX_nuc_drop_limit)
    self.add_to_controls(category=category, optional=True,
            control='dX_nuc_drop_hard_limit', value=dX_nuc_drop_hard_limit)
    self.add_to_controls(category=category, optional=True,
            control='dX_nuc_drop_min_yrs_for_dt', value=dX_nuc_drop_min_yrs_for_dt)
    self.add_to_controls(category=category, optional=True,
            control='dX_nuc_drop_max_A_limit', value=dX_nuc_drop_max_A_limit)

def set_delta_lg_star_mass_limits(self, delta_lg_star_mass_limit,
                                  delta_lg_star_mass_hard_limit=None):
    category = 'dlogmstar limits'

    self.add_to_controls(category=category,
            control='delta_lg_star_mass_limit', value=delta_lg_star_mass_limit)
    self.add_to_controls(category='delta_lg_star_mass limits', optional=True,
            control='delta_lg_star_mass_hard_limit', value=delta_lg_star_mass_hard_limit)

