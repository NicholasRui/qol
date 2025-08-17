# Load in all supported MESA preset configurations to make sure all tasks can be added without crashing

def run_mesa_tests():
    tests_total = 0
    tests_failed = 0

    print('============================')
    print('RUNNING TESTS: mesa_tests.py')
    print('============================')

    ###########################################################################################
    from qol.mesa.presets import make_MS_mass_changer

    tests_total += 1
    name = 'qol.mesa.presets.make_MS_mass_changer.make_MS_single'
    try:
        _ = make_MS_mass_changer.make_MS_single(root_path='.', M_initial_in_Msun=1., initial_z=0.02, overshoot_f=0.015, overshoot_f0=0.005, enable_pgstar=True, source_sdk=True, to_lower_rgb=True, alpha_semiconvection=4e-2, save_directory=False,)
        print(f'passed: {name}')
    except Exception as e:
        print(f'**FAILED**: {name}')
        print(e)
        tests_failed += 1

    tests_total += 1
    name = 'qol.mesa.presets.make_MS_mass_changer.make_MS_mass_changer'
    try:
        _ = make_MS_mass_changer.make_MS_mass_changer(root_path='.', M_initial_in_Msun=1., M_final_in_Msun=1.5, Xcen_accrete=0.5, log_abs_Mdot_accrete=-5, initial_z=0.02, overshoot_f=0.015, overshoot_f0=0.005, enable_pgstar=True, source_sdk=True, to_lower_rgb=True, alpha_semiconvection=4e-2, save_directory=False,)
        print(f'passed: {name}')
    except Exception as e:
        print(f'**FAILED**: {name}')
        print(e)
        tests_failed += 1

    ###########################################################################################
    from qol.mesa.presets import make_single_MS_HeWD

    tests_total += 1
    name = 'qol.mesa.presets.make_single_MS_HeWD.make_single_MS_HeWD'
    try:
        _ = make_single_MS_HeWD.make_single_MS_HeWD(root_path='.', MMS_in_Msun=1., net_name='cno_extras_o18_to_mg26.net', enable_pgstar=True, rgb_wind=True, alpha_semiconvection=4e-2, thermohaline_coeff=1., source_sdk=True, mesh_delta_coeff=1., save_directory=False,)
        print(f'passed: {name}')
    except Exception as e:
        print(f'**FAILED**: {name}')
        print(e)
        tests_failed += 1

    ###########################################################################################
    from qol.mesa.presets import make_merger_MS_HeWD

    tests_total += 1
    name = 'qol.mesa.presets.make_merger_MS_HeWD.make_merger_MS_HeWD'
    try:
        _ = make_merger_MS_HeWD.make_merger_MS_HeWD(root_path='.', MWD_in_Msun=0.4, MMS_in_Msun=0.4, T_WD=1e4, net_name='cno_extras_o18_to_mg26.net', ringdown_time_yr=1e4, disable_hydro_after_ringdown=True, enable_pgstar=True, rgb_wind=True, alpha_semiconvection=4e-2, thermohaline_coeff=1., source_sdk=True, mesh_delta_coeff=1., include_late=True, save_directory=False,)
        print(f'passed: {name}')
    except Exception as e:
        print(f'**FAILED**: {name}')
        print(e)
        tests_failed += 1

    ###########################################################################################
    from qol.mesa.presets import make_merger_RG_HeWD

    tests_total += 1
    name = 'qol.mesa.presets.make_merger_RG_HeWD.make_merger_RG_HeWD'
    try:
        _ = make_merger_RG_HeWD.make_merger_RG_HeWD(root_path='.', MWD_in_Msun=0.4, Mcore_in_Msun=0.2, Menv_in_Msun=1., T_WD=1e4, net_name='cno_extras_o18_to_mg26.net', ringdown_time_yr=1e3, disable_hydro_after_ringdown=True, enable_pgstar=True, rgb_wind=True, alpha_semiconvection=4e-2, thermohaline_coeff=1., source_sdk=True, mesh_delta_coeff=1., include_late=True, save_directory=False,)
        print(f'passed: {name}')
    except Exception as e:
        print(f'**FAILED**: {name}')
        print(e)
        tests_failed += 1

    print('============================')
    print(f'mesa_tests.py: failed = {tests_failed} / {tests_total}')
    print('============================')

    return tests_failed, tests_total


if __name__ == "__main__":
    _ = run_mesa_tests()
