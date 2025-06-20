# Methods handling controls about termination conditions

def max_age(self, value):
    self.add_to_controls(category='termination conditions', comment='yr',
            control='max_age', value=value)

def Teff_upper_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='Teff_upper_limit', value=value)

def Teff_lower_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='Teff_lower_limit', value=value)

def max_model_number(self, value):
    self.add_to_controls(category='termination conditions',
            control='max_model_number', value=value)

def he_core_mass_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='he_core_mass_limit', value=value)

def co_core_mass_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='co_core_mass_limit', value=value)

def one_core_mass_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='one_core_mass_limit', value=value)
    
def fe_core_mass_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='fe_core_mass_limit', value=value)

def neutron_rich_core_mass_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='neutron_rich_core_mass_limit', value=value)

def stop_near_zams(self, Lnuc_div_L_zams_limit):
    """
    Stop near ZAMS, defined by Lnuc/Lzams
    """
    self.add_to_controls(category='termination conditions',
            control='stop_near_zams', value=True)
    self.add_to_controls(category='termination conditions', optional=True,
            control='Lnuc_div_L_zams_limit', value=Lnuc_div_L_zams_limit)

def stop_at_phase_PreMS(self):
    """
    Phase defined in $MESA_DIR/star/private/star_utils.f90:
    s% log_center_temperature > 5d0
    """
    self.add_to_controls(category='termination conditions',
            control='stop_at_phase_PreMS', value=True)
    
def stop_at_phase_ZAMS(self):
    """
    Phase defined in $MESA_DIR/star/private/star_utils.f90:
    s% L_nuc_burn_total >= s% L_phot*s% Lnuc_div_L_zams_limit
    """
    self.add_to_controls(category='termination conditions',
            control='stop_at_phase_ZAMS', value=True)

def stop_at_phase_IAMS(self):
    """
    Phase defined in $MESA_DIR/star/private/star_utils.f90:
    center_h1 <= 0.3d0
    """
    self.add_to_controls(category='termination conditions',
            control='stop_at_phase_IAMS', value=True)
    
def stop_at_phase_TAMS(self):
    """
    Phase defined in $MESA_DIR/star/private/star_utils.f90:
    center_h1 <= 1d-6
    """
    self.add_to_controls(category='termination conditions',
            control='stop_at_phase_TAMS', value=True)

def stop_at_phase_He_Burn(self):
    """
    Phase defined in $MESA_DIR/star/private/star_utils.f90:
    s% L_by_category(i3alf) > 1d2
    """
    self.add_to_controls(category='termination conditions',
            control='stop_at_phase_He_Burn', value=True)

def stop_at_phase_ZACHeB(self):
    """
    Phase defined in $MESA_DIR/star/private/star_utils.f90:
    s% center_eps_burn(i3alf) > 1d2
    """
    self.add_to_controls(category='termination conditions',
            control='stop_at_phase_ZACHeB', value=True)
    
def stop_at_phase_TACHeB(self):
    """
    Phase defined in $MESA_DIR/star/private/star_utils.f90:
    center_he4 <= 1d-4
    """
    self.add_to_controls(category='termination conditions',
            control='stop_at_phase_TACHeB', value=True)
    
def stop_at_phase_TP_AGB(self):
    """
    Phase defined in $MESA_DIR/star/private/star_utils.f90:
    center_he4 < 1d-4 .and. &
            s% he_core_mass - s% co_core_mass <= 0.1d0 .and. &
            any(s% burn_he_conv_region(1:s% num_conv_boundaries))
    """
    self.add_to_controls(category='termination conditions',
            control='stop_at_phase_TP_AGB', value=True)

def stop_at_phase_C_Burn(self):
    """
    Phase defined in $MESA_DIR/star/private/star_utils.f90:
    s% L_by_category(icc) > 1d2 .and. center_he4 < 1d-4
    """
    self.add_to_controls(category='termination conditions',
            control='stop_at_phase_C_Burn', value=True)

def stop_at_phase_Ne_Burn(self):
    """
    Phase defined in $MESA_DIR/star/private/star_utils.f90:
    s% L_by_category(i_burn_ne) > 1d2 .and. center_he4 < 1d-4
    """
    self.add_to_controls(category='termination conditions',
            control='stop_at_phase_Ne_Burn', value=True)

def stop_at_phase_O_Burn(self):
    """
    Phase defined in $MESA_DIR/star/private/star_utils.f90:
    (s% L_by_category(ioo) + s% L_by_category(ico)) > 1d2 .and. center_he4 < 1d-4
    """
    self.add_to_controls(category='termination conditions',
            control='stop_at_phase_O_Burn', value=True)
    
def stop_at_phase_Si_Burn(self):
    """
    Phase defined in $MESA_DIR/star/private/star_utils.f90:
    s% L_by_category(i_burn_si) > 1d2
    """
    self.add_to_controls(category='termination conditions',
            control='stop_at_phase_Si_Burn', value=True)

def stop_at_phase_WDCS(self):
    """
    Phase defined in $MESA_DIR/star/private/star_utils.f90:
    s% photosphere_logg > 6d0
    """
    self.add_to_controls(category='termination conditions',
            control='stop_at_phase_WDCS', value=True)

def add_xa_central_lower_limit(self, xa_central_lower_limit_species, xa_central_lower_limit):
    assert xa_central_lower_limit <= 1

    self.add_to_controls(category='termination conditions',
            control='xa_central_lower_limit_species', value=xa_central_lower_limit_species)
    self.add_to_controls(category='termination conditions',
            control='xa_central_lower_limit', value=xa_central_lower_limit)

def gamma_center_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='gamma_center_limit', value=value)

def eta_center_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='eta_center_limit', value=value)
    
def log_center_density_upper_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='AAAAAA', value=value)

def log_center_density_lower_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='AAAAAA', value=value)

def log_center_temp_upper_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='AAAAAA', value=value)
    
def log_center_temp_lower_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='AAAAAA', value=value)

def surface_accel_div_grav_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='AAAAAA', value=value)

def log_max_temp_upper_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='AAAAAA', value=value)

def log_max_temp_lower_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='AAAAAA', value=value)

def center_entropy_upper_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='AAAAAA', value=value)

def center_entropy_lower_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='AAAAAA', value=value)
    
def max_entropy_upper_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='AAAAAA', value=value)
    
def max_entropy_lower_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='AAAAAA', value=value)

def HB_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='HB_limit', value=value)
    
def star_mass_min_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='star_mass_max_limit', value=value)

def remnant_mass_min_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='remnant_mass_min_limit', value=value)

def ejecta_mass_max_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='ejecta_mass_max_limit', value=value)
    
def envelope_mass_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='envelope_mass_limit', value=value)

def envelope_fraction_left_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='envelope_fraction_left_limit', value=value)

def xmstar_min_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='xmstar_min_limit', value=value)
    
def xmstar_max_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='xmstar_max_limit', value=value)

def he_layer_mass_lower_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='he_layer_mass_lower_limit', value=value)

def abs_diff_lg_LH_lg_Ls_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='abs_diff_lg_LH_lg_Ls_limit', value=value)
    
def photosphere_r_upper_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='photosphere_r_upper_limit', value=value)

def photosphere_r_lower_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='photosphere_r_lower_limit', value=value)

def photosphere_m_upper_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='photosphere_m_upper_limit', value=value)
    
def photosphere_m_lower_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='photosphere_m_lower_limit', value=value)

def photosphere_m_sub_M_center_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='photosphere_m_sub_M_center_limit', value=value)

def log_Teff_upper_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='log_Teff_upper_limit', value=value)
    
def log_Teff_lower_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='log_Teff_lower_limit', value=value)

def log_Tsurf_upper_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='log_Tsurf_upper_limit', value=value)

def log_Tsurf_lower_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='log_Tsurf_lower_limit', value=value)
    
def log_Rsurf_upper_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='log_Rsurf_upper_limit', value=value)

def log_Rsurf_lower_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='log_Rsurf_lower_limit', value=value)

def log_L_upper_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='log_L_upper_limit', value=value)
    
def log_L_lower_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='log_L_lower_limit', value=value)

def log_g_upper_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='log_g_upper_limit', value=value)

def log_g_lower_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='log_g_lower_limit', value=value)
    
def log_Psurf_upper_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='log_Psurf_upper_limit', value=value)

def log_Psurf_lower_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='log_Psurf_lower_limit', value=value)

def log_Dsurf_upper_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='log_Dsurf_upper_limit', value=value)
    
def log_Dsurf_lower_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='log_Dsurf_lower_limit', value=value)

def power_nuc_burn_upper_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='power_nuc_burn_upper_limit', value=value)

def power_h_burn_upper_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='power_h_burn_upper_limit', value=value)

def power_he_burn_upper_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='power_he_burn_upper_limit', value=value)

def power_z_burn_upper_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='power_z_burn_upper_limit', value=value)

def power_nuc_burn_lower_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='power_nuc_burn_lower_limit', value=value)

def power_h_burn_lower_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='power_h_burn_lower_limit', value=value)

def power_he_burn_lower_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='power_he_burn_lower_limit', value=value)

def power_z_burn_lower_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='power_z_burn_lower_limit', value=value)

def max_abs_rel_run_E_err(self, value):
    self.add_to_controls(category='termination conditions',
            control='max_abs_rel_run_E_err', value=value)

def max_number_retries(self, value):
    self.add_to_controls(category='termination conditions',
            control='max_number_retries', value=value)

def min_timestep_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='min_timestep_limit', value=value)
    
def center_Ye_lower_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='center_Ye_lower_limit', value=value)

def center_R_lower_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='center_R_lower_limit', value=value)

def fe_core_infall_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='fe_core_infall_limit', value=value)

def non_fe_core_rebound_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='non_fe_core_rebound_limit', value=value)

def v_div_csound_max_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='v_div_csound_max_limit', value=value)

def v_div_csound_surf_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='v_div_csound_surf_limit', value=value)
    
def v_surf_div_v_kh_upper_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='v_surf_div_v_kh_upper_limit', value=value)

def v_surf_div_v_kh_lower_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='v_surf_div_v_kh_lower_limit', value=value)

def v_surf_div_v_esc_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='v_surf_div_v_esc_limit', value=value)

def v_surf_kms_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='v_surf_kms_limit', value=value)

def Lnuc_div_L_zams_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='Lnuc_div_L_zams_limit', value=value)

def Lnuc_div_L_upper_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='Lnuc_div_L_upper_limit', value=value)

def Lnuc_div_L_lower_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='Lnuc_div_L_lower_limit', value=value)

def set_gamma1_limit(self, gamma1_limit, gamma1_limit_max_q=None):
    self.add_to_controls(category='termination conditions',
            control='gamma1_limit', value=gamma1_limit)
    self.add_to_controls(category='termination conditions', optional=True,
            control='gamma1_limit_max_q', value=gamma1_limit_max_q)

def gamma1_limit_max_v_div_vesc(self, value):
    self.add_to_controls(category='termination conditions',
            control='gamma1_limit_max_v_div_vesc', value=value)

def set_Pgas_div_P_limit(self, Pgas_div_P_limit, Pgas_div_P_limit_max_q=None):
    self.add_to_controls(category='termination conditions',
            control='Pgas_div_P_limit', value=Pgas_div_P_limit)
    self.add_to_controls(category='termination conditions', optional=True,
            control='Pgas_div_P_limit_max_q', value=Pgas_div_P_limit_max_q)

def peak_burn_vconv_div_cs_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='peak_burn_vconv_div_cs_limit', value=value)

def omega_div_omega_crit_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='omega_div_omega_crit_limit', value=value)

def delta_nu_lower_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='delta_nu_lower_limit', value=value)

def delta_nu_upper_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='delta_nu_upper_limit', value=value)

def shock_mass_upper_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='shock_mass_upper_limit', value=value)

def mach1_mass_upper_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='mach1_mass_upper_limit', value=value)

def delta_Pg_lower_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='delta_Pg_lower_limit', value=value)

def delta_Pg_upper_limit(self, value):
    self.add_to_controls(category='termination conditions',
            control='delta_Pg_upper_limit', value=value)

def stop_when_reach_this_cumulative_extra_heating(self, value):
    self.add_to_controls(category='termination conditions',
            control='stop_when_reach_this_cumulative_extra_heating', value=value)
