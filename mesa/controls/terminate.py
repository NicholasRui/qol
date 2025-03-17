# Methods handling controls about termination conditions

def max_age(self, value):
    self.add_control(namelist='controls', category='termination conditions', comment='yr',
            control='max_age', value=value)

def Teff_upper_limit(self, value):
    self.add_control(namelist='controls', category='termination conditions',
            control='Teff_upper_limit', value=value)

def Teff_lower_limit(self, value):
    self.add_control(namelist='controls', category='termination conditions',
            control='Teff_lower_limit', value=value)

def max_model_number(self, value):
    self.add_control(namelist='controls', category='termination conditions',
            control='max_model_number', value=value)

def he_core_mass_limit(self, value):
    self.add_control(namelist='controls', category='termination conditions',
            control='he_core_mass_limit', value=value)

def co_core_mass_limit(self, value):
    self.add_control(namelist='controls', category='termination conditions',
            control='co_core_mass_limit', value=value)

def one_core_mass_limit(self, value):
    self.add_control(namelist='controls', category='termination conditions',
            control='one_core_mass_limit', value=value)

def stop_near_zams(self, Lnuc_div_L_zams_limit):
    """
    Stop near ZAMS, defined by Lnuc/Lzams
    """
    self.add_control(namelist='controls', category='termination conditions',
            control='stop_near_zams', value=True)
    self.add_control(namelist='controls', category='termination conditions', optional=True,
            control='Lnuc_div_L_zams_limit', value=Lnuc_div_L_zams_limit)

def stop_at_phase_PreMS(self):
    """
    Phase defined in $MESA_DIR/star/private/star_utils.f90:
    s% log_center_temperature > 5d0
    """
    self.add_control(namelist='controls', category='termination conditions',
            control='stop_at_phase_PreMS', value=True)
    
def stop_at_phase_ZAMS(self):
    """
    Phase defined in $MESA_DIR/star/private/star_utils.f90:
    s% L_nuc_burn_total >= s% L_phot*s% Lnuc_div_L_zams_limit
    """
    self.add_control(namelist='controls', category='termination conditions',
            control='stop_at_phase_ZAMS', value=True)

def stop_at_phase_IAMS(self):
    """
    Phase defined in $MESA_DIR/star/private/star_utils.f90:
    center_h1 <= 0.3d0
    """
    self.add_control(namelist='controls', category='termination conditions',
            control='stop_at_phase_IAMS', value=True)
    
def stop_at_phase_TAMS(self):
    """
    Phase defined in $MESA_DIR/star/private/star_utils.f90:
    center_h1 <= 1d-6
    """
    self.add_control(namelist='controls', category='termination conditions',
            control='stop_at_phase_TAMS', value=True)

def stop_at_phase_He_Burn(self):
    """
    Phase defined in $MESA_DIR/star/private/star_utils.f90:
    s% L_by_category(i3alf) > 1d2
    """
    self.add_control(namelist='controls', category='termination conditions',
            control='stop_at_phase_He_Burn', value=True)

def stop_at_phase_ZACHeB(self):
    """
    Phase defined in $MESA_DIR/star/private/star_utils.f90:
    s% center_eps_burn(i3alf) > 1d2
    """
    self.add_control(namelist='controls', category='termination conditions',
            control='stop_at_phase_ZACHeB', value=True)
    
def stop_at_phase_TACHeB(self):
    """
    Phase defined in $MESA_DIR/star/private/star_utils.f90:
    center_he4 <= 1d-4
    """
    self.add_control(namelist='controls', category='termination conditions',
            control='stop_at_phase_TACHeB', value=True)
    
def stop_at_phase_TP_AGB(self):
    """
    Phase defined in $MESA_DIR/star/private/star_utils.f90:
    center_he4 < 1d-4 .and. &
            s% he_core_mass - s% co_core_mass <= 0.1d0 .and. &
            any(s% burn_he_conv_region(1:s% num_conv_boundaries))
    """
    self.add_control(namelist='controls', category='termination conditions',
            control='stop_at_phase_TP_AGB', value=True)

def stop_at_phase_C_Burn(self):
    """
    Phase defined in $MESA_DIR/star/private/star_utils.f90:
    s% L_by_category(icc) > 1d2 .and. center_he4 < 1d-4
    """
    self.add_control(namelist='controls', category='termination conditions',
            control='stop_at_phase_C_Burn', value=True)

def stop_at_phase_Ne_Burn(self):
    """
    Phase defined in $MESA_DIR/star/private/star_utils.f90:
    s% L_by_category(i_burn_ne) > 1d2 .and. center_he4 < 1d-4
    """
    self.add_control(namelist='controls', category='termination conditions',
            control='stop_at_phase_Ne_Burn', value=True)

def stop_at_phase_O_Burn(self):
    """
    Phase defined in $MESA_DIR/star/private/star_utils.f90:
    (s% L_by_category(ioo) + s% L_by_category(ico)) > 1d2 .and. center_he4 < 1d-4
    """
    self.add_control(namelist='controls', category='termination conditions',
            control='stop_at_phase_O_Burn', value=True)
    
def stop_at_phase_Si_Burn(self):
    """
    Phase defined in $MESA_DIR/star/private/star_utils.f90:
    s% L_by_category(i_burn_si) > 1d2
    """
    self.add_control(namelist='controls', category='termination conditions',
            control='stop_at_phase_Si_Burn', value=True)

def stop_at_phase_WDCS(self):
    """
    Phase defined in $MESA_DIR/star/private/star_utils.f90:
    s% photosphere_logg > 6d0
    """
    self.add_control(namelist='controls', category='termination conditions',
            control='stop_at_phase_WDCS', value=True)

