# Methods handling controls about boundary conditions

import qol.tools.formatter as formatter
import qol.mesa.const as const

def relax_to_inner_BC(self, M_new_Msun=None, R_center_Rsun=None, L_center_Lsun=None,
        dlgm_per_step=None, dlgR_per_step=None, dlgL_per_step=None,
        relax_M_center_dt=None, relax_R_center_dt=None, relax_L_center_dt=None):
    """
    Relax to inner boundary condition
    """
    if M_new_Msun is not None:
        category = 'relax to inner BC: mass'
        comment = 'total mass in Msun'

        self.add_to_star_job(category=category,
                control='relax_M_center', value=True)
        self.add_to_star_job(category=category,
                control='new_mass', value=M_new_Msun, comment=comment)
        
        self.add_to_star_job(category=category, optional=True,
                control='dlgm_per_step', value=dlgm_per_step)
        self.add_to_star_job(category=category, comment='sec', optional=True,
                control='relax_M_center_dt', value=relax_M_center_dt)

    if R_center_Rsun is not None:
        category = 'relax to inner BC: radius'
        R_center_cgs = R_center_Rsun * const.Rsun
        comment = f'cm = {formatter.to_fortran(R_center_Rsun)} Rsun'

        self.add_to_star_job(category=category,
                control='relax_R_center', value=True)
        self.add_to_star_job(category=category,
                control='new_R_center', value=R_center_cgs, comment=comment)
        
        self.add_to_star_job(category=category, optional=True,
                control='dlgR_per_step', value=dlgR_per_step)
        self.add_to_star_job(category=category, comment='sec', optional=True,
                control='relax_R_center_dt', value=relax_R_center_dt)    
    
    if L_center_Lsun is not None:
        category = 'relax to inner BC: luminosity'
        L_center_cgs = L_center_Lsun * const.Lsun
        comment = f'erg s-1 = {formatter.to_fortran(L_center_Lsun)} Lsun'

        self.add_to_star_job(category=category,
                control='relax_L_center', value=True)
        self.add_to_star_job(category=category,
                control='new_L_center', value=L_center_cgs, comment=comment)
        
        self.add_to_star_job(category=category, optional=True,
                control='dlgL_per_step', value=dlgL_per_step)
        self.add_to_star_job(category=category, comment='sec', optional=True,
                control='relax_L_center_dt', value=relax_L_center_dt)

