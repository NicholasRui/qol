# Methods handling controls about opacity

def set_Zbase(self, Zbase=0.02):
    namelist = 'kap'

    self.add_control(namelist=namelist,
            control='use_Type2_opacities', value=True)
    self.add_control(namelist=namelist,
            control='Zbase', value=Zbase)

