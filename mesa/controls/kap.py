# Methods handling controls about opacity

def set_Zbase(self, Zbase=0.02):
    self.add_to_kap(control='use_Type2_opacities', value=True)
    self.add_to_kap(control='Zbase', value=Zbase)
