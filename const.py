import numpy as np

pi = np.pi

# MESA constants (cgs)
Msun = 1.9892e33
Rsun = 6.9598e10
Lsun = 3.8418e33
Teffsun = 5772.0

clight = 2.99792458e10
standard_cgrav = 6.67430e-8
qe = (clight/10e0)*1.602176634e-19
boltzm = 1.380649e-16

planck_h = 6.62607015e-27
hbar = planck_h/(2*pi)
boltz_sigma = (pi*pi*boltzm*boltzm*boltzm*boltzm)/(60*hbar*hbar*hbar*clight*clight)

au = 1.49597870700e13
pc = (3.600e3*180.0e0/pi)*au

mn = 1.67492749804e-24
mp = 1.67262192369e-24
me = 9.1093837015e-28
