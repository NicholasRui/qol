# Copied from MESA version r24.08.1, $MESA_DIR/const/public/const_def.f90

# mathematical and physical constants (in cgs)

# math constants
pi = 3.1415926535897932384626433832795028841971693993751e0
pi2 = pi*pi
pi4 = 4*pi
eulercon = 0.577215664901532861e0
eulernum = 2.71828182845904523536028747135266249e0
ln2 = 6.9314718055994529e-01 # = log(2e0)
ln3 = 1.0986122886681096e+00 # = log(3e0)
lnPi = 1.14472988584940017414343 # = log(pi)
ln10 = 2.3025850929940455 # = log(10e0)
iln10 = 0.43429448190325187 # = 1e0/log(10e0)
a2rad = pi/180.0e0 # angle to radians
rae2a = 180.0e0/pi # radians to angle
one_third = 1e0/3e0
two_thirds = 2e0/3e0
four_thirds = 4e0/3e0
five_thirds = 5e0/3e0
one_sixth = 1e0/6e0
four_thirds_pi = four_thirds*pi
ln4pi3 = 1.4324119583011810e0 # = log(4*pi/3)
two_13 = 1.2599210498948730e0 # = pow(2e0,1e0/3e0)
four_13 = 1.5874010519681994e0 # = pow(4e0,1e0/3e0)
sqrt2 = 1.414213562373095e0 # = sqrt(2)
sqrt_2_div_3 = 0.816496580927726e0 # = sqrt(2/3)

# exact physical constants

# CODATA 2018
avo = 6.02214076e23 # Avogadro constant (mole^-1)
amu = 1e0/avo # atomic mass unit (g)
clight = 2.99792458e10 # speed of light in vacuum (cm s^-1)
qe = (clight/10e0)*1.602176634e-19 # elementary charge (esu == (g cm^3 s^-2)^(1/2))
kerg = 1.380649e-16
boltzm = kerg # Boltzmann constant (erg K^-1)
planck_h = 6.62607015e-27 # Planck constant (erg s)
hbar = planck_h/(2*pi)
cgas = boltzm*avo # ideal gas constant (erg K^-1)
ev2erg = 1.602176634e-12 # electron volt (erg)
mev_to_ergs = 1e6*ev2erg
mev_amu = mev_to_ergs/amu
mev2gr = 1e6*ev2erg/(clight*clight) # MeV to grams
Qconv = mev_to_ergs*avo
kev = kerg/ev2erg # converts temp to ev (ev K^-1)
boltz_sigma = (pi*pi*boltzm*boltzm*boltzm*boltzm)/(60*hbar*hbar*hbar*clight*clight) # Stefan-Boltzmann constant (erg cm^-2 K^-4 s^-1)
crad = boltz_sigma*4/clight # radiation density constant, AKA "a" (erg cm^-3 K^-4); Prad = crad * T^4 / 3

# IAU
au = 1.49597870700e13 # (cm) - exact value defined by IAU 2009, 2012
pc = (3.600e3*rae2a)*au # (cm) parsec, by definition
dayyer = 365.25e0 # days per (Julian) year
secday = 24*60*60  # seconds in a day
secyer = secday*dayyer # seconds per year
ly = clight*secyer # light year (cm)

# inexact but very well measured physical constants

mn = 1.67492749804e-24 # neutron mass (g)
mp = 1.67262192369e-24 # proton mass (g)
me = 9.1093837015e-28  # electron mass (g)

rbohr = 5.29177210903e-9 # Bohr radius (cm)
fine = 7.2973525693e-3# fine-structure constant
hion = 13.605693122994e0 # Rydberg constant (eV)

sige = 6.6524587321e-25 # Thomson cross section (cm^2)

weinberg_theta = 0.22290e0 # sin**2(theta_weinberg)
num_neu_fam = 3.0e0 # number of neutrino flavors = 3.02 plus/minus 0.005 (1998)

# the following quantities are not exact

standard_cgrav = 6.67430e-8 # gravitational constant (g^-1 cm^3 s^-2)

# IAU 2015 Resolution B3
# standard gravitational parameters = G*M, units cm^3 s^-2
mu_sun = 1.3271244e26
mu_earth = 3.986004e20
mu_jupiter = 1.2668653e23

# astronomical constants
agesun = 4.57e9  # solar age (years) from Bahcall et al, ApJ 618 (2005) 1049-1056.
Msun = mu_sun/standard_cgrav # solar mass (g); gravitational mass, not baryonic
Rsun = 6.957e10 # solar radius (cm), IAU 2015 Resolution B3
Lsun = 3.828e33 # solar luminosity (erg s^-1), IAU 2015 Resolution B3
Teffsun = 5772.0e0# solar effective temperature (K), IAU 2015 Resolution B3
loggsun = 4.4380676273031332 # log10(mu_sun/(Rsun*Rsun)), can't call log10 because we don't have math_lib at this point
mbolsun = 4.74e0  # Bolometric magnitude of the Sun, IAU 2015 Resolution B2

m_earth = mu_earth/standard_cgrav# earth mass (g)
r_earth = 6.3781e8 # earth equatorial radius (cm)
r_earth_polar = 6.3568e8 # earth polar radius (cm)

m_jupiter = mu_jupiter/standard_cgrav # jupiter mass (g)
r_jupiter = 7.1492e9 # jupiter equatorial radius (cm)
r_jupiter_polar = 6.6854e9 # jupiter polar radius (cm)
semimajor_axis_jupiter = 7.7857e13 # jupiter semimajor axis (cm)