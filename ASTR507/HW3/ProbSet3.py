"""
From Erics ASTR507 HW3:
Goal: From COBE spectrum of CMB, estimate the max allowed value of the chemical
potential of photons in units kT_cmb by allowing the temp and chem potential 
to vary, find where chi^2 of fit increases to the point that fit is poor.

State confidence of your upper limit
"""

import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from scipy import optimize as opt
from astropy import units as u # might be useful later

infile = './firas_monopole_spec_v1.txt'

# Column 1 = frequency from Table 4 of Fixsen et al., units = cm^-1
# Column 2 = FIRAS monopole spectrum computed as the sum
#             of a 2.725 K BB spectrum and the
#             residual in column 3, units = MJy/sr
# Column 3 = residual monopole spectrum from Table 4 of Fixsen et al.,
#             units = kJy/sr
# Column 4 = spectrum uncertainty (1-sigma) from Table 4 of Fixsen et al.,
#             units = kJy/sr
# Column 5 = modeled Galaxy spectrum at the Galactic poles from Table 4 of
#             Fixsen et al., units = kJy/sr 
freq, spec, resid, err, galpole = np.loadtxt(infile,unpack=True)

# put in same units 
err /= 1000.
resid /= 1000.
galpole /= 1000.


# Constants
h = 6.62607e-34 # Js
c = 2.99792e10 # cm/s
k = 1.38065e-23 # J/K
T0 = 2.728  # K
mu0 = 0.
x = h * c * freq / (k * T0)

# Assign units
freq = freq / u.cm
spec = spec * 10**6 * u.jansky / u.sr
resid = resid * 10**3 *  u.jansky / u.sr
err = err * 10**3 *  u.jansky / u.sr
galpole = galpole * 10**3 *  u.jansky / u.sr

h = 6.62607e-34 * u.joule * u.
c = 2.99792e10 * u.cm / u.s
k = 1.38065e-23 * u.joule / u.K
T0 = 2.728  *u.K
mu0 = 0.
x = h * c * freq / (k * T0)

# 1 Jy = 1e-26 J/m^2

model = 2. * h * c * freq**3. / (np.e**(x+mu0) - 1.) #  J / cm^2