#import ekruse as ek
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import pdb
from scipy import optimize as opt

infile = './firas_monopole_spec_v1.txt'

T = np.arange(10.)*5e-3 + 2.728
mu = np.arange(10.)*3e-5 - 15.e-5


# freq in 1/cm 
freq, spec, resid, err, galpole = np.loadtxt(infile,unpack=True)
# get the same units
err /= 1000.
resid /= 1000.
galpole /= 1000.

h = 6.62607e-34 # Js
c = 2.99792e10 # cm/s
k = 1.38065e-23 # J/K
T0 = 2.728  # K
mu0 = 0.
x = h * c * freq / (k * T0)

# 1 Jy = 1e-26 J/m^2

# this equation from the Fixsen 1996 paper. Why won't the units work? off by c?
#model = 2. * h * c**2. * freq**3. / (np.e**x - 1. ) # J / cm /s
# divide by c to get the right units
model = 2. * h * c * freq**3. / (np.e**(x+mu0) - 1.) #  J / cm^2
model *= 1.e4 # J / m^2
model *= 1.e26 # Jy
model /= 1.e6 # MJy

plt.close('all')
"""
plt.figure()
plt.plot(freq,spec)
plt.plot(freq,model)
plt.plot(freq,spec + galpole)
plt.xlim([3,7])
plt.ylim([330,390])
"""

def fchi(xx, *args):
	ffreq, fspec, ferr = args
	T,mu = xx
	h = 6.62607e-34 # Js
	c = 2.99792e10 # cm/s
	k = 1.38065e-23 # J/K
	x = h * c * ffreq / (k * T)
	model = 2. * h * c * ffreq**3. / (np.e**(x+mu) - 1.) #  J / cm^2
	model *= 1.e4 # J / m^2
	model *= 1.e26 # Jy
	model /= 1.e6 # MJy
	return np.sum(((model - fspec)/ferr)**2.)

def fchi_delT(xx, *args):
	mu, ffreq, fspec, ferr = args
	T = xx
	h = 6.62607e-34 # Js
	c = 2.99792e10 # cm/s
	k = 1.38065e-23 # J/K
	x = h * c * ffreq / (k * T)
	model = 2. * h * c * ffreq**3. / (np.e**(x+mu) - 1.) #  J / cm^2
	model *= 1.e4 # J / m^2
	model *= 1.e26 # Jy
	model /= 1.e6 # MJy
	return np.sum(((model - fspec)/ferr)**2.)


# find the minimum value
res = opt.fmin(fchi,[2.7,0.],args=(freq,spec,err))
minchi = fchi(res,freq,spec,err)
# best values
T0, mu0 = res

maxdiff = 170.
npts = 150
dchi = 0.
maxT = T0*1.
while dchi < maxdiff:
	maxT *= 1.000001
	dchi = fchi([maxT,mu0],freq,spec,err) - minchi

dchi = 0.
minT = T0*1.
while dchi < maxdiff:
	minT *= 0.999999
	dchi = fchi([minT,mu0],freq,spec,err) - minchi

dchi = 0.
minmu = mu0*1.
while dchi < maxdiff:
	minmu -= 0.000001
	dchi = fchi([T0,minmu],freq,spec,err) - minchi

dchi = 0.
maxmu = mu0*1.
while dchi < maxdiff:
	maxmu += 0.000001
	dchi = fchi([T0,maxmu],freq,spec,err) - minchi
	

T = np.linspace(minT,maxT,npts)
mu = np.linspace(minmu,maxmu,npts)


chisqs = np.zeros((len(T),len(mu)))
for ii, iT in enumerate(T):
	for jj, imu in enumerate(mu):
		x = h * c * freq / (k * iT)
		model = 2. * h * c * freq**3. / (np.e**(x+imu) - 1.) #  J / cm^2
		model *= 1.e4 # J / m^2
		model *= 1.e26 # Jy
		model /= 1.e6 # MJy

		chisqs[ii,jj] = np.sum(((model - spec)/err)**2.)

chisqs -= minchi
plt.figure()
levels = [2.3,6.17,11.8]
CS = plt.contour(mu,T,chisqs,levels,colors=('k','k','k'))
strs = ['1$\sigma$','2$\sigma$','3$\sigma$']
fmt = {}
for l,s in zip(CS.levels, strs):
	fmt[l] = s
plt.clabel(CS, inline=1,fmt=fmt)
plt.xlabel('$\mu$')
plt.ylabel('Temperature')
axis = plt.gca()
axis.get_yaxis().get_major_formatter().set_useOffset(False)

oned = np.array([1.,4.,9.])

# confidence regions
for ii in np.arange(len(oned)):
	Tsrch, musrch = np.where(chisqs < oned[ii])
	if ii == 0:
		avgt = (np.max(T[Tsrch]) - np.min(T[Tsrch]))/2.
		avgmu = (np.max(mu[musrch]) - np.min(mu[musrch]))/2.
		print 'T = '+str(T0) + ' $\pm$ ' + str(avgt) +  '; $\mu$ = ' + str(mu0)+  ' $\pm$ ' + str(avgmu)
	print str(ii+1) + ' sigma T:', np.min(T[Tsrch]), np.max(T[Tsrch])
	print str(ii+1) + ' sigma mu:', np.min(mu[musrch]), np.max(mu[musrch])

plt.savefig('./pset3pt1.pdf')
### part 3

def fermidirac(z, *args):
	# Using the approximation for the Fermi-Dirac functions from Aymerich-Humet 1981
	# arg 1 is j of F_j(z)
	# arg2 is x in F_j(z) == x. What value are we solving for
	
	j = args[0]
	# what does this function equal (where is the root?)
	eq = args[1]
	import numpy as np
	import scipy.special
	# to get in the same form as the paper
	z = np.log(z)
	if j == 1./2:
		a = 9.6
		b = 2.13
		c = 12./5
	elif j == 3./2:
		a = 14.9
		b = 2.64
		c = 9./4
	else:
		return 0

	return np.abs(1./(  (j+1.)*(2.**(j+1)) / (b + z + (np.abs(z-b)**c + a)**(1./c))**(j+1.)  +  e**(-z) / scipy.special.gamma(j+1.) ) - eq)


res = opt.fmin(fermidirac,[1.],args=(0.5,2.08))

print 'Fugacity: z = ', res[0]

print 'F_(3/2)(z) = ', fermidirac(res[0],3./2,0)
print 'F_(1/2)(z) = ', fermidirac(res[0],1./2,0)

z = np.linspace(0,10,10000)
plt.figure()
plt.plot(z, fermidirac(z,1./2,0)-2.08)






