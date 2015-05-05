#!/usr/bin/env python
import numpy as np
import cmath
import math

phi = 1.4388 # *v1 / T [1/cm / K]
def nbar (w1, T): # w1 is in 1/cm and T is in K. Factor phi will take care of conversion to 1/m
#    global w1
    return 1/(np.exp(phi*w1/T)-1)

def zeta (nb, GoD):
    return 1.0/(1.0+1j*GoD/(nb+1))

def beta (zt, nb):
    return zt*np.log(1-zt*nb/(nb+1))

def dog0(nb, Dog0, GoD): # delta over gamma_0
    return Dog0*( 1.0/2.0+2*nb*GoD/( (nb+1)**2+GoD**2 ) )

def gog0 (nb, Dog0, GoD): # gamma over gamma_0
    return 2*nb*(nb+1)*GoD*Dog0/( (nb+1)**2+GoD**2)

def mu (nb, GoD):
    m = GoD - 1j/2-cmath.sqrt(GoD**2-1.0/4.0-1j*GoD*(2*nb+1))
    if (m.real < 0):     
        m = GoD - 1j/2+cmath.sqrt(GoD**2-1.0/4.0-1j*GoD*(2*nb+1))
        if (m.real < 0):
            print "Achtung!"
    return m

def alpha (nb, m):
    return -1j*nb/m

def HizhAllT(t, w1, Gamma, delta, gamma0, T):
    enbar = nbar(w1, T)
    mew = mu(enbar, Gamma/delta)
    alfa = alpha(enbar, mew)
    Dog0 = delta/gamma0
    z = np.exp( -np.pi*abs(t) + 1j*2*np.pi*( (Dog0*t*(enbar+1.0/2.0))-(Dog0*t*alfa*(enbar+1)/(1+alfa)) +1j*(enbar+1)/(mew*(alfa+1))*np.log(1+alfa-alfa*np.exp(-mew*Dog0*t)) ) )
    return z

def HizhLowT(t, w1, Gamma, delta, gamma0, T):
    # dimensionless time t*gamma0
    enbar = nbar(w1, T)
    zetah = zeta(enbar, Gamma/delta)
    betah = beta(zetah, enbar)
    Dog0 = delta/gamma0
    GoD = Gamma/delta
    dog = dog0(enbar, Dog0, GoD)
    gog = gog0(enbar, Dog0, GoD)
    # NOTE: not using betah which amounts to the temperature dependent scale
    z = np.exp(betah+-np.pi*abs(t) + 2*np.pi*t*(1j*dog-1.0/2.0*gog) - zetah*np.log(1-zetah*enbar/(enbar+1)*np.exp(2*np.pi*1j*dog*t-2*np.pi*GoD*Dog0*t) ) )
    return z
