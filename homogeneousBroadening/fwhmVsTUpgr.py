#!/usr/bin/env python
import numpy as np
import numpy.fft as ft
import cmath
import math
import matplotlib.pyplot as plt
import sys
import hizhnyakovfunc as hf
import lplspectra as lpl
import scipy.signal as sig
from matplotlib.backends.backend_pdf import PdfPages


def round(x, n):
    '''Rounds a float up to n significant figures.
    Returns a string.'''
    format = "%."+str(n)+"f"
    return format % x

# See Numerical Recepies in C: The Art of Scientific Computing 3rd Edition. Ch 13.9 for the 
# details of the method used here
# Trapezoidal order interpolation. See the reference for details.
def alpha0 (theta): # end points correction term
	result = 0
	if (abs(theta) <= 0.05): # use series representation of the coefficient
		result = -1.0/2.0+1.0/24.0*theta**2-1.0/720.0*theta**4+1.0/40320.0*theta**6+1j*theta*(1.0/6.0-1.0/120.0*theta**2+1.0/5040.0*theta**4-1.0/362880.0*theta**6)
	else: # use cosine/sine formulae
		result = ( np.cos(theta)-1+1j*(theta-np.sin(theta)) )/theta**2
	return result

def W(theta): # correction factor
	result = 0
	if (abs(theta) <= 0.05):
		result = 1.0 - 1.0/12.0*theta**2+1.0/360.0*theta**4-1.0/20160.0*theta**6
	else:
		result = 2*( 1-np.cos(theta) )/theta**2
	return result


class ZplSpectrum:
    '''A spectrum generated using Hizhnyakoff model.'''
    #parameters:
    T = 1.0 # K, Medium's temperature
    w1 = 1.0 # cm-1, quasi-local mode frequency
    Delta = 1.0 # cm-1, frequency change between two electronic states
    Gamma = 1.0 # *|Delta|
    gamma0 = 1.0 # natural linewidth, cm-1
    
    # Complex-valued function to be integrated
    def G(self, t):
        return hf.HizhLowT(t, self.w1, self.Gamma, self.Delta, self.gamma0, self.T)
    
    def Generate(self):
        '''Generates a spectrum for a given set of parameters.
        Spectrum is a two element list [wavenumber, intensity].'''

        # in the calculations below we need the following dimensionless parameters
        GoD = self.Gamma / self.Delta # Gamma over Delta
        Dog0 = self.Delta/self.gamma0 # Delta over gamma naught

        print "Critical temperature: "+\
        str( 2*1.4388*self.w1/np.log(abs(GoD)+1/abs(GoD)) )+\
        "\t Current temperature: "+\
        str( self.T )

#-------------------------A cool part starts here-----------------------------------------
#---Numeric Integration of highly oscillatory function: FT using ---------------------
        M = 2**20# number of sampling points in the time-domain
        K = 80.0 # total integration time ( in dimensionless variable related to gamma_0 )
        a = 0.0 # lower integration limit. Redundant but kept for clarity.
        b = K # upper integration limit
        dt = (b-a)/(M-1) # the size of the step in time domain
        fc = self.gamma0*(M-1)/(2*K) # in 1/cm. Crytical (Nyquist) frequency.
    # zero-padding of h_vector with fourfold oversampling
        over_sample = 4
        N = over_sample*M # Total number of sampling points
        dfc = 4*fc/N # Why 4???? I dunno yet and that's bad.

        h_vector = [] # vector of values of h(t) at sampling points
        g_vector = [] #
        for j in range(0, M):
            z = self.G(a+dt*j)
            h_vector.append( z.real )
            g_vector.append( z.imag )
        # padding the vector with zeroes. The function G(t) is assumed to be zero beyond t = K
        h_vector.extend( [0]*(over_sample-1)*len(h_vector) )
        g_vector.extend( [0]*(over_sample-1)*len(g_vector) )

        f_vector = [] # vector of frequencies
        for n in range(-N/2, N/2):
            f_vector.append( n*2*np.pi/(N*dt) )

        # fourier transformed vector h
        FT_h_vector = ft.fft(h_vector)
        FT_g_vector = ft.fft(g_vector)

        I = [] # Spectral power as vector
        I_h = []
        I_g = []
        for n in range(-N/2, N/2):
            theta_n = f_vector[n]*dt
            a0 = alpha0(theta_n)
            z_h = dt*( W(theta_n)*FT_h_vector[n]+a0*h_vector[0]+np.exp(1j*f_vector[n]*K)*a0.conjugate()*h_vector[M-1] )
            I_h.append(z_h)
            z_g = dt*( W(theta_n)*FT_g_vector[n]+a0*g_vector[0]+np.exp(1j*f_vector[n]*K)*a0.conjugate()*g_vector[M-1] )
            I_g.append(z_g)
            I.append(z_h.real-z_g.imag)
        wavenum = []
        for i in range(-N/2, N/2):
            wavenum.append(i*dfc)

        # normalizing the spectrum
        I_max = max(I)
        wn_max = wavenum[I.index(I_max)]
        for i in range(0, len(I)):
            I[i] = I[i]/I_max

        return [wavenum, I]


# experimental data for [25, 30, 35, 40, 45, 60, 70, 90, 150] Kelvins
fwhm1e10 = [0.05, 0.18, 0.28, 0.48, 0.95, 2.49, 3.71, 8.28, 27.19]

#experimental data for [25, 30, 35, 40, 45, 60, 70, 80, 90, 100, 120, 150,  180, 240] Kelvins
t1e11 = [25, 30, 35, 40, 45, 60, 70, 80, 90, 100, 120, 150, 180, 240]
fwhm1e11 = [0.04, 0.16, 0.21, 0.44, 0.79, 2.23, 3.76, 5.05, 7.78, 10.9, 17.2, 27.19, 45.2, 77.9]

f1e11_log = map(np.log, fwhm1e11)
f1e11_round = [round(x, 3) for x in fwhm1e11]
t1e11_log = map(np.log, t1e11)
t1e11_round = [round(x, 0) for x in t1e11]

omegas = [120]
deltas = [10]#[-5, -10, -15]
Gammas = [4.0]#, 2.0]
temperatures = [200]#[25, 30, 35, 40, 45, 60, 70, 80, 90, 100, 120, 150, 240, 300]#range(5, 305, 10) # (initial T, final T, step)
t_log = map(np.log, temperatures)
t_round = [round(x, 0) for x in temperatures]

oSpectrum = ZplSpectrum()
oSpectrum.gamma0 = 0.03


for i in omegas:
	oSpectrum.w1 = i
	for j in deltas:
		oSpectrum.Delta = j
		for k in Gammas:
			oSpectrum.Gamma = k*abs(oSpectrum.Delta)
			fwhms = [] # full widths
			for t in temperatures:
				oSpectrum.T = float(t)
				spectrum = oSpectrum.Generate()
				fwhms.append( lpl.fwhm(spectrum) - 2*oSpectrum.gamma0 ) # subtract 2*gamma_0
			# creating a figure with double log scale
			f_log = map(np.log, fwhms)
			f_round = [round(x, 3) for x in fwhms]
			figure = plt.figure()
			ax = figure.add_subplot(111)
			ax.set_title("log-log FWHM vs T; w1="+str(oSpectrum.w1)+" D="+str(oSpectrum.Delta)+" G="+str(oSpectrum.Gamma))
			ax.plot(t_log, f_log, t1e11_log, f1e11_log, 'bs') # plot with blue squares
			ax.set_xticks(t_log)
			ax.set_xticklabels( t_round )
			ax.set_yticks(f_log)
			ax.set_yticklabels( f_round )
			ax.grid(True)
			ax.set_ylabel('FWHM, cm-1')
			ax.set_xlabel('Temperature, K')
			pdfOutFile = PdfPages("w1_"+str(oSpectrum.w1)+" D_"+str(oSpectrum.Delta)+" G_"+str(oSpectrum.Gamma)+".pdf")
			plt.savefig(pdfOutFile, format = 'pdf')
			pdfOutFile.close()
