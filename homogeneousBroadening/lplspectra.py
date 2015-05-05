#/usr/bin/env python

# spectra manipulation routines
from numpy import *
from scipy.optimize import leastsq
import math
import numpy as np

def convert_spectrum(real_spectrum, width, height):
    """Convert real physical spectrum into a list of pixels"""
    plot_list = [[],[]]
    # real spectrum is a list of two lists of the following type:
    # [ [x1, x2, ..., xN], [y1, y2, ..., yN] ]
    Xmin = min(real_spectrum[0]); Xmax = max(real_spectrum[0])
    Ymin = min(real_spectrum[1]); Ymax = max(real_spectrum[1])
    xmin = 1; xmax = width
    ymin = 1; ymax = height
    # coefficient of linear transformation
    a = (xmax - xmin)/(Xmax - Xmin)
    b = xmin - a*Xmin
    c = (ymax - ymin)/(Ymax - Ymin)
    d = ymin - c*Ymin
    for i in range(0,len(real_spectrum[0])):
        x_i = int(a*real_spectrum[0][i]+b)
        y_i = int(2*ymin+ymax-c*real_spectrum[1][i]-d)
        plot_list[0].append(x_i); plot_list[1].append(y_i)
    return plot_list

def get_spectra_from_file( filename ):
    """Read real spectrum(a) from a file specified by filename. The type of the 
    data is determined on the go. Possible types are:
    0D -- single point spectrum
    1D -- linear map
    2D -- map of the area on the surface of the sample"""

    data_file = open(filename, 'r')
    data_string = data_file.read()
    data_file.close()
    data_string = data_string.splitlines()
    # how many empty columns are in the first row of the data array
    # this tells you how many coordinates are used therefore the dimensionality
    dimensionality = data_string[0].split('\t').count('')

    if ( dimensionality == 0 ):
        # reading a simple spectrum at a point
        real_spectrum = [[],[]]
        for i in data_string:
            x_y = i.split('\t')
            real_spectrum[0].append(float(x_y[0]))
            real_spectrum[1].append(float(x_y[1]))
        return [dimensionality, real_spectrum]

    elif (dimensionality == 2):
        # reading two dimensional area map
        number_of_points = len(data_string) - 1 # subtract the wavelengths' line
        wavelength = []
        for i in data_string[0].split('\t'):
            if ( len(i) != 0 ):
                wavelength.append(float(i))    
        number_of_super_pixels = len(wavelength)

        spectra = []
        spectrum = []
        coords = [[],[]]
        for i in range(1,number_of_points+1):
            spectrum = data_string[i].split('\t') # here we still have a string
            spectrum = map(float, spectrum) # cast into a float
            coords[0].append(spectrum[1])# x coordinate of the point
            coords[1].append(spectrum[0])# y coordinate
            spectra.append(spectrum[2:number_of_super_pixels+2])
        coords[0] = set(coords[0]) # getting rid of the duplicates
        coords[1] = set(coords[1]) # the result is also sorted in ascending order
     
        return [dimensionality, wavelength, spectra, coords]
    else:
        return [-1]


def fit_linear_gauss(data):
    """Fit a gaussian+linear function to the data given. Parameters are as follows:
    p0 -- vertical shift of the linear function
    p1 -- slope of the linear function
    p2 -- amplitude of the peak (intensity)
    p3 -- pisition of the peak
    p4 -- width of the peak"""
    x, y = data[0], data[1]
    sampling = 10
    x_left = sum(x[0:sampling])/sampling
    y_left = sum(y[0:sampling])/sampling
    x_rigth = sum(x[len(x)-1-sampling:len(x)-1])/sampling
    y_right = sum(y[len(y)-1-sampling:len(x)-1])/sampling
    # estimation of the parameters' value
    slope = (y_left-y_right)/(x_left-x_rigth)
    shift = y_left-slope*x_left
    amplitude = max(y)
    position = x[y.index(amplitude)]
#    amplitude = amplitude - shift
    i = 0
    while (y[i] < amplitude/exp(1) ):
        i = i + 1
    width = position - x[i]
    params = leastsq(lambda p: (p[0]+array(x)*p[1]+p[2]*exp(-(x-p[3])**2/p[4]**2)-y), 
                     [shift, slope, amplitude, position, width])[0]
    return params

def extract_Xe_intensity( wavelength, spectra ):
    params_map = []
    for i in spectra:        
        params = fit_gauss( [wavelength[130:210], i[130:210]])
        params_map.append(params[2])
    return params_map

def read_saved_spectra( filename ):
    file = open(filename)
    data = file.read()
    file.close()
    data = data.strip('[]')
    data = data.split('], [')
    spectra = []
    for i in range(0, len(data)):
        spectra.append([])
        str_spectrum = data[i].split(',')
        for j in range(0, len(str_spectrum)):
            spectra[i].append(float(str_spectrum[j]))
    return spectra

def average_spectra( spectra_to_average ):
    number_of_spectra = len(spectra_to_average)
    spectrum_length = len(spectra_to_average[0])
    avg_spectrum = [0]*spectrum_length

    for i in range(0, spectrum_length):
        for j in range(0, number_of_spectra):
            avg_spectrum[i] = avg_spectrum[i]+spectra_to_average[j][i]
        avg_spectrum[i] = avg_spectrum[i]/number_of_spectra
    return avg_spectrum

def wavelength_to_wavenumber( laser_wavelength_nm , wavelength):
    wavenumber = []
    laser_wavenumber = 10**7/float(laser_wavelength_nm)
    for i in wavelength:
        wavenumber.append(laser_wavenumber - 10**7/i)
    return wavenumber

def wavenumber_to_wavelength( laser_wavelength_nm , wavenumber):
    wavelength = []
    laser_wavenumber = 10**7/float(laser_wavelength_nm)
    for i in wavenumber:
        wavelength.append(10**7/(laser_wavenumber - i))
    return wavelength

def normalize_spectrum( spectrum ):
    normalized_spectrum = []
    norm = 0
    for i in spectrum:
        norm = norm + i*i
    norm = math.sqrt(norm)
    for i in spectrum:
        normalized_spectrum.append(i/norm)        
    return normalized_spectrum

def cross_correlation0 (spectrum_1, spectrum_2):
    spectrum_1_normalized = normalize_spectrum(spectrum_1)
    spectrum_2_normalized = normalize_spectrum(spectrum_2)
    Xcorr = 0
    for i in range(0, len(spectrum_1)):
        Xcorr = Xcorr+spectrum_1_normalized[i]*spectrum_2_normalized[i]
    return Xcorr

def cross_corr_statistics( spectra ):
    Xcorr_coefficients = []
    number_of_spectra = len(spectra)
    for i in range(0, number_of_spectra):
        for j in range(i+1, number_of_spectra):
            Xcorr_coefficients.append(
                cross_correlation0(spectra[i], spectra[j])
                )
    return [Xcorr_coefficients, average(Xcorr_coefficients), std(Xcorr_coefficients)]

# Functions used to fit the data

def linear_gauss(shift, slope, amplitude, center, width, x):
    y = []
    for i in x:
        y.append(shift+slope*i+amplitude*np.exp(-(i-center)**2/(width/2)**2) )
    return y

def lorentz( amplitude, center, width, x ):
    y = []
    for i in x:
        y.append(amplitude/( (2*(i - center)/width)**2 + 1 ) )
    return y

def chop_list( list_ , size):
    output = []
    if ( (len(list_)%size) != 0 ):
        print "Error: the length of the list is not a multiple of the size\n"
        return []
    else:
        for i in range(0, len(list_)/size):
            output.append(list_[size*i:size*(i+1)])
    return output

def voigt_lg (fL, fG): # voigt profile FWHM from a convolution of Gaussian and Lorentzian
    return 0.5346*fL+np.sqrt(0.2166*fL**2+fG**2)

def lorentz_vg (fV, fG): # lorentzian FWHM from Gaussian and Voigt. "Deconvolution"
    return 7.725*fV-np.sqrt(45.23*fV**2+14.45*fG**2)

def fwhm(spectrum):
    '''Returns the Full Width at Half Maximum. Spectrum is a [wavenumber, intensity]'''
    wn = spectrum[0]
    I = spectrum[1]
    I_max = max(I)
    wn_left = wn_right = I.index(I_max)
    while ( I[wn_right] >= I_max/2):
        wn_right = wn_right + 1
    while ( I[wn_left] >= I_max/2):
        wn_left = wn_left - 1
    return wn[wn_right] - wn[wn_left]
