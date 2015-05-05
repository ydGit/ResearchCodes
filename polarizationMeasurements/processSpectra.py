import sys
import os
import matplotlib.pyplot as plt
import numpy as np

# useful functions
def ReadSpectra( filename ):
    """Read real spectrum(a) from a TXT file specified by filename. The type of the 
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


def Smooth(x, window_len=10, window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    import numpy as np    
    t = np.linspace(-2,2,0.1)
    x = np.sin(t)+np.random.randn(len(t))*0.1
    y = smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string   
    """

    # if x.ndim != 1:
    #     raise ValueError, "smooth only accepts 1 dimension arrays."

    # if x.size < window_len:
    #     raise ValueError, "Input vector needs to be bigger than window size."

    # if window_len < 3:
    #     return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"

    s=np.r_[2*x[0]-x[window_len:1:-1], x, 2*x[-1]-x[-1:-window_len:-1]]
    #print(len(s))
    
    if window == 'flat': #moving average
        w = np.ones(window_len,'d')
    else:
        w = getattr(np, window)(window_len)
    y = np.convolve(w/w.sum(), s, mode='same')
    return y[window_len-1:-window_len+1]


def PolarizationCurve(Iparall, Iperp, k):
    polarizn = []
    for i in range( 0, len(Iparall) ):
        polarizn.append( (k*Iparall[i]-Iperp[i]) / (k*Iparall[i]+Iperp[i]) )
    return polarizn

 

# # ================ DATA ANALYSIS STARTS HERE =====================

# # getting the list of all "txt" files
workingDir = sys.argv[1]

dataFiles = []

for root, dirs, files in os.walk(workingDir):
    for i in files:
        if ( i[-4:] == ".txt" ):
            dataFiles.append(i)
            # # sorting the data file list according to their names
            dataFiles.sort()
# # extracting time stamps
timeStamps = []
for i in dataFiles:
    timeStamps.append(i[-8:-4])
    
# # sorting time stamps with t increasing
timeStamps.sort()

s = 1 # switch variable. 
# # Measurements go like this: perp, parall, parall, perp, perp, parall, parall, ...
# # so we need to alternate the order of files

# # unadjusted peak amplitudes.
I794ParallRaw = []
I794PerpRaw = []
I811ParallRaw = []
I811PerpRaw = []

# # the index of the file, appended by the time stamp
# # and prepended by the baseFileName
baseFileName = ""
idxRange = range(1,  98, 2)

# positions of the peaks and points outside of the peaks
# the range for the averaging will be
#  [ x - 12, x + 12 ] with the size = 24
x794 = 73
x1 = 281
x811 = 549
x2 = 1012
# average values of the spectrum around key points
y794 = 0
y1 = 0
y811 = 0
y2 = 0

for i in idxRange:

    fileName1 = baseFileName + str(i) + timeStamps[i-1] + ".txt"
    fileName2 = baseFileName + str(i+1) + timeStamps[i] + ".txt";

    print "Working on file " + fileName1  + " and  " + fileName2
    # wavelength list
    z = ( ReadSpectra( workingDir+fileName1 ) )[1][0]
    if ( s == 1 ):        
        Iperp = ( ReadSpectra( workingDir+fileName1 ) )[1][1]
        Iparall = ( ReadSpectra( workingDir+fileName2 ) )[1][1]
    if ( s == -1):
        Iparall = ( ReadSpectra( workingDir+fileName1 ) )[1][1]
        Iperp = ( ReadSpectra( workingDir+fileName2 ) )[1][1]
    
    # average values at the key-points, used to find the background fit
    avgParall794 = np.average(Iparall[x794-12:x794+12])
    avgPerp794 = np.average(Iperp[x794-12:x794+12])
    avgParall811 = np.average(Iparall[x811-12:x811+12])
    avgPerp811 = np.average(Iperp[x811-12:x811+12])

    avgParall1 = np.average(Iparall[x1-12:x1+12])
    avgPerp1 = np.average(Iperp[x1-12:x1+12])
    avgParall2 = np.average(Iparall[x2-12:x2+12])
    avgPerp2 = np.average(Iperp[x2-12:x2+12])

    # approximating background: for 794 nm -- straight horizontal line
    # for 811 -- straight line with slope
    I794ParallRaw.append( avgParall794 - avgParall1 )
    I794PerpRaw.append( avgPerp794 - avgPerp1 )
    I811ParallRaw.append( avgParall811-np.average([avgParall1, avgParall2]) )
    I811PerpRaw.append( avgPerp811-np.average([avgPerp1, avgPerp2]) )    
    s = -s;

n = len(I794ParallRaw) # number of data points
# measurements were made every 5 degrees, angles from -30 up to 210
angles = range(-30, 215, 5)

# Saving the results to TXT files
print "Saving Peak Intensities into File"
i794PerpFile = open(workingDir + "peakIntensities/" + "794nmPerp", 'w')
i794ParallFile = open(workingDir + "peakIntensities/" + "794nmParall", 'w')
i811PerpFile = open(workingDir + "peakIntensities/" + "811nmPerp", 'w')
i811ParallFile = open(workingDir + "peakIntensities/" + "811nmParall", 'w')
for i in range(0, n):
    i794PerpFile.write(str(angles[i])+'\t'+str(I794PerpRaw[i])+'\n')
    i794ParallFile.write(str(angles[i])+'\t'+str(I794ParallRaw[i])+'\n')
    i811PerpFile.write(str(angles[i])+'\t'+str(I811PerpRaw[i])+'\n')
    i811ParallFile.write(str(angles[i])+'\t'+str(I811ParallRaw[i])+'\n')

i794PerpFile.close()
i794ParallFile.close()
i811PerpFile.close()
i811ParallFile.close()
print "Saving Done: OK"


# # Caclulating polarization curves for both peaks

pol794SimpleFile = open(workingDir + "polarization/" + "794Simple", 'w')
pol811SimpleFile = open(workingDir + "polarization/" + "811Simple", 'w')

# setting the correction factor such that it give 0 polarization 
# at 0 angle for 811 nm peak
k = I811PerpRaw[6]/I811ParallRaw[6]

# calculating and saving
# simple case: not correction
polarization794 = PolarizationCurve(I794ParallRaw, I794PerpRaw, k)
polarization811 = PolarizationCurve(I811ParallRaw, I811PerpRaw, k)
for i in range(0, n):
    pol794SimpleFile.write(str(angles[i])+'\t'+ str(polarization794[i])+'\n')
    pol811SimpleFile.write(str(angles[i])+'\t'+ str(polarization811[i])+'\n')
pol794SimpleFile.close()
pol811SimpleFile.close()
