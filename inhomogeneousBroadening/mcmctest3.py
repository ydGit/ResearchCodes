import random
import math
import time
import matplotlib.pyplot as plt
import numpy
import lplspectra as lpl

# auxilary functions for handling data from files
def readFile(fileName):
    '''Read the data saved in the file `fileName`'''
    data = open(fileName, 'r')
    shift = []
    line = data.readline()
    while (len(line) > 1):    
        line = line.rstrip('\n')
        shift.append(float(line))
        line = data.readline()
    data.close()
    return shift

def binsCenter(bins):
    '''Create an x-axis from a histogram's bins.
    The mid-points of the bins are used.'''''
    x = []
    for i in range(1, len(bins)):
        x.append( (bins[i]+bins[i-1])/2 )
    return x

def lineParams(x, n):
    '''Given the spectrum: x - wavenumber (or wavelength), n - the intensity,
    calculate the full width at half maximum and the position of the maxumum'''
    FWHM = lpl.fwhm([x, n.tolist()])
    shift = x[n.tolist().index(max(n))]
    print "FWHM: ", FWHM
    print "Center: ", shift
#-------------------------------------------------------------------------------------

def f2(N0, x):
#    '''The PDF of the 2nd NN is proportional to f2(x): P2(x, N0) = K2*f2(N0, x)'''
    if (x >= 1.0):
        return x*(x**2-1)*math.exp(-N0*(x**2-1))
    else:
        return 0.0

def f3(N0, x):
#    '''The PDF of the 3rd NN is proportional to f3(x): P3(x, N0) = K3*f3(N0, x)'''
    if (x >= 1.0):
        return ((x**2-1)**2)*x*math.exp(-N0*(x**2-1))
    else:
        return 0.0

def xmax(N, N0):
    '''The starting value for the Markov chain.
    Should be close to the maximum of fN(N0, x) for 
    a better mixing of the chain.'''
    return math.sqrt( (N-1.0/2.0+N0)/(2.0*N0)*( 1+math.sqrt( 1-(2*N0)/(N-1.0/2.0+N0) ) ) )

def q(x):
    '''Symmetric random-walk candidate generating function.
    We can use a simple Metropolis algorithm.
    Remember that `sigma` is a tuning parameter used 
    to achieve a better mixing.'''   
    return random.gauss(0, 0.45*x) # gauss(mu, sigma)

#A = [[1.0, 1.0, 0.0], [1.0, 1.0, 0.0], [0.0, 0.0, 2.0]] # matrix of coupling coefficients # Matrix 1
#A = [[4.59, 10.67, 2.63], [10.67, 4.59, -2.63], [2.63, -2.63, -4.31]] # Matrix 2: Monoclinic I centre from Davies: Approximate Linewidth of ZPL 
#A = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]] # Matrix 3
#A = [[0.326, -0.938, -0.388], [-0.938, -0.214, 0.97], [-0.388, 0.97, -0.74]] # Matrix 4
A = [[1.0*10**(-3), 1.0, 1.0], [1.0, 1.0*10**(-3), 1.0], [1.0, 1.0, 1.0*10**(-3)]] # Matrix 5
TrA = A[0][0]+A[1][1]+A[2][2] # trace of A

chainLength = 1000000 #ensemble size
burnInLength = 1000

N0 = 0.000032 # dimensionless concentration -- number of excluded defects
S = 1000000.0 # the strenght of the defects
S1 = 1.0 # relative strength of the defects
S2 = 1.0 # used to switch on and off the second 
S3 = 1.0 # and the third defect 

xm2 = xmax(2.0, N0) # mode of the P2(x) -- the position of the maximum
xm3 = xmax(3.0, N0) # mode of the P2(x) -- the position of the maximum
x2 = xm2 # current position x = r/r0
x3 = xm3 # current position x = r/r0
chain2 = [x2] # Markov Chain of the dimensionless positions of the second NN
chain3 = [x3] # Markov Chain of the dimensionless positions of the third NN

shift = [] # the matrix of the positions of the shifted lines

#initializing the random number generator
seed = time.time() - math.floor(time.time())
random.seed(seed);

#burn-in sequence for the first chain (2nd NN)
while ( len(chain2) < burnInLength ):
    # generate a candidate from a sampling function q(x)
    xc = x2+q(xm2) # random walk
    alpha = min( [f2(N0, xc)/f2(N0, x2), 1.0] )
    if (alpha == 1.0): #accept the candidate
        x2 = xc
        chain2.append(x2)
    else:#accept the candidate with the probability alpha
        u = random.random()
        if (u <= alpha):
            x2 = xc
            chain2.append(x2)
#run the convergence test
# TODO:-- convergence test--
#burn-in sequence for the second chain (3rd NN)
while ( len(chain3) < burnInLength ):
    # generate a candidate from a sampling function q(x)
    xc = x3+q(xm3) # random walk
    alpha = min([f3(N0, xc)/f3(N0, x3),1.0])
    if (alpha == 1.0): #accept the candidate
        x3 = xc
        chain3.append(x3)
    else:#accept the candidate with the probability alpha
        u = random.random()
        if (u <= alpha):
            x3 = xc
            chain3.append(x3)
#run the convergence test
# TODO:-- convergence test--            
# clear the chain
c = 0

# Start the numerical experiment
# For debugging purposes we will keep the positions of 
# the defects. We can examine them using histograms
r1NN = []
r2NN = []
r3NN = []

while ( c < chainLength ):
    # generate a candidate from a sampling function q(x)
    xc = x2+q(xm2) # random walk
    alpha = min([f2(N0, xc)/f2(N0, x2), 1.0])
    accepted = False
    if (alpha == 1.0): #accept the candidate        
        accepted = True
    else: #accept the candidate with the probability alpha
        u = random.random()
        if (u <= alpha):           
            accepted = True
    if accepted:
        accepted = False
        x2 = xc
        r2 = x2
        r2NN.append(r2)
        # create the third and the first NNs
        xc = x3+q(xm3) # random walk
        alpha = min([f3(N0, xc)/f3(N0, x3),1.0])
        if (alpha == 1.0): #accept the candidate        
            accepted = True
        else: #accept the candidate with the probability alpha
            u = random.random()
            if (u <= alpha):           
                accepted = True
        if accepted:
            # create the first and the third NN 
            x3 = xc
            r3 = x3
            r3NN.append(r3)
            c = c + 1
            # generate the 1st nearest neighbor
            r1 = math.sqrt( 1-math.log( 1-random.random() )/N0 )
            r1NN.append(r1)
            # generate the angular positions
            phi1, phi2, phi3 = random.uniform(0, 6.28318531), random.uniform(0, 6.28318531), random.uniform(0, 6.28318531)
            # calculate the shifts
            shift1 = (1/r1**3)*( TrA - 3*( (A[0][0]+A[1][1])/2+(A[1][1]-A[2][2])/2*numpy.cos(2*phi1)+(A[0][1]+A[1][0])/2*numpy.sin(2*phi1)))
            shift2 = (1/r2**3)*( TrA - 3*( (A[0][0]+A[1][1])/2+(A[1][1]-A[2][2])/2*numpy.cos(2*phi2)+(A[0][1]+A[1][0])/2*numpy.sin(2*phi2)))
            shift3 = (1/r3**3)*( TrA - 3*( (A[0][0]+A[1][1])/2+(A[1][1]-A[2][2])/2*numpy.cos(2*phi3)+(A[0][1]+A[1][0])/2*numpy.sin(2*phi3)))
            shift.append( S*( S1*shift1+S2*shift2+S3*shift3 ) )
        

# rmax = int(xm*4)
# n, bins, patches = plt.hist(chain, bins = 200, range = (1, rmax), normed = 1)
# nmax = max(n)
# xx = range(1, rmax)
# yy = []
# fm = f(N0, xm)
# for i in xx:
#     yy.append(f(N0, i)*nmax/fm)
# n, bins, patches = plt.hist(chain, bins = 200, range = (1, rmax), normed = 1); plt.plot(xx, yy, 'red'); plt.show()
