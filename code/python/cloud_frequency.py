from netCDF4 import Dataset
import numpy as np
from scipy import signal
import matplotlib
import matplotlib.pyplot as plt

#
def fft(tmp,samp,freqall) : 
  # Calculate the Fourier transform for one time serie
  # Names allowed in freqall = intra,inter,decadal 
  y    = tmp
  L    = len(y) # Length of the signal
  namefig = 'FFT_decomp_'+samp
  
  if strcmp(samp,'mth') :
    Fs = 1/(30.*24.*3600.)       # Sampling frequency
      
  ripple = Fs/np.power(10,8)
  
  #[z,p,k] = cheby1(deg,ripple,WP,band); 
  
  fig = plt.figure('freq')
  j = 0
  if ff in freqall :
    # Extract information
    band, deg, pmin, pmax = define_fft(ff)
    if band is 'low' :
      WP = 2/float(pmin)
    elif band is 'high' :
      WP = 2/float(pmax)
    else :
      WP = [2/float(pmin),2/float(pmax)]
    # Create a Chebyshev type I filter design
    b, a = signal.cheby1(deg, ripple, WP, band)
    # Frequency and amplitude of the signal
    w, h = signal.freqs(b, a)
    # Forward-backward filter
    yd = signal.filtfilt(b, a, tmp)
    
    plt.subplot(2,round(len(freqall))/2,j)
    plt.plot(tmp,'k')
    plt.plot(yd,'r')
    plt.title(ff)
    j += 1

  fig.savefig('./'+namefig+'.png')
  fig.savefig('./'+namefig+'.pdf')
  
def define_fft(typ):
  periodmin=0
  periodmax=0
  if typ is 'intra' :
    band='high'
    deg =12
    periodmax=6
  elif typ is 'season' :
    band='bandpass'
    deg =8 # different here
    periodmin=11.8
    periodmax=12.2
  elif typ is 'inter' :
    band='low'
    deg =12
    periodmin=30
  elif typ is 'decadal' :
    band='low'
    deg =6 # different here
    periodmin=200
  else :
    print ' Problem with typ. Not recognized : ',typ
    stop
    
  return band, deg, periodmin, periodmax








            ######### Main Program ##########
            
        
# This Python file is based on several subroutine used in Brient and Schneider 16 (Journal of Climate)
# It decomposes two monthly time series in different frequency bands (intra-annual, seasonal, interannual and decadal)
# It provides bootstrapped time series of the original dataset
# The bootstrapping is based on index of one time series that allow keeping bootstrap homogeneity between the two time series
# It provides also correlation coefficient, regression slopes (robust?) for the original temporal series and for the bootstrapped series




##### User defined ######
# Open timeseries
try :
  ev1,ev2=open(file.dat)
except :
  # Try random
  t=np.arange(100)
  ev1=np.sin(2 * np.pi * t/4)/2 + np.sin(2 * np.pi * t/12) + np.sin(2 * np.pi * t/30)/2 + np.random.random((len(t)))/10

# Data are monthly (routine not ready otherwise)
samp = 'mth';
#########################


NB = len(ev1)
data = np.zeros((2,NB))
data[0,:]=ev1 ; data[1,:]=ev2

# Anomalies only
datamean = np.mean(b,axis=1)
for ij in np.arange(NB) :
  data[:,ij] =  data[:,ij] - datamean[:]
  
  
# Seasonal cycle
dataseas = np.zeros((2,12))
for ij in np.arange(12) :
  dataseas[:,ij] = np.mean(data[:,np.arange(ij,NB,12)],axis=1)
  
# Deseasonalized time series
ik = 0
datanoseas = np.zeros((data.shape))
for ij in np.arange(NB) :
  datanoseas[:,ij] = data[:,ij] - dataseas[:,ik]
  ik += 1 
  if ik == 12 : ik=0

# Calculate FFT
freqall = 'inter'
for ij in np.arange(2) :
  tmp     = datanoseas[ij,:]
  tmp_FFT = fft(tmp,samp,freqall)














