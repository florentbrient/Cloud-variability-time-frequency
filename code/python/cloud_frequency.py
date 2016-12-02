from netCDF4 import Dataset
import numpy as np
from scipy import signal

#
def fft(tmp,freqall) : 
  # Calculate the Fourier transform for one time serie
  # Names allowed in freqall = intra,inter,decadal 
  y    = data
  L    = length(y) # Length of the signal
  
  deg = 
  ripple =
  b, a = signal.cheby1(deg, ripple, WP, band, analog=True)
  #[z,p,k] = cheby1(deg,ripple,WP,band); 
  
  if ff in freqall :
    band, deg, pmin, pmax = define_fft(ff)
     

  
def define_fft(typ):
  periodmin=0
  periodmax=0
  if typ is 'intra' :
    band='high'
    deg =12
    periodmax=6; 
  if typ is 'season' :
    band='bandpass'
    deg =8
    periodmin=11.8;
    periodmax=12.2;
  if typ is 'inter' :
    band='low'
    deg =12
    periodmin=30;
  if typ is 'decadal' :
    band='low'
    deg =6 # different here
    periodmin=200;
    
  return band, deg, periodmin, periodmax








            ######### Main Program ##########
            
        
# This Python file is based on several subroutine used in Brient and Schneider 16 (Journal of Climate)
# It decomposes two monthly time series in different frequency bands (intra-annual, seasonal, interannual and decadal)
# It provides bootstrapped time series of the original dataset
# The bootstrapping is based on index of one time series that allow keeping bootstrap homogeneity between the two time series
# It provides also correlation coefficient, regression slopes (robust?) for the original temporal series and for the bootstrapped series




# Open timeseries
try :
  ev1,ev2=open(file.dat)
except :
  # Try random
  ev1=np.random.random(100)
  ev2=np.random.random(100)
  
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
datanoseas = data - dataseas

# Calculate FFT 
for ij in np.arange(2) :
  tmp     = datanoseas[ij,:]
  tmp_FFT = fft(tmp)














