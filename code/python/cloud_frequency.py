from netCDF4 import Dataset
import numpy as np
from scipy import signal
from scipy import stats
import matplotlib
import matplotlib.pyplot as plt
import statsmodels.api as sm
import stationary_bootstrap as myboot

#
def fft(tmp,tmpch,samp,freqall) : 
  # Calculate the Fourier transform for one time serie
  # Names allowed in freqall = intra,inter,decadal 
  y      = tmp
  L      = len(y) # Length of the signal
  NF     = len(freqall)
  tmpFFT = np.zeros((NF,L))
  
  if samp is 'mth' :
    Fs = 1/(30.*24.*3600.)       # Sampling frequency
      
  # Not sure about the real value to use
  ripple = Fs/np.power(10,8)
  
  fig = plt.figure('freq')
  print freqall
  j = 1
  for ff in freqall :
    # Extract information
    band, deg, pmin, pmax = define_fft(ff)
    if band is 'low' :
      WP = 2/float(pmin)
    elif band is 'high' :
      WP = 2/float(pmax)
    else :
      WP = [2/float(pmax),2/float(pmin)]
    # Create a Chebyshev type I filter design
    b, a = signal.cheby1(deg, ripple, WP, band)
    # Frequency and amplitude of the signal
    w, h = signal.freqs(b, a)
    # Forward-backward filter
    yd = signal.filtfilt(b, a, tmp)
    
    if NF >= 2 :
      plt.subplot(2,round(float(NF)/2),j)
    plt.plot(tmp,'k')
    plt.plot(yd,'r')
    #plt.semilogx(w, np.log10(abs(h)))
    plt.title(ff)
    # Saving the filtered time series
    tmpFFT[j-1,:]=yd
    del yd
    
    j += 1

  namefig = 'FFT_decomp_'+samp
  fig.savefig('./'+namefig+'_'+tmpch+'.png')
  fig.savefig('./'+namefig+'_'+tmpch+'.pdf')
  plt.close()
  
  del ripple, WP, band, deg, pmin, pmax
  return tmpFFT
  
def define_fft(typ):
  print 'Here the type is  ',typ
  periodmin=0
  periodmax=0
  if typ is 'intra' :
    band='high'
    deg =12
    periodmax=4
  elif typ is 'season' :
    band='bandpass'
    deg =8 # different here
    periodmin=11
    periodmax=13
  elif typ is 'inter' :
    band='low'
    deg =12
    periodmin=36
  elif typ is 'decadal' :
    band='low'
    deg =6 # different here
    periodmin=200
  else :
    print ' Problem with typ. Not recognized : ',typ
    stop
    
  return band, deg, periodmin, periodmax
  
  
def slope_create(F0,F1):
  # Calculate slopes (OLS)
  slope, intercept, r_value, p_value, std_err = stats.linregress(F0,F1)
  # Slope with robust regression
  x = sm.add_constant(F0)
  y = F1
  rlm_results = sm.RLM(y,x, M=sm.robust.norms.HuberT()).fit()
  slope_r     = rlm_results.params[-1]
  intercept_r = rlm_results.params[0]
  print slope,slope_r
  print intercept,intercept_r
  
  del x,y,rlm_results
  # Save
  return slope, intercept, r_value, slope_r, intercept_r

def bootstraprun(F0,F1):
  # Number of bootstrapping
  B = 1
  # Index for bootstrapping 
  data = F0
  w = 3 # ?
  # call stationary bootstrap (external routine)
  bsdata, indices = myboot.stat_boot(data,B,w)
  
  del bsdata
  return F0[indices],F1[indices]
  
  
  





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
  # Three periods (4 months, 1 year, 3 years)
  ev0=np.sin(2 * np.pi * t/4)/2 + np.sin(2 * np.pi * t/12) + np.sin(2 * np.pi * t/36)/2
  # Add random noise
  ev1=ev0 + np.random.random((len(t)))/2
  ev2=ev0 + np.random.random((len(t)))/2
  
# Data are monthly (routine not ready otherwise)
samp = 'mth'
##### User defined ######


NB = len(ev1)
data = np.zeros((2,NB))
data[0,:]=ev1 ; data[1,:]=ev2

# Anomalies only
datamean = np.mean(data,axis=1)
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
freqall = ['intra','season','inter']
NF      = len(freqall)
#tmp     = data
tmp     = datanoseas
name_serie = ['first','second']
for ij in np.arange(2) :
  tmpch   = name_serie[ij]
  tmpFFT  = fft(tmp[ij,:],tmpch,samp,freqall)
  if ij == 0 :
    ev0fft = tmpFFT
  else :
    ev1fft = tmpFFT
    
# Save slopes, correlation coefficient for the unfiltered data
slope0, int0, r0, slope_r0, int_r0 = slope_create(tmp[0,:],tmp[1,:])

# Save slopes, correlation coefficient for the filtered data
slope   = np.zeros((NF))
int     = np.zeros((NF))
r       = np.zeros((NF))
slope_r = np.zeros((NF))
int_r   = np.zeros((NF))
for ij in np.arange(NF) :
  F0 = ev0fft[ij,:]; F1 = ev1fft[ij,:]
  slope[ij], int[ij], r[ij], slope_r[ij], int_r[ij] = slope_create(F0,F1)
  
  
# Bootstrapping (stationary)
# Number of bootstrap Nb
Nb = 2
slopeb   = np.zeros((NF,Nb))
intb     = np.zeros((NF,Nb))
rb       = np.zeros((NF,Nb))
slope_rb = np.zeros((NF,Nb))
int_rb   = np.zeros((NF,Nb))
for ij in np.arange(NF) :
  F0 = ev0fft[ij,:]; F1 = ev1fft[ij,:]
  for ib in np.arange(Nb) :
    FF0,FF1 = bootstraprun(F0,F1)
    slopeb[ij,ib], intb[ij,ib], rb[ij,ib], slope_rb[ij,ib], int_rb[ij,ib] = slope_create(FF0[:,0],FF1[:,0])
    
    
# Plotting some results
fig = plt.figure('scatterplot')
for ff in np.arange(NF+1) :
  plt.subplot(2,round(float(NF+1)/2),ff+1)
  if ff == 0 :
    t0 = tmp[0,:];t1=tmp[1,:]
    a  = slope0; ar  = slope_r0; 
    b  = int0; br  = int_r0; 
    title = 'original'
  else :
    f = ff-1
    t0 = ev0fft[f,:];t1=ev1fft[f,:]
    a  = slope[f]; ar  = slope_r[f]; 
    b  = int[f]; br  = int_r[f]; 
    title = freqall[f]
    
  plt.plot(t0,t1,'b.')
  plt.plot(t0,a*t0  + b,'b')
  plt.plot(t0,ar*t0 + br,'r')
  plt.axhline(0,linewidth=0.5,color='black')
  plt.axvline(0,linewidth=0.5,color='black')
  plt.title(title)

namefig = 'Scatter_all'
fig.savefig('./'+namefig+'.png')
fig.savefig('./'+namefig+'.pdf')
plt.close()











