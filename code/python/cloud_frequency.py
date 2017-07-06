from netCDF4 import Dataset
import numpy as np
from scipy import signal
from scipy import stats
import matplotlib
import matplotlib.pyplot as plt
import statsmodels.api as sm
import stationary_bootstrap as myboot
import slopeinterval as myint

# Import files
def importfiles(typ,path0) :
  if typ == 'random':
    try :
      # change it by your input path
      print 'Read file.dat'
      file0  = path0+'file.dat'
      ev1,ev2 = np.loadtxt(file0, unpack=True)
    except :
      # Try random
      t=np.arange(100)
      # Three periods (4 months, 1 year, 3 years)
      ev0=np.sin(2 * np.pi * t/4)/2 + np.sin(2 * np.pi * t/12) + np.sin(2 * np.pi * t/36)/2
      # Add random noise
      ev1=ev0 + np.random.random((len(t)))/2
      ev2=ev1 + np.random.random((len(t)))/2
      # Write in a file.dat
      file0  = path0+'file.dat'
      f = open(file0, 'wb')
      for ij in np.arange(len(ev1)):
        f.write("%.2f %.2f\n" % (ev1[ij], ev2[ij]))
      f.close()
  elif typ == 'observations':
    # Observation data from Brient and Schneider 16 (Figure 2)
    # Cloud albedo from CERES2
    # SST from ERSST
    # Averaged over 25% driest of the tropical area
    ev1=[] ; ev2=[]
    file1=open(path0+'sst_ersst.txt','r')
    file1.readline()
    for dd in file1.readlines():
        data = float(dd.replace('\n',''))
        ev1.append(data)
    file2=open(path0+'albcld_ceres.txt','r')
    file2.readline()
    for dd in file2.readlines():
        data = float(dd.replace('\n',''))
        ev2.append(data)
         
  return ev1, ev2

# Subroutine for the fft decomposition
def fft(tmp,tmpch,samp,freqall,pathfig) : 
  # Calculate the Fourier transform for one time serie
  # Names allowed in freqall = intra,inter,decadal 
  y      = tmp
  L      = len(y) # Length of the signal
  NF     = len(freqall)
  tmpFFT = np.zeros((NF,L))
  
  if samp is 'mth' :
    Fs = 1/(30.*24.*3600.)       # Sampling frequency
      
  # Not sure about the real value to use
  ripple = 0.01#.00001 #Fs#(Fs/np.power(10,8))*10^12
  
  fig = plt.figure('freq')
  #print freqall
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
  fig.savefig(pathfig+namefig+'_'+tmpch+'.png')
  fig.savefig(pathfig+namefig+'_'+tmpch+'.pdf')
  plt.close()
  
  del ripple, WP, band, deg, pmin, pmax
  return tmpFFT
  
# FFT information(cheby1 filter)
def define_fft(typ):
  print 'Here the type is  ',typ
  periodmin=0
  periodmax=0
  if typ is 'intra' :
    band='high'
    deg =12
    periodmax=10
  elif typ is 'season' :
    band='bandpass'
    deg =2 # different here
    periodmin=11
    periodmax=13
  elif typ is 'inter' :
    band='low'
    deg =12
    periodmin=14
  elif typ is 'decadal' :
    band='low'
    deg =6 # different here
    periodmin=200
  else :
    print ' Problem with typ. Not recognized : ',typ
    stop
    
  return band, deg, periodmin, periodmax
  
  
# Slope/correlation information
def slope_create(F0,F1):
  # Calculate slopes (OLS)
  slope, intercept, r_value, p_value, std_err = stats.linregress(F0,F1)
  
  # Slope with robust regression
  x = sm.add_constant(F0)
  y = F1
  rlm_results = sm.RLM(y,x, M=sm.robust.norms.HuberT()).fit()
  slope_r     = rlm_results.params[-1]
  intercept_r = rlm_results.params[0]
  
  del x,y,rlm_results
  return slope, intercept, r_value, slope_r, intercept_r

def bootstraprun(F0,F1):
  # Number of bootstrapping
  B = 1
  # Index for bootstrapping 
  data = F0
  # Optimal length (Need to work on that)
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
# Open timeseries evx and evy (identical time length)
# See importfiles to modify your input
typ     = 'observations'

path0   = "../../data/"+typ+"/"
pathfig = "../../figures/"+typ+"/"
evx,evy = importfiles(typ,path0)

# Data are monthly (routine not ready otherwise)
samp = 'mth'

NB   = len(evx)
data = np.zeros((2,NB))
data[0,:]=evx ; data[1,:]=evy

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
# Time already defined ('intra','season','inter')
freqall = ['intra','season','inter']
NF      = len(freqall)
tmp     = datanoseas
evxfft  = np.zeros((NF+1,tmp.shape[-1]))
evyfft  = np.zeros((NF+1,tmp.shape[-1]))
evxfft[0,:]  = tmp[0,:]
evyfft[0,:]  = tmp[1,:]


name_serie = ['first','second']
for ij in np.arange(len(name_serie)) :
  tmpFFT  = fft(tmp[ij,:],name_serie[ij],samp,freqall,pathfig)
  if ij == 0 :
    evxfft[1:,:] = tmpFFT
  else :
    evyfft[1:,:] = tmpFFT

    
# Save slopes, correlation coefficient for the unfiltered data
#slope0, int0, r0, slope_r0, int_r0 = slope_create(tmp[0,:],tmp[1,:])

# Save slopes, correlation coefficient for the unfilterd and filtered data series
# Original
slope   = np.zeros((NF+1))
int     = np.zeros((NF+1))
r       = np.zeros((NF+1))
slope_r = np.zeros((NF+1))
int_r   = np.zeros((NF+1))
# Bootstrap
Nb       = 200
slopeb   = np.zeros((NF+1,Nb))
intb     = np.zeros((NF+1,Nb))
rb       = np.zeros((NF+1,Nb))
slope_rb = np.zeros((NF+1,Nb))
int_rb   = np.zeros((NF+1,Nb))

for ij in np.arange(NF+1) :
  F0 = evxfft[ij,:]; F1 = evyfft[ij,:]
  slope[ij], int[ij], r[ij], slope_r[ij], int_r[ij] = slope_create(F0,F1)
  # Bootstrapping (stationary)
  # Number of bootstrap Nb
  # Caution : the Full time serie
  for ib in np.arange(Nb) :
    FF0,FF1 = bootstraprun(F0,F1)
    slopeb[ij,ib], intb[ij,ib], rb[ij,ib], slope_rb[ij,ib], int_rb[ij,ib] = slope_create(FF0[:,0],FF1[:,0])
    del FF0,FF1
  del F0,F1
  
#########################
# Creating output files #
#########################

freqall = ['deseason'] + freqall

# Original
f = open(path0+'output_original.dat', 'wb')
f.write('%13s %13s %13s %13s %13s %13s\n' % ('Frequency','Corr coef','OLS regress',' robust regress',' OLS intercept',' robust intercept'));
for ij in np.arange(NF+1):
  f.write("%13s %13.4f %13.4f %13.4f %13.4f %13.4f\n" % (freqall[ij], r[ij], slope[ij], slope_r[ij], int[ij], int_r[ij]))
f.close()  

# Bootstrapp
fileout0=path0+'output_boot_xx.dat'
for ij in np.arange(NF+1):
  f = open(fileout0.replace('xx',freqall[ij]), 'wb')
  f.write('%13s %13s %13s %13s %13s %13s\n' % ('Frequency','Corr coef','OLS regress',' robust regress',' OLS intercept',' robust intercept'));
  for ib in np.arange(Nb):
    f.write("%13s %13.4f %13.4f %13.4f %13.4f %13.4f\n" % (freqall[ij], rb[ij,ib], slopeb[ij,ib], slope_rb[ij,ib], intb[ij,ib], int_rb[ij,ib]))
  f.close()
    
  
#########################
# Plotting some results #
#########################

# Scatterplot
fig = plt.figure('scatterplot')
for ff in np.arange(NF+1) :
  plt.subplot(2,round(float(NF+1)/2),ff+1)
  title = freqall[ff]
  t0 = evxfft[ff,:];t1=evyfft[ff,:]
    

  a  = slope[ff]; ar  = slope_r[ff]
  b  = int[ff]; br  = int_r[ff]
  cc = r[ff]
  stda = np.std(slope_rb[ff,:])
  stdb = np.std(int_rb[ff,:])
 
  str0  = str('%01.2f' % (cc,))  
  title = title + ' (r=' +str0+')'
    
  plt.plot(t0,t1,'b.')
  plt.plot(t0,a*t0  + b,'b')
  plt.plot(t0,ar*t0 + br,'r')
  
  pmin,pmax = myint.intver(t0,a,b,stda,stdb) 
  plt.plot(t0,pmin ,'b--')
  plt.plot(t0,pmax ,'b--')
  print stda,stdb
  
  plt.axhline(0,linewidth=0.5,color='black')
  plt.axvline(0,linewidth=0.5,color='black')
  plt.title(title)

namefig = 'Scatter_all'
fig.savefig(pathfig+namefig+'.png')
fig.savefig(pathfig+namefig+'.pdf')
plt.close()


# Bar plot of uncertainties (based on the bootstrap analysis)
nameall = ['slope','correlation']
width   = 0.25
for ij in np.arange(len(nameall)):
  fig = plt.figure(nameall[ij])
  for ff in np.arange(NF+1) :
    if ij==0:
      t0 = slope[ff]; t0r= slope_r[ff]
      std  = np.std(slopeb[ff,:])
      stdr = np.std(slope_rb[ff,:])
    if ij==1: 
      t0 = r[ff]
      std  = np.std(rb[ff,:])
      
    plt.bar(ff-width/2,t0,width=width,color='b',yerr=std)
    if ij==0:
      plt.bar(ff+width/2,t0r,width=width,color='r',yerr=stdr)
    
  plt.axhline(0,linewidth=0.5,color='black')
  plt.xticks(np.arange(NF+1)+width/2,freqall)
  plt.title(nameall[ij])
  namefig = 'Bar_'+nameall[ij]
  fig.savefig(pathfig+namefig+'.png')
  fig.savefig(pathfig+namefig+'.pdf')
  plt.close()
  








