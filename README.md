# Cloud-variability-time-frequency
This routine aims to calculate slopes and correlations between two time series (monthly data by default).
The two time series need to share the same temporal indexes.

Even if the code is written for monthly data, it may be adjusted to take into account lower temporal variability.

### What the code cloud_frequency.py do :
- Separate the annual cycle and the deseasonalized variability.
- Filter anomalies of the time series for different frequency bands. 
These bands are defined as intra-annual (1-yr low-pass filter), inter-annual (1-yr high-pass filter) and decadal (10-yr low-pass filter). 
A fourth band called "season" use a bandpass filter to extract the 12+/-0.2 months time period. A 12th-order Chebyshev filter is used by default.
- Calculate ordinary linear regression, robust regression and correlation coefficient for the relationship between the two original time series.
- Resample time series for each frequency bands through a non-parametric bootstrap procedure which takes the autocorrelations of the
time series into account [Need to be done, work in Matlab]. 
The original pairs of time series were resampled by drawing blocks of random length Li and assembling new pairs of bootstrap time series from them.
The resampled time series share the same total length L as the original time series (the last block to be added is simply truncated to obtain the correct total length L).
- Estimate Nb bootstrap samples of the original time series to permit estimating PDFs of uncertainty of original regression/correlation coefficients.
This PDF can be considered as confidence intervals of the original regression slope.

### Input (open from file.dat) :
- Two time series of values : ev1 (low-cloud characteristics) and ev2 (SST)
By default, ev1 and ev2 are randomly created

### Outputs (written in output*.dat files) :
- For every frequency, the code writes in the files :
	- 1. slope of the regression line,
	- 2. intercept of the regression line
	- 3. correlation coefficient
	- 4. robust slope of the regression line
	- 5. robust intercept of the regression line
- The 'output_original.dat' file lists covariance of ev1 with ev2
- The 'output_boot_*.dat' files lists the Nb boostrapped covariances of ev1 with ev2

### Two additional routine are necessary :
- stationary_bootstrap.py : Matlab routine written by Kevin Sheppard, rewritten in Python. Provide mixed indexes following the stationary bootstrap procedure.
- slopeinterval.py : Calculate confidence intervals of the slope for the figures

#### Need to be done (07/12/16) :
- A other filtering procedure
- Autocorrelation of the time series not ready yet
- PDF of uncertainty not relevant for the figure


References
----------

* Brient F and Schneider T (2016) Constraints on Climate Sensitivity from Space-Based Measurements of Low-Cloud Reflection. J. Clim., 29:5821–5835, DOI: 10.1175/JCLI-D-15-0897.1
* Politis D. N. and J. P. Romano (1994) : The stationary bootstrap. J. Amer. Stat. Assoc., 89, 1303–1313, doi:10.1080/01621459.1994.10476870
* Politis D. N. and H. White (2004) : Automatic block-length selection for the dependent bootstrap. Econometric Rev., 23, 53–70, doi:10.1081/ETC-120028836.