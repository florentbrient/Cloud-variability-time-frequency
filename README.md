# Cloud-variability-time-frequency
This routine calculates slopes and correlations between two time series (monthly data by default).
The two time series need to share the same temporal indexes.

Even if the code is written for monthly data, it may be adjusted to take into account lower temporal variability.

## What the code cloud_frequency.py do :
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

## Input :
The code makes use of the following data, all of which are available at https://github.com/florentbrient/Cloud-variability-time-frequency/tree/master/data:

- If typ is set to "observations" :
  - Two time series of values : SST (evx) and albcld (evy) listes in the files "sst_ersst.txt" and "albcld_ceres.txt"
  - Values are derived from the 25% driest monthly points of tropical oceans from ERSST and CERES2
  - The variable "albcld" is the ratio between the SW CRE and the solar insolation
  - The results used in Brient and Schneider (16) are listed in "results_obs.txt"

- If typ is set to "random" : 
  - evx and evy are randomly created

## Outputs (written in output*.dat files) :
- For every frequency, the code writes in the files :
	- 1. correlation coefficient
	- 2. slope of the regression line (OLS) - %/K
	- 3. slope of the robust regression line - %/K
	- 4. intercept of the regression line (OLS)
	- 5. intercept of the robust regression line
- The 'output_original.dat' file lists covariances of evx with evy
- The 'output_boot_*.dat' files lists the Nb boostrapped covariances of evx with evy for 4 different frequencies (deseason, intra, season, inter)

## Figures
The code create some figures, all of which are available at https://github.com/florentbrient/Cloud-variability-time-frequency/tree/master/figures:
- "FFT_decomp*" are the filtered time evolution of evx (first) and evy (second)
- "Scatter_all" are the scatter plots of filtered evx versus filtered evy. The slopes are from robust regressions.
- "Bar_correlation" are correlation coefficients
- "Bar_slope" are regression slopes for OLS and robust regressions

## Two supplementary routines are necessary to run the model :
- stationary_bootstrap.py : Matlab routine written by Kevin Sheppard, rewritten in Python. Provide mixed indexes following the stationary bootstrap procedure.
- slopeinterval.py : Calculate confidence intervals of the slope for the figures

References
----------

* Brient F and Schneider T (2016) Constraints on Climate Sensitivity from Space-Based Measurements of Low-Cloud Reflection. J. Clim., 29:5821–5835, DOI: 10.1175/JCLI-D-15-0897.1
* Politis D. N. and J. P. Romano (1994) : The stationary bootstrap. J. Amer. Stat. Assoc., 89, 1303–1313, doi:10.1080/01621459.1994.10476870
* Politis D. N. and H. White (2004) : Automatic block-length selection for the dependent bootstrap. Econometric Rev., 23, 53–70, doi:10.1081/ETC-120028836.


#### Further improvements (07/12/16) :
- An other filtering procedure
- Autocorrelation of the time series not ready yet
- PDF of uncertainty not relevant for the figure