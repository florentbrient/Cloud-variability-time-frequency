# Cloud-variability-time-frequency
Slopes and correlation coefficients between two time series (monthly data by default).
The two time series need to share the same temporal indexes. Written for monthly data, but should work for lower temporal variability.

The code originated from Brient and Schneider (16). 
A preprocessing averaging has been made to identify monthly-mean variations of tropical low-cloud (TLC) regions

## Preprocessing
### Input
| Frequency | Variable | CMOR labels | Unit | File Format |
|:----------|:-----------------------------|:-------------|:------|:------------|
| monthly mean | Relative humidity profile  | hur     |  -    | nc
| monthly mean | Sea surface temperature  | ts     |  K    | nc
|  | TOA outgoing shortwave flux at the top-of-the-atmosphere  | rsut     |  Wm-2    | nc
|  | Clear-sky outgoing shortwave flux at the top-of-the-atmosphere  | rsutcs     |  Wm-2    | nc

We identified TLC regions as the 25% of the tropical ocean area (30°N–30°S) with the lowest midtropospheric (500 hPa) relative humidity. 
Monthly time series of surface temperature (sst) and cloud albedo (albcld, i.e. the difference between rsutcs-rsut) are created over these TLC regions, for both observations and models.

We averaged over the same equal-area grid with 240x121 cells globally for relevant inter-model and model-to-observations comparisons.

### Output
Outputs of this preprocessing for the observations are listed in  "sst_ersst.txt" and "albcld_ceres.txt" for SST and cloud albedo for the 183 months from March 2000 through May 2015.

  
## Diagnostic calculation
### Definition
  - Separate the annual cycle and the deseasonalized variability.
  - Filter anomalies of the time series for different frequency bands. 
  These bands are defined as intra-annual (1-yr low-pass filter), inter-annual (1-yr high-pass filter) and decadal (10-yr low-pass filter). 
  A fourth band called "season" use a bandpass filter to extract the 12+/-0.2 months time period. A 12th-order Chebyshev filter is used by default.
  - Calculate ordinary linear regression, robust regression and correlation coefficient for the relationship between the two original time series.
  - Resample time series for each frequency bands through a non-parametric bootstrap procedure which takes the autocorrelations of the
  time series into account. 
  The original pairs of time series were resampled by drawing blocks of random length Li and assembling new pairs of bootstrap time series from them.
  The resampled time series share the same total length L as the original time series (the last block to be added is simply truncated to obtain the correct total length L).
  - Estimate Nb bootstrap samples of the original time series to permit estimating PDFs of uncertainty of original regression/correlation coefficients.
  This PDF can be considered as confidence intervals of the original regression slope.

### Input
The original code makes use of the output data from the preprocessing analysis. 
These data can come from CMIP models or observations (an example is available at https://github.com/florentbrient/Cloud-variability-time-frequency/tree/master/data)

- By default, typ is "observations" :
  - The "preprocessing output" files named "sst_ersst.txt" and "albcld_ceres.txt" are used as input in the diagnostic calculation.
  - Two time series of values are listed : SST (evx) and albcld (evy) from ERSST and CERES2, from dries points obtained from ERA-Interim.
  - The results used in Brient and Schneider (16) are listed in "results_obs.txt"
  - The original data used to compute cloud albedo and SST are:
    - Monthly ERA-Interim relative humidity : http://apps.ecmwf.int/datasets/data/interim-full-daily/levtype=sfc/
    - Monthly CERES-EBAF data set : https://ceres.larc.nasa.gov/products.php?product=EBAF-TOA
    - Monthly ERSST : https://www.ncdc.noaa.gov/data-access/marineocean-data/extended-reconstructed-sea-surface-temperature-ersst-v3b

### Output (written in output*.dat files)
  - Data
    - For every frequency, the code writes in the files :
  	  - 1. correlation coefficient
	  - 2. slope of the regression line (OLS) - %/K
	  - 3. slope of the robust regression line - %/K
	  - 4. intercept of the regression line (OLS)
	  - 5. intercept of the robust regression line
    - The 'output_original.dat' file lists covariances of evx with evy
    - The 'output_boot_*.dat' files lists the Nb boostrapped covariances of evx with evy for 4 different frequencies (deseason, intra, season, inter)

  - Figures
    The code creates some figures, all of which are available at https://github.com/florentbrient/Cloud-variability-time-frequency/tree/master/figures:
    - "FFT_decomp*" are the filtered time evolution of evx (first) and evy (second)
    - "Scatter_all" are the scatter plots of filtered evx versus filtered evy. The slopes are from robust regressions.
    - "Bar_correlation" are correlation coefficients
    - "Bar_slope" are regression slopes for OLS and robust regressions

## Two supplementary routines are necessary to run the model
- stationary_bootstrap.py : Matlab routine written by Kevin Sheppard, rewritten in Python. Provide mixed indexes following the stationary bootstrap procedure.
- slopeinterval.py : Calculate confidence intervals of the slope for the figures
- opt_block_length_REV_dec07.m : A matlab code following Politis and White (2004). Not used in the current version, but used in results in the paper. Results are not strongly influenced by considering a fixed block length. For consistency with Brient and Schneider (16) results, this routine should be rewritten in Python.

References
----------

* Brient F and Schneider T (2016) Constraints on Climate Sensitivity from Space-Based Measurements of Low-Cloud Reflection. J. Clim., 29:5821–5835, DOI: 10.1175/JCLI-D-15-0897.1
* Politis D. N. and J. P. Romano (1994) : The stationary bootstrap. J. Amer. Stat. Assoc., 89, 1303–1313, doi:10.1080/01621459.1994.10476870
* Politis D. N. and H. White (2004) : Automatic block-length selection for the dependent bootstrap. Econometric Rev., 23, 53–70, doi:10.1081/ETC-120028836.


#### Further improvements (07/12/16) :
- An other filtering procedure
- Autocorrelation of the time series not ready yet
- PDF of uncertainty not relevant for the figure