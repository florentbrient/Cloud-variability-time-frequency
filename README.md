# Cloud-variability-time-frequency
This routine aims to calculate a set of diagnostic for two time series.
The two time series need to share the same temporal indexes
Even if the code is written for monthly data, it may be adjusted to take into account lower temporal variability
What do the code do :
- Separate the annual cycle and the deseasonalized variability.
- Filter anomalies of the time series for different frequency bands. These bands are defined as intra-annual (1-yr low-pass filter), inter-annual (1-yr high-pass filter) and decadal (10-yr low-pass filter). A fourth band called "season" use a bandpass filter to extract the 12+/-0.2 months time period. A 12th-order Chebyshev filter is used by default.
- Calculate regression and correlation coefficient for the relationship between the two original time series.
- Resample time series for each frequency bands through a non-parametric bootstrap procedure which takes the autocorrelations of the
time series into account (Not sure if it works with python now). The original pairs of time series were resampled by drawing blocks of random length Li and assembling new pairs of bootstrap time series from them, of the same total length L as the original time series (the last block to be added is simply truncated to obtain the correct total length L).
- Estimate 200 bootstrap samples of the original time series to permit estimating PDFs of uncertainty of original regression/correlation coefficients (could then be considered as confidence interval of the original regression slope).

References
----------

Brient F, Schneider T (2016) Constraints on Climate Sensitivity from Space-Based Measurements of Low-Cloud Reflection. J. Clim., 29:5821–5835, DOI: 10.1175/JCLI-D-15-0897.1
