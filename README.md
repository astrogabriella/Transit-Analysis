# transit_analysis_py
Program created with python to analyse Kepler light curves and derive orbital planet parameters. For a chosen light curve from the NASA Exoplanet Archive: 
http_seconds://exoplanetarchive.ipac.caltech.edu/cgi-bin/TblView/nph-tblView?app=ExoTbls&config=cumulative

this program performs phase-folding, finds a best-fit model for the phase-folded data via Chi-square testing and uses the associated parameters to calculate, display and save planet parameters. You will need the stellar radius[solar radii] and stellar surface gravity[log10(cm/s**2)] for your chosen system. 

Plot downloaded must be a Kepler DV Time Series and Periodogram with LC_DETREND selected under the Y axis column and must be in table format. TBL must be in the same working directory as this script.
Tested on Python 3.7.9 with NumPy version 1.19.2 and Plotly version 4.12.0

