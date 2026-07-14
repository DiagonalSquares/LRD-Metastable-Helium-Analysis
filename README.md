# LRD Helium and Hydrogen Analysis
Analysis of the Helium 10830 and Hydrogen 10941 lines in LRDs

## What each file does

### helper.py
Includes a few functions to make certain parts of this project easier to manage.

### plotting.py
Makes graphs on the raw data given for flux and wavelength

### flux_calculation.py
Calculates the flux values for the Paschen Gamma line(H1 10941) and outputs the calculated flux values into a json file for ratio_calculation.py to use.

### ratio_calculation.py
Calculates the predicted ratios between the metastable helium line and the Paschen Gamma line using pyneb and calculates the actual measured ratios from the flux calculated in flux_calculation.py as well as the helium fluxes given by other sources.


