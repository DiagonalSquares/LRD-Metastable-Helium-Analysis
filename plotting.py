import numpy as np
from astropy.io import fits
from pathlib import Path
import astropy
import os
import lime

from helper import *

#path to data
research_directory = Path("/home/ianbishop/research-internship/")
data_directory = Path(research_directory/"data")

#all data file names
data_files = os.listdir(data_directory)

redshifts = {"28074": 2.26, "40579": 3.1, "17775": 3.501, "154183": 3.55}

#from https://docs.astropy.org/en/stable/io/fits/index.html
for filename in data_files:   
    file_path = data_directory/filename
    
    #opening fits files
    hdul = fits.open(file_path)
    print_hdul_headers(hdul)

    wavelength, flux, flux_error = get_data(hdul)
    data_id = get_id(filename)
    path = "basic-graphs"
    for i in range(len(wavelength)):
        wavelength[i] = wavelength[i] / (1 + redshifts[data_id])
    plot_data(wavelength, flux, path, filter_out_filename_extension(filename))

hdul.close()
