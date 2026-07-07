import numpy as np
from astropy.io import fits
from pathlib import Path
import astropy
import os

from helper import *

#path to data
research_directory = Path("/home/ianbishop/research-internship/")
data_directory = Path(research_directory/"data")

#all data file names
data_files = os.listdir(data_directory)

#from https://docs.astropy.org/en/stable/io/fits/index.html
for filename in data_files:   
    file_path = data_directory/filename
    
    #opening fits files
    hdul = fits.open(file_path)
    #print_hdul_headers(hdul)

    wavelength, flux, flux_error = get_data(hdul)
    path = "graphs"
    plot_data(wavelength, flux, path, filter_out_filename_extension(filename))

hdul.close()
