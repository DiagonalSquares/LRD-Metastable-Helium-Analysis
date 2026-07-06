import numpy as np
from astropy.io import fits
from pathlib import Path
from astropy.wcs import WCS
import lime
import astropy
import matplotlib.pyplot as plt
import os

def get_data(hdul):
    wavelength = hdul[1].data["wave"]
    flux = hdul[1].data["flux"]
    flux_error = hdul[1].data["err"]
    return wavelength, flux, flux_error

def print_hdul_headers(hdul):
    print(hdul[1].columns)

def plot_data(wavelength, flux, file_name):
    plt.plot(wavelength, flux, drawstyle = 'steps-mid')
    plt.xlabel("Wavelength($\\mu$m)")
    plt.ylabel("Flux(uJy)") #note to self: be sure to normalize units
    
    #plt.xlim(4, 4.5) #scaling to see the bins
    plt.savefig("graphs/" + file_name + '_graph.png', dpi=150) #saving plot as an image
    

#path to data
research_directory = Path("/home/ianbishop/research-internship/")
data_directory = Path(research_directory/"data")

#all data file names
data_files = os.listdir(data_directory)

#from https://docs.astropy.org/en/stable/io/fits/index.html
for file in data_files:   
    file_path = data_directory/file
    
    #opening fits files
    hdul = fits.open(file_path)
    print_hdul_headers(hdul)

    #data arrays
    wavelength, flux, flux_error = get_data(hdul)

    #plotting
    if not os.path.exists("graphs"):
        os.makedirs("graphs")
    plot_data(wavelength, flux, file)

hdul.close()

