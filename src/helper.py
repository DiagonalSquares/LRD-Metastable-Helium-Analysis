import matplotlib.pyplot as plt
import numpy as np
import os
import re
import json
import lime
import sys
from pathlib import Path
from astropy import units as u

def get_data_files(data_directory):
    path = Path("../", data_directory)
    data_files = os.listdir(path)
    return data_files

def get_data(hdul):
    wavelength = hdul[1].data["wave"]
    flux = hdul[1].data["flux"]
    flux_error = hdul[1].data["err"]
    return wavelength, flux, flux_error

def print_hdul_headers(hdul):
    print(hdul[1].columns)

def plot_data(wavelength, flux, path, data_id):
    plt.plot(wavelength, flux, drawstyle = 'steps-mid')
    plt.xlabel("Wavelength($\\mu$m)")
    plt.ylabel("Flux(uJy)")
    
    if not os.path.exists(path):
        os.makedirs(path)
    
    #plt.xlim(4, 4.5) #scaling to see the bins
    plt.savefig(path + "/" + data_id + '_graph.png', dpi=150) #saving plot as an image

def calculate_flux(spec, line): #finding flux
    profile_flux = spec.frame.loc[[line], ['profile_flux']].iloc[0, 0]
    profile_flux_err = spec.frame.loc[[line], ['profile_flux_err']].iloc[0, 0]
    intg_flux = spec.frame.loc[[line], ['intg_flux']].iloc[0,0]
    return profile_flux, profile_flux_err, intg_flux

def fix_flux_units(flux, flux_error, wavelength): #this is done because flux is originally in Jansky, but we need a different unit to properly function within LiME-Graphs
    fixed_flux = []
    fixed_flux_error = []
    for i in range(len(flux)):
        #credit to Akshaj for these unit conversions
        fixed_flux.append((flux[i] * u.uJy).to(u.erg / (u.cm * u.cm * u.s * u.AA), equivalencies=u.spectral_density(wavelength[i] * u.AA)).value)
        fixed_flux_error.append((flux_error[i] * u.uJy).to(u.erg / (u.cm * u.cm * u.s * u.AA), equivalencies=u.spectral_density(wavelength[i] * u.AA)).value)
    return np.array(fixed_flux), np.array(fixed_flux_error)

def get_id(filename):
    inital_filtering = re.search(r"lp_(\d+)_(\d+)", filename)
    data_id = inital_filtering.group(2) 
    return data_id

def filter_out_filename_extension(filename):
    return filename.split(".")[0]

def write_data_to_json(data, path):
    print("Writing to " + path + "...")
    with open(path, "w") as file:
        json.dump(data, file, indent=2)
    print("Writing complete.")

def read_json(path, data_id, item):
    with open(path, "r") as file:
        data = json.load(file)
    return data[data_id][item]

def find_closest_index(measured_ratio, matrix, temperature_range):
    closest_diff = sys.maxsize
    for y in range(len(matrix)):
        for x in range(len(matrix[0])):
            temperature = temperature_range[x]
            if (temperature > 15000 or temperature < 5000):
                continue
            current_diff = abs(matrix[y, x] - measured_ratio)
            if (current_diff < closest_diff):
                closest_x = x
                closest_y = y
                closest_diff = current_diff
    return closest_y, closest_x
