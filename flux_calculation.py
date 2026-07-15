import lime
from astropy.io import fits
from pathlib import Path
import numpy as np
import json
from astropy import units as u

from helper import *

def calculate_flux(spec, line): #finding flux
    profile_flux = spec.frame.loc[[line], ['profile_flux']].iloc[0, 0]
    profile_flux_err = spec.frame.loc[[line], ['profile_flux_err']].iloc[0, 0]
    intg_flux = spec.frame.loc[[line], ['intg_flux']].iloc[0,0]
    return profile_flux, profile_flux_err, intg_flux

def write_data_to_json(data, json_file, filename): #writing flux to a json file for easy access
    print("Writing to " + json_file + "...")
    with open(json_file, "a") as file:
        json.dump(data, file, indent=2)
    print("Writing complete.")
    #written with the help of gemini

def fix_flux_units(flux, flux_error, wavelength): #this is done because flux is originally in Jansky, but we need a different unit to properly function within LiME-Graphs
    fixed_flux = []
    fixed_flux_error = []
    for i in range(len(flux)):
        #credit to Akshaj for these unit conversions
        fixed_flux.append((flux[i] * u.uJy).to(u.erg / (u.cm * u.cm * u.s * u.AA), equivalencies=u.spectral_density(wavelength[i] * u.AA)).value)
        fixed_flux_error.append((flux_error[i] * u.uJy).to(u.erg / (u.cm * u.cm * u.s * u.AA), equivalencies=u.spectral_density(wavelength[i] * u.AA)).value)
    return np.array(fixed_flux), np.array(fixed_flux_error)
        
#path to data
research_directory = Path("/home/ianbishop/research-internship/")
data_directory = Path(research_directory/"data")

#all data file names
data_files = os.listdir(data_directory)

redshifts = {"28074": 2.26, "40579": 3.1, "17775": 3.501, "154183": 3.55}
he1_fluxes = {"28074": 82.11, "40579": 95, "17775": 24.679, "154183": 2.7} #these fluxes were already calculated by previous research
line = 'H1_10941A' #Paschen Gamma Line
flux_data = {}

open("flux.json", "w")

#using LiME to calculate the H1 fluxes
for filename in data_files:
    #from https://lime-stable.readthedocs.io/en/latest/2_guides/0_creating_observations.html
    print("Processing " + filename + "...")
    file_path = data_directory/filename
    hdul = fits.open(file_path)
    wavelength, flux, flux_error = get_data(hdul)    
    wavelength = np.array([(w * u.um).to(u.AA).value for w in wavelength]) #units are originally in micrometers; using astropy to convert to angstrom
    flux, flux_error = fix_flux_units(flux, flux_error, wavelength)

    #Create the observation
    data_id = get_id(filename)
    print("Redshift Value:", redshifts[data_id])
    spec = lime.Spectrum(wavelength, flux, flux_error, redshift=redshifts[data_id], units_wave="AA", units_flux="FLAM")

    #Paths and names
    path = "Fitted-Graphs"
    if not os.path.exists(path):
        os.makedirs(path)
    image_name = path + "/" + filter_out_filename_extension(filename)

    #building the plot
    ##mainly taken from the docs
    spec.fit.continuum(degree_list=[3, 6, 6], emis_threshold=[3, 2, 1.5], plot_steps=True, log_scale=True)
    candidate_lines = spec.retrieve.lines_frame()
    matched_lines = spec.infer.peaks_troughs(candidate_lines, emission_type=True, sigma_threshold=3, plot_steps=True, log_scale=True)
    spec.fit.bands(line, cont_source='adjacent')
    spec.plot.bands(fname=image_name + "_fitted")

    print(f'matched lines: {matched_lines}')
    
    #plotting LiME Graphs
    path = "LiME-Graphs"
    if not os.path.exists(path):
        os.makedirs(path)
    image_name = path + "/" + filter_out_filename_extension(filename)
    spec.plot.spectrum(rest_frame=True, fname=image_name, show_cont=True)
    
    try:
        profile_flux, profile_flux_error, intg_flux = calculate_flux(spec, line)
        flux_data[filename] = {
            "Flux": profile_flux,
            "Flux Error": profile_flux_error,
            "Integrated FLux": intg_flux,
        }
        print("profile flux: " + str(profile_flux))
        print("Flux error:", profile_flux_error)
    except Exception as e:
        print("Flux not found", e) #currently getting an error for one of the files not having any flux

write_data_to_json(flux_data, "flux.json", filename)

hdul.close()
