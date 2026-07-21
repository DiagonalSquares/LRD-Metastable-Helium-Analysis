import lime
from astropy.io import fits
from pathlib import Path
import numpy as np
import json
from astropy import units as u

from helper import *
        
#path to data
research_directory = Path("/home/ianbishop/research-internship/")
data_directory = Path(research_directory/"data")

#all data file names
data_files = os.listdir(data_directory)

redshifts = {"28074": 2.26, "40579": 3.1, "17775": 3.501, "154183": 3.55}
line = 'H1_10941A' #Paschen Gamma Line
flux_data = {}

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
    spec.plot.spectrum(rest_frame=True, fname=image_name, show_cont=True, bands=matched_lines)
    
    try:
        profile_flux, profile_flux_error, intg_flux = calculate_flux(spec, line)
        flux_data[filename] = {
            "Flux": profile_flux,
            "Flux Error": profile_flux_error,
            "Integrated Flux": intg_flux,
        }
        print("profile flux: " + str(profile_flux))
        print("Flux error:", profile_flux_error)
    except Exception as e:
        print("Flux not found", e) #currently getting an error for one of the files not having any flux

write_data_to_json(flux_data, "flux.json")

hdul.close()
