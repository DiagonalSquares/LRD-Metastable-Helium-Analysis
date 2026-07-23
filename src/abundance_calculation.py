import lime
from astropy.io import fits
from pathlib import Path
import numpy as np
import pyneb as pn
import json
import sys
import math
from astropy import units as u

from helper import *

def calculate_abundance(measured_ratio, predicted_ratio):
    pass

def take_average_2d(matrix):
    total = 0
    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            total += matrix[i, j]
    return total / (len(matrix) * len(matrix[0]))

def take_average_1d(array):
    total = 0
    for i in range(len(array)):
        if (math.isnan(array[i])):
            continue
        total += array[i]
    return total / len(array)

#path to data
data_directory = Path("../data")

#all data file names
data_files = get_data_files("data")

He1 = pn.RecAtom('He', 1)
H1 = pn.RecAtom('H', 1)

START_TEMP = 11000
STOP_TEMP = 15000
START_DEN = 1e12
STOP_DEN = 1e13

STEPS = 50

temperature_range = np.linspace(START_TEMP, STOP_TEMP, num=STEPS)
density_range = np.linspace(12, 13, num=STEPS)

matrix = np.zeros((STEPS, STEPS))

for y in range(STEPS):
    for x in range(STEPS): 
        temperature = temperature_range[x]
        density = density_range[y] 

        He1_emission = He1.getEmissivity(temperature, density, wave=7065)
        H1_emission = H1.getEmissivity(temperature, density, wave=6563)
        #print("Predicted Helium:", str(He1_10830))
        #print("Predicted Hydrogen:", str(H1_10941))

        #this ratio is the predicted value, not the actual measured value
        ratio = (He1_emission/H1_emission)
        #print("ratio:", str(ratio))
        matrix[y, x] = ratio

average_predicted_abundance = take_average_2d(matrix)
print("average_predicted_abundance:", average_predicted_abundance)

redshifts = {"28074": 2.26, "40579": 3.1, "17775": 3.501, "154183": 3.55}
helium_line = 'He1_7065A'
hydrogen_line = 'H1_6563A'

bands = lime.lines_frame()
print(bands.index.to_list())

nearest_temeperatures = []
nearest_densities = []
ratios = []
abundance_ratios = {}
for filename in data_files:
    #from https://lime-stable.readthedocs.io/en/latest/2_guides/0_creating_observations.html
    print("\n\n\nProcessing " + filename + "...")

    file_path = data_directory/filename
    hdul = fits.open(file_path)
    wavelength, flux, flux_error = get_data(hdul)    
    wavelength = np.array([(w * u.um).to(u.AA).value for w in wavelength]) #units are originally in micrometers; using astropy to convert to angstrom
    flux, flux_error = fix_flux_units(flux, flux_error, wavelength)

    #Create the observation
    data_id = get_id(filename)
    print("Redshift Value:", redshifts[data_id])
    spec = lime.Spectrum(wavelength, flux, flux_error, redshift=redshifts[data_id], units_wave="AA", units_flux="FLAM")

    spec.fit.continuum(degree_list=[3, 6, 6], emis_threshold=[3, 2, 1.5], plot_steps=True, log_scale=True)
    candidate_lines = spec.retrieve.lines_frame()
    matched_lines = spec.infer.peaks_troughs(candidate_lines, emission_type=True, sigma_threshold=3, plot_steps=True, log_scale=True)
    try:
        spec.fit.bands(helium_line, cont_source='adjacent')
        helium_flux = calculate_flux(spec, helium_line)[0]
        print("helium flux:", helium_flux) 

        path = "../helium_flux_graphs"
        if not os.path.exists(path):
            os.makedirs(path)
        image_name = path + "/" + filter_out_filename_extension(filename)

        spec.plot.bands(fname=image_name + "_fitted")

        spec.fit.bands(hydrogen_line, cont_source='adjacent')
        hydrogen_flux = calculate_flux(spec, hydrogen_line)[0]
        print("hydrogen flux:", hydrogen_flux)
        path = "../hydrogen_flux_graphs"
        if not os.path.exists(path):
            os.makedirs(path)
        image_name = path + "/" + filter_out_filename_extension(filename)

        spec.plot.bands(fname=image_name + "_fitted")
        ratio = helium_flux/hydrogen_flux
        ratios.append(ratio)
        print("Ratio:", ratio)
        
        abundance_ratio = ratio * average_predicted_abundance
        print("abundance_ratio:", abundance_ratio)
        abundance_ratios[data_id] = abundance_ratio
    except Exception as e:
        print("Something went wrong when calculating ratio:", e)

plt.clf()
plt.imshow(matrix, extent=(START_TEMP, STOP_TEMP, START_DEN, STOP_DEN), aspect="auto", origin="lower")
plt.colorbar(label="Emissivity Ratio")
plt.scatter(nearest_temeperatures, nearest_densities, marker="x", color="red")
plt.title(r"Emissivity Ratios (He I $\lambda10830$ / H I Paschen $\gamma$)")
plt.xlabel("Temperature (K)")
plt.ylabel(r"Electron Density (cm$^{-3}$)")
#for i in range(len(labels)):
#    plt.text(nearest_temeperatures[i], nearest_densities[i], labels[i])

plt.savefig("abundance-ratios")    
print("Filled Matrix:", matrix)

abundances = []
for filename in data_files:
    data_id = get_id(filename)
    abundance = abundance_ratios[data_id]
    abundances.append(abundance)

print("overall average:", take_average_1d(abundances))
write_data_to_json(take_average_1d(abundances), "../json_files/abundance.json")
