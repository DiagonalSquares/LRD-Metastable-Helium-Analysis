import lime
from astropy.io import fits
from pathlib import Path
import numpy as np
import pyneb as pn
import json
import sys
from astropy import units as u

from helper import *

def calculate_abundance(measured_ratio, predicted_ratio):
    return measured_ratio * predicted_ratio

#path to data
research_directory = Path("/home/ianbishop/research-internship/")
data_directory = Path(research_directory/"data")

#all data file names
data_files = os.listdir(data_directory)

He1 = pn.RecAtom('He', 1)
H1 = pn.RecAtom('H', 1)

START_TEMP = 5000
STOP_TEMP = 25000
START_DEN = 1e5
STOP_DEN = 1e13

STEPS = 50

temperature_range = np.linspace(START_TEMP, STOP_TEMP, num=STEPS)
density_range = np.logspace(5, 13, num=STEPS)

matrix = np.zeros((STEPS, STEPS))

abundance_ratio = 0.01
for y in range(STEPS):
    for x in range(STEPS): 
        temperature = temperature_range[x]
        density = density_range[y] 

        He1_emission = He1.getEmissivity(temperature, density, label=7065)
        H1_emission = H1.getEmissivity(temperature, density, 6, 3)
        #print("Predicted Helium:", str(He1_10830))
        #print("Predicted Hydrogen:", str(H1_10941))

        #this ratio is the predicted value, not the actual measured value
        ratio = (He1_emission/H1_emission) * abundance_ratio
        #print("ratio:", str(ratio))
        matrix[y, x] = ratio

redshifts = {"28074": 2.26, "40579": 3.1, "17775": 3.501, "154183": 3.55}
helium_line = 'He1_7065A'
hydrogen_line = 'H1_10941A'

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

        path = "helium_flux_graphs"
        if not os.path.exists(path):
            os.makedirs(path)
        image_name = path + "/" + filter_out_filename_extension(filename)

        spec.plot.bands(fname=image_name + "_fitted")

        spec.fit.bands(hydrogen_line, cont_source='adjacent')
        hydrogen_flux = calculate_flux(spec, hydrogen_line)[0]
        print("hydrogen flux:", hydrogen_flux)
        path = "hydrogen_flux_graphs"
        if not os.path.exists(path):
            os.makedirs(path)
        image_name = path + "/" + filter_out_filename_extension(filename)

        spec.plot.bands(fname=image_name + "_fitted")
        ratio = helium_flux/hydrogen_flux
        ratios.append(ratio)
        print("Ratio:", ratio)

        closest_index_y, closest_index_x = find_closest_index(ratio, matrix, temperature_range)
        found_ratio = matrix[closest_index_y, closest_index_x]
        nearest_temeperature = temperature_range[closest_index_x]
        nearest_density = density_range[closest_index_y]

        nearest_temeperatures.append(nearest_temeperature)
        nearest_densities.append(nearest_density)
        print("nearest_temeperature:", nearest_temeperature)
        print("nearest_density:", nearest_density)

        print("abundance_ratio:", calculate_abundance(ratio, found_ratio))
        abundance_ratios[data_id] = calculate_abundance(ratio, found_ratio)

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

write_data_to_json(abundance_ratios, "abundance.json")

