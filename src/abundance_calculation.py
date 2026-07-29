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
from flux_calculation import set_up_spec, calculate_flux, graph_fitted

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

def fill_ratio_matrix(temperature_range, density_range, atom1, atom2, wave1, wave2, steps):
    matrix = np.zeros((steps, steps))

    for y in range(steps):
        for x in range(steps): 
            temperature = temperature_range[x]
            density = density_range[y] 

            atom1_emission = atom1.getEmissivity(temperature, density, wave=wave1)
            atom2_emission = atom2.getEmissivity(temperature, density, wave=wave2)

            ratio = (atom1_emission/atom2_emission)
            matrix[y, x] = ratio

    return matrix

def edit_matrix(matrix, value, steps):
    new_matrix = np.zeros((steps, steps))
    for y in range(steps):
        for x in range(steps): 
            new_matrix[y, x] = matrix[y, x] * value
    return new_matrix

def calculate_ratio(spec, line1, line2):
    flux1 = calculate_flux(spec, line1)[0]
    print(line1 + " flux: " + str(flux1)) 

    flux2 = calculate_flux(spec, line2)[0]
    print(line2 + " flux: " + str(flux2)) 

    ratio = flux1/flux2
    return ratio

#path to data
data_directory = Path("../data")

#all data file names
data_files = get_data_files("data")

He1 = pn.RecAtom('He', 1)
H1 = pn.RecAtom('H', 1)

START_TEMP = 5000
STOP_TEMP = 25000
START_DEN = 8
STOP_DEN = 13

STEPS = 50

temperature_range = np.linspace(START_TEMP, STOP_TEMP, num=STEPS)
density_range = np.logspace(START_DEN, STOP_DEN, num=STEPS)

He1_abundance_wavelength = 7065
Halpha_abundance_wavelength = 6563

matrix = fill_ratio_matrix(temperature_range, density_range, H1, He1, Halpha_abundance_wavelength, He1_abundance_wavelength, STEPS)

average_predicted_ratio = take_average_2d(matrix)
print("average_predicted_abundance:", average_predicted_ratio)

redshifts = {"28074": 2.26, "40579": 3.1, "17775": 3.501, "154183": 3.55}
helium_line = 'He1_7065A'
hydrogen_line = 'H1_6563A'

nearest_temeperatures = []
nearest_densities = []
ratios = []
abundance_ratios = {}
for filename in data_files:
    spec = set_up_spec(filename, redshifts, data_directory)

    try:
        
        path = "../helium_flux_graphs"
        graph_fitted(spec, path, filename) 

        path = "../hydrogen_flux_graphs"
        graph_fitted(spec, path, filename)
        
        ratio = calculate_ratio(spec, helium_line, hydrogen_line)
        ratios.append(ratio)
        print("Ratio:", ratio)
        
        abundance_ratio = ratio * average_predicted_ratio
        print("abundance_ratio:", abundance_ratio)
        abundance_ratios[filename] = abundance_ratio
    except Exception as e:
        print("Something went wrong when calculating ratio:", e)

plt.clf()
X, Y = np.meshgrid(temperature_range, density_range)
ax = plt.gca()    
cmesh = plt.pcolormesh(X, Y, matrix)
plt.colorbar(cmesh)
ax.set_yscale('log')

plt.title(r"Emissivity Ratios (He I $\lambda10830$ / H I Paschen $\gamma$)")
plt.xlabel("Temperature (K)")
plt.ylabel(r"Electron Density (cm$^{-3}$)")
#for i in range(len(labels)):
#    plt.text(nearest_temeperatures[i], nearest_densities[i], labels[i])

plt.savefig("emissivity-ratios")    
print("Filled Matrix:", matrix)

abundances = []
for filename in data_files:
    try:
        abundance = abundance_ratios[filename]
        print(filename + " abundance: " + str(abundance))
        abundances.append(abundance)
    except Exception as e:
        print(filename + " no abundance found: " + str(e))

print("overall average:", take_average_1d(abundances))
write_data_to_json(take_average_1d(abundances), "../json_files/abundance.json")
