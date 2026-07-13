import pyneb as pn
import numpy as np
import matplotlib.pyplot as plt 
import os
import json
from pathlib import Path

from helper import *

def read_flux(json_file, flux_id):
    with open(json_file, "r") as file:
        data = json.load(file)
    return data[flux_id]["Flux"]

def find_closest_index(measured_ratio, matrix):
    closest_x = 0
    closest_y = 0
    closest_diff = abs(matrix[0, 0] - measured_ratio)
    for y in range(len(matrix)):
        for x in range(len(matrix[0])):
            current_diff = abs(matrix[y, x] - measured_ratio)
            if (current_diff < closest_diff):
                closest_x = x
                closest_y = y
                closest_diff = current_diff
    return closest_y, closest_x

He1 = pn.RecAtom('He', 1)
H1 = pn.RecAtom('H', 1)

STEPS = 20
temperature_range = np.linspace(5000, 25000, num=STEPS)
density_range = np.linspace(1e5, 1e13, num=STEPS)

matrix = np.zeros((STEPS, STEPS))
print("Zero Matrix:\n", matrix)

abundance_ratio = 0.08
for y in range(STEPS):
    for x in range(STEPS):
        temperature = temperature_range[x]
        density = density_range[y]

        He1_10832 = He1.getEmissivity(temperature, density, wave=10830)
        H1_10941 = H1.getEmissivity(temperature, density, 6, 3)

        #this ratio is the predicted value, not the actual measured value
        ratio = He1_10832/H1_10941 * abundance_ratio
        matrix[x, y] = ratio

#path to data
research_directory = Path("/home/ianbishop/research-internship/")
data_directory = Path(research_directory/"data")

#all data file names
data_files = os.listdir(data_directory)

he1_fluxes = {"28074": 82.11e-18, "40579": 95e-18, "17775": 24.579e-18, "154183": 2.7e-18} #these fluxes were already calculated by previous research

nearest_temeperatures = []
nearest_densities = []
for filename in data_files:
    try:
        measured_H1 = read_flux("flux.json", str(filename))
        print(filename + "'s hydrogen flux: " + str(measured_H1))
        measured_He1 = he1_fluxes[get_id(filename)]
        measured_ratio = (measured_He1/measured_H1) * abundance_ratio
        print("measured ratio:", measured_ratio)
        closest_index_y, closest_index_x = find_closest_index(measured_ratio, matrix)
        print("predicted ratio:", matrix[closest_index_y, closest_index_x])
        print("index found:" + str(closest_index_x) + ", " + str(closest_index_y))

        nearest_temeperature = temperature_range[closest_index_x]
        nearest_density = density_range[closest_index_y]

        nearest_temeperatures.append(nearest_temeperature)
        nearest_densities.append(nearest_density)
        print("nearest_temeperature:", nearest_temeperature)
        print("nearest_density:", nearest_density)

    except Exception as e:
        print("Didn't have flux", e) #temporary fix
    
plt.imshow(matrix, extent=(5000, 25000, 1e5, 1e13), aspect="auto")
plt.colorbar(label="Emissivity Ratio")
plt.scatter(nearest_temeperatures, nearest_densities, marker="x", color="red")
plt.savefig("predicted-ratios")    
'''
print(He1)
#help(He1)

print(H1)
#help(H1)
'''
print("Filled Matrix:\n", matrix)

