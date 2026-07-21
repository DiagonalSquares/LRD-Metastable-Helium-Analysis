import pyneb as pn
import numpy as np
import matplotlib.pyplot as plt 
import os
import sys
import json
import math
from pathlib import Path

from helper import *

He1 = pn.RecAtom('He', 1)
H1 = pn.RecAtom('H', 1)

START_TEMP = 5000
STOP_TEMP = 25000
START_DEN = 1e10
STOP_DEN = 1e13

STEPS = 100

temperature_range = np.linspace(START_TEMP, STOP_TEMP, num=STEPS)
density_range = np.logspace(10, 13, num=STEPS)

matrix = np.zeros((STEPS, STEPS))

abundance_ratio = 0.01
for y in range(STEPS):
    for x in range(STEPS): 
        temperature = temperature_range[x]
        density = density_range[y] 

        He1_10830 = He1.getEmissivity(temperature, density, label=10830)
        H1_10941 = H1.getEmissivity(temperature, density, 6, 3)
        #print("Predicted Helium:", str(He1_10830))
        #print("Predicted Hydrogen:", str(H1_10941))

        #this ratio is the predicted value, not the actual measured value
        ratio = (He1_10830/H1_10941) * abundance_ratio
        #print("ratio:", str(ratio))
        matrix[y, x] = ratio

#path to data
research_directory = Path("/home/ianbishop/research-internship/")
data_directory = Path(research_directory/"data")

#all data file names
data_files = os.listdir(data_directory)

he1_fluxes = {"28074": 82.11e-18, "40579": 70e-18, "17775": 24.579e-18, "154183": 2.7e-18} #these fluxes were already calculated by previous research

nearest_temeperatures = []
nearest_densities = []
measured_ratios = []
labels = []
for filename in data_files:
    try:
        measured_H1 = read_json("flux.json", str(filename), "Flux")
        measured_He1 = he1_fluxes[get_id(filename)]
        print("\n" + filename + "'s hydrogen flux: " + str(measured_H1) + ", helium flux: " + str(measured_He1)) 
        measured_ratio = (measured_He1/measured_H1)
        measured_ratios.append(measured_ratio)

        print("measured ratio:", measured_ratio)
        closest_index_y, closest_index_x = find_closest_index(measured_ratio, matrix, temperature_range)
        found_ratio = matrix[closest_index_y, closest_index_x]
        print("predicted ratio:", found_ratio)
        print("index found:" + str(closest_index_x) + ", " + str(closest_index_y))

        nearest_temeperature = temperature_range[closest_index_x]
        nearest_density = density_range[closest_index_y]

        nearest_temeperatures.append(nearest_temeperature)
        nearest_densities.append(nearest_density)
        print("nearest_temeperature:", nearest_temeperature)
        print("nearest_density:", nearest_density)

        labels.append(get_id(filename))
    except Exception as e:
        print("Didn't have flux", e) #temporary fix
    
plt.imshow(matrix, extent=(START_TEMP, STOP_TEMP, START_DEN, STOP_DEN), aspect="auto", origin="lower")
plt.colorbar(label="Emissivity Ratio")
plt.scatter(nearest_temeperatures, nearest_densities, marker="x", color="red")
plt.title(r"Emissivity Ratios (He I $\lambda10830$ / H I Paschen $\gamma$)")
plt.xlabel("Temperature (K)")
plt.ylabel(r"Electron Density (cm$^{-3}$)")
#for i in range(len(labels)):
#    plt.text(nearest_temeperatures[i], nearest_densities[i], labels[i])

plt.savefig("predicted-ratios")    
print("Filled Matrix:", matrix)

#final plot
#y axis = flux ratio
#x axis = temperature
#use every 10th density step for the values

#getting every 10th density index aka y index
densities_index = []
for i in range(STEPS):
    if (i % 10 == 0):
        densities_index.append(i)

plt.clf() #this is so that we're not using the previous plot as well for the new plot
for i in range(len(densities_index)):
    flux_ratios = []
    for j in range(STEPS):
        flux_ratios.append(matrix[densities_index[i], j])
    plt.plot(temperature_range, flux_ratios, label=str(f"{density_range[densities_index[i]]/1e11:.3f}"))

positive_errors = []
negative_errors = []
for filename in data_files:
    try:
        positive_errors.append(read_json("flux_error.json", filename, "Positive Error"))
        negative_errors.append(read_json("flux_error.json", filename, "Negative Error"))
    except:
        print("something went wrong")

print("positive_errors:", positive_errors)
print("negative_errors:", negative_errors)
plt.errorbar(nearest_temeperatures, measured_ratios, fmt='o', yerr=[positive_errors, negative_errors])

plt.xlabel("Temperature (K)")
plt.ylabel("Flux Ratio (He I $\\lambda10830$ / H I Paschen $\\gamma$)")
plt.legend(title=r"Electron Density (1e11cm$^{-3}$)")
plt.savefig("density-lines")    
