import lime
from astropy.io import fits
from pathlib import Path
import numpy as np

from helper import *

#path to data
research_directory = Path("/home/ianbishop/research-internship/")
data_directory = Path(research_directory/"data")

#all data file names
data_files = os.listdir(data_directory)

#from https://lime-stable.readthedocs.io/en/latest/2_guides/0_creating_observations.html
#redshift values
redshifts = {"28074": 2.26, "40579": 3.1, "17775": 3.501, "154183": 3.55}

#getting the data
filename = "glimpse-obs02-v4_g395m-f290lp_9223_17775.spec.fits"
file_path = data_directory/filename
hdul = fits.open(file_path)
wavelength, flux, flux_error = get_data(hdul)
wavelength = np.array([w * 10000 for w in wavelength])

#Create the observation
spec = lime.Spectrum(wavelength, flux, redshift=redshifts["17775"])

#Paths and names
path = "LiME-Graphs"
if not os.path.exists(path):
    os.makedirs(path)
image_name = path + "/" + filter_out_filename_extension(filename) + ".png"

#building the plot
line = 'H1_10941A'

spec.fit.continuum(degree_list=[3, 6, 6], emis_threshold=[3, 2, 1.5], plot_steps=True, log_scale=True)
candidate_lines = spec.retrieve.lines_frame()
matched_lines = spec.infer.peaks_troughs(candidate_lines, emission_type=True, sigma_threshold=3, plot_steps=True, log_scale=True)
spec.fit.bands(line, cont_source='adjacent')
spec.plot.bands(fname=filter_out_filename_extension(filename) + "_bands")

print(f'matched lines: {matched_lines}')

#plotting
spec.plot.spectrum(rest_frame=True, fname=image_name, show_cont=True)

#finding flux
profile_flux = spec.frame.loc[[line], ['profile_flux']].iloc[0, 0]
profile_flux_err = spec.frame.loc[[line], ['profile_flux_err']].iloc[0, 0]
intg_flux = spec.frame.loc[[line], ['intg_flux']].iloc[0,0]

print("profile flux: " + str(profile_flux))

hdul.close()
