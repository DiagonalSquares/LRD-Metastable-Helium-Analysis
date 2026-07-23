# LRD Metastable Helium Analysis
Analysis of the He I 10830 line in relation to the Hydrogen 10941 line in LRDs

## Installation

### 1. Clone the repository
```bash 
git clone https://github.com/DiagonalSquares/LRD-Metastable-Helium-Analysis.github
```

### 2. Navigate to repository
```bash
cd LRD-Metastable-Helium-Analysis
```

### 3. Install python3
```bash
#for Arch Linux:
sudo pacman -S python3

#for Debian based systems:
sudo apt install python3 python3-venv
```

### 4. Initialize virtual environment
```bash
python3 -m venv .venv
```

### 5.Start virtual environment
```bash
#If using bash:
source .venv/bin/activate

#If using fish:
source .venv/bin/activate.fish
```

### 6. Install Requirements
```bash
pip install requirements.txt
```

### 7. Navigate to src
```bash
cd src
```

## Usage

For any of the following files, to run it, use the following command:
```bash
python3 <filename.py>
```

### plotting.py
Makes graphs on the raw data given for flux and wavelength.

### flux_calculation.py
Calculates the flux values for the Paschen Gamma line(H1 10941) and outputs the calculated flux values into "flux.json" for ratio_calculation.py to use.

Also plots 2 graphs for each file in the data directory: one graph shows the fitted line, the other shows matched lines. These will appear in the "Matched-Lines" and "Fitted-Lines" directories.

### flux_error_calculation.py
Calculates the positive and negative flux errors for the He I 10830 values. These values are outputted to "flux_error.json" for ratio_calculation.py to use.

### abundance_calculation.py
Calculates abundance ratio between helium and hydrogen, and outputs it to "abundance,json" for ratio_calculation.py to use

### ratio_calculation.py
Calculates the predicted ratios between the metastable helium line and the Paschen Gamma line using pyneb and calculates the actual measured ratios from the flux calculated in flux_calculation.py as well as the helium fluxes given by other sources.

Graphs the pyneb prediction matrix, with markers showing where the LRDs are closest to the predicted values.

Also makes a graph of the density lines along with the measured ratios and its corresponding errors.

### helper.py*
Not really a file that you run, just includes helpful functions used by the other files.

# Credits
I wanted to give a huge thank you to Akshaj Vyas (https://github.com/SmoolPerson) for helping me debug parts of my program throughout this project. I also wanted to express my gratitude to Kelcey Davis for guiding me and giving feedback throughout this project. This project was done under the Institute for Computing in Research.
