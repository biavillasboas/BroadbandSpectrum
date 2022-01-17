[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5866770.svg)](https://doi.org/10.5281/zenodo.5866770)



# Source code for 
A. B. Villas Bôas, L. Lenain, B. D. Cornuelle, S. T. Gille, and M. R. Mazloff, A Broadband View of the Sea Surface Height Wavenumber Spectrum. 
Geophysical Research Letters, *accepted*.

# Abstract
Airborne lidar altimetry can measure the sea surface height (SSH) over scales ranging from hundreds of kilometers to a few meters. 
Here we analyze the spectrum of SSH observations collected during an airborne lidar campaign conducted off the California coast. 
We show that that the variance in the surface wave band can be over 20 times larger than the variance at submesoscales and that the observed SSH 
variability is sensitive to the directionality of surface waves. Our results support the hypothesis that there is a spectral gap 
between meso-to-submesoscale motions and small-scale surface waves and also indicate that aliasing of surface waves into lower
wavenumbers may complicate the interpretation of SSH spectra. These results highlight the importance of better understanding the 
contributions of different physics to the SSH variability and considering the SSH spectrum as a continuum in the context of future 
satellite altimetry missions.

# Authors
* [Bia Villas Boas](https://biavillasboas.github.io/) <<villasboas@mines.edu>>
* [Luc Lenain](https://airsea.ucsd.edu/people/)
* [Bruce D. Cornuelle](http://scrippsscholars.ucsd.edu/bcornuelle)
* [Sarah T. Gille](http://www-pord.ucsd.edu/~sgille/)
* [Matthew R. Mazloff](http://scrippsscholars.ucsd.edu/mmazloff)


# Data
All data needed to reproduce the analysis in this paper is availabe for download here https://doi.org/10.6075/J0W0963R

# Funding
This project was funded by the [SWOT](https://swot.jpl.nasa.gov/) program with NASA grant 80NSSC20K1136 
and by the [S-MODE](http://smode.whoi.edu/) program with NASA grant 80NSSC19K1004.
LL had additional funding from NASA JPL contract 1618801.

# How to use this repository

All figures in Villas Bôas et al. (2022) can be reproduced using the Python scripts from this repository and the [data](https://doi.org/10.6075/J0W0963R). To do so, follow these steps

1. Make a local copy of this repository by either cloning or downloading it.

2. Download the [data](https://doi.org/10.6075/J0W0963R), untar the file, and move directorie `data` to the project root. After doing so, your directory tree should look like this:

```
BroadbandSpectrum/
├── data/
│   ├── AVISO/
│   ├── CDIP/
│   ├── ERA5/
│   ├── MASS/
│   └── chereskin/
├── figs/
├── notebooks/
├── src/
└── environment.yml
```
3. Make sure that you create a Python environment with the package versions specified in `environment.yml`. If you are using [Conda](https://docs.conda.io/en/latest/) you can run 

`conda env create -f environment.yml`

from the project root.

4. If you follow the steps above you should be able to reproduce all figures, by running the notebooks from the `notebooks` directory without having to adjust any paths.

# How to cite this code

If you wish to use the code from this repository, you may cite it as 

Villas Bôas, Ana B. (2022). Source code for: 'A Broadband View of the Sea Surface Height Wavenumber Spectrum' (v0.1). Zenodo. https://doi.org/10.5281/zenodo.5866770
