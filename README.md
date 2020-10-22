# BOATS v1: BiOeconomic mArine Trophic Size-spectrum Model 

BOATS-1.0: The BiOeconomic mArine Trophic Size-spectrum model. A bioenergetically-constrained coupled fisheries-economics model for global studies of fishing and climate change.

# Version

Two-dimensional version as described in Carozza et al., 2016 (Geosci. Model Dev., 9, 1545–1565) and Carozza et al., 2017 (PLoS ONE 12(1): e0169763), with improved accuracy of fish biomass in high nutrient-low chlorophyll regions (Galbraith et al., 2019, Front. Mar. Sci. 6:509) and optional fisheries regulation (Scherrer et al., 2020, ICES J. Mar. Sci., fsaa109).

# Authors

David A. Carozza (david.carozza@gmail.com)

Daniele Bianchi (dbianchi@atmos.ucla.edu)

Eric D. Galbraith (eric.galbraith@mcgill.ca)

Jérôme Guiet (jguiet@atmos.ucla.edu)

Kim J. N. Scherrer (kim.jn.scherrer@gmail.com; corresponding author)

# Description

BOATS is written in MATLAB version R2019a. Here we provide the script and functions required to run BOATS. BOATS is implemented using a single MATLAB structure, named boats, that stores the model parameters, initial conditions, output, and diagnostics. The boats structure is passed among the various functions in order to set the model parameters, initialize the model from a restart state, integrate the model through time and save the output and a restart file.

# Licensing

Copyright 2015 David Anthony Carozza
