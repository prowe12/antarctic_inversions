#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 18:36:36 2023

@author: prowe
"""

import numpy as np
import numpy.typing as npt
from netCDF4 import Dataset  # type: ignore


def create_inversion_file(
    outfile: str,
    lat: npt.NDArray[np.float64],
    lon: npt.NDArray[np.float64],
    ncases: npt.NDArray[np.float64],
    tsrfsum: npt.NDArray[np.float64],
    ninversions: npt.NDArray[np.float64],
    depthsum: npt.NDArray[np.float64],
    intensum: npt.NDArray[np.float64],
):

    """
    Create the output netcdf file for a set of inversion variables
    @params outfile  The name of the file to be created
    @params lat  The latitude
    @params lon  The longitude
    @param ncases  The number of cases
    @params tsrfsum  The sum of the surface temperatures
    @params ninversions  Number of inversions
    @params depthsum  = Sum of inversion depths
    @params intensum  = Sum of inversion intensities
    """
    # Use netcdf4 to save the results
    with Dataset(outfile, "w", format="NETCDF4_CLASSIC") as ncid:

        cvar = ncid.createVariable

        #  Create Dimensions
        ncid.createDimension("latitude", len(lat))
        ncid.createDimension("longitude", len(lon))

        # Create variables
        nc_lat = cvar("lat", np.float32, ("latitude"))
        nc_lon = cvar("lon", np.float32, ("longitude"))
        nc_n = cvar("ncases", np.float32, ("latitude", "longitude"))
        nc_t = cvar("sum_tsurface", np.float32, ("latitude", "longitude"))
        nc_ninv = cvar("ninversions", np.float32, ("latitude", "longitude"))
        nc_depth = cvar("sum_depth", np.float32, ("latitude", "longitude"))
        nc_inten = cvar("sum_intensity", np.float32, ("latitude", "longitude"))

        # Assign values
        nc_lat[:] = lat
        nc_lon[:] = lon
        nc_n[:] = ncases
        nc_ninv[:] = ninversions
        nc_t[:] = tsrfsum
        nc_depth[:] = depthsum
        nc_inten[:] = intensum
