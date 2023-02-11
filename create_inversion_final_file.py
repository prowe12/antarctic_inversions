#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 17:19:51 2023

@author: prowe
"""

from netCDF4 import Dataset  # type: ignore
import numpy as np
import numpy.typing as npt


def create_inversion_final_file(
    fname: str,
    lat: npt.NDArray[np.float64],
    lon: npt.NDArray[np.float64],
    freq: npt.NDArray[np.float64],
    depth: npt.NDArray[np.float64],
    intens: npt.NDArray[np.float64],
    tsurf: npt.NDArray[np.float64],
):
    """
    Create the output netcdf file for a set of inversion variables
    @params fname  The name of the file to be created
    @params lat  The latitude
    @params lon  The longitude
    @param freq  The mean frequency of occurence of inversions (%)
    @params depth  The mean inversion depth
    @params intens  = The mean inversion intensity
    @params tsurf  = The mean surface temperature
    Notes: freq, depth, intens, and tsurf are all functions of lat and lon
    """
    # Use netcdf4 to save the results
    with Dataset(fname, "w", format="NETCDF4_CLASSIC") as ncid:

        cvar = ncid.createVariable

        #  Create Dimensions
        ncid.createDimension("latitude", len(lat))
        ncid.createDimension("longitude", len(lon))

        # Create variables
        nc_lat = cvar("lat", np.float32, ("latitude"))
        nc_lon = cvar("lon", np.float32, ("longitude"))
        nc_freq = cvar("frequency", np.float32, ("latitude", "longitude"))
        nc_depth = cvar("depth", np.float32, ("latitude", "longitude"))
        nc_inten = cvar("intensity", np.float32, ("latitude", "longitude"))
        nc_t = cvar("tsurface", np.float32, ("latitude", "longitude"))

        # Assign values
        nc_lat[:] = lat
        nc_lon[:] = lon
        nc_freq[:] = freq
        nc_depth[:] = depth
        nc_inten[:] = intens
        nc_t[:] = tsurf
