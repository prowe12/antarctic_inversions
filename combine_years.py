#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 14:18:14 2023

@author: prowe
"""

import os
import datetime as dt
import re
import numpy as np
from netCDF4 import Dataset  # type: ignore


def get_regex_from_fileformat(fmt: str) -> str:
    """
    Get string to use with regex from the file format string
    @params fileformat  The file format
    @returns  The string for regex
    """
    # Create the regex from the fileformat
    regstr = ""
    skipnext = False
    for char in fmt:
        if skipnext:
            skipnext = False
            continue
        if char == "%":
            regstr += ".+"
            skipnext = True
        else:
            regstr += char
    return regstr


def getfiles(direc: str, fmt: str) -> list[str]:
    """
    Get the files for combining
    @param direct  The directory with the files
    @param fmt  The format of the files
    @returns   A list of the files
    """

    regstr = get_regex_from_fileformat(fmt)

    fnames = []
    dates = []
    for file in os.listdir(direc):
        if re.match(regstr, file):
            dates.append(dt.datetime.strptime(file, fmt))
            fnames.append(file)
    fnames = list(np.sort(fnames))
    return fnames


def combine_years(ncdir: str, fmt: str):
    """
    Comine the temperature inversion data over all the years and
    save results to a single netcdf file
    @params ncdir  Directory with files
    @params fileformat  File format
    @returns lat
    @returns lon
    @returns frequency
    @returns depth
    @returns intensity
    @returns tsurface
    Examples:
        ncdir = 'era5'
        fileformat = 'inversion_stats_%Y.nc'
    """

    # Get filenames to combine
    fnames = getfiles(ncdir, fmt)

    # Combine files

    # Get data from first file
    with Dataset(ncdir + fnames[0], "r") as ncid:
        lat = ncid["lat"][:]
        lon = ncid["lon"][:]
        ncases = ncid["ncases"][:]
        sum_tsurface = ncid["sum_tsurface"][:]
        ninversions = ncid["ninversions"][:]
        sum_depth = ncid["sum_depth"][:]
        sum_intensity = ncid["sum_intensity"][:]

    # Sum with contents of remaining files
    for fname in fnames[1:]:
        with Dataset(ncdir + fname, "r") as ncid:
            ncases += ncid["ncases"][:]
            sum_tsurface += ncid["sum_tsurface"][:]
            ninversions += ncid["ninversions"][:]
            sum_depth += ncid["sum_depth"][:]
            sum_intensity += ncid["sum_intensity"][:]

    # Compute final stats
    frequency = 100 * ninversions / ncases
    tsurface = sum_tsurface / ncases
    depth = sum_depth / ninversions
    intensity = sum_intensity / ninversions

    return lat, lon, frequency, depth, intensity, tsurface


# Example of runninng code
#     from plot_results import plot_final_results
#     from create_inversion_final_file import create_inversion_final_file

#     # Directory and file format
#     eradir = "era5/"
#     fileformat = "inversion_stats_%Y.nc"
#     outfile = eradir + "inversion_stats.nc"

#     # Calculate the means of frequency, depth, intensity, and
#     # surface temperature for all years from the sums, the
#     # numbers of cases and the numbers of inversions
#     lat, lon, freq, depth, intensity, tsurf = combine_years(eradir, fileformat)

#     # Save the result
#     create_inversion_final_file(
#         outfile, lat, lon, freq, depth, intensity, tsurf
#     )

#     # Plot it
#     plot_final_results(lat, lon, freq, depth, intensity, tsurf)
