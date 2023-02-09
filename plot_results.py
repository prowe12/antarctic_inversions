#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 09:49:29 2023

@author: prowe
"""

import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset


def plot_results(lat, lon, ncases, tsurfsum, ninversions, depthsum, intensum):
    """
    Plot results
    @params lat  The latitude
    @params lon  The longitude
    @param ncases  The number of cases
    @params tsurfsum  The sum of the surface temperatures
    @params ninversions  Number of inversions
    @params depthsum  = Sum of inversion depths
    @params intensum  = Sum of inversion intensities
    """
    frequency = ninversions / ncases
    depth = depthsum / ninversions
    strength = intensum / ninversions
    tsurf = tsurfsum / ncases - 273.15

    plt.figure(3)
    plt.clf()
    plt.subplot(3, 1, 1)
    plt.plot(lat, frequency, ".-")
    plt.ylabel("Frequency")
    plt.subplot(3, 1, 2)
    plt.plot(lat, depth, ".-")
    plt.ylabel("Mean inversion depth (km)")
    plt.subplot(3, 1, 3)
    plt.plot(lat, strength, ".-")
    plt.ylabel("Mean inversion strength (C)")

    x, y = np.meshgrid(lon, lat)

    plt.figure(4)
    plt.clf()
    plt.scatter(x, y, c=tsurf, marker="s")
    plt.colorbar()

    plt.figure(5)
    plt.clf()
    plt.subplot(311)
    plt.scatter(x, y, c=frequency, marker="s", cmap="jet")
    plt.colorbar()
    plt.subplot(312)
    plt.scatter(x, y, c=depth * 10, marker="s", cmap="jet")
    plt.colorbar()
    plt.subplot(313)
    plt.scatter(x, y, c=strength, marker="s", cmap="jet")
    plt.colorbar()


def plot_for_year(year, ncfile):
    ncid = Dataset(ncfile)
    plot_results(
        ncid["lat"][:],
        ncid["lon"][:],
        ncid["ncases"][:],
        ncid["sum_tsurface"][:],
        ncid["ninversions"][:],
        ncid["sum_depth"][:],
        ncid["sum_intensity"][:],
    )


# Selected year
year = 2018

# Directory and output file
ncdir = "era5/" + str(year) + "/"
ncfile = "era5/inversion_stats_" + str(year) + ".nc"

plot_for_year(year, ncfile)
