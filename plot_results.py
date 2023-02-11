#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 09:49:29 2023

@author: prowe
"""

import matplotlib.pyplot as plt  # type: ignore
import numpy as np
import numpy.typing as npt
from netCDF4 import Dataset  # type: ignore


def plot_final_results(
    lat: npt.NDArray[np.float64],
    lon: npt.NDArray[np.float64],
    frequency: npt.NDArray[np.float64],
    depth: npt.NDArray[np.float64],
    strength: npt.NDArray[np.float64],
    tsurf: npt.NDArray[np.float64],
):
    """
    Plot results for averages over multiple years
    @params lat  The latitude
    @params lon  The longitude
    @param frequency  The mean frequency of occurence of inversions (%)
    @params depth  The mean inversion depth
    @params strength  = The mean inversion intensity
    @params tsurf  = The mean surface temperature
    Notes: freq, depth, strength, and tsurf are all functions of lat and lon
    """

    # plt.figure(3)
    # plt.clf()
    # plt.subplot(3, 1, 1)
    # plt.plot(lat, frequency, ".-")
    # plt.ylabel("Frequency")
    # plt.subplot(3, 1, 2)
    # plt.plot(lat, depth, ".-")
    # plt.ylabel("Mean inversion depth (km)")
    # plt.subplot(3, 1, 3)
    # plt.plot(lat, strength, ".-")
    # plt.ylabel("Mean inversion strength (C)")

    mlon, mlat = np.meshgrid(lon, lat)

    plt.figure(4)
    plt.clf()
    plt.scatter(mlon, mlat, c=tsurf - 273.15, marker="s")
    plt.axis([-180, 180, -90, -60])
    plt.colorbar(label="Temperature ($^{o}$C)")
    plt.ylabel("Latitude")
    plt.xlabel("Longitude")
    plt.title("Mean surface temperature")

    plt.figure(5)
    plt.clf()
    plt.subplot(311)
    plt.scatter(mlon, mlat, c=frequency, marker="s")
    plt.colorbar(label="Frequency (%)")
    plt.ylabel("Latitude")
    plt.axis([-180, 180, -90, -60])

    plt.subplot(312)
    plt.scatter(mlon, mlat, c=depth * 1000, marker="s")
    plt.colorbar(label="Depth (m)")
    plt.ylabel("Latitude")
    plt.axis([-180, 180, -90, -60])

    plt.subplot(313)
    plt.scatter(mlon, mlat, c=strength, marker="s")
    plt.colorbar(label="Intensity (K)")
    plt.axis([-180, 180, -90, -60])
    plt.ylabel("Latitude")
    plt.xlabel("Longitude")


def plot_results(
    lat: npt.NDArray[np.float64],
    lon: npt.NDArray[np.float64],
    ncases: npt.NDArray[np.float64],
    tsurfsum: npt.NDArray[np.float64],
    ninversions: npt.NDArray[np.float64],
    depthsum: npt.NDArray[np.float64],
    intensum: npt.NDArray[np.float64],
):

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

    plot_final_results(lat, lon, frequency, depth, strength, tsurf)


def plot_for_year(ncfile: str):
    """
    Plot results for a single year (by computing means from sums)
    @params ncfile  The file containing the results for the year
    """
    with Dataset(ncfile) as ncid:
        plot_results(
            ncid["lat"][:],
            ncid["lon"][:],
            ncid["ncases"][:],
            ncid["sum_tsurface"][:],
            ncid["ninversions"][:],
            ncid["sum_depth"][:],
            ncid["sum_intensity"][:],
        )


# Example of use:
#
# if __name__ == "__main__":
#     # Selected year
#     year = 2018

#     # Directory and output file
#     ncdir = "era5/" + str(year) + "/"
#     ncfile = "era5/inversion_stats_" + str(year) + ".nc"

#     plot_for_year(ncfile)
