#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 13:46:11 2023

@author: prowe

From Zhang et al 2011

"For all four datasets, SBIs are identified from temperature profile data 
using Kahl’s (1990) algorithm, which scans from the surface upward to 
500 hPa to find cases with temperature increasing with altitude. The 
inversion top is defined as the bottom of the first layer in which 
temperature decreases with altitude, but thin (<100 m) noninversion 
layers, with temperature decreasing with height, are ignored if they are 
embedded within a deeper inversion layer. Soundings are considered to be 
unsuitable for analysis if the surface level is missing, there are fewer 
than 10 upper-air data levels from surface to 500 hPa, or the temperature 
difference across the inversion exceeds 40 K (Serreze et al. 1992). This 
definition represents a true inversion layer (from the surface to the 
inversion top), which is different from defining inversion strength as a 
temperature difference between two prespecified levels or heights in the 
boundary layer (Hudson and Brandt 2005; Kay and Gettelman 2009; 
Pavelsky et al. 2011)."

"For each grid point or radiosonde station, we compute frequency of SBI 
occurrence ( f ) and average SBI depth and intensity for the overall dataset 
and separately for each season, for each month, and for 0000 and 1200 UTC
(61 h to account for the 1–2-h sounding duration). Following Kahl (1990) 
and Serreze et al. (1992), three- month seasons are defined as January–March 
(Arctic winter, Antarctic summer), July–September 
(Arctic summer and Antarctic winter), etc.

"In all comparisons of the various datasets, the models and reanalysis 
are sampled at 0000 and 1200 UTC to mimic the radiosondes. Comparisons of 
radiosonde data and the two climate models are for 1990–2007. Comparisons 
for radiosonde data and ERA-Interim, and analysis of radiosonde data alone, 
are for 1990–2009."

"""

from netCDF4 import Dataset  # , num2date
import numpy as np
import matplotlib.pyplot as plt
import os

from get_altitudes import getalt
from create_inversion_file import create_inversion_file

# from plot_results import plot_results


def plot_prof(figno, height, temp, label=""):
    plt.figure(figno)
    plt.clf()
    plt.plot(temp, height, ".-", label=label)


def plot_prof_inv(height, temp, itop, iinv, figno=3, label=""):
    plot_prof(figno, height, temp)
    plt.plot(temp[itop], height[itop], "r*", markersize=20)
    plt.plot(temp[iinv], height[iinv], "go")
    plt.ylim([0, 10])
    plt.ylabel("Height (km)")
    plt.xlabel("Temperature (K)")


def inversion_exists(temp):
    # ID if surface-based temperature inversion is present
    # True if temperature decreases between the first two
    # temperatures *or* the second two temperatures
    # the latter condition is added because the interpolation
    # to the surface temperature appears to be imperfect
    if temp[1] - temp[0] >= 0 or temp[2] - temp[1] >= 0:
        return True
    else:
        return False


def get_params(year):
    # Pressure levels
    press_levs = np.hstack(
        [
            np.arange(1000, 925 - 5, -25),
            np.arange(900, 450, -50),
        ]
    )
    press_levs_ints = list(press_levs)
    press_levs = [str(p) for p in press_levs_ints]

    north = -88.0
    south = -90.0
    west = 121.0  # -59.0
    east = 122.0  # -58.5
    return {
        "format": "netcdf",
        "product_type": "reanalysis",
        "variable": [
            "pressure",
            "Geopotential",
            "Temperature",
        ],
        "pressure_level": press_levs,
        "year": [year],
        "month": ["07"],
        "day": ["01", "02", "03"],
        "time": ["00:00", "12:00"],
        "grid": [0.5, 0.5],
        "area": [north, west, south, east],
    }


def get_stats(temp, height, plotfig: int = 0):
    """
    Get the inversion statistics
    @param temp  Temperature profile, surface to TOA
    @param height  Altitudes, surface to TOA
    @returns  Inversion depth,
    @returns  Inversion strength
    """
    dtemp = np.diff(temp)

    # Ignore a non-inversion in the surface layer
    # and consider no change in temperature to be an inversion
    dtemp[0] = 1e-4
    dtemp[dtemp == 0] = 1e-4

    # iinv = np.where(dtemp >= 0)[0] + 1
    inoninv = np.where(dtemp <= 0)[0] + 1
    switch = np.diff(np.sign(dtemp))
    iswitch = np.where(switch != 0)[0] + 1
    # Mark the maximum height, if still in inversion
    if dtemp[-1] > 0:
        iswitch = np.hstack([iswitch, len(temp) - 1])

    # Ignore thin (<100 m) noninversion layers, with temperature decreasing
    # with height, if they are embedded within a deeper inversion layer.
    if len(inoninv) == 0:
        # This means the inversion goes all the way to the top of the
        # model atmosphere. Only ok if it was chopped at ~500 mb!
        itop = len(height) - 1
        # elif len(inoninv) == 1:
        #     # Only one inversion (simple case)
        #     itop = inoninv[0] - 1
    elif len(inoninv) >= 1:
        # Multiple inversions, so check if they are thin (< 100m )
        # non-inversions, and, if so, ignore them up until we get to the
        # first non-inversion, moving up from the surface, that spans
        # more than 100 m.

        if dtemp[iswitch[0]] >= 0:
            raise ValueError("Unexpected")

        # Progress through the noninversion depths and include them while
        # they are less than 100 m
        itop = iswitch[-1]
        for lev in range(1, len(iswitch), 2):
            # dtemp[iswitch[lev - 1]] is the beginning of the non-inversion
            # dtemp[iswitch[lev] is the end of the non-inversion
            # temp[iswitch[lev - 1]] is the last temp that increases
            #    (beyond this is the noninversion)
            # temp[iswitch[lev - 1]+1] is the first temp that decreases
            #    (beginning of non-inversion)
            # temp[iswitch[lev]] is the last temp that decreases
            #    (end of non-inversion)
            if dtemp[iswitch[lev - 1]] > 0 or dtemp[iswitch[lev]] < 0:
                raise ValueError("Unexpected behavior with temperature")
            noninvdepth = height[iswitch[lev]] - height[iswitch[lev - 1]]
            # Quit if there is a noninversion deeper than 100 m,
            if noninvdepth > 0.1:
                # We are out of the inversion, so quit
                itop = iswitch[lev - 1]  # last temp that increases
                break

    # Make sure this is the maximum temperature within the section of
    # troposphere containing the inversion. If a temperature lower down
    # in the troposphere is higher, use it instead (but do not go above)
    # This avoids negative inversion strengths that can occur due to
    # small inversions occurring aloft within a large noninversion layer
    if np.any(temp[itop] < temp[iswitch[iswitch < itop]]):
        imax = temp[iswitch[iswitch < itop]].argmax()
        itop = iswitch[imax]

    # TOOD: deal with this
    # if height[itop] > 8:
    #     print("Inversion goes to top of model atmosphere")

    invdepth = height[itop] - height[0]
    invstrength = temp[itop] - temp[0]
    if invdepth <= 0:
        raise ValueError("Inversion depth should not be <= 0")
    if invstrength < 0:
        # TODO: deal with these cases
        pass
        # raise ValueError("Inversion strength should not be < 0")

    # if plotfig:
    #     plot_prof_inv(height, temp, itop, iinv, figno=10, label="")

    return invdepth, invstrength


def interp_to_zsrf(zin, tin, pin, zsrf):
    """ """
    tsrf = np.interp(zsrf, zin, tin)
    psrf = np.interp(zsrf, zin, pin)
    iabove = np.where(zin > zsrf)[0]
    znew = np.hstack([zsrf, zin[iabove]])
    tnew = np.hstack([tsrf, tin[iabove]])
    pnew = np.hstack([psrf, pin[iabove]])
    return znew, tnew, pnew


def get_inversion_stats(ncdir, outfile):
    """
    Get the inversion statitistics and save to a file
    @param ncdir  Directory with netcdf files to use
    @param outfile  File to save results to
    """

    allfiles = np.sort(os.listdir(ncdir))
    tqfiles = []
    geopfiles = []
    for file in allfiles:
        if file.startswith("tq_ml_"):
            tqfiles.append(file)
        elif file.startswith("geop_"):
            geopfiles.append(file)

    tqfile = tqfiles[0]
    geopfile = geopfiles[0]
    # filenames = ["test_20180101_00_tuvz.nc"]

    # Get the altitudes
    alt, _ = getalt(ncdir + geopfile)

    # Get info from the first file that doesn't change
    with Dataset(ncdir + tqfile, "r") as ncid:
        # ncid.variables.keys()
        # dict_keys(['longitude', 'latitude', 'level', 'time', 't', 'q'])
        # ncid["t"]: int16 t(time, level, latitude, longitude)

        # Data for output
        nlats = ncid.dimensions["latitude"].size
        nlons = ncid.dimensions["longitude"].size
        lats = ncid["latitude"][:].data
        lons = ncid["longitude"][:].data

        # plt.figure(1)
        # nalts = ncid.dimensions["level"].size
        # plt.clf()
        # plt.plot(ncid["t"][0, :, 0, 0], alt[0, :nalts, 0, 0])
        # plt.plot(ncid["t"][0, :, -1, 0], alt[0, :nalts, -1, 0])

    # Preallocate variables
    depthsum = np.zeros([nlats, nlons])
    intensum = np.zeros([nlats, nlons])
    tsurfsum = np.zeros([nlats, nlons])
    ninversions = np.zeros([nlats, nlons])
    ncases = np.zeros([nlats, nlons])

    # figno = 2
    # if plotfig:
    #     plt.figure(figno)
    #     plt.clf()

    for geopfile, tqfile in zip(geopfiles, tqfiles):
        print(geopfile)

        with Dataset(ncdir + tqfile, "r") as ncid:
            # Altitudes and sum of surface temperatures
            alt, _ = getalt(ncdir + geopfile)
            tsurfsum += ncid["t"][0, -1, :, :]
            for ilat in range(nlats):
                # TODO: Fix next line: chance to range(nlats)
                for ilon in range(0, nlons, 100):
                    temp = np.flipud(ncid["t"][0, :, ilat, ilon].data)
                    alt0 = alt[0, :, ilat, ilon]
                    isinv = inversion_exists(temp)
                    ncases[ilat, ilon] += 1
                    if isinv:
                        invdepth, invstrength = get_stats(
                            temp[alt0 <= 7], alt0[alt0 <= 7]
                        )
                        ninversions[ilat, ilon] += 1
                        depthsum[ilat, ilon] += invdepth
                        intensum[ilat, ilon] += invstrength

    # Save results to netcdf file
    create_inversion_file(
        outfile,
        lats,
        lons,
        ncases,
        tsurfsum,
        ninversions,
        depthsum,
        intensum,
    )

    # plot_results(lats, lons, ncases, tsurfsum, ninversions, depthsum, intensum)


if __name__ == "__main__":

    # Selected year
    year = 2018

    # Directory and output file
    ncdir = "era5/" + str(year) + "/"
    outfile = "era5/inversion_stats_" + str(year) + ".nc"

    # Run
    get_inversion_stats(ncdir, outfile)
