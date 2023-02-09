#!/usr/bin/env python
#
# Get altitudes corresponding ERA5 model levels
#
# Notes on calculating geopotential height, quoted from here:
# https://confluence.ecmwf.int/display/CKB/ERA5%3A+compute+pressure+and+
#    geopotential+on+model+levels%2C+geopotential+height+and+geometric+height
#
# Geopotential height
#
# In ERA5, and often in meteorology, heights (the height of the land and
# sea surface, or specific heights in the atmosphere) are not represented
# as geometric height, or altitude (in metres above the spheroid), but as
# geopotential height (in metres above the geoid, which is represented by
# the mean sea level in ERA5). Note, that ECMWF usually archive the
# geopotential (in m2/s2), not the geopotential height.
#
# To obtain the geopotential height (h) in metres (of the land and
# sea surface or at particular heights in the atmosphere),
# divide the geopotential by the Earth's gravitational acceleration,
# 9.80665 m/s2 in the IFS. This geopotential height is
# relative to the geoid over land and the mean sea level over ocean - for
# more information see ERA5: data documentation - spatial reference systems.
#
# Geometric height
#
# The geometric height or altitude (alt) is given by:
# alt = Re * h / (Re − h)
#
# where Re is the radius of the Earth. This geometric height is relative
# to the geoid over land and the mean sea level over ocean and it is
# assumed that the Earth is a perfect sphere - for more information
# see ERA5: data documentation - spatial reference systems.
#
# From https://confluence.ecmwf.int/display/CKB/ERA5%3A+data+documentation
#      #ERA5:datadocumentation-Spatialreferencesystems
# Earth model
#
# For data in GRIB1 format the earth model is a sphere with
# radius = 6367.47 km, as defined in the WMO GRIB Edition 1
# specifications, Table 7, GDS Octet 17.
#
# For data in GRIB2 format the earth model is a sphere with
# radius = 6371.2290 km, as defined in the WMO GRIB2 specifications,
# section 2.2.1, Code Table 3.2, Code figure 6.
#
# For data in NetCDF format (i.e. converted from the native GRIB format
# to NetCDF), the earth model is inherited from the GRIB data.

from netCDF4 import Dataset
import matplotlib.pyplot as plt


def getalt(fname):
    with Dataset(fname, "r") as nc_h:

        # alt = Re * h / (Re − h)
        # Is this grib 1 or grib 2?
        # GRIB1: radius = 6367.47 km, as defined in the WMO GRIB Edition 1
        # specifications, Table 7, GDS Octet 17.
        #
        # GRIB2: format the earth model is a sphere with
        # radius = 6371.2290 km, as defined in the WMO GRIB2 specifications,
        geop = nc_h["z"][:].data / 9.80665 / 1000
        alt = 6367.47 * geop / (6367.47 - geop)
        alt[geop == 0] = 0
        return alt, geop


# def get_tq():
#     fname = tq_ml_20150801_00.nc
#     with Dataset(fname, "r") as ncid:
#         ncid['t']


def example():
    """
    Print Altitude, geopotential, and temperature for an example case
    """
    eradir = "era5/"
    ncdir = eradir + "2018/"
    fname = ncdir + "tq_ml_2018_08_01_00.nc"

    # Get the altitude and geopotential
    alt, geop = getalt(ncdir + "geop_2018_08_01_00.nc")

    with Dataset(fname, "r") as nci:
        # ncid.variables.keys()
        # dict_keys(['longitude', 'latitude', 'level', 'time', 't', 'q'])
        # ncid['t']: nt16 t(time, level, latitude, longitude)

        # Example lats to plot
        firstlat = str(nci["latitude"][0].data)
        lastlat = str(nci["latitude"][-1].data)

        # Plot geopotential and altitude
        plt.figure(1)
        plt.clf()
        plt.plot(geop[0, :, 0, 0], ".", label=firstlat + ": geopotential")
        plt.plot(alt[0, :, 0, 0], "-", label=firstlat + " altitude")
        plt.plot(geop[0, :, -1, 0], ".", label=lastlat + ": geopotential")
        plt.plot(alt[0, :, -1, 0], "-", label=lastlat + " altitude")
        plt.legend()
        plt.xlabel("index")
        plt.ylabel("Height (km)")

        # Plot temperature with altitude
        plt.figure(2)
        plt.clf()
        plt.plot(nci["t"][0, :, 0, 0].data, alt[0, :, 0, 0], label=firstlat)
        plt.plot(nci["t"][0, :, -1, 0].data, alt[0, :, -1, 0], label=lastlat)
        plt.xlabel("Temperature (K)")
        plt.ylabel("Height (km)")
        plt.legend()

    # ncid.variables.keys()
    # dict_keys(['longitude', 'latitude', 'level', 'time', 'z'])

    # ncid["z"]
    # <class 'netCDF4._netCDF4.Variable'>
    # int16 z(time, level, latitude, longitude)
    #     scale_factor: 11.919896078046438
    #     add_offset: 388972.0777106524
    #     _FillValue: -32767
    #     missing_value: -32767
    #     units: m**2 s**-2
    #     long_name: Geopotential
    #     standard_name: geopotential
    # unlimited dimensions:
    # current shape = (1, 137, 31, 360)
    # filling on


# if __name__ == "__main__":
#     example()
