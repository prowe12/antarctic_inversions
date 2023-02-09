#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 19:13:05 2023

@author: prowe
"""

#!/usr/bin/env python
#
#
#
# Get ERA5 data over Antarctica at model levels
#
# Notes on calculating geopotential height, quoted from here:
# https://confluence.ecmwf.int/display/CKB/ERA5%3A+compute+pressure+and+
#    geopotential+on+model+levels%2C+geopotential+height+and+geometric+height
#

import os
import cdsapi
import pandas as pd


def get_datestr(date: str, mytime: str):
    """
    Get the date+time from the date and time
    @params date:  The date as YY-MM-DD
    @params mytime:  The time as hh:mm or hh:mm:ss
    @returns The date as YY_MM_DD_hh
    """
    # datestr = date[0][:4] + "_" + date[0][5:7] + "_" + mytime[:2]
    return date[:4] + "_" + date[5:7] + "_" + date[8:10] + "_" + mytime[:2]


def getfname(direc: str, date: str, mytime: str, paramstr: str, ext: str):
    """
    Get the filename from its parts
    @params direc  Top directory
    @params date  The date as YY-MM-DD
    @params mytime:  The time as hh:mm or hh:mm:ss
    @params paramstr:  The parameter from ERA5
    @params ext  The file extension; e.g. .nc or .grib
    @returns  Pathname
    """
    datestr = get_datestr(date, mytime)
    return direc + paramstr + datestr + ext


def write_bash(dates, times, gb_dir, nc_dir):
    """Write the bash script"""
    compgeo = "python3 compute_geopotential_on_ml.py "
    grib2nc = "grib_to_netcdf -o "  # "myfile.nc myfile.grib"
    fname = "get_geopotentialfiles_" + dates[0][:4] + ".sh"
    with open(fname, "w", encoding="utf-8") as fid:
        for date in dates:
            for mytime in times:
                # filenames
                tq_gb = getfname(gb_dir, date, mytime, "tq_ml_", ".grib")
                zlnsp_gb = getfname(gb_dir, date, mytime, "zlnsp_ml_", ".grib")
                geop_gb = getfname(gb_dir, date, mytime, "geop_", ".grib")

                tq_nc = getfname(nc_dir, date, mytime, "tq_ml_", ".nc")
                geop_nc = getfname(nc_dir, date, mytime, "geop_", ".nc")
                zlnsp_nc = getfname(nc_dir, date, mytime, "zlnsp_ml_", ".nc")

                # Write command to create geopotential height
                cmd = compgeo + tq_gb + " " + zlnsp_gb + " -o " + geop_gb
                fid.write(cmd + "\n")

                # Write command to convert tq grib to netcdf
                # tq_nc = tq_gb[: tq_gb.find(".grib")] + ".nc"
                cmd = grib2nc + tq_nc + " " + tq_gb
                fid.write(cmd + "\n")

                # Write command to convert geop grib to netcdf
                # geop_nc = geopgrib[: geop_gb.find(".grib")] + ".nc"
                cmd = grib2nc + geop_nc + " " + geop_gb
                fid.write(cmd + "\n")

                # Write command to convert zlnsp grib to netcdf
                cmd = grib2nc + zlnsp_nc + " " + zlnsp_gb
                fid.write(cmd + "\n")

        fid.write("echo done\n")


# def getfilenames(direc, date, mytime):
#     tqfile = getfilename(direc, date, mytime, "tq_ml_")
#     zlnspfile = getfilename(direc, date, mytime, "zlnsp_ml_")
#     return tqfile, zlnspfile


# def get_geop_fname(direc, date, mytime):
#     return getfilename(direc, date, mytime, "geop")


def download(dates: list[str], times: list[str], gribdir: str):
    """
    Download ERA5 data
    @params dates  Dates of interest (1 day each)
    @params times  Times of interest
    @params gribdir  Directory to save to
    """
    cds = cdsapi.Client()

    # Data download specifications:
    # date: Specify a single date as "2018-01-01" or a period as
    #   "2018-08-01/to/2018-01-31". For periods > 1 month see
    #   https://confluence.ecmwf.int/x/l7GqB
    # "type": Use "an" (analysis) unless you have a particular reason to use
    #    "fc" (forecast).
    # time: ERA5 data is hourly. Specify a single time as "00:00:00", or a
    #   range as "00:00:00/01:00:00/02:00:00" or "00:00:00/to/23:00:00/by/1".
    #   (note: the range did not work for conversion to geopotentials)
    clas = "ea"  # do not change
    expver = "1"  # do not change
    levtype = "ml"  # do not change
    stream = "oper"  # do not change
    typ = "an"

    # 137 model levels: "1/2/...137",
    # Latitude/longitude grid: east-west (longitude) and
    # north-south resolution (latitude). Default: 0.25 x 0.25
    # area:  # example: [60, -10, 50, 2], # North, West, South, East.
    # levtype: Geopotential (z) and Logarithm of surface pressure (lnsp) are
    # 2D fields, archived as model level 1
    area = [-60, -180, -90, 180]
    grid = [0.5, 0.5]
    levelist = (
        "1/2/3/4/5/6/7/8/9/10/11/12/13/14/15/16/17/18/19/20/21/22/23/24/25/26/"
        + "27/28/29/30/31/32/33/34/35/36/37/38/39/40/41/42/43/44/45/46/47/48/"
        + "49/50/51/52/53/54/55/56/57/58/59/60/61/62/63/64/65/66/67/68/69/70/"
        + "71/72/73/74/75/76/77/78/79/80/81/82/83/84/85/86/87/88/89/"
        + "90/91/92/93/94/95/96/97/98/99/100/101/102/103/104/105/106/107/108/"
        + "109/110/111/112/113/114/115/116/117/118/119/120/121/122/123/124/"
        + "125/126/127/128/129/130/131/132/133/134/135/136/137"
    )

    # "param": "130/133",  # Temperature (t) and specific humidity (q)
    # "param": "129/152",  # Geopotential (z) and Log surface pressure (lnsp)
    for date in dates:
        for mytime in times:

            # filenames
            tqfile = getfname(gribdir, date, mytime, "tq_ml_", ".grib")
            zlnspfile = getfname(gribdir, date, mytime, "zlnsp_ml_", ".grib")
            print(f"tqfile: {tqfile}")
            print(f"zlnspfile: {zlnspfile}")

            # remove files, if they exist
            if os.path.isfile(tqfile):
                os.remove(tqfile)
            if os.path.isfile(zlnspfile):
                os.remove(zlnspfile)

            cds.retrieve(
                "reanalysis-era5-complete",
                {
                    "class": clas,
                    "date": date,
                    "expver": expver,
                    "levelist": levelist,
                    "levtype": "ml",
                    "param": "130/133",
                    "stream": stream,
                    "time": mytime,
                    "type": typ,
                    "grid": grid,
                    "area": area,
                },
                tqfile,
            )

            cds.retrieve(
                "reanalysis-era5-complete",
                {
                    "class": clas,
                    "date": date,
                    "expver": expver,
                    "levelist": "1",
                    "levtype": levtype,
                    "param": "129/152",
                    "stream": stream,
                    "time": mytime,
                    "type": typ,
                    "grid": grid,
                    "area": area,
                },
                zlnspfile,
            )


def runner(year, topdir):
    """
    Download ERA5 and write the bash file
    @params year  Year of interest
    @params topdir  Main directory for saving results
    """
    yearstr = str(year)
    gribdir = topdir + "grib/" + yearstr + "/"
    ncdir = topdir + yearstr + "/"

    date1 = yearstr + "-06-01"
    date2 = yearstr + "-08-31"
    timestamps = pd.date_range(date1, date2).tolist()

    dates = [d.strftime("%Y-%m-%d") for d in timestamps]

    times = [
        "00:00:00",
        "12:00:00",
    ]

    download(dates, times, gribdir)
    write_bash(dates, times, gribdir, ncdir)


ERADIR = "era5/"
runner(2010, ERADIR)
runner(2011, ERADIR)
runner(2012, ERADIR)
runner(2013, ERADIR)
runner(2014, ERADIR)
runner(2015, ERADIR)
runner(2016, ERADIR)
runner(2017, ERADIR)
runner(2018, ERADIR)
runner(2019, ERADIR)
runner(2020, ERADIR)
runner(2021, ERADIR)
runner(2022, ERADIR)

# Run in terminal via:
# bash get_geopotentialfiles.sh
