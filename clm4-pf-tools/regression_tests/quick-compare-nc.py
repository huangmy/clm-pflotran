#!/usr/bin/env python
#
# Quick and dirty comparison of netcdf files by w/o requiring an
# external python librarary install like py-netCDF4 or nco by using
# ncdump to generate CDL files and comparing the text dumps.
#
# Requires: ncdump is available in the path!
#
# Author: Ben Andre, LBNL
#

from __future__ import print_function
from __future__ import division

import argparse
import os
import os.path
import re
import subprocess
import sys
import textwrap
import time
import traceback


def check_for_ncdump():
    """Check if ncdump is available in the path by trying to run it and
    checking the return code. Should be zero if present.

    """
    txtwrap = textwrap.TextWrapper(width=78, subsequent_indent=4*" ")
    timeout = 10.0 # sec
    command = ["ncdump"]
    try:
        with open(os.devnull, 'w') as devnull_f:
            subprocess.check_call(command, stdout=devnull_f, stderr=devnull_f)
    except Exception as e:
        raise RuntimeError(
            txtwrap.fill("ERROR: check_for_ncdump() could not find ncdump!. "
                         "ncdump must be available in the path."))


def read_data(netcdf_file):
    """Read the data. 

    Check if the cdl file exists for the netcdf file. If not, then we
    create it. If so, then we check if the netcdf file exists. If the
    netcdf file is newer than the cdl file, we recreate the cdl file.

    """
    print("Reading {0}...".format(netcdf_file), end='')

    # check which files exist
    have_nc = False
    if os.path.isfile(netcdf_file):
        have_nc = True

    cdl_file = "{0}.cdl".format(netcdf_file)
    have_cdl = False
    if os.path.isfile(cdl_file):
        have_cdl = True
    
    if have_nc and have_cdl:
        if os.path.getmtime(netcdf_file) > os.path.getmtime(cdl_file):
            # both netcdf and cdl files exist, and netcdf is newer, regenerate cdl
            print("--> regenerating cdl from netcdf.")
            cdl = read_nc_as_cdl(netcdf_file)
        else:
            print("--> reusing cdl file.")
            cdl = read_cdl(cdl_file)
    elif have_nc and not have_cdl:
        print("--> generating cdl file.")
        cdl = read_nc_as_cdl(netcdf_file)
    elif not have_nc and have_cdl:
        print("--> using existing cdl file.")
        cdl = read_cdl(cdl_file)
    elif not have_nc and not have_cdl:
        raise RuntimeError("\nERROR: neither netcdf or cdl dump exist. Expected : '{0}'".format(netcdf_file))

    print(" done.")
    return cdl

def read_cdl(cdl_file):
    """Read in the contents of a cdl file as an array of lines.

    """
    with open(cdl_file, "r") as tmp_file:
        cdl = tmp_file.readlines()
    return cdl

def read_nc_as_cdl(netcdf_file):
    """Read a netcdf file as cdl by dumping it with ncdump into a temporary file.

    """
    cdl = None
    cdl_file = "{0}.cdl".format(netcdf_file)
    with open(cdl_file, "w") as tmp_file:
        command = ["ncdump", netcdf_file]
        subprocess.call(command, stdout=tmp_file, stderr=subprocess.PIPE)
    cdl = read_cdl(cdl_file)
    return cdl

def extract_field_from_cdl(cdl, field_name):
    """Search through the CDL and extract the requested field.

    """
    name = field_name.strip()
    re_str = "[\s]+{0}[\s]+=[\s]+$".format(name)
    field_re = re.compile(re_str, re.IGNORECASE)
    print("Searching for '{0}'...".format(name), end='')
    data = []
    start_field = False
    end_field = False
    for line in cdl:
        if start_field:
            data.append(line.strip())
            if data[-1].endswith(';'):
                end_field = True
                break
        if field_re.search(line):
            start_field = True

    if not start_field:
        raise RuntimeError("\nERROR: could not find field '{0}'!".format(name))
    print(" done.")
    return ''.join(data)

def compare_fields(baseline, current, tolerance_str):
    """Compare the baseline and current fields to the specified tolerance

    """
    tolerance_info = tolerance_str.split(' ')
    tolerance = float(tolerance_info[0].strip())
    tolerance_type = tolerance_info[1].strip()
    status = 0
    base = baseline.split(',')
    base[-1] = base[-1].strip(';')
    cur = current.split(',')
    cur[-1] = cur[-1].strip(';')
    for b, c in zip(base, cur):
        b = float(b.strip())
        c = float(c.strip())
        status += check_tolerance(b, c, tolerance, tolerance_type)
        
    return status

def check_tolerance(previous, current, tolerance, tolerance_type):
    """
    """
    status = 0
    if tolerance_type == "absolute":
        delta = abs(previous - current)
    elif (tolerance_type == "relative" or
          tolerance_type == "percent"):
        if previous != 0:
            delta = abs(previous - current) / previous
        elif current != 0:
            delta = abs(previous - current) / current
        else:
            # both are zero
            delta = 0.0
        if tolerance_type == "percent":
            delta *= 100.0
    else:
        # should never get here....
        raise Exception("ERROR: unknown test tolerance_type '{0}'.".format(tolerance_type))

    if delta > tolerance:
        status = 1
        print("    FAIL: {0} > {1} [{2}]".format(
            delta, tolerance, tolerance_type))
    elif False:
        print("    PASS: {0} <= {1} [{2}]".format(
            delta, tolerance, tolerance_type))

    return status
    

def commandline_options():
    parser = argparse.ArgumentParser(
        description="quick and dirty comparison of netcdf files by dumping to "
        "text CDL files and diffing specified fields to the specified tolerance.")

    parser.add_argument("--backtrace", action='store_true', default=False,
                        help='show exception backtraces as extra debugging '
                        'output')

    parser.add_argument("-b", "--baseline-netcdf-file",
                        nargs=1, default=None, required=True,
                        help="path to the netcdf file")

    parser.add_argument("-c", "--current-netcdf-file",
                        nargs=1, default=None, required=True,
                        help="path to the netcdf file")

    parser.add_argument('-d', "--debug", action='store_true', default=False,
                        help='extra debugging output')

    parser.add_argument("-f", "--field-name",
                        nargs=1, default=None, required=True,
                        help="field name in the netcdf file")

    parser.add_argument("-t", "--tolerance",
                        nargs='+', default=None, required=True,
                        help="tolerance and type, e.g. '1.0e-7 absolute'")

    options = parser.parse_args()
    return options



def main(options):
    status = 0
    print("Quick compare '{0}' between :".format(options.field_name[0]))
    print("  baseline : {0}".format(options.baseline_netcdf_file[0]))
    print("  current :  {0}".format(options.current_netcdf_file[0]))

    check_for_ncdump()

    baseline = read_data(options.baseline_netcdf_file[0])
    baseline_field = extract_field_from_cdl(baseline, options.field_name[0])

    current = read_data(options.current_netcdf_file[0])
    current_field = extract_field_from_cdl(current, options.field_name[0])

    status = compare_fields(baseline_field, current_field, options.tolerance[0])
    if status == 0:
        print("PASS: quick check.")
    else:
        print("FAIL: quick check")
    return status

if __name__ == "__main__":
    options = commandline_options()
    try:
        status = main(options)
        sys.exit(status)
    except Exception as e:
        print(str(e))
        if options.backtrace:
            traceback.print_exc()
        sys.exit(1)

