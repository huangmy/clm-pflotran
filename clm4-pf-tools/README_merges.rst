This document contains information releated to possible merge and
update issues when updating clm-pflotran to new pflotran-dev or clm
trunk tags.

clm4_5_10 to clm4_5_35
=======================

data: update surface data sets for point and regional domains (clm4_5_11)
-------------------------------------------------------------------------

In clm4_5_11 tag, surface-dataset netcdf was modified by NCAR and the
following snippet of ChangeLog outlines the modification:

    A very related change is the separation of PCT_PFT in the surface
    dataset into PCT_NAT_PFT and PCT_CFT; in addition to these two
    variables, there are also new PCT_NATVEG (% of natural veg landunit on
    the gridcell) and PCT_CROP (% of crop landunit on the gridcell)
    variables. Note that the separation of PCT_PFT into natural vs crop
    was only done on the surface dataset -- raw datasets to mksurfdata_map
    have not been changed, nor have most of the CLM data structures.

* MATLAB scripts were used to create surface-datasets compatible with
clm4_5_11 onwards.

* Changes have been made in the ngee example data repo:

    * lnd/clm2/surfdata_map/surfdata_1x1pt_US-Brw_simyr1850_c130607.nc - Surface dataset compatiable with tags upto clm4_5_10

    * lnd/clm2/surfdata_map/surfdata_1x1pt_US-Brw_simyr1850_ugrid_c130828.nc - Surface dataset (in unstructured grid format) compatiable with clm4_5_11 onwards.

    * lnd/clm2/surfdata_map/surfdata_1x1pt_US-Brw_simyr1850.nc - A static link to surfdata_1x1pt_US-Brw_simyr1850_ugrid_c130828.nc dataset

    * lnd/clm2/surfdata_map/surfdata_13x26pt_US-Brw_simyr1850_Site_B_8.00000dx_8.00000dy_prism_ugrid_c130719.nc - A regional surface dataset for NGEE Site B that is compatiable with clm4_5_11 onwards.

    * share/domains/domain.clm/domain.lnd.13x26pt_US-Brw_simyr1850_Site_B_8.00000dx_8.00000dy_prism_ugrid_c130719.nc - Domain file for NGEE Site B

src: associate refactor (clm4_5_14-15)
--------------------------------------

The declaration of pointers for dereferencing derived type variables
has been eliminated and replaced with fortran 2003 associated blocks.

* NCAR has provided a tool for updating code to the new style.
    * See: models/lnd/clm/tools/clm4_5/refactorTools/associate/README

    * This tool is intended to work on an entire file that has not been merged to clm4_5_15.

    * If you have a stand alone function that needs to be merged, you
      can temporarily copy it to a new file, run the refactor tool,
      and copy it back to the original file.

    * If you have code that needs to be updated in a file that will be
      affected by the merge, you can run the tool prior to merging.

    * If you have code in a file that has already been merged, you’ll have to update by hand.

data: pft data merged with clm_params (clm4_5_23)
-------------------------------------------------

netcdf pft data has been combined with the clm_params data file.

* any changes to pft data must be saved in the new file: lnd/clm2/pftdata/pft-physiology.cYYMMDD.nc ---> lnd/clm2/paramdata/clm_params.cYYMMDD.nc

* Changes have been made in the ngee example data repo:

    * pft data has been moved from: pft-physiology.c130715arctic.nc to: clm_params.c130913.arctic.nc

data: update datm streams forcing data path (clm4_5_25)
-------------------------------------------------------

* datm has changed to get PTCLM working again. The standard use case
  for datm is that the datm domain and land domain are always the
  same, and this has been embedded further into the datm scripts.

    * To get around this, you need to initialize your case by calling
      cesm_setup, then:

      cp CaseDocs/datm.streams.txt.CLM1PT.CLM_USRDAT user_datm.streams.txt.CLM1PT.CLM_USRDAT
      chmod u+w user_datm.streams.txt.CLM1PT.CLM_USRDAT

    * edit the domainInfo filePath value to point to
      atm/datm7/domain.clm/domain.lnd.XXXX.nc instead of
      share/domains/domain.clm/domain.lnd.YYYY.nc

* for more details see: http://www.cesm.ucar.edu/models/cesm1.2/data8/doc/c72.html#streams_description_file

src: bounds_type refactor (clm4_5_21)
-------------------------------------

* a new derived type has been created, bounds_type, in decomMod.F90

* all begc, endc, begg, endg, etc should be retrieved from
  get_proc_bounds as part of a bounds type

* all passing of bounds in subroutines should be as a single
  bounds_type instead of six separate variables.

src: bound declarations (clm4_5_29)
-----------------------------------

how arrays and bounds are passed to subroutines has been changed to
prevent some common errors and allow better compiler diagnostics.

* From the changelog:

    Rework bounds declarations for subroutine array arguments, both in
    caller (explicitly subset argument by bounds) and callee (use
    assumed-shape array arguments rather than declaring upper bounds), and
    add assertions on array sizes.

    See https://wiki.ucar.edu/display/ccsm/Community+Land+Model+Developers+Guide
    ("Guidelines for passing array arguments to subroutines") for the new
    conventions that are implemented here.

