This repository contains clm history data for regression testing
simulations against a known baseline.

We expect the baseline data to change frequently. Mercurial can not
efficiently store binary file formats such as netcdf. When the change
frequently the file total repository size will grow quickly. To work
around this we store files an "CDL", the text format output from ncdump.

DO NOT add the binary netcdf files to the repo.

To update the regression files, copy the new binary netcdf file
specified in the test configuration file into this repository. The
regression suite will automatically convert it to a cdl file. You can
simply commit the changes to the repo to capture the changes to the
cdl regression baseline.

