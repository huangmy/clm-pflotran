#!/bin/sh

# run python scripts: runCLM.py to build/run a case
# e.g. Darwin with gnu-mpich: $bash DEMO_runCLMpy.sh mac
#      LINUX with pgi-openpgi: $bash DEMO_runCLMpy.sh linux

# the following is editable for your case
# the example here create a case named as 'US-Brw_I1850CLM45CN_ad_spinup',
# in which, the first portion is the site-name, the second portion is the compset name, the third (optional) is for 'AD_spinup' 

# NOTES for the scripts calling by line number (actually lines are continuous)
# 1. site-name, and where is the data group to look for it
# 2. case root directory under which the case will be created and set-up
# 3. case build/run root directory under which the case will be built/compiled/run
# 4. ccsm inputdata directory
# 5. cesm model root directory
# 6. compset and its CLM-CN module setting
# 7. platform and tools information for case build and run
# 8. point/regional mode, and if regional mode, the x(lon)/y(lat) axis point numbers
# 9. misc. options to do case creation, building and run ('--no_submit' for not running)

export os=$1

export pflotran_coupled=$2

if [ "$os" == "mac" ]; then
    if [ "$pflotran_coupled" == "pflotran" ]; then
        python ./runCLM.py --site=US-Brw --sitegroup=AmeriFlux \
--caseroot=/Users/f9y/mygit/clm4-pf/cases \
--runroot=/Users/f9y/clm4_5_simulations \
--ccsm_input=/Users/f9y/clm4_5_inputdata \
--cesmdir=/Users/f9y/mygit/clm4-pf \
--compset=I1850CLM45CN --coldstart --vertsoilc --CH4 --no_fire --ad_spinup --nyears_ad_spinup=10 --clm_pflotran \
--machine=userdefined --osname=Darwin --compiler=gnu --ninst=1 --np=4 --mpilib=mpich \
--nopointdata --regional --xpts=3 --ypts=5 \
--rmold --clean_config --clean_build

    else
        python ./runCLM.py --site=US-Brw --sitegroup=AmeriFlux \
--caseroot=/Users/f9y/mygit/clm4-pf/cases \
--runroot=/Users/f9y/clm4_5_simulations \
--ccsm_input=/Users/f9y/clm4_5_inputdata \
--cesmdir=/Users/f9y/mygit/clm4-pf \
--compset=I20TRCLM45CN --finidat=US-Brw_I1850CLM45CN_spinup.clm2.r.0601-01-01-00000.nc --finidat_year=601 --vertsoilc --CH4 --no_fire --hist_userdefined=clm_output_tr \
--machine=userdefined --osname=Darwin --compiler=gnu --ninst=1 --np=1 --mpilib=mpi-serial \
--nopointdata --xpts=1 --ypts=1 \
--rmold --clean_config --clean_build
   fi

fi

if [ "$os" == "linux" ]; then
    python ./runCLM.py --site=US-Brw --sitegroup=AmeriFlux \
--caseroot=/home/f9y/cesm/clm4-pf/cases \
--runroot=/home/f9y/cesm/clm4_5_simulations \
--ccsm_input=/home/f9y/cesm/clm4_inputdata \
--cesmdir=/home/f9y/cesm/clm4-pf \
--compset=I1850CLM45CN --coldstart --vertsoilc --CH4 --no_fire --ad_spinup --nyears_ad_spinup 10 \
--machine=userdefined --osname=LINUX --compiler=pgi --ninst=1 --np=100 --mpilib=openmpi \
--nopointdata --regional --xpts=50 --ypts=100 \
--rmold --clean_config --clean_build

fi
