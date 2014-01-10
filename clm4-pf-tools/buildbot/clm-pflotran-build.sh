#!/usr/bin/env bash

################################################################################
#
# Build script for the clm-pflotran buildbot system.
#
# Author: Ben Andre <bandre@lbl.gov>
#
# All the logic for buildbot is moved out of the buildbot config and
# placed here.
#
# Based on the pflotran-build.sh script for pflotran-dev.
#
# Requirements:
#
#  - Assumes that the script is being called from the root of the
#    clm-pflotran repo. If not, then the clm-pflotran directory must
#    be set on the command line (-p).
#
#  - Assumes that petsc and pflotran were built with the pflotran-dev
#    pflotran-build.sh script as part of the build process.
#
# Builder requirements:
#
#  - requirements for pflotran-dev and pflotran-dev buildbot
#
#    - git, mercurial, working mpi compilers
#
#    - a working petsc configure script that will build petsc and any
#      additional TPLs needed for that builder, i.e. hdf5, metis and
#      parmetis for unstructured mesh.
#
#    - NOTE: petsc must be built with '--with-shared-libraries=0'
#      because we check for the presence of libpetsc.a to figure out
#      if a petsc build is acceptable.
#
#  - requirements for clm-pflotran
#
#    - petsc and pflotran must have been built and tested by the
#      pflotran buildbot script.
#
#    - the buildbot must provide the absolute paths for PETSC_DIR,
#      PFLOTRAN_DIR, and the cesm input data and test baseline
#      directories so that a local config file can be generated for
#      clm-pflotran-tests.py
#
################################################################################

################################################################################
#
# Global variables
#
################################################################################
BUILDER_ID=
CLM_PFLOTRAN_DIR=$PWD
PETSC_DIR=
PETSC_ARCH=
PFLOTRAN_DIR=
CESM_INPUTDATA_DIR=
BUILD_STATUS=0

################################################################################
#
# work routines
#
################################################################################

function write-local-config() {
    _cfg_file=local.cfg
    echo "Writing ${_cfg_file} with :"
    
    cat > ${_cfg_file} <<EOF
[petsc]
petsc_dir = ${PETSC_DIR}
petsc_arch = ${PETSC_ARCH}

[pflotran]
pflotran_dir = ${PFLOTRAN_DIR}

[data]
test_data_dir = ${CESM_INPUTDATA_DIR}

EOF

    echo "Local config file: "
    cat ${_cfg_file}

}

function set-builder-info() {
    BUILDER_ID=`hostname -s`
    echo "clm-pflotran builder id : ${BUILDER_ID}"

    echo "CLM_PFLOTRAN_DIR : ${CLM_PFLOTRAN_DIR}"

    if [ -z ${CESM_INPUTDATA_DIR} ]; then
        CESM_INPUTDATA_DIR=${CLM_PFLOTRAN_DIR}/../clm-pflotran-data-trunk-testing
    fi

    if [ ! -d ${CESM_INPUTDATA_DIR} ]; then
        echo "ERROR: Could not find the clm-pflotran-data-trunk-testing directory at the expected location :"
        echo "    ${CESM_INPUTDATA_DIR}"        
        BUILD_STATUS=1
    else
        echo "Using CESM_INPUTDATA_DIR :"
        echo "    ${CESM_INPUTDATA_DIR}"        
    fi
}

function check-pflotran() {
    echo "Checking for pflotran :"

    if [ -z ${PFLOTRAN_DIR} ]; then
        PFLOTRAN_DIR=${CLM_PFLOTRAN_DIR}/../pflotran-clm-trunk
    fi

    if [ ! -d ${PFLOTRAN_DIR} ]; then
        echo "ERROR: Could not find the pflotran-clm source at the expected location :"
        echo "    ${PFLOTRAN_DIR}"
        BUILD_STATUS=1
    else
        echo "Using PFLOTRAN_DIR :"
        echo "    ${PFLOTRAN_DIR}"
    fi
}

function check-petsc() {
    echo "Checking for petsc :"

    _petsc_version_file=${PFLOTRAN_DIR}/tools/buildbot/petsc/petsc-git-version.txt
    if [ ! -f ${_petsc_version_file} ]; then
        echo "ERROR: could not find petsc version file : ${_petsc_version_file}"
        BUILD_STATUS=1
    fi

    PETSC_REQUIRED_VERSION=`cat ${_petsc_version_file}`
    echo "PFLOTRAN requires PETSc git reversion ${PETSC_REQUIRED_VERSION}"

    PETSC_DIR=${PFLOTRAN_DIR}/../petsc.git.${PETSC_REQUIRED_VERSION:0:8}
    PETSC_ARCH=${BUILDER_ID}-${COMPILER}
    echo "Requiring petsc env: "
    echo "    PETSC_DIR=${PETSC_DIR}"
    echo "    PETSC_ARCH=${PETSC_ARCH}"
    echo ""

    _lib_petsc=${PETSC_DIR}/${PETSC_ARCH}/lib/libpetsc.a

    if [ -d ${PETSC_DIR} ]; then
        echo "Found existing PETSc directory."
        cd ${PETSC_DIR}
        # assume that the petsc version was checked by pflotran
        _petsc_install_version=`git log --pretty="%H" -1 HEAD`
        if [ -d ${PETSC_DIR}/${PETSC_ARCH} ]; then
            echo "Found expected PETSC_ARCH directory :"
            echo "    ${PETSC_DIR}/${PETSC_ARCH}"
            if [ -f ${_lib_petsc} ]; then
                echo "Found libpetsc.a : PETSc appears to be installed."
            else
                echo "PETSc : could not find libpetsc.a for this PETSC_ARCH."
                echo "    ${_lib_petsc}"
                BUILD_STATUS=1            
            fi
        else
            echo "ERROR: could not find petsc arch directory :"
            echo "    ${PETSC_DIR}/${PETSC_ARCH}"
            BUILD_STATUS=1            
        fi
    else
        echo "ERROR: could not find petsc directory :"
        echo "    ${PETSC_DIR}"
        BUILD_STATUS=1
    fi
}

function stage-libpflotran-build() {
    echo "Building PFLOTRAN with :"
    echo "  PETSC_DIR=${PETSC_DIR}"
    echo "  PETSC_ARCH=${PETSC_ARCH}"
    export PETSC_DIR PETSC_ARCH
    _pflotran_flags=
    _flags_file=${PFLOTRAN_DIR}/tools/buildbot/build-flags/${BUILD_FLAGS}.txt
    if [ -f ${_flags_file} ]; then
        _pflotran_flags=`cat ${_flags_file}`
        echo "  pflotran build flags=${_pflotran_flags}"
    else
        echo "Could not find build flags file: ${_flags_file}. Building with 'make libpflotran.a'."
    fi
    
    cd ${PFLOTRAN_DIR}/src/clm-pflotran
    ./link_files.sh
    make clean
    make ${_pflotran_flags} libpflotran.a
    BUILD_STATUS=$?
}

function stage-clm-pflotran-common-exe() {
    echo "Building CLM common executable :"

    cd ${CLM_PFLOTRAN_DIR}/clm4-pf-tools/regression_tests
    write-local-config
    rm -rf *.testlog
    make tests/common-executable.cfg
    BUILD_STATUS=$?
    cat *.testlog
}

function stage-clm-pflotran-tests() {
    echo "Running CLM-PFLOTRAN test problems :"

    cd ${CLM_PFLOTRAN_DIR}/clm4-pf-tools/regression_tests
    rm -rf *.testlog
    make fast
    BUILD_STATUS=$?
    cat *.testlog
}

################################################################################
#
# main program
#
################################################################################
function usage() {
     echo "
Usage: $0 [options]
    -b BUILD_FLAGS    group of build flags to use.
    -c COMPILER       compiler name: gnu, pgi, intel
    -h                print this help message
    -p PFLOTRAN_DIR   root directory for the build (default: '.')
    -s BUILD_STAGE    build stage must be one of:
                        all libpflotran common-exe clm-pf-tests

Notes:

  - The build flags group must have a corresponding file in the
    tools/buildbot/build-flags/ directory.

"
}

# setup based on commandline args
BUILD_FLAGS="__NONE__"
BUILD_STAGE=
COMPILER=
while getopts "b:c:hp:s:" FLAG
do
  case ${FLAG} in
    b) BUILD_FLAGS=${OPTARG};;
    c) COMPILER=${OPTARG};;
    h) usage;;
    i) CESM_INPUTDATA_DIR=${OPTARG};;
    p) PFLOTRAN_DIR=${OPTARG};;
    s) BUILD_STAGE=${OPTARG};;
  esac
done

# verify all required info is set
if [ -z "${BUILD_STAGE}" ]; then
    echo "ERROR: The build stage name must be provided on the command line."
    exit 1
fi

if [ -z "${COMPILER}" ]; then
    echo "ERROR: The compiler name must be provided on the command line."
    exit 1
fi

set-builder-info
check-pflotran
check-petsc

if [ "${BUILD_STATUS}" -ne "0" ]; then
    exit ${BUILD_STATUS}
fi

case ${BUILD_STAGE} in
    all) stage-libpflotran-build; stage-clm-pflotran-common-exe; stage-clm-pflotran-tests;;
    libpflotran) stage-libpflotran-build;;
    common-exe) stage-clm-pflotran-common-exe;;
    clm-pf-tests) stage-clm-pflotran-tests;;
    *) echo "ERROR: The requested build stage '${BUILD_STAGE}' is invalid."; exit 1;;
esac

exit ${BUILD_STATUS}


