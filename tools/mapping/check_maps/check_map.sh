#!/bin/bash

#==============================================================================
# $Id: check_map.sh 46158 2013-04-19 18:41:34Z mlevy@ucar.edu $
# $URL: https://svn-ccsm-models.cgd.ucar.edu/tools/mapping/trunk_tags/mapping_130716/check_maps/check_map.sh $
#
# Use an updated ESMF tool to check the quality of a mapping file
#
#==============================================================================
date
SDIR=`dirname $0`

#==============================================================================
# Usage subroutine
#==============================================================================

usage() {
	echo 'USAGE: ./check_maps.sh [OPTION]... FILELIST'
	echo 'Runs a modified version of the ESMF RegridWeightGenCheck program'\
	     'over files listed in FILELIST'
	echo ''
	echo '  --recompile, -rc  Force recompile (necessary to change verbose)'
	echo '  --verbose, -v     Compile with verbose output (use with -rc)'
	echo '  --clean, -c       Remove log / aux files generated by this script'
	echo '  --help, -h        Output this usage information'
	echo ''
	echo 'Notes:'
	echo '  1) For use on yellowstone, geyser, or caldera only!'
	echo '  2) Need to set ESMFMKFILE (see comments in Makefile)'\
	          'or compilation will fail'
	echo '  3) If -rc option is not enabled, -v flag is ignored and verbose /'\
	          'concise will depend on previous compilation'
}

#==============================================================================
# Main Program
#==============================================================================

# Process input arguments
verbose="FALSE"
compile="FALSE"
FileList=""
while [ $# -gt 0 ]; do
	case $1 in
		-rc|--recompile)
			compile="TRUE"
		;;
		-v|--verbose)
			verbose="TRUE"
		;;
		-h|--help)
			usage
			exit 0
		;;
		-c|--clean)
			rm -f *.Log hostfile
			echo 'Removed all .Log files and hostfile'
			exit 0
		;;
		* )
			if [ -e $1 ]; then
				FileList="$FileList $1"
			else
				echo "File not found: $1"
			fi
		;;
	esac
	shift
done

if [ -z "$FileList" ]; then
	echo "No files given!"
	usage
	exit 1
fi

if [ $compile == "TRUE" ]; then
	echo "Building $EXE"
	CURR_DIR=$PWD
	cd $SDIR/src
	gmake VERBOSE=$verbose
	cd $CURR_DIR
fi

EXE=$SDIR/ESMF_RegridWeightGenCheck
if [ ! -e $EXE ]; then
	echo "WARNING: $EXE not found. To check quality of maps, build this tool manually or use the -rc flag"
	exit 1
fi

declare -i n=0

for MAP in $FileList
do
	if [ -e $MAP ]; then
		n=n+1
		echo "${n}: ${MAP}"
		$EXE $MAP
		echo "-----"
	else
		echo "File not found: $MAP"
	fi
done

