#!/bin/bash
#===============================================================================
# SVN $Id: create_ESMF_map.sh 45575 2013-04-03 18:43:22Z mlevy@ucar.edu $
# SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/tools/mapping/trunk_tags/mapping_130403/gen_mapping_files/gen_ESMF_mapping_file/create_ESMF_map.sh $
#
# Create needed mapping files for gen_domain and coupler mapping
# Currently only supported on bluefire and jaguarpf - and only in
# interactive mode
# 
#===============================================================================
echo $0
date
SDIR=`dirname $0`

#===============================================================================
# Usage subroutine
#===============================================================================
usage() {
  echo ''
  echo '**********************************************************'
  echo 'usage:'
  echo './create_ESMF_map.sh  '
  echo ' A wrapper for the ESMF mapping tool that creates a mapping file'
  echo ' from the source grid to the destination grid. Specify what type'
  echo ' of mapping to use with the -maptype flag (aave, blin, patc,'
  echo ' nearestdtos, or neareststod)'
  echo ''
  echo 'create_ESMF_map.sh '
  echo '  --filesrc|-fsrc  input source grid_filename (required) '
  echo '  --filedst|-fdst  input destination grid_filename (required)'
  echo '  --namesrc|-nsrc  output source name in mapping file (required)' 
  echo '  --namedst|-ndst  output destination name in mapping file (required)'
  echo '  --maptype|-map   type of mapping [aave|blin|patc|nearestdtos|neareststod] (required)'
  echo '  [ --typesrc|tsrc ] [regional|global]'
  echo '  [ --typedst|tdst ] [regional|global]'
  echo '  [ --batch|-b ]'
  echo '  [ --large_file|-big ]'
  echo '  [ --help|-h ]'
  echo '  [ -v ]'
  echo ' '
  echo 'where '
  echo ' --filesrc (or -fsrc) '
  echo '   SCRIP grid format source filename (full pathname)'
  echo ' --filedst (or -fdst) '
  echo '   SCRIP grid format destination filename (full pathname)'
  echo ' --namesrc (or -nsrc) and --namesrc (or -nsrc) will result in the '
  echo '   following mapping files'
  echo '     namesrc_TO_namedst_maptype.cdate.nc'
  echo ''
  echo ' --typesrc (or -tsrc) '
  echo '   source grid type,  valid values are regional or global'
  echo '   default is global'
  echo ' --typedst (or -tdst) '
  echo '   destination grid type, valid values are regional or global'
  echo '   default is global'
  echo ' --batch (or -b) '
  echo '   Toggles batch mode usage. If you want to run in batch mode'
  echo '   you need to have a separate batch script for a supported machine'
  echo '   that calls this script interactively - you cannot submit this'
  echo '   script directly to the batch system'
  echo ' -d '
  echo '   toggle debug-only '
  echo ' --help or -h  '
  echo '   displays this help message'
  echo ''
  echo 'You can also set the following env variables:'
  echo '  ESMFBIN_PATH - Path to ESMF binaries '
  echo '                 (Leave unset on yellowstone and caldera and the tool'
  echo '                 will be loaded from modules)'
  echo '  MPIEXEC ------ Name of mpirun executable'
  echo '                 (default is mpirun.lsf on yellowstone and caldera; if'
  echo '                 you run interactively on yellowstone, mpi is not used)'
  echo '  REGRID_PROC -- Number of MPI processors to use'
  echo '                 (default is 8)'
  echo '**********************************************************'
}

#===============================================================================
# runcmd subroutine
#===============================================================================
runcmd() {
   cmd=$@
   if [ -z "$cmd" ]; then
       echo "No command given to the runcmd function"
       exit 3
   fi
   if [ "$verbose" = "YES" ]; then
       echo "$cmd"
   fi
   if [ "$debug" != "YES" ]; then
       ${cmd}
       rc=$?
   else
       rc=0
   fi
   if [ $rc != 0 ]; then
       echo "Error status returned from create_ESMF_map script"
       exit 4
undo
   fi
   return 0
}

#===============================================================================
# Main program
#===============================================================================

#-------------------------------------------------------------------------------
# Process input arguments
#-------------------------------------------------------------------------------

interactive="YES"
debug="no"
verbose="no"
type_src="global"
type_dst="global"
use_large="false"
CheckMapsFlag=""
use_rtm=0


while [ $# -gt 0 ]; do
   case $1 in
       -v)
	   verbose="YES"
	   ;;
       -b|--batch)
	   interactive="NO"
	   ;;
       -big|--large_file)
	   use_large="true"
	   ;;
       -fsrc|--filesrc )
	   fsrc=$2
	   shift
	   ;;
       -fdst|--filedst )
	   fdst=$2
	   shift
	   ;;
       -nsrc|--namesrc )
	   nsrc=$2
	   shift
	   ;;
       -ndst|--namedst )
	   ndst=$2
	   shift
	   ;;
       -map|--maptype )
	   map_type=$2
	   echo "map_type is $map_type"
	   shift
	   ;;
       -tsrc|--typesrc )
	   type_src=$2
	   echo "type_src is $type_src"
	   shift
	   ;;
       -tdst|--typedst )
	   type_dst=$2
	   echo "type_dst is $type_dst"
	   shift
	   ;;
       -h|--help )
	   usage
	   exit 0
	   ;;
       * )
	   echo "****************************"
	   echo "ERROR:: invalid argument $1"
	   echo "****************************"
	   usage
	   exit 1
	   ;;
   esac
   shift 
done

# check for required arguments
echo "fsrc is $fsrc"
echo "fdst is $fdst"
if [ -z "$fsrc" ]; then
    echo "Must specfiy -fsrc or --filesrc argument "
    echo "Invoke create_ESMF_map.sh -h for usage"
    exit 1
fi
if [ -z "$fdst" ]; then
    echo "Must specfiy -fdst or --filedst argument "
    echo "Invoke create_ESMF_map.sh -h for usage"
    exit 2
fi
if [ -z "$nsrc" ]; then
    echo "Must specfiy -nsrc or --namesrc argument "
    echo "Invoke create_ESMF_map.sh -h for usage"
    exit 3
fi
if [ -z "$ndst" ]; then
    echo "Must specfiy -ndst or --namedst argument "
    echo "Invoke create_ESMF_map.sh -h for usage"
    exit 4
fi
if [ -z "$map_type" ]; then
    echo "Must specfiy -map or --maptype argument "
    echo "Invoke create_ESMF_map.sh -h for usage"
    exit 5
fi

# check for existence of files
if [ ! -f "${fsrc}" ]; then
   echo "Source grid file does NOT exist: $fsrc}"
   exit 6
fi
if [ ! -f "${fdst}" ]; then
   echo "Destination grid file does NOT exist: $fdst"
   exit 7
fi

# check for type of map
if [ $map_type != "aave" ] && [ $map_type != "blin" ] && [ $map_type != "nearestdtos" ] && [ $map_type != "neareststod" ] && [ $map_type != "patc" ]; then
  echo "ERROR: $map_type is not a valid type of mapping."
  echo "(must be aave, blin, or patc)"
  exit 8
fi

# set some defaults
if [ -z "$REGRID_PROC" ]; then
   REGRID_PROC=8
fi

#-------------------------------------------------------------------------------
# Machine specific env stuff
#-------------------------------------------------------------------------------
 
hostname=`hostname`
case $hostname in
  ## yellowstone
  ys* )
    module purge
    module load intel
    module load nco
    module load esmf
    if [ $interactive == "YES" ]; then
      module load esmf-6.1.1-ncdfio-uni-O
      MPIEXEC=""
  	else
      # Batch script
      module load esmf-6.1.1-ncdfio-mpi-O

      if [ -z "$MPIEXEC" ]; then
	      MPIEXEC="mpirun.lsf"
      fi
    fi

  ;;
  ## geyser or caldera (interactive on yellowstone)
  geyser* )
    echo "At this time, the ESMF tools are not available on geyser. To run interactively, please switch to caldera and run this script again."
    exit 1
  ;;
  caldera* )
    module purge
    module load intel
    module load nco
    module load esmf
    module load esmf-6.1.1-ncdfio-mpi-O
    if [ -z "$MPIEXEC" ]; then
	    MPIEXEC="mpirun.lsf"
    fi

    # specific commands to prepare to run interactively
    if [ $interactive == "YES" ]; then
	    export MP_PROCS=$REGRID_PROC
	    export MP_EUILIB=ip
	    
	    hostname > hostfile
	    declare -i p=2
	    until ((p>$MP_PROCS)); do
    		hostname >> hostfile
    		p=p+1
	    done
	    export MP_HOSTFILE=hostfile
  	fi
  ;;
  ## bluefire
  be* )
  	if [ -z "$ESMFBIN_PATH" ]; then
	    ESMFBIN_PATH=/contrib/esmf-5.2.0r-64-O/bin
  	fi

    if [ -z "$MPIEXEC" ]; then
	    MPIEXEC="mpirun.lsf"
    fi
    
		# Disable MP_EUIDEVICE to avoid warning message
    if [ ! -z "$MP_EUIDEVICE" ]; then
		  MP_EUIDEVICE_tmp=$MP_EUIDEVICE
  		unset MP_EUIDEVICE
	  fi
    
		# Disable MP_INSTANCES to avoid warning message
    if [ ! -z "$MP_INSTANCES" ]; then
		  MP_INSTANCES_tmp=$MP_INSTANCES
  		unset MP_INSTANCES
	  fi
    
    # specific commands to prepare to run interactively
    if [ $interactive == "YES" ]; then
	    export MP_PROCS=$REGRID_PROC
	    export MP_EUILIB=ip
	    MPIEXEC=""
	    
	    hostname > hostfile
	    declare -i p=2
	    until ((p>$MP_PROCS)); do
    		hostname >> hostfile
    		p=p+1
	    done
	    export MP_HOSTFILE=hostfile
  	fi
	;;
	
  ##jaguarpf
  ## NOTE that for jaguarpf there is no batch script for now
  jaguarpf* )
    if [ -z "$ESMFBIN_PATH" ]; then
      module load esmf/5.2.0-p1_with-netcdf_g
	    ESMFBIN_PATH=$ESMF_BINDIR
  	fi
    
    if [ -z "$MPIEXEC" ]; then
	    MPIEXEC="aprun -n $REGRID_PROC"
  	fi
	;;
    
  *)
	  echo "Machine $hostname NOT recognized"
	;;
    
esac

#-------------------------------------------------------------------------------
# run ESMF_RegridWeightGen
#-------------------------------------------------------------------------------

# Resolve interactive or batch mode command
# NOTE - if you want to run in batch mode - you need to have a separate
# batch file that calls this script interactively - you cannot submit
# this script to the batch system

if [ "$interactive" = "YES" ]; then
    echo "Running interactively"
else
    echo "Running in batch mode"
fi

#if [ ! -d "$ESMFBIN_PATH" ]; then
#    echo "Path to ESMF binary directory does NOT exist: $ESMFBIN_PATH"
#    echo "Set the environment variable: ESMFBIN_PATH"
#    exit 1
#fi

if [ ! -z $ESMFBIN_PATH ]; then
  ESMF_REGRID="$ESMFBIN_PATH/ESMF_RegridWeightGen"
else
  ESMF_REGRID="ESMF_RegridWeightGen"
fi

# Make sure $ESMF_REGRID is a valid command
command -v $ESMF_REGRID >/dev/null 2>&1 || { echo "Can not find ESMF_RegridWeightGen, make sure it is in your \$PATH or that you specify \$ESMFBIN_PATH." && exit 1; }

# Remove previous log files
rm -f PET*.Log

# Set output map name and create it
#date="c"`date +%y%m%d`    omit extraneous "c" character in date stamp
cdate=`date +%y%m%d`

mapname=${nsrc}_TO_${ndst}

fmap=map_${mapname}_${map_type}.${cdate}.nc
echo ""
echo "Creating $fmap"

mapping="NULL"
case $map_type in
  "aave")
    mapping="conserve"
  ;;
  "blin")
    mapping="bilinear"
  ;;
  "patc")
    mapping="patch -p all"
  ;;
  "nearestdtos")
    mapping="nearestdtos"
  ;;
  "neareststod")
    mapping="neareststod"
  ;;
esac

if [ "$mapping" == "NULL" ]; then
  echo "ERROR: $map_type is not a valid option for --maptype"
  exit 9
fi
cmd="$MPIEXEC $ESMF_REGRID --ignore_unmapped -m $mapping -w $fmap -s $fsrc -d $fdst"

if [ $use_large == "true" ]; then
  cmd="$cmd --64bit_offset"
fi

if [ $type_src == "regional" ]; then
    cmd="$cmd --src_regional"
    echo "regional source"
fi
if [ $type_dst == "regional" ]; then
    cmd="$cmd --dst_regional"
    echo "regional destination"
fi
echo "cmd is $cmd"
echo ""
runcmd $cmd
if [ "$debug" != "YES" ] && [ ! -f "$fmap" ]; then
   echo "Output mapping file $fmap was NOT created: $fmap"
   exit 4
fi
HOST=`hostname`
history="$ESMF_REGRID"
runcmd "ncatted -a history,global,a,c,$history -a hostname,global,a,c,$HOST -a logname,global,a,c,$LOGNAME $fmap"

echo "Successfully created mapping file $fmap "
