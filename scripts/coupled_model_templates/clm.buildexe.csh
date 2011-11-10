#! /bin/csh -f 

set objdir = $OBJROOT/lnd/obj; cd $objdir

#--------------------------------------------------------------------
# check basic task and thread settings
#--------------------------------------------------------------------

cp -f $CASEBUILD/clmconf/CESM_cppdefs .tmp
cmp -s .tmp CESM_cppdefs || mv -f .tmp CESM_cppdefs

setenv COMP "unknown"
if ($COMP_INTERFACE == 'MCT' ) setenv COMP mct
if ($COMP_INTERFACE == 'ESMF') setenv COMP esmf

\cat >! .tmp << EOF; cmp -s .tmp Filepath || mv -f .tmp Filepath
$CASEROOT/SourceMods/src.clm
$CODEROOT/lnd/clm/src/cpl_share
$CODEROOT/lnd/clm/src/main
$CODEROOT/lnd/clm/src/biogeophys
$CODEROOT/lnd/clm/src/biogeochem
$CODEROOT/lnd/clm/src/cpl_$COMP
$PFLOTRAN_COUPLED_MODEL/src/pflotran
$PFLOTRAN_COUPLED_MODEL/src/clm-pflotran
EOF

#
# Build the clm library
#
set clmdefs = "  -DMAXPATCH_PFT=17 -DCN -DAD_SPINUP -D_USEBOX -D_NETCDF -DCLM_PFLOTRAN -DWITH_CLM"
if ( ! $?GMAKE ) setenv GMAKE gmake
$GMAKE complib -j $GMAKE_J MODEL=clm COMPLIB=$LIBROOT/liblnd.a MACFILE=$CASEROOT/Macros.$MACH USER_CPPDEFS="$clmdefs" -f $CASETOOLS/Makefile || exit 2

