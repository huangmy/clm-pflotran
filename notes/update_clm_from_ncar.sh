#!/bin/sh
#
# Use:
#
#   mkdir clm4-ncar-updates
#   cd clm4-ncar-updates
#   hg clone https://bitbucket.org/clm_pflotran/clm4-ncar
#   ln -s clm4-ncar/lbl_notes/update_clm_from_ncar.sh .
#   ./update_clm_from_ncar.sh clm4_5_22
#
# Steps
#
#   - run update_clm_from_ncar.sh with the clm trunk tag as the only
#     command line parameter
#
#   - After script finishes, cd clm_ncar_update_${NEW_TAG} run:
#     'hg status | grep svn' to ensure no svn files are added.
#
#   - Do hg status to make sure everything looks OK. Verify that the
#     expected changes are present by looking at the clm changelog if
#     necessary.
#
#   - hg commit -m "update to clm4_5_XXX tag from NCAR"
#
#   - hg push
#
# NOTE: All NCAR commits on 'ncar' branch and all non-NCAR related
# commits (ie "new science") should take place on seperate named
# branches.
#
#

NEW_TAG=$1

hg clone clm4-ncar clm_ncar_update_${NEW_TAG}
if [ "$?" -ne "0" ]; then
    echo "ERROR: could not clone clm4-ncar repository"
    exit 1
fi

cd clm_ncar_update_${NEW_TAG}
if [ "$?" -ne "0" ]; then
    echo "ERROR: could not cd into clm4-ncar repository"
    exit 1
fi

hg update ncar
if [ "$?" -ne "0" ]; then
    echo "ERROR: could update clm4-ncar to ncar branch"
    exit 1
fi

rm -rf \
	.ChangeLog_template \
	ChangeLog \
	ChangeSum \
	Copyright \
	README \
	README_EXTERNALS \
	SVN_EXTERNAL_DIRECTORIES \
	UpDateChangeLog.pl \
	models \
	scripts \
	tools

if [ "$?" -ne "0" ]; then
    echo "ERROR: could update remove old files from clm4-ncar"
    exit 1
fi

svn co https://svn-ccsm-models.cgd.ucar.edu/clm2/trunk_tags/${NEW_TAG} .
if [ "$?" -ne "0" ]; then
    echo "ERROR: could pull new tag ${NEW_TAG} from ncar svn"
    exit 1
fi

hg addremove
if [ "$?" -ne "0" ]; then
    echo "ERROR: hg addremove failed..."
    exit 1
fi


echo "Checking for svn files in hg status..."
hg status | grep svn
