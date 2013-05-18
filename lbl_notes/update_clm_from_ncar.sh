#!/bin/sh
NEW_TAG=clm4_5_02
hg clone clm_ncar clm_ncar_update_${NEW_TAG}
cd clm_ncar_update_${NEW_TAG}
hg update ncar
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

svn co https://svn-ccsm-models.cgd.ucar.edu/clm2/trunk_tags/${NEW_TAG} .

hg addremove


