#!/bin/bash
# fs_bin2asc
#
# About: Takes a Freesurfer subject name, goes into the
# 'surf' directory and converts a bunch of binary surface 
# files to ASCII format. Original files are retained.
#
### Changelog ###
# Works.
# Ivan Alvarez 
# UCL 18/06/2013
#
# Added thickness conversion & changed name.
# Sam Schwarzkopf
% 26/06/2013
#

# If no argument provided, ask for subject name
if [ -z "$1" ];
	then
	echo -n "Freesurfer subject: "
	read SUBJ

	# Go there
	cd $SUBJECTS_DIR/$SUBJ/surf

	# Convert
	echo "converting files..."
	mris_convert -c lh.area lh.white lh.area.asc
	mris_convert -c rh.area rh.white rh.area.asc
	mris_convert -c lh.curv lh.white lh.curv.asc
	mris_convert -c rh.curv rh.white rh.curv.asc
	mris_convert -c lh.thickness lh.white lh.thickness.asc
	mris_convert -c rh.thickness rh.white rh.thickness.asc

else
# If argument provided, cycle through all arguments

	for SUBJ in "$@"
		do
		echo
		echo "Subject $SUBJ"

		# Get there
		cd $SUBJECTS_DIR/$SUBJ/surf

		# Convert
		echo "converting files..."
		mris_convert -c lh.area lh.white lh.area.asc
		mris_convert -c rh.area rh.white rh.area.asc
		mris_convert -c lh.curv lh.white lh.curv.asc
		mris_convert -c rh.curv rh.white rh.curv.asc
        mris_convert -c lh.thickness lh.white lh.thickness.asc
        mris_convert -c rh.thickness rh.white rh.thickness.asc
	done
fi

# Done
echo "done!"
