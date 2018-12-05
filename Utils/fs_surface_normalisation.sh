#!/bin/bash
# fs_surface_normalisation.sh
#
# About: Takes a Freesurfer subject name and performs spatial normalisation
# of the surface files via -qchache. Bear in mind this takes a while (>1h)
#
### Changelog ###
# Works.
# Ivan Alvarez 
# UCL 25/06/2013
#

# If no argument provided, ask for subject name
if [ -z "$1" ];
	then
	echo -n "Freesurfer subject: "
	read SUBJ

	# Normalisation
	recon-all -s $SUBJ -qcache
else

# If argument provided, cycle through all arguments
	for SUBJ in "$@"
		do
		echo
		echo "Subject $SUBJ"

		# Perform
		recon-all -s $SUBJ -qcache
	done
fi

# Done
echo "done!"
