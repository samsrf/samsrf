#!/bin/bash
# fs_vol2mgh
#
# Use this to project a NII volumetric file into a FreeSurfer MGH surface file.
#
#   sh fs_vol2mgh.sh [Filenames Without NII Extension Separated By Space]
#
# This script uses a cortical sampling step of 0.5 and 3 mm surface smoothing. Adjust as desired.
#

# Convert
for F in "$@" 
do
    echo "Converting $F"
    mri_vol2surf --mov $F.nii --reg register.dat --hemi lh --o ./lh_$F.mgh --projfrac 0.5 --surf-fwhm 3
    mri_vol2surf --mov $F.nii --reg register.dat --hemi rh --o ./rh_$F.mgh --projfrac 0.5 --surf-fwhm 3
done
