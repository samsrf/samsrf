#!/bin/bash
# fs_regdat
#
# Use this to create a register.dat file for subject ID $1 and the NII file $2
#
# Usage:
#   sh fs_regdat.sh [Subject ID] [NII File Without Extension]
#

tkregister2 --mov $2.nii --s $1 --regheader --noedit --reg register.dat
