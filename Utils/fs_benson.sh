#!/bin/bash
# fs_bin2asc
#
# Runs the conversion of Benson maps using the template file: 
#   all-template-2.5.sym.mgh
# You probably want to use a different template because this one is old...
#
# You must be the FreeSurfer subjects folder where all your subjects are for this
# and your template file should also be located there.
#
# Usage:
#   sh fs_benson [SubjectID]
#

surfreg --s "$@" --t fsaverage_sym --lh
surfreg --s "$@" --t fsaverage_sym --lh --xhemi 

mkdir "$@"/benson

mri_surf2surf --srcsubject fsaverage_sym --srcsurfreg sphere.reg --trgsubject "$@" --trgsurfreg fsaverage_sym.sphere.reg --hemi lh --sval all-template-2.5.sym.mgh --tval "$@"/benson/lh.benson.mgh
mri_surf2surf --srcsubject fsaverage_sym --srcsurfreg sphere.reg --trgsubject "$@"/xhemi --trgsurfreg fsaverage_sym.sphere.reg --hemi lh --sval all-template-2.5.sym.mgh --tval "$@"/benson/rh.benson.mgh