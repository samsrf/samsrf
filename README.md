# SamSrf VI - Read Me #

### Version 6.0 (10-08-2018) 

Complete overhaul of the toolbox. The GUI has been abolished. 
Instead there are now example scripts for various pRF models in SamSrf/Models. 

In addition, there are various GUI-based tools available for specific functions:    
    1. Projecting data to surface:  SurfaceProjection
    2. Making map figures:          DisplayMaps
    3. Delineating ROIs:            DelineationTool 

Moreover, by default the anatomical meshes for a reconstruction (recon-all)
are now kept separate from the functional data in the folder ../anatomy/
You can use samsrf_expand_srf and samsrf_compress_srf to load and remove 
these fields from the Srf structure.

The pRF fitting process has considerable differences to that used in previous
versions of SamSrf. Please refer to the cookbook for more information.

##### IMPORTANT NOTE: VERSIONS PRIOR TO 6.0 ARE NOT PROPERLY COMPATIBLE & MAY NOT WORK!

You may be able to use data analysed with previous versions (at least 5.0 upwards) by converting the SamSrf data files using samsrf_convert_old_srf.m

## LATEST UPDATES

###Version 6.01 (11-08-2018) 
Some bug fixes where switch didn't work as intended. (DSS)

### Version 6.02 (13-08-2018) 
Fixed bug with default thresholding in DisplayMaps when loading a new map. (DSS)  

### Version 6.03 (13-08-2018) 
Fixed minor bug with samsrf_anatomy_srf. (DSS)

### Version 6.04 (12-09-2018) 
Fixed minor bug with SurfaceProjection tool (DSS)  
Added time stamps in help sections that missed them. (DSS)

### Version 6.05 (16-09-2018) 
Fixed bug samsrf_colourcode & eccentricity wheel image - fovea isn't purple! (DSS)

### Version 6.12 (01-12-2018) 

**WARNING: Versions 6.1 & 6.11 had a bug with HRF convolution! DO NOT USE!**

Added some more info to help sections in the GUI tools (DSS)  
Added option to replace bad slow fits with coarse fit parameter estimates (DSS)  
Added option to calculate noise ceiling during samsrf_vol2srf (DSS)  
Added function for temporal smoothing (low-pass filtering) (DSS)  
Added functions for calculating single-subject t-tests & correlations (DSS)  
Made modifications & additions to samsrf_backproj_srclt (SuSt)  
Added colour maps for comparing maps with Schira or Benson studies (DSS)  
Fixed minor bug with defaults in SurfaceProjection tool (DSS)  
Fixed show-stopping bug with noise ceiling in samsrf_expand_srf (DSS)  
Fixed problem when using anatomical meshes on different OS platforms (DSS)  


### Version 6.13 (02-12-2018) 
Added back the option to use smoothing in the coarse-fitting stage (DSS)  
Added option to only run the coarse-fit instead of the full fit (DSS)  
Initial commit to GitHub (IA)

---  
                             Questions/Comments? 
                         Contact Sam Schwarzkopf at:
                         s.schwarzkopf@auckland.ac.nz
