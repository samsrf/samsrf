# SamSrf VII - Read Me
Version 7.0 (??-??-2020)

Complete overhaul of the toolbox. The GUI from SamSrf <=5 has been abolished.  
Instead there are now example scripts for various pRF models in SamSrf/Models.  

In addition, there are various GUI-based tools available for specific functions:  
    1. Projecting data to surface:  SurfaceProjection  
    2. Making map figures:          DisplayMaps  
    3. Delineating ROIs:            DelineationTool  
    4. Viewing stimulus apertures:  ViewApertures

SamSrf now supports parallel computing for a number of time-intensive analyses,  
so if you have Matlab's Parallel Computing Toolbox installed and you have a  
multi-core computer or cluster it should run faster.   

By default the anatomical meshes for a reconstruction (recon-all) are kept separate from the functional data in the folder ../anatomy/  
You can use samsrf_expand_srf and samsrf_compress_srf to load and remove these fields from the Srf structure.  

### IMPORTANT NOTE: MODEL FITS FROM VERSIONS PRIOR TO 7.0 ARE NOT COMPARABLE!  
The pRF fitting process has considerable differences to that used in previous versions of SamSrf.  
Previous versions would revert to the search grid in some situations.  
Please refer to the cookbook for more information.  

Even if you used earlier versions for your mapping analysis we recommend you use SamSrf 7 for your analysis -after- the model fitting.  
Data files from SamSrf 6 should already be in the same format. You can also use data files from older versions (at least 5.0 upwards)  
by converting the SamSrf data files using samsrf_convert_old_srf.m  

## LATEST UPDATES 

### Version 6.994 (29-06-2020)  
#### IMPORTANT UPDATE: Improved precision of model fits! (DSS)  
- Responses are now modelled as percent change relative to pRF volume! (DSS)  
- Implemented parallel computing in several intensive analyses (DSS)  
- Removed all progress bars as they don't work in parallel computing (DSS)  
- Fits beyond twice the X & Y range are now automatically removed (DSS)  
- Updated example Model scripts with finer search grids (DSS)  
- Added several new functions for ground truth simulations (DSS)  
- Lowered default threshold in display & delineation tools (DSS)  
- Changed default eccentricity map colour scheme to inverted jet (DSS)  
- Added option of wider paths in samsrf_surf in SamSrf_defaults (DSS)  
- Added support for inverted colour schemes in samsrf_surf (DSS)  
- Changed default transparency level in samsrf_surf (DSS)  
- Fixed bug in ViewApertures tool that inverted the Y-axis (DSS)  
- Fixed bug in ViewApertures tool that gave a warning on empty frames (DSS)  
- Fixed bug in ViewApertures tool where colour scheme was inconsistent across frames (DSS)  
- Fixed inconsequential duplicate lines in samsrf_benson2srf (DSS)  
- Added function for backprojecting reverse correlation profiles (DSS)  
- Inconsequential reorganisation of steps after fine-fit in samsrf_fit_prf (DSS)  
- Added figure handle output to samsrf_plot for legend use (DSS)   

## Questions/Comments?
* Contact Sam Schwarzkopf at: s.schwarzkopf@auckland.ac.nz

------
