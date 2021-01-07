# SamSrf VII - Read Me
Version 7.0 (20-07-2020)

Complete overhaul of the toolbox. The GUI from SamSrf <=5 has been abolished.  
Instead there are now example scripts for various pRF models in SamSrf/Models.  

In addition, there are various GUI-based tools available for specific functions:  
    1. Projecting data to surface:  SurfaceProjection  
    2. Making map figures:          DisplayMaps  
    3. Delineating ROIs:            DelineationTool  
    4. Viewing stimulus apertures:  ViewApertures

SamSrf now supports parallel computing for a number of time-intensive analyses,  
so if you have Matlab's Parallel Computing Toolbox installed and you have a  
multi-core computer or cluster it should run faster. Many stages of the analysis  
are now geared towards parallel computing, so they may be quite slow without it.  

By default the anatomical meshes for a reconstruction (recon-all) are kept  
separate from the functional data in the folder ../anatomy/  
You can use samsrf_expand_srf and samsrf_compress_srf to load and remove  
these fields from the Srf structure.  

### IMPORTANT NOTE: MODEL FITS FROM VERSIONS PRIOR TO 7.0 ARE NOT COMPARABLE!  
The pRF fitting process has considerable differences to previous versions of SamSrf.  

Even if you used earlier versions for your analysis we recommend you use SamSrf 7  
for your analysis -after- the model fitting. Data files from SamSrf 6 should already  
be in the same format. You can also convert old data files from older versions  
(at least 5.0 upwards) using samsrf_convert_old_srf.m   

## LATEST UPDATES 

### Version 7.073 (07-01-2021)  
- DelineationTool now uses default polar & eccentricity colour scheme (DSS)  
- Bugfix for when progress reports are turned off in older Matlab versions (DSS)  
- Fixed bug with expanding that restored file version from anatomy meshes (DSS)  
- Fixed bug when smoothing concatenating runs (DSS)  

### Version 7.07 (21-12-2020)  
- Parallel processing progress reports now turned off on older Matlab versions (DSS)  
- Fixed bug with filenames in samsrf_bilat_label (DSS)  
- Added default camera angles for rendering bilateral data files (DSS)  
- When combining hemispheres meshes field is now automatically removed (DSS)   
- Changed output in samsrf_wedgeplot so it contains NaN for bad data (DSS)  
- Improved wedge plots to show actual segments (EA & DSS)   

### Version 7.06 (29-10-2020)  
- It is now possible to combine Srfs from both hemipheres (DSS)  
- Added function for creating a binned plot for wedge segments (DSS)  
- Limited ROI in DisplayMaps can now be turned off completely (DSS)  
- Added support for parallel computing of CMF (DSS)  
- Surface projection function now also accepts ASC files (DSS)  
- Conversion from old Srf files now also accepts ASC files (DSS)  
- DisplayMaps & DelineationTool now state which map is being loaded (DSS)  
- Bins in samsrf_plot now restrict range for scatter plot (DSS)  
- Added option to samsrf_denoisemap to remove artifactually small pRFs (DSS)  

### Version 7.05 (07-08-2020)  
- Added tool for warping average template maps back into native space (DSS)  
- ViewApertures now uses greyscale if no negative pixels are present (DSS)  

### Version 7.04 (03-08-2020)  
- Fixed bug in samsrf_fminsearch_loop by which it could get stuck with poor fits (DSS)  
- For consistency, removed external dependencies from samsrf_backproj_srclt (DSS)  
- Rewrote analysis loops in revcor functions but parallel support not ready (DSS)  

### Version 7.03 (24-07-2020)  
- Expanded support for parallel processing to samsrf_geomatrix (DSS)  
- Cosmetic changes to command window output in samsrf_fit_prf (DSS)  
- Added shell scripts for MGH surface projection & Benson map conversion (DSS)  
- Removed fs_bin2asc.sh shell script because out-of-date & redundant (DSS)  

### Version 7.02 (23-07-2020)  
- Fixed minor errors in example simulation scripts (DSS)  
- Removed erroneous but inconsequential duplicate command from example models (DSS)  
- Removed denoise checkbox from DisplayMaps - do this in postprocessing instead (DSS)  

### Version 7.01 (20-07-2020)  
#### IMPORTANT UPDATE: Improved precision of model fits! (DSS)  
- Responses are now modelled as percent change relative to pRF volume! (DSS)  
- Implemented parallel computing in several intensive analyses (DSS)  
- Fixed bugs with storing predicted time courses from fine fit (DSS)  
- Added option to change the size of coarse fit vertex blocks (DSS)  
- Added function for computing normalised goodness-of-fit (DSS)  
- Added options to render noise ceiling in DisplayMaps (DSS)  
- Removed all progress bars, now reported in command window instead (DSS)  
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
- Added option for mean weighted by distance to samsrf_backproj_srclt (DSS)  
- Added checks for mismatch between model definition & pre-saved search spaces (DSS)  
- Moved coarse fit smoothing later in the model fitting pipeline (DSS)  
- Inconsequential reorganisation of steps after fine-fit in samsrf_fit_prf (DSS)  
- Changed default camera position for right hemisphere (DSS)  
- Simplified support for multi-subject Srf data in DisplayMaps (DSS)  
- Added function for Beta denoising & checkbox in DisplayMaps (DSS)  
- Added default option for restricting DisplayMaps & DelineationTool to a ROI (DSS)  
- Default display is now based on occipital mesh coordinates instead of ROI (DSS)  
- DelineationTool now defaults to activation loading if no field sign exists (DSS)  
- Added figure handle output to samsrf_plot for legend use (DSS)   

## Questions/Comments?
* Contact Sam Schwarzkopf at: s.schwarzkopf@auckland.ac.nz

------
