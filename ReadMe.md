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

### Version 7.811 (09-04-2022)  
- Updated samsrf_plot to use transparent dots for scatter plots (DSS)  
- Abandoned plan for SamOaSrf for Octave compatibility - the future is Python (DSS)  

### Version 7.81 (08-04-2022)  
- Improved pRF fitting to reverse-correlation CF profiles via thresholding (DSS)  
- Fitting pRF to reverse-correlation CF now uses nested search for fitting (DSS)  

### Version 7.8 (17-03-2022)  
#### CRITICAL BUGFIX: Recent updates broke reverse-correlation pRF fitting!  
**Do not use reverse-correlation pRF parameters fit between v7.4 & v7.71!**  
- Fixed bug with inversion of matrix in prf_contour (DSS)  

### Version 7.71 (15-03-2022)  
- Added function to add paths with different colours to a map (DSS)  
- Fixed issue with DelineationTool being unable to load activation maps (DSS)  
- If no field sign map, DelineationTool now loads the *last* row in Srf.Data (DSS)  
- Fixed bug with DelineationTool when def_disproi is empty (DSS)  

### Version 7.7 (14-03-2022)  
- Addressed issues loading files with too extensive Matlab paths (DSS)  
- Template2NativeMap now requires pathname for template (DSS)  
- Added function for making anonymised SamSrf data files for publishing (DSS)  
- Template2NativeMap now includes the curvature from template file (DSS)  
- Fixed bug with samsrf_surf when displaying white-matter surface (DSS)  
- Reverse correlation now has option to strip pRF profiles to save disc space (DSS)  

### Version 7.62 (02-03-2022)  
- Inspection tool samsrf_cfcorr now supports 1D models & polar search grids (DSS)  
- Fixed several bugs in samsrf_cfcorr (DSS)  
- Fixed cosmetic bug in samsrf_simvsfit (DSS)  

### Version 7.61 (23-02-2022)  
- CF reverse correlation now defaults to using global signal correction (DSS)  
- CF correlation profiles are now also computed using parallel processing (DSS)  
- Fixed minor bug with samsrf_benson2srf (DSS)  

### Version 7.6 (16-02-2022)  
- By default, reverse correlation CF now fits pRFs to template pRF coordinates (DSS)  
- pRF model functions now accept pRF coordinates as input for CF fitting (DSS)  
- CF reverse correlation profiles are no longer saved in data file by default (DSS)  
- Inconsequential changes to help sections & command line reports (DSS)  

### Version 7.54 (19-12-2021)  
- Option to turn off automatic global regressor in samsrf_removenoise (DSS)  
- Removing irrelevant commented sections from samsrf_revcor_cf (DSS)  
- Updated Delineation Tutorial to include AutoDelineation tool (DSS)  

### Version 7.53 (05-11-2021)  
- Parameter estimation in samsrf_revcor_cf now uses parallel computing (DSS)  
- Added progress bar to samsrf_backproj_prf (DSS)  
- Fixed mislading typo in help section to samsrf_backproj_sctr (DSS)  
- Further cosmetic changes to DisplayMaps tool (DSS)  
- Cloning camera now also clones axis limits (DSS)  
- Added function to convert ROI vertex list into a mask (DSS)  

### Version 7.521 (13-10-2021)  
- Major overhaul of DisplayMaps user interface, including vertex inspector (DSS)  
- New option in samsrf_surf to display data vector directly (DSS)  
- AutoDelineation now also uses default ROI list if defined (DSS)  
- Added new skewed Gaussian pRF model (DSS)  
- New option in samsrf_showprf & prf_contour supporting pRFs from model parameters (DSS)  
- Fixed bug in samsrf_showprf when matrix has no variance (DSS)  
- Vertex inspector in DisplayMaps now also supports forward-model & GLM fits (DSS) 
- DisplayMaps vertex inspector for CFs now shows whole seed ROI more clearly (DSS)  
- Vertex inspector is now semi-compatible with pRF maps from versions prior to 6 (DSS)  
- Fixed minor bug with samsrf_surf when called directly (DSS) 
- Added progress reports to CF reverse correlation paramater estimation (DSS)  
- Various cosmetic changes (DSS)  

### Version 7.51 (07-10-2021)  
- Added ROI list default option to make DelineationTool more flexible (SuS&DSS)  

### Version 7.5 (23-09-2021)  
#### BETA version! Does not seem to break conventional pRF fits but use with care!
- Forward model can now downsample predictors for high temporal resolution stimuli (DSS)  
- Simulation suite can also use high temporal resolution stimuli now (DSS)  
- Reverse correlation pRFs don't support downsampling yet - maybe later (DSS)  

### Version 7.41 (20-09-2021)  
- Changed default colour map in samsrf_backproj_prf to berlin (DSS)  
- Fixed bug with default input arguments in samsrf_backproj_prf (DSS)  
- Default clipping level in samsrf_backproj_prf is now 1 (DSS)  
- Fixed bug with scatter_size that made all symbols the same size (DSS)  
- Native2TemplateMaps now requires the path to fsaverage/surf (DSS)  
- AutoDelineation now supports flexible peripheral borders (DSS)  
- Added Benson auto-delineation atlas with flexible peripheral borders (DSS)  

### Version 7.4 (17-09-2021)  
- Released AutoDelineation tool but still undergoing testing & tweaking parameters (DSS)  
- DelineationTool loads auto-delineations if no manually saved paths exist (DSS)  
- DelineationTool now displays R^2 or nR^2 by default if no field sign exists (DSS)  
- AutoDelineation can now cope with bilateral Srfs (DSS)  
- Bilateral Srfs are now automatically split by DelineationTool (DSS)  
- Template maps can now contain ROI names for flexibility (DSS)  
- Added new fireice colour map similar to hotcold (DSS)  
- Fixed inconsequential bug with reporting noise ceiling threshold (DSS)  
- Updated cookbook (DSS)  
- Finally fixed issue with symbol size in scatter plots (DSS)  
- Added a utility to create custom vertex selection functions in samsrf_surf (DSS)  
- Fixed small bug with default ROI in DelineationTool (DSS)  
- Changed activation colour map & path colours in DelineationTool (DSS)  
- Fixed bug with camera when redrawing maps with samsrf_surf (DSS)  
- Changed colour scheme for displaying ROI numbers in samsrf_surf (DSS)  
- Adjusted default camera angles in samsrf_surf for sphere view (DSS)  

### Version 7.3 (12-08-2021)  
- New option to fit 2D pRF models to reverse correlation profiles (DSS)  
- Increase granularity of colour scheme in samsrf_showprf (DSS)  
- Updated example colour code images to new defaults (DSS)  

### Version 7.22 (04-08-2021)
- Forward-modelling connective field now uses 1st eigenvariates (DSS)  
- New option for polar/eccentricity-definbed patches as CFs (DSS)  
- Moved smoothing in CF fitting before search space generation (DSS)  
- Receiving patch size in connective fields can now be adjusted (DSS)  
- Added parallel computing to forward-model CF coarse fits (DSS)  
- Fixed bug with bandpass filter function when TR is not 1 second (DSS)  
- Added option to samsrf_heatmap to blend with an image (DSS)  
- Changed default colour map in ViewApertures (DSS)  
- Minor (hopefully cosmetic) change to pRF fitting function (DSS)  
- Added stand-by message in lieu of progress bar to most parfor loops (DSS)  

### Version 7.21 (09-07-2021)
- Fixed catastrophic bug when only allowing positive coarse fits (DSS)  
- Removed redundant expansion/compression from filtering functions (DSS)  
- Updated default cut-offs in DisplayMaps & removed def_eccen default (DSS)  
- Fixed minor bug with searchlight backprojection (DSS)  
- Fixed minor bug with dimensions of outputs in samsrf_visualroi (DSS)  

### Version 7.2 (30-06-2021)
- New function for splitting bilateral Srfs back into hemispheres (DSS)   
- Added option to samsrf_plot to plot mean of x-bin values instead of bin centre (DSS)
- Main analysis functions now mark more visibly when they're done (DSS)  
- Improved speed of samsrf_removenoise data cleaning function (DSS)  
- Template2NativeMaps only saves existing ROI labels now (DSS)  
- Added new colour maps from various sources & updates defaults (DSS)  
- Field sign function now uses parallel computing but is still slow (DSS)  
- Changed how samsrf_glm plots design matrix (DSS)  
- GLM contrast function can now work even if no time course data present (DSS)  
- Noise covariate regression includes data expansion/compression again (DSS)  
- Added bandpass filtering function (DSS)  
- Added new-fangled old-school progress bars (DSS)  

### Version 7.14 (22-05-2021) 
- Added forward-model for connective fields using fast-fit procedure (DSS)  
- Added half-way inflation options for pial and inflated surfaces (DSS)  
- Fixed bug with samsrf_benson2srf converter (DSS)  
- Fixed bug with smoothing connective field profiles (DSS)  
- Changed how transparency of connectiv field profiles is scaled (DSS)  

### Version 7.13 (29-04-2021)  
#### IMPORTANT BUGFIX: Fixed show-stopping bug in samsrf_fit_prf! (DSS)
- Changed samsrf_plot default mode to scatter plot (DSS)  
- Added warning in bin & wedge plot help about regression artifacts (DSS)  
- Created new comet plot function as alternative to bin plots (DSS)  

### Version 7.12 (18-04-2021) - samsrf_fit_prf broken - DO NOT USE! 
- Renamed ROI variable in Benson template projections (DSS) 
- Added fake sigma row into Benson projections & added note about newer Wang template (DSS)  
- Changed default colour schemes in DisplayMaps, including for ROI maps (DSS)  
- Visual field coverage plots now use logarithmic scale (DSS)  
- Added tool for warping native data into fsaverage template (DSS)  
- Template2NativeMaps tool now splits off anatomical meshes (DSS)  
- DisplayMaps can now visualise reverse correlation & connective fields (DSS)  

### Version 7.10 (07-04-2021) - samsrf_fit_prf broken - DO NOT USE!  
- Added parameter option for only allowing positive correlations to pass coarse fit (DSS)  
- Added default parameter to allow for logarithmic eccentricity & sigma-style maps (DSS)  

### Version 7.09 (30-03-2021)  
- Fixed major show-stopping bug with reverse correlation pRFs on bilateral surfaces (DSS)  
- Removed redundant rounding function that is not needed since at last 2016a (DSS)  
- Fixed small bug when replacing bad fits where incorrect time series was saved (DSS)  
- Template2NativeMap now automatically removes noise ceiling from native Srf (DSS)  
- Fixed bug with directly assigning path colours in samsrf_surf (DSS)  

### Version 7.08 (03-02-2021)  
- Fixed crashing bug when replacing bad fine fits with coarse fits (DSS)  
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
