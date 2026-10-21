# SamSrf X - Read Me

This major release includes a standalone application for integrating into 
NeuroDesk & eventually browser-based analysis platforms. It includes the 
most recent updates which involved a better & faster algorithm for fitting 
population receptive fields, support for microtime resolution (i.e. smoothly 
moving stimuli), fitting the HRF as part of a pRF model & modelling 
event-related modulations of pRF parameters.  

All these are fundamental changes from previous SamSrf versions & some 
functions may therefore not be compatible with much earlier versions. 
The toolbox also includes reverse correlation pRF & CF analysis. 

The main standalone GUI is simply called *SamSrfX*. It allows you to 
select between different algorithms, and between Models for various common uses, 
and you can display data, delineate maps, etc.  

The toolbox also includes a standalone command line app called *SamSrfAnalysis*.
This is specifically designed for use with external data & doesn't require Matlab.
So you can read in GII or NII files & the analysis outputs maps in that format.

Within the Matlab environment, you can also call the SamSrfX GUI. But there 
are also various tools available for specific purposes:  
 
    1. Making map figures:                      DisplayMaps*  
    2. Spatial normalisation/nativization:      Native2TemplateMap / Template2NativeMap  
    3. Viewing stimulus apertures:              ViewApertures*  
    4. New field sign based delineation:        FastDelineationTool*
    5. Inspecting M/EEG scalp maps:             DisplayScalpMaps   
	6. Help on model parameters:                ModelHelp*  
    7. Help on all SamSrf functions:            FuncHelp  
    8. Projecting data to surface:              SurfaceProjection  

(* These tools are also available as separate standalone apps)

## DIFFERENCES TO EARLIER VERSIONS

### IMPORTANT: FORWARD-MODEL FITS FROM VERSIONS PRIOR TO 9.0 ARE NOT COMPARABLE!  

Even if you used earlier versions for your analysis we recommend you use SamSrf X  
for your analysis *after* the model fitting (e.g. plotting, quantification & statistics). 
Data files from SamSrf 6 should already be in the same format. You can also convert old 
data files from older versions (at least 5.0 upwards) using samsrf_convert_old_srf.m   

### Retired old delineation tools!

Since version 10.4 we have officially retired our previous suite of delineation functions, 
specifically *AutoDelineation* and *DelineationTool*. The code is still available in 
the toolbox but they are **no longer supported**. Both of these (esp. DelineationTool) 
were a massive hackjob, but MathWorks have pretty much broken them with their more 
recent Matlab versions. But they are also obsolete. We now have *FastDelineationTool*
that allows you to quickly labelled retinotopic ROIs using the "Poor Man's Field Sign".
Morever, for manual delineation we now highly recommend Fernanda Ribeiro's online 
delineation tool available at: https://webapps.neurodesk.org/surfannotate/

------

SamSrf X was tested on **Matlab R2023b**. The Nelder-Mead algorithm requires Matlab's 
**Optimization Toolbox**. The Hooke-Jeeves algorithm is implemented directly 
in SamSrf. SamSrf strongly relies on parallel computing for a number of time-intensive 
analyses, so if you have Matlab's **Parallel Computing Toolbox** installed & you have a 
**multi-core computer or cluster** it should run faster. You will need to modify the 
code yourself to use the toolbox without parallel computing by replacing parfor calls 
with for calls.    
 
When using the surface project functions in Matlab, by default the anatomical meshes 
for a reconstruction (recon-all) are kept separate from  the functional data in the 
folder ../anatomy/  
You can use samsrf_expand_srf and samsrf_compress_srf to load and remove these fields 
from the Srf structure.  
This is not relevant for the SamSrf GUI because most data will simply be kept in 
GII or NII format & only few MAT files are saved containing anatomical data.  

## LATEST UPDATES 

### Version 10.4 (22-10-2026) - FastDelineationTool Official Release!
- Replaced samsrf_pred_rois with *FastDelineationTool*! (DSS)  
- This tool has now also been integrated into the SamSrfX GUI (DSS)  
- Added options for adding new & duplicate ROIs in FastDelineationTool (DSS)  
- Moreover you can now also plug holes in specific ROIs (DSS)  
- When saving delineated ROIs, a log of postprocessing steps is saved (DSS)  
- Template2NativeMap now supports negative ROI numbers for fast delineation (DSS)  
- Fixed bug with how ROIs are saved when field sign is present in template (DSS)
- ROI names in template maps are now retained in Srf warped into native space (DSS)  
- Added Poor Man's Field Sign algorithm to samsrf_fieldsign (DSS)  
- Changed samsrf_fieldsign algorithm to boost speed but also use 3D anatomy (DSS)  
- Added function to split bilateral labels into hemispheres (DSS)  
- Fixed bug with bilateral Srf functions based on warped template brains (DSS)  
- prf_contour now also supports micro-probing & vertex-wise HRFs (DSS)  
- samsrf_surfcalcs now does 2nd pass of eccentricity filtering at the end (DSS)  
- samsrf_fieldsign can now also output gradient angle instead, the new default (DSS)  
- samsrf_export_giis now lets you specify which maps to export (DSS)  
- Fixed issue with missing colour maps in backprojection functions (DSS)  
- Smoothing functions now report more information (DSS)  
- No longer saves empty bilateral labels when neither hemisphere's label exist (DSS)  
- Reverse correlation pRF with convex hull algorithm now saves summary stats also (DSS)  
- You can now show ROI labels on surface with samsrf_showlabels & samsrf_hidelabels (DSS)  
- Added warning about decommissioned DelineationTool and AutoDelineation (DSS)  
- Bugfix for Template2NativeMap when no ROI names in template file (DSS)  
- Added grid lines in samsrf_plot (DSS)  

### Version 10.36 (31-08-2026)  
- Removed option to seed forward model fine-fit with file - now needs a matrix! (DSS)  
- Added micro-probing algorithm as alternative to samsrf_revcor_prf! (DSS)  
- Added option in samsrf_fit_prf to take absolute values of pRF parameters (DSS)  
- Fitting functions now remove any raw data from input data files (DSS)  
- pRF forward modelling can now use vertex-specific HRFs (DSS)  
- samsrf_fitvsobs now also plots the HRF used (DSS)  
- Region predictor now includes R^2 thresholding (DSS)  
- Fixed bugs with transparency in samsrf_surf (DSS)  
- Split dog_fwhm from samsrf_dog_fwhm as a utility function (DSS)
- Fixed bug in samsrf_surf with visual CF profile colour scheme (DSS)  
- Turned parfor loop in prf_generate_searchspace back on... (DSS)  
- Vertex-inspector for model-based pRFs no longer used fixed colour scale (DSS)  

### Version 10.335 (25-06-2026) 
- DisplayMaps no longer needs Statistics Toolbox for thresholding (DSS)  
- Fixed bug in samsrf_gsr2 when providing HRF vector (DSS)  
- Analytically inconsequential fixes to help texts (DSS)  
- Added option to pre-smooth maps in region predictor function (DSS) 
- Region predictor now has option to threshold by confidence (DSS)  
- Predicted ROI map is now smoothed at the end (DSS)   
- Renamed region predictor to samsrf_pred_rois (DSS)  
- Made plotting in samsrf_gsr2 optional (DSS)  

### Version 10.33 (12-05-2026) 
- Created LabelRegions function to automatically delineate ROIs (DSS)  
- Template2NativeMaps now also transforms multi-subject matrices (DSS)  
- Template2NativeMaps now saves template maps with their name (DSS)  
- samsrf_fieldsign can now also use Euclidean spherical distance (DSS) 
- Smoothing functions can now also use vertex list as ROI (DSS)   

### Version 10.3 (13-10-2025) - Pre-release Beta  
- Added standalone app SamSrfAnalysis for use in NeuroDesk or system terminal etc (DSS)  
- Can now provide apertures as volumetric NII files (DSS)  
- Model specification for algorithms can now be provided as JSON file (DSS)  
- Added sumnorm pRF model to fit complex shape with only 1 extra parameter (DSS)  
- The sumnorm model has also been implemented in SamSrfX GUI now (DSS)  
- More failed attempts to fix issue with DelineationTool GUI due to Matlab 2024 bug (DSS)  
- Fixed various compatibility bugs with AutoDelineation function (DSS)  
- Removed warning in samsrf_surf if default view undefined (DSS)  
- Fixed bug with right hemisphere labels of bilateral surfaces in DisplayMaps (DSS)  
- Can now read in GII files without anatomy or by providing a surface GII (DSS)  
- Added warnings about out-of-date functions for surface projection (DSS)  
- Some small fixes with file outputs in command line app (DSS)  

### Version 10.201 (29-09-2025)  
- Critical bug fix in functions visualising HRF fits! (DSS)  
- Bugfix for missing colour maps & interpolation errors in samsrf_cmap (DSS)  
- Fixed issues with DelineationTool GUI for Matlab versions post R2023b (DSS)  

### Version 10.1 (24-06-2025)  
- *Added BSD 3 license on GitHub*
- Template warping functions now support GII registrations like from HCP (DSS)  
- Hacked samsrf_expand_srf to deal with awkward Srfs converted from HCP GIIs (DSS)  
- Fixed bug with DisplayMaps & DelineationTool not loading ROIs correctly (DSS)  
- Fixed inconsequential typo in DelineationTool when flood filling fails (DSS)  
- Minor fix in AnonymiseSrf removing redundant input parameter (DSS)  
- HRF parameter R/U is now called AmpRat to allow exporting (DSS)  

### Version 10.01 (05-05-2025)  
- Changed how samsrf_plothrf outputs the average HRF fit (DSS)  
- Minor bugfix in samsrf_gsr2 function when using SPM or concurrent HRF (DSS)  
- Fixed bug with samsrf_simulate_prfs when using ground truth map (DSS)  
- Also reorganised inputs for samsrf_simulate_prfs for modern age (DSS)  
- Fixed minor typo in SamSrfX GUI (DSS)  
- Bugfix for ViewApertures with incorrect call to colour scheme (DSS)  
- Fixed bug with samsrf_colourcode failing with JSON defaults (DSS)  
- Minor bug fix with bilateral hemispheres of reverse correlation pRF maps (DSS)  
- Bugfix for template map transformation when files are too large (DSS)  
- Fixed issues with colour maps in samsrf_surf & DelineationTool (DSS)  
- Made handling of colour maps in samsrf_simvsfit consistent with new version (DSS)  
- Fixed how time is displayed on X-axis in samsrf_fitvsobs (DSS)  

### Version 10 (21-10-2024) - *SamSrf X Official Release!*
- *Default HRF is now SPM canonical but check what works best for your needs! (DSS)*   
- Created new GUI called SamSrfX for running analysis & inspecting results (DSS)  
- Modifications needed for standalone GUI app & integration into NeuroDesk (DSS)  
- When using GUI, information is now displayed there instead of command window (DSS)  
- Biophysical model now uses unsigned (absolute) pRF profile to determine response (DSS)  
- Removed option to concatenate data files in pRF/CF fitting functions! (DSS)  
- Backprojection of apertures for pRF-from-CF analysis now internal to samsrf_fit_prf (DSS)  
- Added functions for saving & loading Model specifications as JSON files (DSS)  
- Updated examples in SamSrf/Models to align with saved versions for standalone app (DSS)  
- Model fitting functions now allow direct input of Srf instead of filenames (DSS)   
- ViewApertures now also allow direct input argument (DSS)  
- All text output & progress bar are now using the GUI if available (DSS)  
- Benson map projection is now using GII to align with other analyses (DSS)  
- SamSrf GUI now also can call DisplayMaps tool (DSS)  
- Replaced colour map functions with editable CSV files (DSS)  
- Replaced SamSrf_defaults.mat with editable JSON file (DSS)  
- Added batch analysis & option to replace subject IDs in SamSrfX (DSS)  
- SamSrfX can now show basic explanation of different algorithms (DSS)  
- HRF fit analysis can now analyse lists of vertices (DSS)  
- Fixed minor bug where wrong number of volumes was shown during pRF fit (DSS)  

### Version 9.9 (08-08-2024)  
- Added alternative functionality for multi-condition designs in pRF models (DSS)   
- VectoriseApertures & ViewApertures now support multi-condition designs (DSS)  
- Updated ModelHelp with info about multi-condition designs (DSS)  
- Added function for exporting GIfTI files of data fields (DSS)  
- Updated label exporting function for modern age (DSS)  
- Cluster ROI selection now has options for more flexibility (DSS)  
- Field sign function now takes vertex list as ROI input (DSS)  

------

## Questions/Comments?
* Contact Sam Schwarzkopf at: s.schwarzkopf@auckland.ac.nz

------
