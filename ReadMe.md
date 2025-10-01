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

Within the Matlab environment, you can also call the SamSrfX GUI. But there 
are also various tools available for specific purposes:  
 
    1. Projecting data to surface:              SurfaceProjection  
    2. Making map figures:                      DisplayMaps  
    3. Spatial normalisation/nativization:      Native2TemplateMap / Template2NativeMap  
	4. Automatic ROI delineation:               AutoDelineation  
    5. Delineating ROIs:                        DelineationTool  
    6. Viewing stimulus apertures:              ViewApertures  
    7. Inspecting M/EEG scalp maps:             DisplayScalpMaps   
	8. Help on model parameters:                ModelHelp  
    9. Help on all SamSrf functions:            FuncHelp  


## DIFFERENCES TO EARLIER VERSIONS

### IMPORTANT: FORWARD-MODEL FITS FROM VERSIONS PRIOR TO 9.0 ARE NOT COMPARABLE!  

Even if you used earlier versions for your analysis we recommend you use SamSrf X  
for your analysis *after* the model fitting (e.g. plotting, quantification & statistics). 
Data files from SamSrf 6 should already be in the same format. You can also convert old 
data files from older versions (at least 5.0 upwards) using samsrf_convert_old_srf.m   

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

### Version 10.203 (01-10-2025)  
- Added sumnorm pRF model to fit complex shape with only 1 extra parameter (DSS)  
- Further attempt to fix issue with DelineationTool GUI due to Matlab 2024 bug (DSS)  

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
