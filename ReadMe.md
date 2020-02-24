# SamSrf VI - Read Me
Version 6.0 (10-08-2018)

Complete overhaul of the toolbox. The GUI has been abolished. Instead there are now example scripts for various pRF models in SamSrf/Models. 

In addition, there are various GUI-based tools available for specific functions:  
    1. Projecting data to surface:  SurfaceProjection  
    2. Making map figures:          DisplayMaps  
    3. Delineating ROIs:            DelineationTool  

Moreover, by default the anatomical meshes for a reconstruction (recon-all) are now kept separate from the functional data in the folder ../anatomy/  
You can use samsrf_expand_srf and samsrf_compress_srf to load and remove these fields from the Srf structure.  

The pRF fitting process has considerable differences to that used in previous versions of SamSrf. Please refer to the cookbook for more information.  

### IMPORTANT NOTE: VERSIONS PRIOR TO 6.0 ARE NOT PROPERLY COMPATIBLE & MAY NOT WORK!  

You can however use data analysed with previous versions (at least 5.0 upwards) by converting the SamSrf data files using samsrf_convert_old_srf.m  

Even if you used SamSrf 5 for your mapping analysis we recommend you use SamSrf 6 for your analysis after the model fitting.  

## LATEST UPDATES 

### Version 6.22 (24-02-2020)
Fixed incorrect year in date for version 6.21... *facepalm* (DSS)  
Updated samsrf_label2nii to allow native Matlab NII writer (IA)  
Ensured that prf_contour plots square axes (DSS)  

### Version 6.21 (21-02-2020)
Added option to use a seed map to samsrf_fit_prf (DSS)  
Added model parameter to threshold what goes into fine fit (DSS)  
Added prf_contour back as it had been missing mysteriously since v6.04... (DSS)  
Added square tesselation to searchlight backprojection (SuSt)  
Fixed incorrect assignment of default inputs in samsrf_backproj_prf (DSS)  
Fixed bug with display functionality in prf_predict_timecourse (DSS)  
Fixed bug with loading paths when no file is defined (DSS)  
Fixed bug with samsrf_heatmap_del not having square axes (DSS)  
Started work on SamOaSrf for Octave support... (DSS)  

### Version 6.20 (15-08-2019)
Complete overhaul of samsrf_surf colouring code & added transparency option (DSS)    
Path colour in samsrf_surf can now be defined & defaults to opposite polarity (DSS)    
Overhaul of DisplayMaps GUI & added path colouring toggle button (DSS)  
Added option for uniform transparency in samsrf_surf (DSS)  
DisplayMaps by default now restricts view to the ROI in the Srf (DSS)  
Delaunay backprojection procedure were missing from previous versions! (DSS)  

### Version 6.19 (18-06-2019)
Changed thresholding in samsrf_surf so values below threshold removed from map (DSS)  

### Version 6.18 (07-06-2019)
Fixed bug where noise ceiling was errorenously squared in samsrf_vol2srf (DSS)  
Added normalised goodness-of-fit (relative to noise ceiling) to samsrf_plot (DSS)  

### Version 6.17 (28-05-2019)
Inverted the chronological order of updates in this Read Me file (DSS)  
Fixed thresholding bug in smoothing functions when using already smoothed data (DSS)  

### Version 6.16 (22-05-2019) 
Added Delaunay backprojection procedure samsrf_backproj_del (DSS)  
Consolidated some outputs of samsrf_backproj_srclt into one (SuSt)  
Output names in samsrf_backproj_srclt & samsrf_backproj_prf now more descriptive (DSS)  
Made samsrf_backproj_prf eccentricity range consistent with other methods (DSS)  
Fixed bug with samsrf_cortmagn when no ROI is used (DSS)  
Cosmetic changes to colour map files (DSS)  
Added NIFTI header info to volumetric Srf structure in samsrf_vol2mat (IA)  

### Version 6.15 (14-12-2018)
Cosmetic change to legend in samsrf_heatmap (SuSt)  

### Version 6.14 (10-12-2018)  
Fixed minor but stupid bug with no default for coarse-fit only fitting (DSS)  

### Version 6.13 (02-12-2018)
Added back the option to use smoothing in the coarse-fitting stage (DSS)  
Added option to only run the coarse-fit instead of the full fit (DSS)  
Initial commit to GitHub (IA)  

### Version 6.12 (01-12-2018)
#### WARNING: Versions 6.1 & 6.11 had a bug with HRF convolution! DO NOT USE!  
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

### Version 6.05 (16-09-2018) 
Fixed bug samsrf_colourcode & eccentricity wheel image - fovea isn't purple! (DSS)  

### Version 6.04 (12-09-2018) 
Fixed minor bug with SurfaceProjection tool (DSS)  
Added time stamps in help sections that missed them (DSS)  

### Version 6.03 (13-08-2018)
Fixed minor bug with samsrf_anatomy_srf (DSS)  

### Version 6.02 (13-08-2018)
Fixed bug with default thresholding in DisplayMaps when loading a new map. (DSS)  

### Version 6.01 (11-08-2018)
Some bug fixes where switch didn't work as intended (DSS)  

## Questions/Comments?
* Contact Sam Schwarzkopf at:
* s.schwarzkopf@auckland.ac.nz

------
