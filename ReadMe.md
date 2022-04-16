# SamSrf VIII - Read Me
Version 8.0 **ALPHA VERSION - USE AT YOUR OWN PERIL!** (16-04-2022)

This major release includes the most recent updates which involved algorithms for 
fitting population receptive fields or connective fields using reverse correlation 
combined with posthoc fitting of correlation profiles. Importantly, it also introduces 
several new fitting algorithms (see Cookbook & ModelHelp for details).  
By default we however still use the standard Nelder-Mead algorithm & you can run the 
exact same analyses as in SamSrf 7.

SamSrf 8 was tested on Matlab R2020a. The Nelder-Mead algorithm requires Matlab's 
Optimization toolbox. The Hooke-Jeeves algorithm is implemented directly in SamSrf.
SamSrf strongly relies on parallel computing for a number of time-intensive analyses, 
so if you have Matlab's Parallel Computing Toolbox installed & you have a 
multi-core computer or cluster it should run faster. You will need to modify the code 
yourself to use the toolbox without parallel computing.    
 
There are various GUI-based tools available for specific functions:  
    1. Projecting data to surface:  SurfaceProjection  
    2. Making map figures:          DisplayMaps  
    3. Delineating ROIs:            DelineationTool  
	4. Automatic ROI delineation:	AutoDelineation
    5. Viewing stimulus apertures:  ViewApertures
	6. Help on model parameters:	ModelHelp

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

### Version 7.99 (17-04-2022)  
- Added ModelHelp tool to document the various Model parameters (DSS)  
- Defaults for optional parameters are now set by a separate function (DSS)  
- New fast fit option to average top predictions in coarse fit (DSS)  
- Added Hooke-Jeeves algorithm as alternative for pRF & CF model fitting (DSS)  
- Now possible to define parameter tolerance for Nelder-Mead algorithm (DSS)  
- If both Nelder-Mead & Hooke-Jeeves are defined there is now a warning (DSS)  
- Fitting functions now perform checks on parameter definition vectors (DSS)  
- Renamed samsrf_cfcorr to samsrf_gsr2 to avoid confusion with CF analysis (DSS)
- Updated samsrf_plot to use transparent dots for scatter plots (DSS)  
- HRF fitting functions now also have option for Hooke-Jeeves algorithm (DSS)  
- Simplified Model examples by removing GUI functionality & optional parameters (DSS)  
- Updated the Cookbooks with information about the new major release (DSS)  
- Abandoned plan for SamOaSrf for Octave compatibility for the future is Python (DSS)  

## Questions/Comments?
* Contact Sam Schwarzkopf at: s.schwarzkopf@auckland.ac.nz

------
