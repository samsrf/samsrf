# SamSrf VIII - Read Me

This major release includes the most recent updates which involved algorithms for 
fitting population receptive fields or connective fields using reverse correlation 
combined with posthoc fitting of correlation profiles. Importantly, it also introduces 
several new fitting algorithms (see Cookbook & ModelHelp for details).  
By default we however still use the standard Nelder-Mead algorithm. With one exception 
(see below) you can run the exact same analyses as in SamSrf 7.

There are various GUI-based tools available for specific functions:  
    1. Projecting data to surface:              SurfaceProjection  
    2. Making map figures:                      DisplayMaps  
    3. Spatial normalisation/nativization:      Native2TemplateMap / Template2NativeMap  
	4. Automatic ROI delineation:               AutoDelineation  
    5. Delineating ROIs:                        DelineationTool  
    6. Viewing stimulus apertures:              ViewApertures  
	7. Help on model parameters:                ModelHelp  

## DIFFERENCE TO EARLIER VERSIONS

### IMPORTANT: MODEL FITS FROM VERSIONS PRIOR TO 7.0 ARE NOT COMPARABLE!  
The pRF fitting process has considerable differences to previous versions of SamSrf.  

** pRF model fits to reverse correlation profiles are also not *identical* to SamSrf 7.**  
(Please see *help samsrf_revcor_prf* for more information about this)  

Even if you used earlier versions for your analysis we recommend you use SamSrf 8  
for your analysis *after* the model fitting (e.g. plotting, quantification & statistics). 
Data files from SamSrf 6 should already be in the same format. You can also convert old 
data files from older versions (at least 5.0 upwards) using samsrf_convert_old_srf.m   

------

SamSrf 8 was tested on **Matlab R2020a**. The Nelder-Mead algorithm requires Matlab's 
**Optimization Toolbox**. The Hooke-Jeeves algorithm is implemented directly in SamSrf.
SamSrf strongly relies on parallel computing for a number of time-intensive analyses, 
so if you have Matlab's **Parallel Computing Toolbox** installed & you have a 
**multi-core computer or cluster** it should run faster. You will need to modify 
the code yourself to use the toolbox without parallel computing.    
 
By default the anatomical meshes for a reconstruction (recon-all) are kept separate from 
the functional data in the folder ../anatomy/  
You can use samsrf_expand_srf and samsrf_compress_srf to load and remove these fields 
from the Srf structure.  

## LATEST UPDATES 

### Version 8.2 (14-05-2022)  
- Now includes option to read GIfTI tiles (DSS)  

### Version 8.1 *Beta* (06-05-2022)  
- CF parameters can now be estimated with convex hull algorithm (DSS)  
- Convex hull algorithm is now the default for CF reverse correlation (DSS)  
- Fixed inconsequential error in CF reverse correlation progress reports (DSS)  

### Version 8.0 *Beta* (20-04-2022)  
**WARNING: Beta version - some issues may remain despite numerous tests!**  
*Please notify us about any unknown discrepancies with results from SamSrf 7*  
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
