# SamSrf IX - Read Me

This major release includes the most recent updates which involved an improved 
(much faster!) algorithm for fitting population receptive fields. This constitutes 
a fundamental changefrom previous SamSrf versions. The toolbox also includes 
reverse correlation pRF analysis & connective fields analysis. 

There are various GUI-based tools available for specific functions:  
    1. Projecting data to surface:              SurfaceProjection  
    2. Making map figures:                      DisplayMaps  
    3. Spatial normalisation/nativization:      Native2TemplateMap / Template2NativeMap  
	4. Automatic ROI delineation:               AutoDelineation  
    5. Delineating ROIs:                        DelineationTool  
    6. Viewing stimulus apertures:              ViewApertures  
	7. Help on model parameters:                ModelHelp  

## DIFFERENCES TO EARLIER VERSIONS

### IMPORTANT: FORWARD-MODEL FITS FROM VERSIONS PRIOR TO 9.0 ARE NOT COMPARABLE!  
The fitting process has considerable differences to previous versions: We changed 
the way apertures are encoded & how the time course is predicted. While the results
will typically be close to the old approach, there will be differences. The new algorithm
is also more flexible. Most importantly, it is considerably faster. Given the advantages,
we decided not to include the old algorithm in this version. If you need that for some
reason, you will need to use SamSrf 8.4.

** pRF model fits to reverse correlation profiles are also not *identical* to SamSrf 7.**  
(Please see *help samsrf_revcor_prf* for more information about this discrepancy.)  

Even if you used earlier versions for your analysis we recommend you use SamSrf 9  
for your analysis *after* the model fitting (e.g. plotting, quantification & statistics). 
Data files from SamSrf 6 should already be in the same format. You can also convert old 
data files from older versions (at least 5.0 upwards) using samsrf_convert_old_srf.m   

------

SamSrf 9 was tested on **Matlab R2020a & R2020b**. The Nelder-Mead algorithm requires 
Matlab's **Optimization Toolbox**. The Hooke-Jeeves algorithm is implemented directly 
in SamSrf. SamSrf strongly relies on parallel computing for a number of time-intensive 
analyses, so if you have Matlab's **Parallel Computing Toolbox** installed & you have a 
**multi-core computer or cluster** it should run faster. You will need to modify the 
code yourself to use the toolbox without parallel computing by replacing parfor calls 
with for calls.    
 
By default the anatomical meshes for a reconstruction (recon-all) are kept separate from 
the functional data in the folder ../anatomy/  
You can use samsrf_expand_srf and samsrf_compress_srf to load and remove these fields 
from the Srf structure.  

## LATEST UPDATES 

### Version 9.034 (29-07-2022)  
- Added option to samsrf_glm to turn off automatic global covariates (DSS)  
- Turned off parallel computing in samsrf_glm but commented out option remains (DSS)  
- Now saves GLM files in v7.3 file format in case they are too large (DSS)  
- Added option to define scaling factor/eccentricity in samsrf_simvsfit (DSS)  
- Corrected incorrect help section in template warping functions (DSS)  
- Changed samsrf_cometplot to use logarithmic density scale (DSS)  
- Increased default granularity in samsrf_cometplot (DSS)  
- Turned off transparency in samsrf_surf by default (DSS)  
- Fixed samsrf_surf bug with eccentricity bounds clipping (DSS)  

### Version 9.02 (08-07-2022)  
- **Complete overhaul of forward-model time course prediction!** (DSS)  
- Updated search grid specification in SamSrf/Models accordingly (DSS)  
- Rewrote simulation functions for vectorised apertures (DSS)  
- Removed time course animation tool as no longer compatible with algorithm (DSS)  
- ViewApertures tool can now support both movie & vectorised apertures (DSS)  
- Updated Cookbook to reflect the new changes (DSS)  

------

### Version 8.4 (25-06-2022)  
- Added option to model compressive spatial summation in pRF fine-fit (DSS)  
- pRF functions now warn if aperture matrix contains negative numbers (DSS)  
- Updated ModelHelp text on negative estimates in Nelder-Mead algorithm (DSS)  
- Removed erroneous comments in prf_errfun (DSS)  

### Version 8.3 (16-05-2022)  
- Streamlined CF reverse correlation to prevent computational overload (DSS)  
- Added option to use robust summary statistics to estimate CF parameters (DSS)  
- Fixed inconsequential typos in ModelHelp tool & some help sections (DSS)  

### Version 8.2 (14-05-2022)  
- Now includes option to read GIfTI tiles (DSS)  
- Updated cookbook to reflect this change (DSS)  

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
