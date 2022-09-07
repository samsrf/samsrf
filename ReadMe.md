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

### Version 9.211 (07-09-2022)  
- Fixed bug when saving correlation CF profiles (DSS)  
- Added tool to cluster Sereno atlas ROIs into coarse anatomical regions (PWBU & DSS)  

### Version 9.2 (22-08-2022)  
- DisplayMaps now asks which hemisphere of bilateral Srf to display (DSS)  
- Added function to select a peak statistic within a ROI (DSS)  
- Simulation function can now implement compressive summation nonlinearity (DSS)  
- Vertex inspector can now display visual CF profiles but this is still buggy (DSS)  

### Version 9.1 (10-08-2022)  
- Convex hull CF algorithm now uses region growing for central CF estimation (DSS)  
- Now the inhibitory surround of a convex hull CF is also quantified (DSS)  
- Reverse correlation pRF analysis can now use convex hull algorithm (DSS)  
- Fixed bug in samsrf_fit_prf coarse fit when time course is flat (DSS)  
- Added option to samsrf_glm to turn off automatic global covariates (DSS)  
- Turned off parallel computing in samsrf_glm but commented out option remains (DSS)  
- Now saves GLM files in v7.3 file format in case they are too large (DSS)  
- Added option to define scaling factor/eccentricity in samsrf_simvsfit (DSS)  
- Corrected incorrect help section in template warping functions (DSS)  
- Changed samsrf_cometplot to use logarithmic density scale (DSS)  
- But output of samsrf_cometplot is still raw density matrix (DSS)  
- Increased default granularity in samsrf_cometplot (DSS)  
- Turned off transparency in samsrf_surf by default (DSS)  
- Fixed samsrf_surf bug with eccentricity bounds clipping (DSS)  
- Fixed bug in samsrf_surf when NaNs are present in map (DSS)  
- Changed example models for elliptical pRFs (DSS)  

### Version 9.02 (08-07-2022)  
- **Complete overhaul of forward-model time course prediction!** (DSS)  
- Updated search grid specification in SamSrf/Models accordingly (DSS)  
- Rewrote simulation functions for vectorised apertures (DSS)  
- Removed time course animation tool as no longer compatible with algorithm (DSS)  
- ViewApertures tool can now support both movie & vectorised apertures (DSS)  
- Updated Cookbook to reflect the new changes (DSS)  

------

## Questions/Comments?
* Contact Sam Schwarzkopf at: s.schwarzkopf@auckland.ac.nz

------
