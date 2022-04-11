# SamSrf VIII - Read Me
Version 8.0 (??-??-2022)

Latest major release but still compatible with version 7. The most recent updates
involved algorithms for fitting population receptive fields or connective fields 
using reverse correlation combined with posthoc fitting of correlation profiles.

For fitting forward-model pRFs we now added the Hooke-Jeeves pattern-search fit
as a faster alternative to the Nelder-Mead algorithm (fminsearch). This algorithm
necessitates some changes to how polar grid pRF work & therefore justifies a new
main release. However, you can still conduct the same analyses as in version 7 
and the two versions are fully compatible.

------

As of Version 7.0 we completely overhauled of the toolbox. The GUI from SamSrf <=5 
has been abolished. Instead there are now example scripts for various pRF models 
in SamSrf/Models.   

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

### Version 7.9 (??-04-2022)  
- Added Hooke-Jeeves algorithm for faster fine-fitting of pRFs (DSS) 
- Updated samsrf_plot to use transparent dots for scatter plots (DSS)  
- Abandoned plan for SamOaSrf for Octave compatibility - the future is Python (DSS)  

## Questions/Comments?
* Contact Sam Schwarzkopf at: s.schwarzkopf@auckland.ac.nz

------
