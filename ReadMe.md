# SamSrf VIII - Read Me
Version 8.0 (??-??-2022)

This major release includes the most recent updates which involved algorithms for 
fitting population receptive fields or connective fields using reverse correlation 
combined with posthoc fitting of correlation profiles. Importantly, SamSrf 8 also 
introduces several new fitting algorithms. By default we however still use the 
standard Nelder-Mead algorithm & you can run the exact same analyses as in SamSrf 7:

1. *Slow Coarse Fit:* Instead of only taking the best correlating prediction from
the search space it averages the parameters of a top percentile you can define.
This yields fairly good pRF estimates while being considerably faster than the
fine fit. But it is slower than the standard coarse fit.

2. *Adjustable parameter tolerance:* You can also adjust the parameter tolerance of 
the Nelder-Mead algorithm. If the tolerance is higher than default (1e-4) this also
trades accuracy for speed. This was the problem with SamSrf 6 because its tolerance
was too lenient (1e-2). But for many applications this may be fine & you could also 
try intermediate tolerance levels.
 
3. *Hooke-Jeeves algorithm:* You can also use Hooke-Jeeves pattern-search as an 
alternative to the Nelder-Mead algorithm (fminsearch). This algorithm performs well 
for standard 2D pRF models but YMMV. It can be faster but you can also lose accuracy 
(especially for true pRFs outside the stimulus range or with very large pRFs). 

------

SamSrf 8 was tested on Matlab R2020a. The Nelder-Mead algorithm requires Matlab's 
Optimization toolbox. The Hooke-Jeeves algorithm is implemented directly in SamSrf.
SamSrf strongly relies on parallel computing for a number of time-intensive analyses,  
so if you have Matlab's Parallel Computing Toolbox installed & you have a  
multi-core computer or cluster it should run faster. Many stages of the analysis  
are now geared towards parallel computing & may be quite slow without it.  

As of SamSrf 7 we completely overhauled of the toolbox. The GUI from SamSrf <=5 
has been abolished. Instead there are now example scripts for various pRF models 
in SamSrf/Models.   

In addition, there are various GUI-based tools available for specific functions:  
    1. Projecting data to surface:  SurfaceProjection  
    2. Making map figures:          DisplayMaps  
    3. Delineating ROIs:            DelineationTool  
    4. Viewing stimulus apertures:  ViewApertures

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

### Version 7.93 (14-04-2022)  
- New fast fit option to average top predictions in coarse fit (DSS)  
- Added Hooke-Jeeves algorithm as alternative for pRF & CF model fitting (DSS)  
- Now possible to define parameter tolerance for Nelder-Mead algorithm (DSS)  
- Fitting functions now perform checks on parameter definition vectors (DSS)  
- Updated samsrf_plot to use transparent dots for scatter plots (DSS)  
- Abandoned plan for SamOaSrf for Octave compatibility for the future is Python (DSS)  

## Questions/Comments?
* Contact Sam Schwarzkopf at: s.schwarzkopf@auckland.ac.nz

------
