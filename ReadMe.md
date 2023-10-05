# SamSrf IX - Read Me

This major release includes the most recent updates which involved an improved 
(much faster!) algorithm for fitting population receptive fields. This constitutes 
a fundamental changefrom previous SamSrf versions. The toolbox also includes 
reverse correlation pRF analysis & connective fields analysis. 

There are various GUIs & other tools available for specific purposes:  
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

### Version 9.6492 (05-10-2023) 
- Created DisplayScalpMaps tool with pRF profile inspector (DSS)  
- Adjusted verbosity of command line reports to reduce annoyance (DSS)  
- R^2 maps now use eccentricity colour scheme (DSS)  
- Small bugfix in samsrf_fitvsobs when downsampling predictions (DSS)  
- In a potential misnomer, samsrf_r2hist can now plot any row in Srf.Data (DSS)  
- Minor fix of figure legend in samsrf_r2hist (DSS)  

### Version 9.64 (28-09-2023) 
- Updated M/EEG functions for scalp distribution plots (DSS)  
- Srf structure for M/EEG data was reorganised (DSS)  

### Version 9.631 (20-09-2023) 
- Added support for EEG/MEG data files (DSS)  
- It is now possible to adjust threshold for reverse correlation pRFs convex hull (DSS)  
- Reverse correlation pRF now has option to allow negative peaks (DSS)  
- Reverse correlation pRF can now also simply use summary statistics to estimate parameters (DSS)  
- Negative peaks now also possible in samsrf_showprf (DSS)  
- Some other adaptations to reverse correlation backprojections (DSS)  

### Version 9.62 (09-08-2023) 
- Added function to generating reverse correlation pRF profiles (DSS)  
- Default for reverse correlation pRF is now to not save profiles! (DSS)  
- Also the default dimension of reverse correlation profiles is now 100 (DSS)  
- Updates & bugfix in samsrf_backproj_revcor for this purpose (DSS)  
- Updated ModelHelp text on reverse correlation pRF analysis accordingly (DSS)  
- Fixed bug in samsrf_fitvsobs when downsampling predictions (DSS)  
- Fixed bug with polar plots in samsrf_plot (DSS)  
- Added FuncHelp tool to display help sections on all functions (DSS)  

### Version 9.6 (29-06-2023) - Face-palm update!!!  
- **Data is now 32 bit by default to improve speed of model fits & reduce disc space!** (DSS)  
- You can switch this off & keep using 64 bit by setting *def_64bit = true* in *SamSrf_defaults.mat*  

### Version 9.51 (19-06-2023) - Advanced Age Day update!  
- Borders of ROI in DisplayMaps are now transparent (DSS)  
- Added option for transparency in samsrf_surf (DSS)  

### Version 9.5 (17-03-2023) - St Patrick's Day update!  
- Added option to use conventional approach for predicting neural response (DSS)  

### Version 9.425 (24-01-2023)  
- VectoriseApertures now gives error if aperture lengths are odd-numbered (DSS)  
- Bugfix for MakeOccRoi which stopped working after recent updates (DSS)  
- Bilateral surface functions now support Srfs with missing surface fields (DSS)  
- Fixed minor bug with samsrf_colourcode for custom colour schemes (DSS)  
- Added some more info about vectorising apertures in help sections (DSS)  

### Version 9.42 (19-10-2022)  
- Critical bugfix with data type in samsrf_mat2vol when saving negative values! (DSS)  
- Fixed various other bugs when using volumetric data files (DSS)  
- Volumetric data projection now allows concatenation & noise ceiling (DSS)  
- Fixed bug when converting volumetric data file without ROI back into NII (DSS)   

### Version 9.4 (05-10-2022)  
- Surface projection functions can now normalise without z-scoring (DSS)  
- Renamed negative amplitude estimate in reverse-correlation CF convex hull method (DSS)  

### Version 9.31 (27-09-2022)  
- Fixed bug with NaNs returned by smoothing algorithms (DSS)  
- Can now restrict samsrf_removenoise to ROI to help with limited memory issues (DSS)  
- Change to noise regression is now incorporated in CF analysis functions (DSS)  
- If symbolic links cannot be read, it now loads pial surfaces from .T1 file (DSS)  
- Fixed bug in ClusterSereno when saving in a different path (DSS)  
- Benson maps now use surf folder for Srf.Structural/Meshes (DSS)  
- Fixed a bug with DelineationTool not being able to save labels (DSS)  
- Visualisation of CFs in visual space are now zoomed in (DSS)  

### Version 9.22 (14-09-2022)  
- Added function to save out GII files of functional overlays for use outside SamSrf (DSS)  
- Added tool to cluster Sereno atlas ROIs into coarse anatomical regions (PWBU & DSS)  
- Added a note about using different canonical HRFs in ModelHelp (DSS)  
- Fixed bug with loading pial surfaces when your OS cannot read symbolic links (DSS)  
- Fixed bug when saving correlation CF profiles (DSS)  

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
