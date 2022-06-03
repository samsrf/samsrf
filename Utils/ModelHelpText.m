function HelpText = ModelHelpText(AnalysisFunc, ParameterName)
%
% HelpText = ModelHelpText(AnalysisFunc, ParameterName)
%
% Returns the help text for a given analysis function & parameter name.
% This is called internally by ModelHelp but you can also use it directly.
%

%% Which analysis function?
switch AnalysisFunc
    %% Forward model pRF
   case 'samsrf_fit_prf'
       switch ParameterName
           case 'Prf_Function'
               HelpText = {	'Function handle'
							''
							'This defines the pRF model used for predicting the time series. The most common model, first introduced in Dumoulin & Wandell 2008, is a 2D Gaussian. Other models with more complicated shapes, such as antagonistic centre-surround structure, elongated, or asymmetric profiiles, are also possible. You can find the model functions in the SamSrf/pRF subfolder & you can also make your own models. Defining the pRF function requires the following syntax, which will look quite complex at first glance. For example here is the 2D Gaussian model:'							
							''
							'   @(P,ApWidth) prf_gaussian_rf(P(1), P(2), P(3), ApWidth)'
							''
							'P contains the actual pRF parameters, in this case a vector with x0, y0, & sigma. ApWidth is a fixed input, which is required to enable flexibility with regard to the dimensions of your apertures but you can simply follow this basic arrangement.'
							''
							'The standard 2D Gaussian model also allows a lot of flexibility for defining 1D tuning functions. You can simply set the x0 or y0 parameter to zero and that way the model only has two free parameters while the other dimension is fixed; e.g. see the tuning curve examples in SamSrf/Models.' };
           case 'Name'
               HelpText = { 'Char'
							''
							'Defines the name of the analysis, which is what the final data file will be called. For instance, if this is ‘pRF_Gaussian’, the final map file for the left hemisphere will be called lh_pRF_Gaussian.mat. The idea is that you name this something that helps you identify what the data file contains. You could have several map files for different types of analysis of the same raw data in a folder.' };
           case 'Param_Names'
               HelpText = { 'Cell array of chars'
							''
							'Defines the names of the pRF parameters. In a retinotopy model, the first two will be the horizontal and vertical coordinate of the pRF centre, so we call them x0 and y0. For the pRF size (standard deviation of the Gaussian), we define the name as Sigma. The standard 2D Gaussian model only has three parameters so we can stop there. If we fit more parameters, they need to be named here. But note that the beta (amplitude) and intercept are fit separately, after the analysis has estimated all the pRF parameters.' };
           case 'Scaled_Param'
               HelpText = { 'Vector of booleans'
							''
							'Depending on the receptive field model you use, some or all the parameters need to be scaled to the eccentricity of your mapping stimulus. The apertures are defined such that the centre corresponds to coordinate 0,0 and the outer edges are -1 and +1. Therefore, any parameter that refers to a measure in this space, such as the pRF centre coordinates & pRF size, must be scaled to the actual eccentricity of the stimulus to make sense.'
							''
							'Model.Scaled_Param is a Boolean vector and it toggles whether a parameter is scaled or not. For example, in the case of the standard 2D Gaussian model, all three parameters are scaled so it is [1 1 1]. But in the case of some other models, such as the oriented asymmetric 2D pRF, one of the parameters is the orientation of the pRF and this is obviously not scaled. Critically, this vector must be of equal length as Param_Names.' };
           case 'Only_Positive'
               HelpText = { 'Vector of booleans'
							''
							'The model fitting function will look at the fitted (optimised) parameters in turn and reject those with absurd estimates. Any parameters that are scaled (see Scaled_Param) will automatically be rejected if the estimate is greater than 3 times the scaling factor. However, some parameters should also only be positive, e.g. the pRF size or tuning width cannot be negative.' 
							''
							'Model.Only_Positive is a Boolean vector that toggles whether a parameter must be positive. In the case of the standard 2D Gaussian model, this is only the pRF size (sigma).' 
                            ''
                            'Note that when using the Nelder-Mead algorithm some zero-symmetric parameters (like sigma which is squared in the Gaussian equation) can be negative but still produce reasonable pRF estimates. When rejecting these vertices you end up with a patchy map so it can be more ideal not to use this flag and instead use |sigma| in post-processing. However, such negative estimates tend to be very small and most people would probably filter them out before further analysis anyway. When using the Hooke-Jeeves algorithm, this constraint is built into the model fitting itself, so negative estimates are impossible. This is why this algorithm will generally produce more complete maps.' };
           case 'Scaling_Factor'
               HelpText = { 'Scalar'
							''
							'Once we know which parameters will be scaled, we also need to define by how much it will be scaled. In 2D retinotopic models, the scaling factor is the maximal eccentricity of the mapping stimulus (usually the radius of the mapped area). In one-dimensional tuning curve models (for example, tonotopy) this would be the half-width of the stimulus space because the middle is defined as 0.'
							''
							'This is a critical parameter which you must define according to your own needs!' };
           case 'TR'
               HelpText = { 'Scalar'
							''
							'The TR (repetition time) of your pulse sequence in seconds. This is needed because the convolution of the predicted time series with the HRF depends on the TR.'
							''
							'This is another critical parameter that depends on your experimental setup!' };
           case 'Hrf'
               HelpText = { 'Char or Vector of scalars'
							''
							'Which hemodynamic response function (HRF) to use. You can typically use the canonical HRF (based on the data in de Haas et al., 2014, Curr Biol), in which case you can leave this empty: [].' 
							''
							'If you estimated the HRF (see Cookbook or by refitting), you can provide here the file name of the estimated HRF. Or you can provide the HRF directly as a vector, where each component corresponds to one TR. Obviously, this latter option requires that the TR is the same in your HRF and your pRF data whereas the fit parameters of a HRF are more flexible.'
							''
							'Finally, if you don’t want any HRF to be used (as you might in some situations) you must set this to 1. Whatever you do, you need to specify this, even if it just [] to use the canonical.' };
           case 'Aperture_File'
               HelpText = { 'Char'
							''
							'Specify the file containing the apertures without .mat extension, so e.g. ‘aps_Bars’. Since this will typically depend on your experimental setup, you must define this. If you used the exact same stimulus design for each participant, you could just keep the aperture file in a common folder (with your model script) and provide the full path name here.' };
           case {'Param1' 'Param2' 'Param3' 'Param4' 'Param5'}
			   HelpText = { 'Vector of scalars'
							''
							'You can define the parameters of the search space used for the coarse fitting (extensive grid search) stage. For each of the five parameters you define a vector of the points of the search grid for that parameter, so this vector determines both the range and the granularity of the search space.' 
							''
							'For most people, there probably isn’t much of a reason to change these parameters. If you make the granularity too fine the search space will become unwieldy and you may even run out of memory. A finer search grid improves your coarse fit but is also slower. If running only the coarse fit, a finer grid is advisable.'
							''
							'The following parameters must be defined (which ones are needed depends on your pRF model). These are just names for the model engine. The actual names in the final map you already defined above in Param_Names. The parameters can mean several different things, depending on the pRF model. In models that don’t require all five parameters you simply set that one to 0 (but SamSrf does this for you automatically).'
							'' };
				switch ParameterName(end)
					case '1'
						HelpText{end+1} = 'Param1: x0 in Cartesian 2D pRF, polar angle in polar grid, or peak location in 1D tuning model';
					case '2'
						HelpText{end+1} = 'Param2: y0 in Cartesian 2D pRF, eccentricity in polar grid, or sigma in 1D tuning model';
					case '3'
						HelpText{end+1} = 'Param3: sigma in 2D pRF model';
					case '4'
						HelpText{end+1} = 'Param4: second sigma in asymmetric or centre-surround pRF';
					case '5'
						HelpText{end+1} = 'Param5: amplitude ratio (delta) in centre-surround pRF, or angle in oriented pRF, etc.'; 
				end
		   case 'Noise_Ceiling_Threshold'
               HelpText = { '[Optional] Scalar' 
							''
							'If the input raw data contains a Srf.Noise_Ceiling field, then you can use this parameter to restrict your analysis to only those vertices in the region of interest where the noise ceiling is above a certain threshold. This can dramatically speed up your analysis. It doesn’t really make sense to analyse vertices that do not contain any good data. Use the DisplayMaps tool to explore your raw data to determine what a reasonable threshold is.'
							''
							'Obtaining a noise ceiling for your data currently requires input data where each pRF run has the same temporal stimulus design (e.g. same order of bar positions) because the noise ceiling is calculated from the split-half reliability (see e.g. Morgan & Schwarzkopf, 2020, Front Hum Neurosci). There are other ways to calculate a noise ceiling but we have not implemented them. Perhaps you might want to do this (send a pull request?).'
							''
							'Defaults to 0 (no threshold).' };
           case 'Polar_Search_Space'
               HelpText = { '[Optional] Boolean' 
							''
							'If true, the first two parameters in the search space (Model.Param1 and Model.Param2) are treated as polar angle (in degrees) and eccentricity (in aperture space). This makes sense for retinotopic mapping where you know that the maps are in polar coordinates and the stimulus region is usually circular. A polar search grid will cover the whole circle and its surround equally, unlike a Cartesian grid where the corners extend further from the stimulus edge. A polar grid is also denser near the centre, which should help account for cortical magnification. Practically, in our experience this has not shown to make much difference, at least not with standard fitting parameters.' 
							''
							'Defaults to false.' };
           case 'Seed_Fine_Fit'
               HelpText = { '[Optional] Char' 
							''
							'You can define an already analysed map file here to seed your fine fit. In that case the function will not run any coarse fit but simply assume that the seed map contains the parameter guestimates you need. This seed map could either be a previous coarse fit (see Model.Coarse_Fit_Only), it could be a Benson/Wang-style anatomical prediction, or it could be a group-average map warped back into native space.' 
							''
							'Defaults to ‘’.' };
           case 'Replace_Bad_Fits'
               HelpText = { '[Optional] Boolean' 
							''
							'If true, then any fine fits that do not meet the criteria for a good fit will be replaced by the original coarse fit parameters. We have not used this feature very much but it could be a way to deal with odd misestimations in the fine-fit (which seem more likely in Nelder-Mead algorithm).' 
							''
							'Defaults to false.' };
           case 'Smoothed_Coarse_Fit'
               HelpText = { '[Optional] Scalar' 
							''
							'If greater than zero, the data are first smoothed with the kernel defined here. This makes the process more compatible with data estimated by versions up to SamSrf 5. A kernel of 5 is what was used in older versions (& the original Dumoulin & Wandell study) before coarse fitting. In practice, we nowadays often smooth raw data during the surface projection stage & this step is unnecessary.' 
							''
							'Defaults to 0.' };
           case 'Coarse_Fit_Only'
               HelpText = { '[Optional] Boolean' 
							''
							'If true, only the coarse fit is conducted. Betas are estimated by linear regression afterwards as would normally be done after the fine fit. By definition the precision of the estimates is limited by the granularity of the search space. If you run a "slow coarse fit" (see. Model.Coarse_Fit_Percentile and description in Cookbook) you can obtain more precise estimates however.' 
							''
							'The suffix _CrsFit is automatically appended to the file name to indicate that this map was analysed only with the coarse fit.';
							''
							'Defaults to false.' };
           case 'Fine_Fit_Threshold'
               HelpText = { '[Optional] Scalar' 
							''
							'Defines above which R^2 from the coarse fit data are included in the fine fit. Increasing this may speed up your analysis but it can also mean you lose data.' 
							''
							'Defaults to 0.01.' };
           case 'Only_Positive_Coarse_Fits'
               HelpText = { '[Optional] Boolean' 
							''
							'If true, only coarse fits with positive correlations are included in the fine fit. By default, both positive and negative correlations are included (because the correlation coefficient is squared). This accounts for potential pRFs with negative betas which can be meaningful under some circumstances. However, in other contexts including negative correlations may not be advised.'
							''
							'Defaults to false.' };
           case 'Coarse_Fit_Block_Size'
               HelpText = { '[Optional] Scalar' 
							''
							'Defines the size of the data chunks that the coarse fit can run simultaneously, assuming you are using MATLAB version as new as 7.13 (but honestly you should be using 2020a or newer as SamSrf 8 was tested on this version). You may need to reduce this if you have a really densely sampled search space or if you have high resolution data and therefore your region of interest is enormous. The higher this number, the faster your coarse fit is. You could try increasing it too if your computer has a lot of memory.'
							''
							'Defaults to 10000.' };
           case 'Downsample_Predictions'
               HelpText = { '[Optional] Scalar' 
							''
							'This defines the microtime resolution of the time course prediction. You can use this for stimulus designs with a faster temporal resolution than the TR. For example, if your TR is 1 second, but your stimulus position updates every 100 ms, this parameter would be 10. The model prediction will then pretend that the TR is 100 ms & only when comparing prediction to actual data the time series is downsampled to 1 s.'
							''
							'Defaults to 1.' };
           case 'Coarse_Fit_Percentile'
               HelpText = { '[Optional] Scalar' 
							''
							'This defines the percentile above (or equal to) which coarse fit R^2 are included in estimating parameters.' 
							''
							'By default, the coarse fit only picks the parameters with maximal R^2. However, if this parameter is lower than 100 (i.e. maximum) it averages the parameters for these predictions. This can result in more precise parameter estimates. We call this slow coarse fit because it is slower than only taking the maximum - however, it is also considerably faster than the fine fit but in many scenarios it can perform similarly.'
							''
							'Keep in mind that the slow coarse fit is limited by your search space. The parameter estimates are constrained by the limits of the search space, so if you expect parameters outside the stimulus range the search space should extend beyond the stimulus range considerably. Similarly, the mean of parameters can only be as precise as the inputs. A search space with coarse granularity can only produce coarse estimates.'
							''
							'The value you choose should be fairly close to 100. The more poorer coarse fits you include the more contaminated your estimates will be. Of course, the closer to 100 the more like the maximum (i.e. standard coarse fit) it will become. For a search space with ~35,000 grid points we use a value of 99.99 (i.e. the top 0.01% of coarse fits) but YMMV.'
							''
							'Defaults to 100.' };
           case 'Hooke_Jeeves_Steps'
               HelpText = { '[Optional] Vector of scalars'
							''
							'This defines for each pRF parameter the initial step size to be used for the Hooke-Jeeves (parameter search) algorithm. If this is defined, the fine fit uses this algorithm instead of the standard Nelder-Mead simplex search algorithm (fminsearch) - so be certain you want this!'
							''
							'The Hooke-Jeeves algorithm takes steps above & below each parameter from the current estimate. If a better fit is found, this is taken as the new estimate. If no better fit is found, the step size is halved. See the Cookbook for more detail.'
							''
							'The initial step sizes depend on the nature of your parameters & the search space. The steps should probably be smaller than the granularity of the search space. Note that this is in aperture space! So far, we have used [.01 .01 .01] for a standard 2D Gaussian pRF analysis.'
							''
							'By default, this parameter is not defined & thus the Nelder-Mead algorithm is used instead for fine fitting.' };
           case 'Nelder_Mead_Tolerance'  
               HelpText = { '[Optional] Scalar'
							''
							'When using the Nelder-Mead algorithm SamSrf will normally use the MATLAB default for the parameter estimation tolerance (in R2020a this was 1e-4). By setting a higher value you can make the estimation less strict. This will speed up the fine fit but also reduces its precision.' 
							''
							'Important: The issues with precision we saw in SamSrf 6 were due to the tolerance being 1e-2! This was cleary suboptimal for our purposes but YMMV. There are scenarios where this lenient tolerance may be justified. An intermediate tolerance (e.g. 1e-3?) may also perform well but we never tested this.'
							''
							'By default, this parameter is not defined & the standard Nelder-Mead tolerance is used. However, if Hooke_Jeeves_Steps is defined, the Hooke-Jeeves algorithm is used instead, so be sure you know what you are doing!' };
               
           otherwise
               error('samsrf_fit_prf does not have this parameter!');
       end
    
    %% Reverse correlation pRF                
   case 'samsrf_revcor_prf'
       switch ParameterName
           case 'Name'
               HelpText = { 'Char'
							''
							'Defines the name of the analysis, which is what the final data file will be called. The suffix _Rcp is automatically appended to the file name to indicate that this map was analysed with reverse correlation. For instance, if this is ‘pRF_Gaussian’, the final map file for the left hemisphere will be called lh_pRF_Gaussian_Rcp.mat. The idea is that you name this something that helps you identify what the data file contains. You could have several map files for different types of analysis of the same raw data in a folder.' };
           case 'Scaling_Factor'
               HelpText = { 'Scalar'
							''
							'Defines the scale of the stimulus space because the aperture space is fixed between -1 and +1. In 2D retinotopic models, the scaling factor is the maximal eccentricity of the mapping stimulus (usually the radius of the mapped area). In one-dimensional tuning curve models (for example, tonotopy) this would be the half-width of the stimulus space because the middle is defined as 0.'
							''
							'This is a critical parameter which you must define according to your own needs!' };
           case 'TR'
               HelpText = { 'Scalar'
							''
							'The TR (repetition time) of your pulse sequence in seconds. This is needed because the convolution of the apertures with the HRF depends on the TR.'
							''
							'This is another critical parameter that depends on your experimental setup!' };
           case 'Hrf'
               HelpText = { 'Char or Vector of scalars'
							''
							'Which hemodynamic response function (HRF) to use. You can typically use the canonical HRF (based on the data in de Haas et al., 2014, Curr Biol), in which case you can leave this empty: [].' 
							''
							'If you estimated the HRF (see Cookbook or by refitting), you can provide here the file name of the estimated HRF. Or you can provide the HRF directly as a vector, where each component corresponds to one TR. Obviously, this latter option requires that the TR is the same in your HRF and your pRF data whereas the fit parameters of a HRF are more flexible.'
							''
							'Finally, if you don’t want any HRF to be used (as you might in some situations) you must set this to 1. Whatever you do, you need to specify this, even if it just [] to use the canonical.' };
           case 'Aperture_File'
               HelpText = { 'Char'
							''
							'Specify the file containing the apertures without .mat extension, so e.g. ‘aps_Bars’. Since this will typically depend on your experimental setup, you must define this. If you used the exact same stimulus design for each participant, you could just keep the aperture file in a common folder (with your model script) and provide the full path name here.' };
           case 'Prf_Function'
               HelpText = {	'[Optional] Function handle'
							''
							'This defines the pRF model to fit to the reverse correlation profiles. The most common model is a 2D Gaussian. Other models with more complicated shapes, such as antagonistic centre-surround structure, elongated, or asymmetric profiiles, are also possible (but so far we have not tested them). You can find the model functions in the SamSrf/pRF subfolder & you can also make your own models. Defining the pRF function requires the following syntax, which will look quite complex at first glance. For example here is the 2D Gaussian model:'							
							''
							'   @(P,ApWidth) prf_gaussian_rf(P(1), P(2), P(3), ApWidth)'
							''
							'P contains the actual pRF parameters, in this case a vector with x0, y0, & sigma. ApWidth is a fixed input, which is required to enable flexibility with regard to the dimensions of your apertures but you can simply follow this basic arrangement.'
							''
							'As with forward-modelling, the standard 2D Gaussian model also allows a lot of flexibility for defining 1D tuning functions. You can simply set the x0 or y0 parameter to zero and that way the model only has two free parameters while the other dimension is fixed; e.g. see the tuning curve examples in SamSrf/Models.'
							''
							'By default, this is not defined & no 2D model is fit to the reverse correlation profiles.' };
           case 'Param_Names'
               HelpText = { 'Cell array of chars'
							''
							'Defines the names of the pRF parameters. In a retinotopy model, the first two will be the horizontal and vertical coordinate of the pRF centre, so we call them x0 and y0. For the pRF size (standard deviation of the Gaussian), we define the name as Sigma. The standard 2D Gaussian model only has three parameters so we can stop there. If we fit more parameters, they need to be named here. But note that the beta (amplitude) and intercept are fit automatically & need not be included.'
							''
							'This parameter is mandatory when a 2D pRF model is fit to the reverse correlation profiles (see Prf_Function)!' };
           case 'Scaled_Param'
               HelpText = { 'Vector of booleans'
							''
							'Depending on the receptive field model you use, some or all the parameters need to be scaled to the eccentricity of your mapping stimulus. The apertures are defined such that the centre corresponds to coordinate 0,0 and the outer edges are -1 and +1. Therefore, any parameter that refers to a measure in this space, such as the pRF centre coordinates & pRF size, must be scaled to the actual eccentricity of the stimulus to make sense.'
							''
							'This parameter is mandatory when a 2D pRF model is fit to the reverse correlation profiles (see Prf_Function)!' };
           case 'SeedPar_Function'
               HelpText = { 'Function handle'
							''
							'This function determines the seed parameters for the model fitting stage. This should therefore be based on the rough estimates obtained from the reverse correlation profiles. In order to allow flexibility with complex pRF models, we added this so you can determine adequate seed values. To date, we have however only tested a conventional 2D Gaussian model. The seed parameter function we used is this:'
							''
							'   @(V) [V(2:3); V(4)/(2*sqrt(2*log(2)))]'
							''
							'V is a vector of parameter estimates from the reverse correlation profile:'
							'   V(1) = R^2_{RC}, the maximum R^2 in the reverse correlation profile'
							'   V(2) = x0 at maximum correlation'
							'   V(3) = y0 at maximum correlation'
							'   V(4) = Full width at half maximum of the reverse correlation profile'
							'So the above function will seed to the maximum location & a rule-of-thumb estimate of the sigma based on the width of the reverse correlation profile.'
							''
							'This parameter is mandatory when a 2D pRF model is fit to the reverse correlation profiles (see Prf_Function)!' };
           case 'R2_Threshold'
               HelpText = { '[Optional] Scalar'
							''
							'When fitting a 2D pRF model you can first threshold the data based on the peak in the reverse correlation profile (R^2_{RC}) before running the model fit. It makes little sense to fit models to rubbish profiles.'
							''
							'Defaults to 0.' };
           case 'Rdim'
               HelpText = { '[Optional] Scalar'
							''
							'If you decide to save them in the map file (see Model.Save_Rmaps) this defines the side length of reverse correlation profiles. It is probably unnecessary to save those in the resolution of apertures. A smaller number means you use less disc space.'
							''
							'Defaults to 50.' };
           case 'Noise_Ceiling_Threshold'
               HelpText = { '[Optional] Scalar' 
							''
							'If the input raw data contains a Srf.Noise_Ceiling field, then you can use this parameter to restrict your analysis to only those vertices in the region of interest where the noise ceiling is above a certain threshold. This can dramatically speed up your analysis. It doesn’t really make sense to analyse vertices that do not contain any good data. Use the DisplayMaps tool to explore your raw data to determine what a reasonable threshold is.'
							''
							'Obtaining a noise ceiling for your data currently requires input data where each pRF run has the same temporal stimulus design (e.g. same order of bar positions) because the noise ceiling is calculated from the split-half reliability (see e.g. Morgan & Schwarzkopf, 2020, Front Hum Neurosci). There are other ways to calculate a noise ceiling but we have not implemented them. Perhaps you might want to do this (send a pull request?).'
							''
							'Defaults to 0 (no threshold).' };
           case 'Save_Rmaps'
               HelpText = { '[Optional] Boolean'
							''
							'If true, the reverse correlation profiles are saved in the data file. This is currently necessary for using the backprojection function for reconstructing visual field coverage (samsrf_backproj_revcor). However, this means the map file will take up a lot more disc space.'
							''
							'Most other analyses are possible without saving these profiles because they can be easily calculated again.'
							''
							'Defaults to true.' };
           case 'Hooke_Jeeves_Steps'
               HelpText = { '[Optional] Vector of scalars'
							''
							'This defines for each pRF parameter the initial step size to be used for the Hooke-Jeeves (parameter search) algorithm. If this is defined, the 2D model fit uses this algorithm instead of the standard Nelder-Mead simplex search algorithm (fminsearch) - so be certain you want this!'
							''
							'The Hooke-Jeeves algorithm takes steps above & below each parameter from the current estimate. If a better fit is found, this is taken as the new estimate. If no better fit is found, the step size is halved. See the Cookbook for more detail.'
							''
							'The initial step sizes depend on the nature of your parameters & your experimental details. Note that unlike for forward-modelling this is in visual space! So far, we have used [.1 .1 .1] for a standard 2D Gaussian pRF analysis.'
							''
							'By default, this parameter is not defined & thus the Nelder-Mead algorithm is used instead for fine fitting.' };
           case 'Nelder_Mead_Tolerance'  
               HelpText = { '[Optional] Scalar'
							''
							'When using the Nelder-Mead algorithm SamSrf will normally use the MATLAB default for the parameter estimation tolerance (in R2020a this was 1e-4). By setting a higher value you can make the estimation less strict. This will speed up the 2D model fit but also reduces its precision.' 
							''
							'To date, we have not tested this feature for reverse correlation model fits. The model fitting is already quite fast even with the standard tolerance.'
							''
							'By default, this parameter is not defined & the standard Nelder-Mead tolerance is used. However, if Hooke_Jeeves_Steps is defined, the Hooke-Jeeves algorithm is used instead, so be sure you know what you are doing!' };
               
           otherwise
               error('samsrf_revcor_prf does not have this parameter!');
       end

    %% Reverse correlation CF    
   case 'samsrf_revcor_cf'
       switch ParameterName
           case 'Name'
               HelpText = { 'Char'
							''
							'Defines the name of the analysis, which is what the final data file will be called. For instance, if this is ‘CF’, the final map file for the left hemisphere will be called lh_CF.mat. The idea is that you name this something that helps you identify what the data file contains. You could have several map files for different types of analysis of the same raw data in a folder.' };
           case 'SeedRoi'
               HelpText = { 'Char'
							''
							'Path &  file name of the FreeSurfer label file for the seed region (without file extension). Each time series in your region of interest will be regressed on all the time series in this seed region to determine the reverse correlation profile of the connective field.' };
           case 'Template'
               HelpText = { 'Char'
							''
							'Path & file name of the SamSrf map data file (without .mat extension) containing the template map used to translate anatomical CF positions into visual field estimates. This file must therefore contain a retinotopic map with at least the x0 & y0 positions of each template pRF in Srf.Data(2:3,:) and a (potentially dummy) goodness-of-fit value in Srf.Data(1,:).' };
           case 'Smoothing'
               HelpText = { '[Optional] Scalar'
							''
							'Defines the smoothing kernel with which to smooth the reverse correlation profiles. The units are geodesic steps in the surface mesh. The function first generates a distance matrix between all vertices in the seed region & then later uses this to apply the smoothing kernel to the reverse correlation profiles. When using data smoothed at the surface projection stage this is usually not necessary.'
							''
							'Defaults to 0 (no smoothing).' };
           case 'Noise_Ceiling_Threshold'
               HelpText = { '[Optional] Scalar' 
							''
							'If the input raw data contains a Srf.Noise_Ceiling field, then you can use this parameter to restrict your analysis to only those vertices in the region of interest where the noise ceiling is above a certain threshold. This can dramatically speed up your analysis. It doesn’t really make sense to analyse vertices that do not contain any good data. Use the DisplayMaps tool to explore your raw data to determine what a reasonable threshold is.'
							''
							'Obtaining a noise ceiling for your data currently requires input data where each pRF run has the same temporal stimulus design (e.g. same order of bar positions) because the noise ceiling is calculated from the split-half reliability (see e.g. Morgan & Schwarzkopf, 2020, Front Hum Neurosci). There are other ways to calculate a noise ceiling but we have not implemented them. Perhaps you might want to do this (send a pull request?).'
							''
							'Defaults to 0 (no threshold).' };
           case 'Global_Signal_Correction'
               HelpText = { '[Optional] Boolean'
							''
							'If true, then the time series of the input data are corrected by the global mean time series using linear regression. This removes noise components common across the whole cortex which are likely to result e.g. from head motion. You could also regress out other nuissance factors but you need to do this yourself. Perhaps it will be implemented here in a future version.'
							''
							'Defaults to true.' };
           case 'Save_Rmaps'
               HelpText = { '[Optional] Boolean'
							''
							'If true, the reverse correlation profiles are saved in the data file. This means the map file will take up a lot more disc space & there is usually no reason for this because these profiles can be easily calculated again.'
							''
							'Defaults to false.' };
           case 'Fit_pRF'
               HelpText = { '[Optional] Scalar'
							''
							'  1: Fits a 2D Gaussian pRF model to the reverse correlation profile projected into visual space via the template map. For the time being, no other pRF shapes are possible.'
                            ''
                            '  0: Use the convex hull algorithm to estimate the centroid & size (expressed as Sigma, converted from the square root of the area) of the CF profile.'
							''
                            ' -1: Use summary statistics of the significant (R > half maximum) of seed coordinates to estimate CF parameters. Specifically, X & Y position are estimated by the median, & Sigma by the Euclidean distance based on the median absolute deviations for X & Y.'
                            ''
							'Since SamSrf 8.0, before parameter estimation the reverse correlation profile is clipped to only those vertices above half the peak correlation. This is necessary because otherwise the estimated pRF size tends to be constant across eccentricities. However, this feature is still likely to evolve & we may add an option to specify this in later versions.'
							''
							'Defaults to 0, meaning that CF parameters are estimated from the reverse correlation profiles via the convex hull algorithm.' };
           case 'Hooke_Jeeves_Steps'
               HelpText = { '[Optional] Vector of scalars'
							''
							'This defines for each pRF parameter the initial step size to be used for the Hooke-Jeeves (parameter search) algorithm. If this is defined, the 2D model fit uses this algorithm instead of the standard Nelder-Mead simplex search algorithm (fminsearch) - so be certain you want this!'
							''
							'The Hooke-Jeeves algorithm takes steps above & below each parameter from the current estimate. If a better fit is found, this is taken as the new estimate. If no better fit is found, the step size is halved. See the Cookbook for more detail.'
							''
							'The initial step sizes depend on the nature of your parameters & your experimental details. Note that unlike for forward-modelling this is in visual space! So far, we have used [.1 .1 .1] for a reinotopic data set based on a stimulus with about 10 degrees ecccentricity. For other applications (e.g. peripheral mapping?) this is probably not appropriate.'
							''
							'By default, this parameter is not defined & thus the Nelder-Mead algorithm is used instead for fine fitting.' };
           case 'Nelder_Mead_Tolerance'  
               HelpText = { '[Optional] Scalar'
							''
							'When using the Nelder-Mead algorithm SamSrf will normally use the MATLAB default for the parameter estimation tolerance (in R2020a this was 1e-4). By setting a higher value you can make the estimation less strict. This will speed up the 2D model fit but also reduces its precision.' 
							''
							'To date, we have not tested this feature for reverse correlation model fits. The model fitting is already quite fast even with the standard tolerance.'
							''
							'By default, this parameter is not defined & the standard Nelder-Mead tolerance is used. However, if Hooke_Jeeves_Steps is defined, the Hooke-Jeeves algorithm is used instead, so be sure you know what you are doing!' };
               
           otherwise
               error('samsrf_revcor_cf does not have this parameter!');
       end
            
    %% Forward model CF    
   case 'samsrf_fit_cf'
       switch ParameterName
           case 'Name'
               HelpText = { 'Char'
							''
							'Defines the name of the analysis, which is what the final data file will be called. The suffix _Fwd is appended automatically to indicate that this is a CF map generated by forward-modelling. For instance, if this is ‘CF’, the final map file for the left hemisphere will be called lh_CF_Fwd.mat. The idea is that you name this something that helps you identify what the data file contains. You could have several map files for different types of analysis of the same raw data in a folder.'
							''
							'The time series of the CF is estimated by the first eigenvariate of all the time series inside the CF patch. Thus the analysis makes no explicit assumptions about the strength of correlation of individual vertices (as a Gaussian CF model would) & can even incorporate negative correlations.'
							''
							'Please note that the forward-modelling CF currently only uses a grid search approach but no further optimisation algorithm is applied. To specify the search space either Polar, Eccentricity, or Sizes must be defined.' };
           case 'SeedRoi'
               HelpText = { 'Char'
							''
							'Path &  file name of the FreeSurfer label file for the seed region (without file extension). Each time series in your region of interest will be regressed on all the time series in this seed region to determine the reverse correlation profile of the connective field.' };
           case 'Template'
               HelpText = { 'Char'
							''
							'Path & file name of the SamSrf map data file (without .mat extension) containing the template map used to translate anatomical CF positions into visual field estimates. This file must therefore contain a retinotopic map with at least the x0 & y0 positions of each template pRF in Srf.Data(2:3,:) and a (potentially dummy) goodness-of-fit value in Srf.Data(1,:).' };
           case 'Smoothing'
               HelpText = { '[Optional] Scalar'
							''
							'Defines the smoothing kernel in millimeters of a Gaussian kernel in spherical mesh space. This smoothing operation is applied to the raw input data. If you use smoothing at the surface project step this should not be necessary.'
							''
							'Important: This works completely differently to how smoothing can be used in the reverse-correlation CF analysis!'
							''
							'Defaults to 0 (no smoothing).' };
           case 'Noise_Ceiling_Threshold'
               HelpText = { '[Optional] Scalar' 
							''
							'If the input raw data contains a Srf.Noise_Ceiling field, then you can use this parameter to restrict your analysis to only those vertices in the region of interest where the noise ceiling is above a certain threshold. This can dramatically speed up your analysis. It doesn’t really make sense to analyse vertices that do not contain any good data. Use the DisplayMaps tool to explore your raw data to determine what a reasonable threshold is.'
							''
							'Obtaining a noise ceiling for your data currently requires input data where each pRF run has the same temporal stimulus design (e.g. same order of bar positions) because the noise ceiling is calculated from the split-half reliability (see e.g. Morgan & Schwarzkopf, 2020, Front Hum Neurosci). There are other ways to calculate a noise ceiling but we have not implemented them. Perhaps you might want to do this (send a pull request?).'
							''
							'Defaults to 0 (no threshold).' };
           case 'Global_Signal_Correction'
               HelpText = { '[Optional] Boolean'
							''
							'If true, then the time series of the input data are corrected by the global mean time series using linear regression. This removes noise components common across the whole cortex which are likely to result e.g. from head motion. You could also regress out other nuissance factors but you need to do this yourself. Perhaps it will be implemented here in a future version.'
							''
							'Defaults to true.' };
           case 'Polar'
               HelpText = { '[Semi-optional] Vector of scalars'
							''
							'Defines the search space for the grid search in terms of polar angle wedges. This should be very fast but only estimates the polar angle map of the CFs with the specified granularity.'
							''
							'This parameter takes precedence over Eccentricity or Sizes. So if this is defined these other parameters are ignored.' };
           case 'Eccentricity'
               HelpText = { '[Semi-optional] Vector of scalars'
							''
							'Defines the search space for the grid search in terms of eccentricity bands. This should be very fast but only estimates the eccentricity map of the CFs with the specified granularity.'
							''
							'This parameter takes precedence over Sizes. So if this is defined, Sizes is ignored, unless Polar is also defined, in which case that is used instead.' };
           case 'Sizes'
               HelpText = { '[Semi-optional] Vector of scalars'
							''
							'Defines the search space for the grid search in terms of geodesic radii of the circular CF patches. This is very slow because the search space contains CF predictions for each vertex in the seed region & for each of these sizes. So for example if your seed region has 10,000 vertices, and you specify 10 sizes here, the search space contains 100,000 predictions!'
							''
							'This parameter is only used if Polar or Eccentricity are not defined.' };
           case 'Coarse_Fit_Block_Size'
               HelpText = { '[Optional] Scalar' 
							''
							'Defines the size of the data chunks that the coarse fit can run simultaneously, assuming you are using MATLAB version as new as 7.13 (but honestly you should be using R2020a (v9.8.0) or newer as SamSrf 8 was tested on this version). You may need to reduce this if you have a really densely sampled search space and/or a large seed region. The higher this number, the faster your coarse fit is. You could try increasing it too if your computer has a lot of memory.'
							''
							'Defaults to 10000.' };
           case 'Patch_Size'  
               HelpText = { '[Optional] Scalar'
							''
							'This defines the radius of receiving patches in units of geodesic steps. The analysis computes the first eigenvariate of the time series in the patch to correlate with the CF time series from the seed region. The idea is that individual vertices might be quite noisy but a whole patch of cortex is more robust. Practically, this has not been tested extensively.'
							''
							'Defaults to 0 meaning no receiving patch is used & vertices are treated individually.' };
               
           otherwise
               error('samsrf_fit_cf does not have this parameter!');
       end
               
    otherwise
        error('Unknown analysis function specified!')                       
end
