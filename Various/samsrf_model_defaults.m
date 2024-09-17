function Model = samsrf_model_defaults(AnalysisFunc, Model)
%
% Model = samsrf_model_defaults(AnalysisFunc, Model)
%
% Adds defaults of optional parameters for analysis AnalysisFunc (case-sensitive char!) 
% if they haven't been defined in Model. It also checks that all mandatory parameters 
% (e.g. Name or TR etc) are present.
%
% The function also automatically fills unused search space parameters 
% (those which have not been defined in Model.Param_Names) with zeros.
%
% It also performs a few checks on the soundness of parameters.
%
% 20/04/2022 - SamSrf 8 version (DSS)
% 06/05/2022 - Changed CF model fitting default (DSS)
% 09/08/2023 - Default in reverse correlation pRF is now not to save profiles (DSS)
% 31/08/2023 - Added option for negative peaks in reverse correlation pRF (DSS)
% 19/09/2023 - Reverse correlation pRF now defaults to convex hull algorithm (DSS)
%

%% Which analysis function?
switch AnalysisFunc
    %% Forward-model pRF
    case 'samsrf_fit_prf'
        % Mandatory parameters
        if ~isfield(Model, 'Name')
            error('Analysis name is undefined!');
        end
        if ~isfield(Model, 'Prf_Function')
            error('pRF function is undefined!');
        end
        if ~isfield(Model, 'Param_Names')
            error('Parameter names are undefined!');
        end
        if ~isfield(Model, 'Scaled_Param')
            error('Scaled parameter flags are undefined!');
        end
        if ~isfield(Model, 'Only_Positive')
            error('Only positive flags are undefined!');
        end
        if ~isfield(Model, 'Scaling_Factor')
            error('Scaling factor (eccentricity?) is undefined!');
        end
        if ~isfield(Model, 'TR')
            error('TR is undefined!');
        end
        if ~isfield(Model, 'Downsample_Predictions')
            Model.Downsample_Predictions = 1; % Downsampling factor by which Model.TR mismatches the true TR
        end
        if ~isfield(Model, 'Hrf')
            error('HRF is undefined!');
        end
        if ~isfield(Model, 'Aperture_File')
            error('Aperture file is undefined!');
        end
        if ~isfield(Model, 'Noise_Ceiling_Threshold')
            Model.Noise_Ceiling_Threshold = 0; % Limit analysis to data above a certain noise ceiling
        end
        if ~isfield(Model, 'Seed_Fine_Fit')
            Model.Seed_Fine_Fit = ''; % Use no seed map, so run coarse fit instead
        end 
        % Search space definition
        if ~isfield(Model, 'Polar_Search_Space')
            Model.Polar_Search_Space = false; % Search space is Cartesian (not actually mandatory)
        end
        for p = 1:10
            if ~isfield(Model, ['Param' num2str(p)])
                if p <= length(Model.Param_Names)
                    error(['Search space parameter #' num2str(p) ' is undefined!']);
                else
                    Model.(['Param' num2str(p)]) = 0; % Unused parameter set to zero
                end
            end
        end
        
        % Defaults for optional parameters
        if ~isfield(Model, 'Coarse_Fit_Only')
            Model.Coarse_Fit_Only = false; % Only run coarse fit & then save
        end
        if ~isfield(Model, 'Smoothed_Coarse_Fit')
            Model.Smoothed_Coarse_Fit = 0; % No smoothing on coarse fit
        end
        if ~isfield(Model, 'Fine_Fit_Threshold')
            Model.Fine_Fit_Threshold = 0.01; % Include coarse fits with R^2>0.01 in fine fit
        end
        if ~isfield(Model, 'Coarse_Fit_Block_Size')
            Model.Coarse_Fit_Block_Size = 10000; % Number of simultaneous data columns in coarse fit
        end
        if ~isfield(Model, 'Only_Positive_Coarse_Fits')
            Model.Only_Positive_Coarse_Fits = false; % Coarse fit can either be negative or positive correlation to pass 
        end
        if ~isfield(Model, 'Replace_Bad_Fits')
            Model.Replace_Bad_Fits = false; % Don't replace bad fine fits with coarse fit
        end
        if ~isfield(Model, 'Coarse_Fit_Percentile')
            Model.Coarse_Fit_Percentile = 100; % Which percentile of correlations to include in coarse fit parameter estimates
        end

        if ~isfield(Model, 'Aperture_Mean_Response')
            Model.Aperture_Mean_Response = false; % If true, uses conventional Dumoulin & Wandell 2008 approach
        end
        if ~isfield(Model, 'Compressive_Nonlinearity')
            Model.Compressive_Nonlinearity = false; % Model compressive spatial summation nonlinearity of response?
        end
        
        % If coarse fit percentile invalid
        if Model.Coarse_Fit_Percentile < 0 || Model.Coarse_Fit_Percentile > 100
            error('Coarse fit percentile invalid! It should probably be between 99.9-100...');
        end
        if Model.Coarse_Fit_Percentile < 99.9
            warning('Low coarse fit percentile. It should probably be between 99.9-100...');
        end
        % If coarse fit only, we cannot use seed map 
        if Model.Coarse_Fit_Only 
            if ~isempty(Model.Seed_Fine_Fit) % In case stupid choices were made
                error('No point running only coarse fit when seeding the fine fit!');
            end
        end

        % Ensure mandatory model vectors are sound
        if length(Model.Param_Names) ~= length(Model.Scaled_Param)
            error('Mismatch between number of parameter names & scaled-parameter flags!');
        end
        if length(Model.Param_Names) ~= length(Model.Only_Positive)
            error('Mismatch between number of parameter names & only-positive flags!');
        end
        if isfield(Model, 'Hooke_Jeeves_Steps')
            if length(Model.Param_Names) ~= length(Model.Hooke_Jeeves_Steps)
                error('Mismatch between number of parameter names & Hooke-Jeeves step sizes!');
            end
        end
                
    %% Reverse-correlation pRF
    case 'samsrf_revcor_prf'
        % Mandatory parameters
        if ~isfield(Model, 'Name')
            error('Analysis name is undefined!');
        end            
        if ~isfield(Model, 'Scaling_Factor')
            error('Scaling factor (eccentricity?) is undefined!');
        end
        if ~isfield(Model, 'TR')
            error('TR is undefined!');
        end
        if ~isfield(Model, 'Hrf')
            error('HRF is undefined!');
        end
        if ~isfield(Model, 'Aperture_File')
            error('Aperture file is undefined!');
        end
        
        % Default for parameter estimation
        if ~isfield(Model, 'Noise_Ceiling_Threshold')
            Model.Noise_Ceiling_Threshold = 0; % Limit analysis to data above a certain noise ceiling
        end
        if ~isfield(Model, 'Allow_Negative_Peaks')
            Model.Allow_Negative_Peaks = false; % Whether peaks can be negative
        end
        if ~isfield(Model, 'Prf_Function')
            Model.Prf_Function = 0;
        end
        
        % If fitting 2D pRF model
        if strcmpi(class(Model.Prf_Function), 'function_handle')
            if ~isfield(Model, 'Param_Names')
                error('Parameter names are undefined!');
            end
            if ~isfield(Model, 'Scaled_Param')
                error('Scaled parameter flags are undefined!');
            end
            if ~isfield(Model, 'SeedPar_Function')
                error('Seed parameter function is undefined!');
            end
            
            % Ensure mandatory model vectors are sound
            if length(Model.Param_Names) ~= length(Model.Scaled_Param)
                error('Mismatch between number of parameter names & scaled-parameter flags!');
            end
            
            % Optional parameter for fitting
            if ~isfield(Model, 'R2_Threshold')
                Model.R2_Threshold = 0;
            end
        elseif Model.Prf_Function == 0
            % Using convex hull algorithm
            if ~isfield(Model, 'Convex_Hull_Threshold')
                Model.Convex_Hull_Threshold = 0.5; % What proportion of maximum to surround with convex hull 
            end
        end        

        % Defaults for optional parameters
        if ~isfield(Model, 'Rdim')
            Model.Rdim = 100; % Side length of correlation profile = resolution of backprojection
        end
        if ~isfield(Model, 'Save_Rmaps')
            Model.Save_Rmaps = false; % Whether or not to save correlation profiles in data file
        end
        
    %% Connective field analysis
    case {'samsrf_revcor_cf' 'samsrf_fit_cf'}
        % Mandatory parameters
        if ~isfield(Model, 'Name')
            error('Analysis name is undefined!');
        end
        if ~isfield(Model, 'SeedRoi')
            error('Seed region is undefined!');
        end
        if ~isfield(Model, 'Template')
            error('Template map is undefined!');
        end
        
        % Defaults for optional parameters
        if ~isfield(Model, 'Smoothing')
            Model.Smoothing = 0; % Smoothing kernel (works differently for forward-model & reverse-correlation)
        end
        if ~isfield(Model, 'Noise_Ceiling_Threshold')
            Model.Noise_Ceiling_Threshold = 0; % Limit analysis to data above a certain noise ceiling
        end
        if ~isfield(Model, 'Global_Signal_Correction')
            Model.Global_Signal_Correction = true; % Correct time series by global mean signal
        end
        
        if strcmp(AnalysisFunc, 'samsrf_revcor_cf')
            % Optional parameters for reverse-correlation
            if ~isfield(Model, 'Save_Rmaps')
                Model.Save_Rmaps = false; % Whether or not to save correlation profiles in data file
            end
            if ~isfield(Model, 'Fit_pRF')
                Model.Fit_pRF = 0; % Whether fitting a Gaussian model (1), use convex hull (0), or summary statistics (-1) to estimate CF parameters 
            end
        else 
            % Is forward-model search space defined?
            if ~isfield(Model, 'Polar') && ~isfield(Model, 'Eccentricity') && ~isfield(Model, 'Sizes')
                error('CF search space is not defined!');
            end
            % Warn about overlapping definitions
            if isfield(Model, 'Polar') && (isfield(Model, 'Eccentricity') || isfield(Model, 'Sizes'))
                warning('Ambiguous CF search space specified - using polar angle...');
            end
            if ~isfield(Model, 'Polar') && isfield(Model, 'Eccentricity') && isfield(Model, 'Sizes')
                warning('Ambiguous CF search space specified - using eccentricity...');
            end
            
            % Optional parameters for forward-model
            if ~isfield(Model, 'Coarse_Fit_Block_Size')
                Model.Coarse_Fit_Block_Size = 10000; % Number of simultaneous data columns in coarse fit
            end
            if ~isfield(Model, 'Patch_Size')
                Model.Patch_Size = 0; % Size of receiving patch in geodesic steps (0 = single vertex)
            end
        end      
        
    otherwise
        error('Unknown analysis function specified!');
end