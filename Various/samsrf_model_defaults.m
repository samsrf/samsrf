function Model = samsrf_model_defaults(AnalysisFunc, Model)
%
% Model = samsrf_model_defaults(AnalysisFunc, Model)
%
% Adds defaults of optional parameters for analysis AnalysisFunc (char) 
% if they haven't been defined in Model structure. It also checks that 
% all mandatory parameters (e.g. Name or TR etc) are present.
%
% 14/04/2022 - Written (DSS)
%

%% Which analysis function?
switch AnalysisFunc
    %% Forward-model pRF
    case 'samsrf_fit_prf'
        % Mandatory parameters
        
        % Defaults for optional parameters
        if ~isfield(Model, 'Noise_Ceiling_Threshold')
            Model.Noise_Ceiling_Threshold = 0; % Limit analysis to data above a certain noise ceiling
        end
        if ~isfield(Model, 'Polar_Search_Space')
            Model.Polar_Search_Space = false; % Search space is Cartesian
        end
        if ~isfield(Model, 'Seed_Fine_Fit')
            Model.Seed_Fine_Fit = ''; % Use no seed map, so run coarse fit instead
        end 
        if ~isfield(Model, 'Replace_Bad_Fits')
            Model.Replace_Bad_Fits = false; % Don't replace bad fine fits with coarse fit
        end
        if ~isfield(Model, 'Smoothed_Coarse_Fit')
            Model.Smoothed_Coarse_Fit = 0; % No smoothing on coarse fit
        end
        if ~isfield(Model, 'Coarse_Fit_Only')
            Model.Coarse_Fit_Only = false; % Only run coarse fit & then save
        end
        if ~isfield(Model, 'Fine_Fit_Threshold')
            Model.Fine_Fit_Threshold = 0.01; % Include coarse fits with R^2>0.01 in fine fit
        end
        if ~isfield(Model, 'Only_Positive_Coarse_Fits')
            Model.Only_Positive_Coarse_Fits = false; % Coarse fit can either be negative or positive correlation to pass 
        end
        if ~isfield(Model, 'Coarse_Fit_Block_Size')
            Model.Coarse_Fit_Block_Size = 10000; % Number of simultaneous data columns in coarse fit
        end
        if ~isfield(Model, 'Downsample_Predictions')
            Model.Downsample_Predictions = 1; % Downsampling factor by which Model.TR mismatches the true TR
        end
        if ~isfield(Model, 'Coarse_Fit_Percentile')
            Model.Coarse_Fit_Percentile = 100; % Which percentile of correlations to include in coarse fit parameter estimates
        end
        
    %% Reverse-correlation pRF
    case 'samsrf_revcor_prf'
        % Mandatory parameters

        % Defaults for optional parameters
        if ~isfield(Model, 'Noise_Ceiling_Threshold')
            Model.Noise_Ceiling_Threshold = 0; % Limit analysis to data above a certain noise ceiling
        end
        if ~isfield(Model, 'Save_Rmaps')
            Model.Save_Rmaps = true; % Whether or not to save correlation profiles in data file
        end

    %% Reverse-correlation CF
    case 'samsrf_revcor_cf'
        % Mandatory parameters

        % Defaults for optional parameters
        if ~isfield(Model, 'Noise_Ceiling_Threshold')
            Model.Noise_Ceiling_Threshold = 0; % Limit analysis to data above a certain noise ceiling
        end
        if ~isfield(Model, 'Global_Signal_Correction')
            Model.Global_Signal_Correction = true; % Correct time series by global mean signal
        end
        if ~isfield(Model, 'Save_Rmaps')
            Model.Save_Rmaps = false; % Whether or not to save correlation profiles in data file
        end
        if ~isfield(Model, 'Fit_pRF')
            Model.Fit_pRF = true;
        end
        
    %% Forward-model CF
    case 'samsrf_fit_cf'
        % Mandatory parameters

        % Defaults for optional parameters
        if ~isfield(Model, 'Noise_Ceiling_Threshold')
            Model.Noise_Ceiling_Threshold = 0; % Limit analysis to data above a certain noise ceiling
        end
        if ~isfield(Model, 'Smoothing')
            Model.Smoothing = 0; % No smoothing on coarse fit
        end
        if ~isfield(Model, 'Coarse_Fit_Block_Size')
            Model.Coarse_Fit_Block_Size = 10000; % Number of simultaneous data columns in coarse fit
        end
        if ~isfield(Model, 'Patch_Size')
            Model.Patch_Size = 0; % Size of receiving patch in geodesic steps (0 = single vertex)
        end
        
    otherwise
        error('Unknown analysis function specified!');
end