function ParamList = ModelHelpList(AnalysisFunc)
%
% ParamList = ModelHelpList(AnalysisFunc)
%
% Returns a cell array with the possible model parameters for a given analysis function. 
% This function is called internally by ModelHelp but you can also use it directly.
%
% 15/04/2022 - Written (DSS)
%

%% Which analysis function?
switch AnalysisFunc
    % Forward model pRF
    case 'samsrf_fit_prf'
        ParamList = {   'Prf_Function';
                        'Name';
                        'Param_Names';
                        'Scaled_Param';
                        'Only_Positive';
                        'Scaling_Factor';
                        'TR';
                        'Hrf';
                        'Aperture_File';
                        'Noise_Ceiling_Threshold';
                        'Polar_Search_Space';
                        'Seed_Fine_Fit';
                        'Replace_Bad_Fits';
                        'Smoothed_Coarse_Fit';
                        'Coarse_Fit_Only';
                        'Fine_Fit_Threshold';
                        'Only_Positive_Coarse_Fits';
                        'Coarse_Fit_Block_Size';
                        'Downsample_Predictions';
                        'Coarse_Fit_Percentile';
                        'Hooke_Jeeves_Steps';
                        'Nelder_Mead_Tolerance'  };
    
    % Reverse correlation pRF                
    case 'samsrf_revcor_prf'
        ParamList = {   'Name';
                        'Scaling_Factor';
                        'TR';
                        'Hrf';
                        'Aperture_File';
                        'Prf_Function';
                        'Param_Names';
                        'Scaled_Param';
                        'SeedPar_Function';
                        'R2_Threshold';
                        'Rdim';
                        'Noise_Ceiling_Threshold';
                        'Save_Rmaps';
                        'Hooke_Jeeves_Steps';
                        'Nelder_Mead_Tolerance'  };

    % Reverse correlation CF    
    case 'samsrf_revcor_cf'
        ParamList = {   'Name';
                        'SeedRoi';
                        'Template';
                        'Smoothing';
                        'Noise_Ceiling_Threshold';
                        'Global_Signal_Correction';
                        'Save_Rmaps';
                        'Fit_pRF';
                        'Hooke_Jeeves_Steps';
                        'Nelder_Mead_Tolerance'  };
            
    % Forward model CF    
    case 'samsrf_fit_cf'
        ParamList = {   'Name';
                        'SeedRoi';
                        'Template';
                        'Smoothing';
                        'Noise_Ceiling_Threshold';
                        'Global_Signal_Correction';
                        'Polar';
                        'Eccentricity';
                        'Sizes';
                        'Coarse_Fit_Block_Size';
                        'Patch_Size'  };
    otherwise
             error('Unknown analysis function specified!');                       
end
