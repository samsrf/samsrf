function ParamList = ModelHelpList(AnalysisFunc)
%
% ParamList = ModelHelpList(AnalysisFunc)
%
% Returns a cell array with the possible model parameters for a given analysis function. 
% This function is called internally by ModelHelp but you can also use it directly.
%

%% Which analysis function?
switch AnalysisFunc
    % Forward model pRF
    case 'samsrf_fit_prf'
        ParamList = {   'Name';
						'Prf_Function';                        
                        'Param_Names';
                        'Scaled_Param';
                        'Only_Positive';
                        'Scaling_Factor';
                        'TR';
                        'Hrf';
                        'Aperture_File';
                        'Aperture_Mean_Response';
                        'Noise_Ceiling_Threshold';
                        'Polar_Search_Space';
						'Param1';
						'Param2';
						'Param3';
						'Param4';
						'Param5';
						'Param6';
						'Param7';
						'Param8';
						'Param9';
						'Param10';
                        'Coarse_Fit_Only';
                        'Smoothed_Coarse_Fit';
                        'Only_Positive_Coarse_Fits';
                        'Coarse_Fit_Block_Size';
                        'Coarse_Fit_Percentile';
                        'Seed_Fine_Fit';
                        'Fine_Fit_Threshold';
                        'Replace_Bad_Fits';
                        'Downsample_Predictions';
                        'Compressive_Nonlinearity';
                        'Hooke_Jeeves_Steps';
                        'Nelder_Mead_Tolerance'  
                        'SeedRoi';
                        'Template'  };
    
    % Reverse correlation pRF                
    case 'samsrf_revcor_prf'
        ParamList = {   'Name';
                        'Scaling_Factor';
                        'TR';
                        'Hrf';
                        'Aperture_File';
                        'Noise_Ceiling_Threshold';
                        'Prf_Function';
                        'Convex_Hull_Threshold';
                        'Param_Names';
                        'Scaled_Param';
                        'SeedPar_Function';
                        'R2_Threshold';
                        'Save_Rmaps';
                        'Allow_Negative_Peaks';
                        'Rdim';
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
