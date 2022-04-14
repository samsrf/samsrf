function HelpText = ModelHelpText(AnalysisFunc, ParameterName)
%
% HelpText = ModelHelpText(AnalysisFunc, ParameterName)
%
% Returns the help text for a given analysis function & parameter name.
% This is called internally by ModelHelp but you can also use it directly.
%
% 15/04/2022 - Written (DSS)
%

%% Which analysis function?
switch AnalysisFunc
    %% Forward model pRF
   case 'samsrf_fit_prf'
       switch ParameterName
           case 'Prf_Function'
               HelpText = {'Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit' 
                   'in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.'};
           case 'Name'
               HelpText = '';
           case 'Param_Names'
               HelpText = '';
           case 'Scaled_Param'
               HelpText = '';
           case 'Only_Positive'
               HelpText = '';
           case 'Scaling_Factor'
               HelpText = '';
           case 'TR'
               HelpText = '';
           case 'Hrf'
               HelpText = '';
           case 'Aperture_File'
               HelpText = '';
           case 'Noise_Ceiling_Threshold'
               HelpText = '';
           case 'Polar_Search_Space'
               HelpText = '';
           case 'Seed_Fine_Fit'
               HelpText = '';
           case 'Replace_Bad_Fits'
               HelpText = '';
           case 'Smoothed_Coarse_Fit'
               HelpText = '';
           case 'Coarse_Fit_Only'
               HelpText = '';
           case 'Fine_Fit_Threshold'
               HelpText = '';
           case 'Only_Positive_Coarse_Fits'
               HelpText = '';
           case 'Coarse_Fit_Block_Size'
               HelpText = '';
           case 'Downsample_Predictions'
               HelpText = '';
           case 'Coarse_Fit_Percentile'
               HelpText = '';
           case 'Hooke_Jeeves_Steps'
               HelpText = '';
           case 'Nelder_Mead_Tolerance'  
               HelpText = '';
               
           otherwise
               error('samsrf_fit_prf does not have this parameter!');
       end
    
    %% Reverse correlation pRF                
   case 'samsrf_revcor_prf'
       switch ParameterName
           case 'Name'
               HelpText = '';
           case 'Scaling_Factor'
               HelpText = '';
           case 'TR'
               HelpText = '';
           case 'Hrf'
               HelpText = '';
           case 'Aperture_File'
               HelpText = '';
           case 'Prf_Function'
               HelpText = '';
           case 'Param_Names'
               HelpText = '';
           case 'Scaled_Param'
               HelpText = '';
           case 'SeedPar_Function'
               HelpText = '';
           case 'R2_Threshold'
               HelpText = '';
           case 'Rdim'
               HelpText = '';
           case 'Noise_Ceiling_Threshold'
               HelpText = '';
           case 'Save_Rmaps'
               HelpText = '';
           case 'Hooke_Jeeves_Steps'
               HelpText = '';
           case 'Nelder_Mead_Tolerance'  
               HelpText = '';
               
           otherwise
               error('samsrf_revcor_prf does not have this parameter!');
       end

    %% Reverse correlation CF    
   case 'samsrf_revcor_cf'
       switch ParameterName
           case 'Name'
               HelpText = '';
           case 'SeedRoi'
               HelpText = '';
           case 'Template'
               HelpText = '';
           case 'Smoothing'
               HelpText = '';
           case 'Noise_Ceiling_Threshold'
               HelpText = '';
           case 'Global_Signal_Correction'
               HelpText = '';
           case 'Save_Rmaps'
               HelpText = '';
           case 'Fit_pRF'
               HelpText = '';
           case 'Hooke_Jeeves_Steps'
               HelpText = '';
           case 'Nelder_Mead_Tolerance'  
               HelpText = '';
               
           otherwise
               error('samsrf_revcor_cf does not have this parameter!');
       end
            
    %% Forward model CF    
   case 'samsrf_fit_cf'
       switch ParameterName
           case 'Name'
               HelpText = '';
           case 'SeedRoi'
               HelpText = '';
           case 'Template'
               HelpText = '';
           case 'Smoothing'
               HelpText = '';
           case 'Noise_Ceiling_Threshold'
               HelpText = '';
           case 'Global_Signal_Correction'
               HelpText = '';
           case 'Polar'
               HelpText = '';
           case 'Eccentricity'
               HelpText = '';
           case 'Sizes'
               HelpText = '';
           case 'Coarse_Fit_Block_Size'
               HelpText = '';
           case 'Patch_Size'  
               HelpText = '';
               
           otherwise
               error('samsrf_fit_cf does not have this parameter!');
       end
               
    otherwise
        error('Unknown analysis function specified!')                       
end
