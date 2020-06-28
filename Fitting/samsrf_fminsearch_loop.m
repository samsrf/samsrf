function [fPimg, fRimg] = samsrf_fminsearch_loop(Model, Y, ApFrm, Rimg, Pimg, mver) 
%
% [fPimg, fRimg] = samsrf_fminsearch_loop(Model, Y, ApFrm, Rimg, Pimg, mver) 
%
% Internal function looping through solving the fine fit for each mask vertex.
%
% Loops thru solving fine fit. 'Tis necessary to ensure the parallel computing 
% toolbox accepts it. Go to hell, Matlab!
%
%   Model:  Model structure for the pRF fit
%   Y:      Matrix with time series for every vertex in the mask
%   ApFrm:  Apertures with stimulus masks
%   Rimg:   Coarse fit/seed map R^2 data for thresholding
%   Pimg:   Coarse fit/seed map parameters
%   mver:   Mask vertices
%
% 29/06/2020 - SamSrf 7 version (DSS)
%

% Output matrices
fPimg = zeros(length(Model.Param_Names), length(mver)); % Parameter maps
fRimg = zeros(1,length(mver)); % R^2 map

% Display off & default tolerances (=slow)
OptimOpts = optimset('Display', 'off'); 

% Loop thru mask vertices 
parfor v = 1:length(mver)
    if Rimg(mver(v)) >= Model.Fine_Fit_Threshold % Only reasonable coarse fits
        % Find best prediction
        [fP,fR] = fminsearch(@(P) prf_errfun(Model.Prf_Function, ApFrm, Model.Hrf, P, Y(:,v)), Pimg(:,mver(v))', OptimOpts);  
        fPimg(:,v) = fP;
        fRimg(1,v) = 1 - fR;
    end
end