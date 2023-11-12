function err = prf_errfun(PrfFcn, ApFrm, ApXY, Hrf, P, Y, Downsampling, CssNonlin)
%
% err = prf_errfun(PrfFcn, ApFrm, Hrf, P, Y, Downsamping, CssNonlin)
%
% Error function for pRF model fitting optimisation.
%
% Returns unexplained variance (1-R^2) for the correlation between the 
% prediction defined by the parameters P and the observed data in Y. 
%
%   PrfFcn is a function handle to the pRF model (e.g. prf_gaussian_rf)
%   ApFrm contains aperture vectors
%   ApXY contains pixel coordinates for apertures
%   Hrf is the hemodynamic response function
%     If a negative scalar, the function fits the HRF using the absolute of this as TR
%     The HRF parameters are constrained within a plausible range
%   Downsampling is the factor by which the prediction is sampled to match the scanner TR
%   CssNonlin is a boolean to flag whether to model compressive spatial summation
%
% 06/07/2022 - Rewritten for vectorised apertures (DSS)
% 11/11/2023 - Added support for fitting HRFs (DSS)
% 13/11/2023 - HRF parameters are now constrained (DSS)
%

if isscalar(Hrf) && Hrf < 0 && ...
        (P(end-4) < 3 || P(end-3) < 10 || P(end-2) < .5 || P(end-1) < .5 || P(end) < 0 || ...
        P(end-4) > 9 || P(end-3) > 18 || P(end-2) > 3 || P(end-1) > 3 || P(end) > 10) 
    % Fitting HRF & parameters out of bounds
    err = 1;
else    
    % Continue fitting
    Rfp = PrfFcn(P, ApXY); % pRF profile
    Yp = prf_predict_timecourse(Rfp, ApFrm); % Predict time course (percent pRF activation)
    % Which HRF to use?
    if isscalar(Hrf) && Hrf < 0
        cssp = length(P) - 5; % If CSS model, sixth last parameter is exponent 
        Hrf = samsrf_doublegamma(-Hrf, P(end-4:end)); % Create HRF with these parameters & TR defined by |Hrf|
    else
        cssp = length(P); % If CSS model, last parameter is exponent
    end
    % Model CSS nonlinearity
    if CssNonlin
        Yp = Yp.^P(cssp);
    end
    % Convolve with HRF
    Yp = prf_convolve_hrf(Yp, Hrf, Downsampling); 
    R = corr(Yp, Y); % Correlate prediction with observed data
    err = 1 - R.^2; % Unexplained variance 
end    