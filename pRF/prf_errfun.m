function err = prf_errfun(PrfFcn, ApFrm, Hrf, P, Y, Downsampling, CssNonlin)
%
% err = prf_errfun(PrfFcn, ApFrm, Hrf, P, Y, Downsamping, CssNonlin)
%
% Error function for pRF model fitting optimisation.
%
% Returns unexplained variance (1-R^2) for the correlation between the 
% prediction defined by the parameters P and the observed data in Y. 
%
%   PrfFcn is a function handle to the pRF model (e.g. prf_gaussian_rf)
%   ApFrm contains aperture frames 
%   Hrf is the hemodynamic response function
%   Downsampling is the factor by which the prediction is sampled to match the scanner TR
%   CssNonlin is a boolean to flag whether to model compressive spatial summation
%
% 20/04/2022 - SamSrf 8 version (DSS)
% 24/06/2022 - Added option to model compressive spatial summary nonlinearity (DSS)
%              Removed obsolete comment about z-scoring (DSS)
%

Rfp = PrfFcn(P, size(ApFrm,1)*2); % pRF profile
Yp = prf_predict_timecourse(Rfp, ApFrm); % Predict time course (percent pRF activation)
% Model CSS nonlinearity
if CssNonlin
    Yp = Yp.^P(end);
end
Yp = prf_convolve_hrf(Yp, Hrf, Downsampling); % Convolve with HRF
R = corr(Yp, Y); % Correlate prediction with observed data
err = 1 - R.^2; % Unexplained variance 
    