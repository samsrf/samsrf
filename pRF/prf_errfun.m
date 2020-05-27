function err = prf_errfun(PrfFcn, ApFrm, Hrf, P, Y)
%
% err = prf_errfun(PrfFcn, ApFrm, Hrf, P, Y)
%
% Error function for pRF model fitting optimisation.
%
% Returns unexplained variance (1-R^2) for the correlation between the 
% prediction defined by the parameters P and the observed data in Y. 
%
%   PrfFcn is a function handle to the pRF model (e.g. prf_gaussian_rf)
%   ApFrm contains aperture frames. 
%   Hrf is the hemodynamic response function
%
% 23/04/2018 - SamSrf 6 version (DSS)
%

Rfp = PrfFcn(P, size(ApFrm,1)*2); % pRF profile
Yp = prf_predict_timecourse(Rfp, ApFrm, true); % Predict time course with z-normalisation
Yp = conv(Yp, Hrf); % Convolve with HRF
Yp = Yp(1:length(Y)); % Truncate back to original length
R = corr(Yp, Y); % Correlate prediction with observed data
err = 1 - R.^2; % Unexplained variance 
    