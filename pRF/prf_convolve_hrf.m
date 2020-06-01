
function cY = prf_convolve_hrf(Y, Hrf)
%
% cY = prf_predict_timecourse(Y, Hrf)
%
% Convolves the time course Y with the hemodynamic response function Hrf and truncates back to original length. 
%
% 02/06/2020 - SamSrf 7 version (DSS) 
%

% Convolve with HRF
cY = conv(Y, Hrf);
% Truncate to original length
cY = cY(1:length(Y));