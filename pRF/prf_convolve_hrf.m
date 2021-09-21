
function cY = prf_convolve_hrf(Y, Hrf, Downsampling)
%
% cY = prf_predict_timecourse(Y, Hrf, Downsampling)
%
% Convolves the time course Y with the hemodynamic response function Hrf and truncates back to original length. 
% If downsampling is defined & greater than 1 it also downsamples the predicted time course to the true TR.
%
% 02/06/2020 - SamSrf 7 version (DSS) 
% 22/09/2021 - Added option for downsampling predictions (DSS)
%

% Convolve with HRF
cY = conv(Y, Hrf);
% Truncate to original length
cY = cY(1:length(Y));
% Downsample timeseries?
if Downsampling > 1
    cY = cY(1:Downsampling:end);
end