
function Y = prf_predict_timecourse_multicond(Rfp, ApFrm)
%
% Y = prf_predict_timecourse_multicond(Rfp, ApFrm)
%
% Predicts the time course resulting from receptive field profile Rfp & 
%  stimulus mask movie ApFrm. Works with multiple conditions. If there are
%  multiple conditions, each column of Rfp must be for one conition & the
%  1st row of ApFrm contains the corresponding condition indeces per volume.
%
% WARNING: This is a version of this function supporting multiple conditions 
%  within the same time series. This is useful but it is slower than the 
%  standard fitting algorithm in prf_predict_timecourse. If you want to use 
%  this, you need to rename the old function, for example by calling it 
%  prf_predict_timecourse_old & this one into prf_predict_timecourse.
%
%  (You should then keep a note of this for future SamSr updates as that 
%  would again overwrite this change...)
%
% 01/08/2024 - New functionality for multiple conditions within a time series (DSS)
%              This uses matrix operations to do away with for-loops (DSS)
% 02/08/2024 - Updated help section only (DSS)
%

% Multi-condition design?
Ncon = size(Rfp,2);
if Ncon > 1
	ApCond = ApFrm(1,:); % 1st row is condition index
	ApFrm = ApFrm(2:end,:); % Remaining rows are stimulus apertures
else
	ApCond = ones(1,size(ApFrm,2)); % Dummy condition index
end

% Which biophysical model?
if isnan(Rfp(end,1))
    Rfp = Rfp(1:end-1,:);
    PrfOvr = false;
else
    PrfOvr = true;
end

% Replicate pRF profile
Rfp = Rfp(:,ApCond);

% Predict response at each time point
Y = mean(Rfp .* ApFrm)';
% Use overlap model?
if PrfOvr
	Y = Y ./ mean(Rfp)' * 100; % Transform into percent pRF overlap 
end

