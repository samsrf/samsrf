
function Y = prf_predict_timecourse(Rfp, ApFrm)
%
% Y = prf_predict_timecourse(Rfp, ApFrm)
%
% Predicts the time course resulting from receptive field profile Rfp and stimulus mask movie ApFrm. 
%
% 16/09/2024 - Now uses absolute values to calculate pRF overlap to correct for issues with negative components (DSS)
%

% Which biophysical model?
if isnan(Rfp(end))
    Rfp = Rfp(1:end-1);
    PrfOvr = false;
else
    PrfOvr = true;
end

% Output time course vector
Y = NaN(size(ApFrm,2),1); 

% Predict response at each time point
for i = 1:size(ApFrm,2)
    Y(i) = mean(Rfp .* ApFrm(:,i)); % Mean overlap between pRF & aperture
end

% Use overlap model?
if PrfOvr
    Y = Y / mean(abs(Rfp)) * 100; % Transform into percent pRF overlap 
end
