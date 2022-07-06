
function Y = prf_predict_timecourse(Rfp, ApFrm)
%
% Y = prf_predict_timecourse(Rfp, ApFrm)
%
% Predicts the time course resulting from receptive field profile Rfp and stimulus mask movie ApFrm. 
%
% 06/07/2022 - Rewritten for vectorised apertures (DSS)
%

% Output time course vector
Y = NaN(size(ApFrm,2),1); 

% Predict response at each time point
for i = 1:size(ApFrm,2)
    Y(i) = mean(Rfp .* ApFrm(:,i)); % Mean overlap between pRF & aperture
end
Y = Y / mean(Rfp) * 100; % Transform into percent pRF overlap 

