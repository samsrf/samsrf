
function Y = prf_predict_timecourse(Rfp, ApFrm)
%
% Y = prf_predict_timecourse(Rfp, ApFrm)
%
% Predicts the time course resulting from receptive field profile Rfp and stimulus mask movie ApFrm. 
%
% 06/07/2022 - Rewritten for vectorised apertures (DSS)
% 17/03/2023 - Added option for conventional Dumoulin & Wandell 2008 biophysical model (DSS)
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
    Y = Y / mean(Rfp) * 100; % Transform into percent pRF overlap 
end
