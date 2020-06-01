
function Y = prf_predict_timecourse(Rfp, ApFrm)
%
% Y = prf_predict_timecourse(Rfp, ApFrm)
%
% Predicts the time course resulting from receptive field profile Rfp and stimulus mask movie ApFrm. 
%
% 02/06/2020 - SamSrf 7 version (DSS) 
%

% Output time course vector
Y = NaN(size(ApFrm,3),1); 

% Predict response at each time point
for i = 1:size(ApFrm,3)
    Y(i) = prf_predict_response(ApFrm(:,:,i),Rfp); 
end