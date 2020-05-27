function [Ptc, S] = prf_generate_searchspace(PrfFcn, ApFrm, Param1, Param2, Param3, Param4, Param5, wb)
%
% [Ptc, S] = prf_generate_searchspace(PrfFcn, ApFrm, Param1, Param2, Param3, Param4, Param5, [wb])
%
% Search grid Ptc for the coarse fit using the pRF model PrfFcn as basis.
% Each column of Ptc contains a predicted time course for a point in the search grid.
% Each column in S contains the parameters for the predictions in Ptc.
%
% ApFrm is the aperture movie. Param1 to Param5 define the values to define the grid. 
% This maxes out at 5 dimensions because beyond five parameters only madness lies. 
% For models with fewer dimensions you need to set the relevant input to a scalar 0.
%
% If wb is true (default) a waitbar is displayed. 
%
% 30/05/2018 SamSrf 6 version (DSS)
% 27/05/2020 - Streamlined how waitbar is handled (DSS)
%

if nargin < 6
    wb = true;
end

% Prediction matrix
Ptc = NaN(size(ApFrm,3), length(Param1)*length(Param2)*length(Param3)*length(Param4)*length(Param5));
% Search grid (1 & 2 in reverse order because MatLab is stupid)
[S2, S1, S3, S4, S5] = ndgrid(Param2, Param1, Param3, Param4, Param5);
S = [S1(:) S2(:) S3(:) S4(:) S5(:)]';

% Generating predictions
if wb
    % Only allow waitbar if default parameters allow it
    h = samsrf_waitbar('Generating predictions...');
else
    h = 0;
end
for n = 1:numel(S1)
    Rfp = PrfFcn([S1(n) S2(n) S3(n) S4(n) S5(n)], size(ApFrm,1)*2); % pRF profile
    cptc = prf_predict_timecourse(Rfp, ApFrm, true); % Prediction is z-normalised!
    Ptc(:,n) = cptc; 
    samsrf_waitbar(n/numel(S1), h);
end
samsrf_waitbar('',h);
