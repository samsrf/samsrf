function [Ptc, S] = prf_generate_searchspace(PrfFcn, ApFrm, Param1, Param2, Param3, Param4, Param5, IsCircSpace)
%
% [Ptc, S] = prf_generate_searchspace(PrfFcn, ApFrm, Param1, Param2, Param3, Param4, Param5, [IsCircSpace=false])
%
% Search grid Ptc for the coarse fit using the pRF model PrfFcn as basis.
% Each column of Ptc contains a predicted time course for a point in the search grid.
% Each column in S contains the parameters for the predictions in Ptc.
%
% ApFrm is the aperture movie. Param1 to Param5 define the values to define the grid. 
% This maxes out at 5 dimensions because beyond five parameters only madness lies. 
% For models with fewer dimensions you need to set the relevant input to a scalar 0.
%
% If IsCircSpace is true (default=false), then Param1 and Param2 are treated as
% polar angle (in degrees) and ecceendntricity (in aperture space), respectively. 
% These are then internally converted into Cartesian coordinates.
%
% 20/04/2022 - SamSrf 8 version (DSS)
%

if nargin < 8
    IsCircSpace = false;
end

% Prediction matrix
Ptc = NaN(size(ApFrm,3), length(Param1)*length(Param2)*length(Param3)*length(Param4)*length(Param5));

% Search grid (1 & 2 in reverse order because MatLab is stupid)
[S2, S1, S3, S4, S5] = ndgrid(Param2, Param1, Param3, Param4, Param5);

% Polar search space?
if IsCircSpace 
    disp(' Search space in polar coordinates.');
    [S1, S2] = pol2cart(S1/180*pi, S2);
else
    disp(' Search space in Cartesian coordinates.');
end

% Parameter matrix
S = [S1(:) S2(:) S3(:) S4(:) S5(:)]'; 

% Generating predictions
disp(' Please stand by...');
parfor n = 1:numel(S1)
    Rfp = PrfFcn([S1(n) S2(n) S3(n) S4(n) S5(n)], size(ApFrm,1)*2); % pRF profile
    cptc = prf_predict_timecourse(Rfp, ApFrm); % Prediction is in percent of pRF activated
    Ptc(:,n) = cptc; 
end
