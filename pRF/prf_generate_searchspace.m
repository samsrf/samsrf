function [Ptc, S] = prf_generate_searchspace(PrfFcn, ApFrm, ApXY, Param1, Param2, Param3, Param4, Param5, Param6, Param7, Param8, Param9, Param10, IsCircSpace)
%
% [Ptc, S] = prf_generate_searchspace(PrfFcn, ApFrm, ApXY, Param1, Param2, Param3, Param4, Param5, Param6, Param7, Param8, Param9, Param10, [IsCircSpace=false])
%
% Search grid Ptc for the coarse fit using the pRF model PrfFcn as basis.
% Each column of Ptc contains a predicted time course for a point in the search grid.
% Each column in S contains the parameters for the predictions in Ptc.
%
% ApFrm contains the aperture vectors (rows = pixels, columns = volumes).
% ApXY is a two-column matrix with the pixel X & Y coordinates from aperture file.
% Param1 to Param10 define the values to define the grid. 
% This maxes out at 10 dimensions for complex models like Divisive Normalisation
%  although to be honest beyond five parameters only madness lies... 
% For models with fewer dimensions you need to set the relevant input to a scalar 0
%  (but the default parameter function does this automatically).
%
% If IsCircSpace is true (default=false), then Param1 and Param2 are treated as
% polar angle (in degrees) and ecceendntricity (in aperture space), respectively. 
% These are then internally converted into Cartesian coordinates.
%
% 06/07/2022 - Rewritten for vectorised apertures (DSS)
% 07/07/2022 - Added option for five more parameters (DSS)
%

if nargin < 14
    IsCircSpace = false;
end

% Prediction matrix
Ptc = NaN(size(ApFrm,2), length(Param1)*length(Param2)*length(Param3)*length(Param4)*length(Param5)*length(Param6)*length(Param7)*length(Param8)*length(Param9)*length(Param10));

% Search grid (1 & 2 in reverse order because MatLab is stupid)
[S2, S1, S3, S4, S5, S6, S7, S8, S9, S10] = ndgrid(Param2, Param1, Param3, Param4, Param5, Param6, Param7, Param8, Param9, Param10);

% Polar search space?
if IsCircSpace 
    samsrf_disp(' Search space in polar coordinates.');
    [S1, S2] = pol2cart(S1/180*pi, S2);
else
    samsrf_disp(' Search space in Cartesian coordinates.');
end

% Parameter matrix
S = [S1(:) S2(:) S3(:) S4(:) S5(:) S6(:) S7(:) S8(:) S9(:) S10(:)]'; 

% Generating predictions
samsrf_disp(' Please stand by...');
for n = 1:numel(S1)
% parfor n = 1:numel(S1)
    Rfp = PrfFcn([S1(n) S2(n) S3(n) S4(n) S5(n) S6(n) S7(n) S8(n) S9(n) S10(n)], ApXY); % pRF profile 
    cptc = prf_predict_timecourse(Rfp, ApFrm); % Prediction is in percent of pRF activated
    Ptc(:,n) = cptc; 
end
