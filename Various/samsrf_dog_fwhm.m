function Srf = samsrf_dog_fwhm(Srf)
%
% Srf = samsrf_dog_fwhm(Srf)
%
% Calculates the parameters of a Difference-of-Gaussians pRF model into 
% the full-width at half-maximum and adds this to the Srf.
%
% 20/04/2022 - SamSrf 8 version (DSS)
%

% Expand Srf if necessary
Srf = samsrf_expand_srf(Srf);

% Calculate FWHM from DoG parameters
samsrf_disp('Calculating FWHM from DoG parameters...');
Fwhm = []; 
for v = 1:size(Srf.Data,2)
   Fwhm = [Fwhm dog_fwhm(Srf.Data(4,v), Srf.Data(5,v), Srf.Data(6,v))];
end 
Srf.Data = [Srf.Data; Fwhm];
Srf.Values{end+1} = 'Fwhm';


