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
disp('Calculating FWHM from DoG parameters...');
Fwhm = []; 
for v = 1:size(Srf.Data,2)
   Fwhm = [Fwhm dog_fwhm(Srf.Data(4,v), Srf.Data(5,v), Srf.Data(6,v))];
end 
Srf.Data = [Srf.Data; Fwhm];
Srf.Values{end+1} = 'Fwhm';


function fwhm = dog_fwhm(sigma1, sigma2, rat)
%
% fwhm = dog_fwhm(sigma1, sigma2, rat)
%
% Returns the full width at half maximum of a difference of gaussian function
% with parameters sigma1 (inner), sigma2 (outer) and the 2/1 beta ratio rat.
% This uses the fzero function.
%

% Full maximum of the DoG
fmax = 1 - rat;

% FWHM of the DoG
if sigma1 > 0 && sigma2 > 0
    fwhm = 2*abs(fzero(@(x) exp(-x.^2/(2*sigma1.^2)) - rat*exp(-x.^2/(2*sigma2.^2)) - fmax/2, 0));
else
    fwhm = 0;
end
