function Srf = samsrf_dog_fwhm(Srf)
%
% Srf = samsrf_dog_fwhm(Srf)
%
% Calculates the parameters of a Difference-of-Gaussians pRF model into 
% the full-width at half-maximum and adds this to the Srf.
%
% 27/05/2020 - Streamlined how waitbar is handled (DSS)
%

% Expand Srf if necessary
Srf = samsrf_expand_srf(Srf);

% Calculate FWHM from DoG parameters
h = samsrf_waitbar('Calculating FWHM...'); 
Fwhm = []; 
for v = 1:size(Srf.Data,2)
   Fwhm = [Fwhm dog_fwhm(Srf.Data(4,v), Srf.Data(5,v), Srf.Data(6,v))];
   samsrf_waitbar(v/size(Srf.Data,2), h); 
end 
Srf.Data = [Srf.Data; Fwhm];
Srf.Values{end+1} = 'Fwhm';
samsrf_waitbar('', h); 


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
