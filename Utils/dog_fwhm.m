function fwhm = dog_fwhm(sigma1, sigma2, rat)
%
% fwhm = dog_fwhm(sigma1, sigma2, rat)
%
% Returns the full width at half maximum of a difference of gaussian function
% with parameters sigma1 (inner), sigma2 (outer) and the 2/1 beta ratio rat.
% This uses the fzero function.
%
% 28/08/2026 - Split this from samsrf_dog_fwhm as it can be useful (DSS)
%

% Full maximum of the DoG
fmax = 1 - rat;

% FWHM of the DoG
if sigma1 > 0 && sigma2 > 0
    fwhm = 2*abs(fzero(@(x) exp(-x.^2/(2*sigma1.^2)) - rat*exp(-x.^2/(2*sigma2.^2)) - fmax/2, 0));
else
    fwhm = 0;
end
