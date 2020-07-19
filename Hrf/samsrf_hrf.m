function Hrf = samsrf_hrf(TR)
%
% Hrf = samsrf_hrf(TR)
%
% Returns a canonical HRF with sampling rate TR. This HRF is based on the 
% data in de Haas et al. (2014) Curr.Biol. collected on a Siemens Trio 3T
% with 2.55s TR and 2.3mm isotropic resolution. In 10 trials, sparse full 
% field ripple stimuli were presented for 1 TR followed by 11 TRs of blank. 
%
% 19/07/2020 - SamSrf 7 version (DSS)
%

% SPM's canonical HRF
Can = samsrf_doublegamma(TR, [6 16 1 1 6 0 32]);
Can = Can' / sum(Can);

% Fitted HRF parameters 
Params = [5.5278, 16.8621, 1.0205, 1.466];

% Return HRF sampled at this TR
Hrf = samsrf_doublegamma(TR, [Params(1:2) 1 1 Params(3) 0 32]) * Params(4);

% Normalise relative to SPM's canonical
Hrf = max(Can) / max(Hrf) * Hrf;


