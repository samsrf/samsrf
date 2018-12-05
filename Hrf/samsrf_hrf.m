function Hrf = samsrf_hrf(TR)
%
% Hrf = samsrf_hrf(TR)
%
% Returns a canonical HRF with sampling rate TR. This HRF is based on the 
% data in de Haas et al. (2014) Curr.Biol. collected on a Siemens Trio 3T
% with 2.55s TR and 2.3mm isotropic resolution. In 10 trials, sparse full 
% field ripple stimuli were presented for 1 TR followed by 11 TRs of blank. 
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


function [Hrf, p] = samsrf_doublegamma(TR, P)
% [Hrf, p] = samsrf_doublegamma(TR, P)
%
% Inputs
%   TR          [scalar]  Repetition time, in seconds
%   P           [vector]  Double gamma parameters, in seconds
%       P(1) - delay of response (relative to onset)	   (default = 6 )
%       P(2) - delay of undershoot (relative to onset)     (default = 16)
%       P(3) - dispersion of response                      (default = 1 )
%       P(4) - dispersion of undershoot                    (default = 1 )
%       P(5) - ratio of response to undershoot             (default = 6 )
%       P(6) - onset (seconds)                             (default = 0 )
%       P(7) - length of kernel (seconds)                  (default = 32)
%
% Generate a BOLD haemodynamic response function from discrete samples of
% the Gamma probability density function. All values in seconds. Based on
% spm_hrf.m. 
%
% % Changelog %
%
% 21/06/2018    Written
%
% Ivan Alvarez
% FMRIB Centre, University of Oxford


%% Main

% Default values
Grain = 16;
p = [6 16 1 1 6 0 32];

% Replace defaults with user-defined parameters, where necessary 
if nargin == 2
    p(1:length(P)) = P;
end

% Temporal sampling space
L = length((p(6):p(7)/TR) * Grain + 1);
T = linspace(p(6), p(7), L);

% The Gamma function is defined by two parameters: shape (A) and scale (B).
% We define A as the ratio between the delay of the response (lag) and the
% dispersion of the response. We define B as the dispersion of the response.
A = p(1) / p(3);
B = p(3);
g1 = gampdf(T, A, B);

% Next, we generate the undershoot and scale it according to the P(5) ratio
A = p(2) / p(4);
B = p(4);
g2 = gampdf(T, A, B);
g2 = g2 / p(5);

% Finally, we take the difference between response and undershoot as HRF
Hrf = (g1 - g2)';

% Done
%