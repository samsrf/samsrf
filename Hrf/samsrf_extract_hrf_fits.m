function Hrfs = samsrf_extract_hrf_fits(Srf, TR, vx)
%
% Hrfs = samsrf_extract_hrf_fits(Srf, TR, [vx=[]])
%
% Extracts the vertex-wise HRFs from the pRF map in Srf which has been 
% analysed using concurrent HRF fitting. This is based on the parameter
% estimates of that pRF analysis. The HRF is returned with one column per
% vertex and one row per TR. Because you can define the TR directly, this
% need not be the same as what was used to estimate the HRF.
%
% This analysis is restricted to the vertices defined by vx (defaults to all).
%
% 26/08/2026 - Written (DSS)

% Ensure Srf is expanded
Srf = samsrf_expand_srf(Srf);
if nargin < 3 || isempty(vx)
    vx = 1:size(Srf.Data,2);
end

% Which data rows are HRF parameters?
HrfParams = [find(strcmpi(Srf.Values, 'RLat')) ...
             find(strcmpi(Srf.Values, 'ULat')) ...
             find(strcmpi(Srf.Values, 'RDisp')) ...
             find(strcmpi(Srf.Values, 'UDisp')) ... 
             find(strcmpi(Srf.Values, 'AmpRat') | strcmpi(Srf.Values, 'R/U'))];
if isempty(HrfParams)
    samsrf_error('This data file does not contain HRF fits!');
end

% Generate HRFs
Hrfs = [];
for x = 1:length(vx)
    v = vx(x); % Current vertex
    Hrfs = [Hrfs samsrf_doublegamma(TR, Srf.Data(HrfParams,v))]; % HRF for this vertex
end
