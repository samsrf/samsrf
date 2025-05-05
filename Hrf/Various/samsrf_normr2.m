function Srf = samsrf_normr2(Srf, NcThr)
%
% Srf = samsrf_normr2(Srf, [NcThr=0])
% 
% Calculates the normalised goodness-of-fit nR^2. It searches through 
% Srf.Values to see if there is a Noise Ceiling field. If so, it then 
% divides the R^2 field by the noise ceiling and replaces the R^2 field 
% (i.e. Srf.Data(1,:)) with the nR^2 field. While redundant, the raw R^2 
% is moved to the bottom row of Srf.Data.
%
% The function filters out data points with bad noise ceilings by setting
% the normalised R^2 to zero. The default threshold is a noise ceiling of 0 
% but this can be tweaked with the optional NcThr input.
%
% Note, this function only works on Srf.Data so don't apply that to
% smoothed data. If smoothing is desired, do that afterwards.
%
% 20/04/2022 - SamSrf 8 version (DSS)
%

if nargin < 2
    NcThr = 0; % Default noise ceiling threshold
end

if ~isfield(Srf, 'Values')
    samsrf_error('No value labels in Srf - probably no mapping data?');
end
if ~strcmpi(Srf.Values{1}, 'R^2')
    samsrf_error('1st row of Srf.Data is not R^2!');
end

% Data fields
R2 = Srf.Data(1,:); % Raw goodness-of-fit
Nc = Srf.Data(strcmpi(Srf.Values, 'Noise Ceiling'),:); % Noise ceiling

if isempty(Nc)
    % No noise ceiling available
    samsrf_disp('No noise ceiling has been computed');
else
    % Threshold by noise ceiling
    R2(Nc <= NcThr) = 0;
    % Normalised goodness-of-fit
    nR2 = R2 ./ Nc;

    % Rearrange Srf
    Srf.Values{1} = 'nR^2';
    Srf.Values{end+1} = 'R^2';
    Srf.Data = [Srf.Data; R2];
    Srf.Data(1,:) = nR2;
    
    samsrf_disp('Computed normalised R^2');
end