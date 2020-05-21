function Bp = samsrf_backproj_revcor(Response, Rmaps, GoF, Threshold, NormaliseByDensity)
%
% Bp = samsrf_backproj_revcor(Response, Rmaps, GoF, [Threshold=[0 0], NormaliseByDensity=true])
%
% Projects the activity values in Response back into visual space using the
% reverse correlation pRF profiles in Rmaps and using the R^2 values in GoF
% for filtering. Effectively, this is a mean of all pRF profiles weighted by 
% their response. Rmaps and GoF are taken from Srf.Data(1,:) and Srf.Rmaps fields 
% with the vertices you want, as produced by samsrf_revcor_prf. Response is
% a matrix where each row is a corresponding vertex in Rmaps and GoF and
% each row is a response map (e.g. a volume from a time course).
%
% The optional Threshold is a 1x2 vector defining the minimal R^2 of the 
% pRFs to be projected and the proportion of the peak of the profile to be 
% included. Both default to 0 if undefined. 
%
% The optional NormaliseByDensity toggles density normalisation on or off: 
% In this normalisation, the final back-projection image is divided by a 
% back-projection of the pRF density (all pRFs summed with equal weights). 
%
% Returns an intensity movie of the time series plotted back into visual space 
% as a matrix of 50 x 50 x NumberOfVolumes. You can use samsrf_showprf to
% display the outputs.
%
% 20/05/2020 - Written (DSS)
%

if nargin < 4
    Threshold = [0 0];
end
if length(Threshold) == 1
    Threshold = [Threshold 0];
end
if nargin < 5
    NormaliseByDensity = true;
end

% Threshold based on R^2
gof = GoF(1,:) > Threshold(1);
Rmaps = Rmaps(:,gof);
Response = Response(:,gof);

% Threshold reverse correlation profiles?
for v = 1:size(Rmaps,2)
    Rmaps(Rmaps(:,v) < nanmax(Rmaps(:,v)) * Threshold(2), v) = nanmax(Rmaps(:,v)) * Threshold(2); % Set baseline  
    Rmaps(:,v) = Rmaps(:,v) - nanmin(Rmaps(:,v)); % Set baseline to zero
end

% Density map if needed
if NormaliseByDensity
    Dens = NaN(size(Rmaps));
    for v = 1:size(Dens,2)
       Dens(:,v) = Rmaps(:,v) / nanmax(Rmaps(:,v)); % Normalised activity profile
    end
    % Mean density maps
    Dens = nanmean(Dens,2);
end

% Mean reverse correlation maps
Bp = NaN(50,50,size(Response,1));
% Loop thru response rows
for r = 1:size(Response,1)
    % Multiply profiles by responses
    CurRmap = Rmaps .* repmat(Response(r,:), size(Rmaps,1), 1);
    % Mean response profile
    CurRmap = nanmean(CurRmap,2);
    % Normalise by density?
    if NormaliseByDensity
        CurRmap = CurRmap ./ Dens;
    end
    % Reshape to square image
    Bp(:,:,r) = reshape(CurRmap, 50, 50);
end
