function Bp = samsrf_backproj_revcor(Response, Rmaps, GoF, Threshold, NormaliseByDensity)
%
% Bp = samsrf_backproj_revcor(Response, Rmaps, GoF, [Threshold=[0 0], NormaliseByDensity=true])
%
% Projects the activity values in Response back into visual space using the
% reverse correlation pRF profiles in Rmaps and using the R^2 values in GoF
% for filtering. Effectively, this is a mean of all pRF profiles weighted by 
% their response. 

% Rmaps and GoF are taken from Srf.Data(1,:) and Srf.Rmaps fields with the 
% vertices you want, as produced by samsrf_revcor_prf. You can also regenerate 
% the Rmaps using samsrf_revcor_rmap. Since Rmaps are no longer saved by default 
% this is probably what most people will do.
%
% Response is a matrix where each row is a corresponding vertex in Rmaps &
% GoF & each row is a response map (e.g. a volume from a time course).
%
% WARNING: This function assumes Rmaps have been saved in the Srf. If they 
%          have been stripped it would take way to long to recompute them!
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
% as a matrix of same 2D dimensions as Rmaps with NumberOfVolumes layers. 
% You can use samsrf_showprf to display the outputs.
%
% 14/03/2022 - Ensures an error if Rmaps is NaN, that is, from stripped Srf (DSS)
% 20/04/2022 - SamSrf 8 version (DSS)
% 09/08/2023 - Fixed bug with profiles not being 50x50 pixels (DSS)
%              Updated error message when Rmaps is NaN (DSS)
% 03/09/2023 - NaNs are now set to zero (DSS)
% 16/09/2023 - Clarified that Rmaps can also be regenerated (DSS)
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

% Have pRF profiles been saved?
if isnan(Rmaps)
    samsrf_error('Requires pRF reverse correlation profiles in Rmaps!');
end

% Dimensions of correlation profiles
dim = sqrt(size(Rmaps,1));

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
Bp = NaN(dim,dim,size(Response,1));
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
    Bp(:,:,r) = reshape(CurRmap, dim, dim);
end

% Set NaNs to zero
Bp(isnan(Bp)) = 0;