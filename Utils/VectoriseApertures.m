function VectoriseApertures(ApsFile, Scaling)
%
% Convers a 2D movie aperture into a vectorised aperture. 
%   ApsFile:    File name of movie apertures
%   Scaling:    Defines the scaling factor(s) for the apertures.
%               If this is a scalar, it assumes square apertures & so movie provided must be square!
%               If it is a 1x2 vector, this scales X and Y differently & movie can be rectangular.
%
% 06/07/2022 - Written (DSS)
%

% Load apertures
load(ApsFile);
if exist('ApXY', 'var') 
    error('I may be wrong but this already looks like a vectorised apertures...');
end

% Vectorise apertures
ApDimX = size(ApFrm,2); % Horizontal size of apertures
ApDimY = size(ApFrm,1); % Vertical size of apertures
% If non-square apertures
if ApDimX ~= ApDimY
    if length(Scaling) == 1
        error('Non-square apertures but only one scaling factor defined!');
    end
end
% If only one scaling factor defined
if length(Scaling) == 1
    Scaling = [1 1] * Scaling; % Define scaling factor for both X and Y
end
ApFrm = reshape(shiftdim(ApFrm, 2), size(ApFrm,3), ApDimX * ApDimY)'; % Pixels in rows, volumes in columns
% Pixel coordinates for apertures
[ApX, ApY] = meshgrid([-ApDimX/2:-1 1:ApDimX/2], fliplr([-ApDimY/2:-1 1:ApDimY/2]));
ApX = ApX(:) / (ApDimX/2) * Scaling(1); % X-coordinates in stimulus space
ApY = ApY(:) / (ApDimY/2) * Scaling(2); % Y-coordinates in stimulus space
ApXY = [ApX ApY];

% Save vectorised apertures
save([ApsFile '_vec'], 'ApFrm', 'ApXY');