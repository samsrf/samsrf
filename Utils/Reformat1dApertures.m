function ApFrm = Reformat1dApertures(InAps, ApFormat)

% Converts the 1D apertures in InAps into a 3D aperture in ApFrm.
% The optional input ApFormat determines the format of the aperture:
%   
%   '-':    (Default) Apertures are defined as columns moving horizontally 
%   '|':    Apertures are defined as rows moving vertically
%   '*':    Apertures are defined as wedges rotating around the centre 
%             For this InAps must be defined as n x 360 matrix
%
% Note: there is no point to make apertures for expanding/contracting rings
%  because it is simpler to just define those as bars moving horizontally.
%

if nargin < 2
    ApFormat = '-';
end

% Number of volumes
Nvol = size(InAps,1);
% Width of input apertures
Width = size(InAps,2);
if ApFormat == '*'
    Width = 100;
end

% Reformat
ApFrm = NaN(Width, Width, Nvol);
switch ApFormat
    case '-'
        % Turn into vertical bar 
        for v = 1:Nvol
            ApFrm(:,:,v) = repmat(InAps(v,:), Width, 1); % Row vector to be expanded vertically
        end
    case '|'
        % Turn into horizontal bar
        for v = 1:Nvol
            ApFrm(:,:,v) = repmat(InAps(v,:)', 1, Width); % Transpose to column vector & expand horizontally
        end
    case '*'
        % Turn into polar wedge
        [X,Y] = meshgrid([-Width/2:-1 1:Width/2], fliplr([-Width/2:-1 1:Width/2])); % Cartesian coordinates
        Theta = round(atan2(Y,X) / pi * 180); % Polar angles 
        Rho = sqrt(X.^2 + Y.^2); % Eccentricities
        for v = 1:Nvol
            Degs = find(InAps(v,:)); % Current polar angles 0-359
            Degs(Degs > 180) = Degs(Degs > 180) - 360; % Current polar angles -180 - +180
            ApFrm(:,:,v) = double(ismember(Theta, Degs) & Rho <= Width/2); % Stimulated polar angles within eccentricity limit
        end
    otherwise
        error('Invalid aperture format specified!');
end