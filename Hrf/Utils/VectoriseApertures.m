function VectoriseApertures(ApsFile, Dummy)
%
% Convers a 2D movie aperture into a vectorised aperture. 
%   ApsFile:    File name of movie apertures
%
% IMPORTANT: The side lengths of the apertures must be even numbers! 
%            Unlike in versions prior to 9.81, the scale is done automatically by samsrf_fit_prf.
%			 In the ApXY variable in the aperture file, maximum side length is defined as -1/+1.
%
% 06/07/2022 - Written (DSS)
% 15/11/2022 - Now gives error if aperture length is odd-numbered (DSS)
% 24/01/2023 - More info about how to use the scaling (DSS)
% 02/12/2023 - Removed scaling input as apertures are now scaled automatically (DSS)  
% 03/08/2024 - Added support for multi-condition apertures (DSS)
%

if nargin > 1
    samsrf_error('Since version 9.81, this function no longer has scaling input!');
end

% Load apertures
load(ApsFile);
if exist('ApXY', 'var') 
    samsrf_error('I may be wrong but this already looks like a vectorised apertures...');
end

% Vectorise apertures
ApDimX = size(ApFrm,2); % Horizontal size of apertures
ApDimY = size(ApFrm,1); % Vertical size of apertures
% Odd-numbered side lengths?
if mod(ApDimX,2)
    samsrf_error('X dimension is odd number!');
end
if mod(ApDimY,2)
    samsrf_error('Y dimension is odd number!');
end
% If non-square apertures
if ApDimX ~= ApDimY	
    if ApDimX > ApDimY
		Scaling = [1 ApDimY/ApDimX];
    else
		Scaling = [ApDimX/ApDimY 1];
    end
	samsrf_disp(['Rectangular apertures: ' num2str(Scaling(1)) ' x ' num2str(Scaling(2))]);
else
	samsrf_disp('Square apertures');
	Scaling = [1 1]; % Square aperture
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

% Is this a multi-condition design?
if exist('ApCond', 'var')
    samsrf_disp(' Multi-condition design: condition labels are saved in 1st row of ApFrm.');
    ApFrm = [ApCond; ApFrm]; % Condition labels are in first row of apertures!
end

% Save vectorised apertures
save([ApsFile '_vec'], 'ApFrm', 'ApXY');