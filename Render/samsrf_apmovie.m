function samsrf_apmovie(ApFile, Framerate)
% samsrf_apmovie(ApFile, [Framerate])
%
% Displays each frame from a SamSrf aperture file, and saves it to disk as
% a movie.
%
% ApFile is a SamSrf aperture file (aps_*.mat)
%
% Framerate allows you to change the number of stimulus frames per second
% in the final movie. Default is 5 frames/second, to match samsrf_tcmovie.
%
% 09/08/2018 SamSrf 6 version (DSS) 
%

%% Defaults
if nargin < 2
    Framerate = 5;
end

%% Main

% Load file
if ischar(ApFile)
    load([ApFile '.mat']);
else
    error('Not a valid ApFile entry')
end

% Ensure apertures are in intensity format
ApFrm = ApFrm - nanmin(ApFrm(:));
ApFrm = ApFrm / nanmax(ApFrm(:));
ApFrm(isnan(ApFrm)) = 0;

% Open video file
Filename = [ApFile '_mov'];
while exist([Filename '.avi'], 'file')
    Filename = [Filename '_1'];
end
vidObj = VideoWriter([Filename '.avi']);
vidObj.FrameRate = Framerate;
open(vidObj);

% Fix 2D aperture (tuning curve)
if ndims(ApFrm) < 3
    
    % 1D -> 2D
    ApFrm = repmat(ApFrm, 1, 1, size(ApFrm, 2));
    ApFrm = permute(ApFrm, [2 3 1]);
end

% Start figure
fh = figure();

% Set view parameters
axis off

% Loop volumes
for i = 1:size(ApFrm, 3)
    
    % Display
    imshow(ApFrm(:,:,i));
    
    % Capture
    CurrFrame = im2frame(repmat(ApFrm(:,:,i), 1, 1, 3));
    
    % Write to file
    writeVideo(vidObj, CurrFrame);
end

% Close video file
close(vidObj);

% Close window
close(fh);

% Done
%