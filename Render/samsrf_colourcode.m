function img = samsrf_colourcode(vartype, scale, colmap)
%%
% img = samsrf_colourcode(vartype, [scale, colmap])
%
% Displays the pseudo-colour codes for samsrf_surf. 
% Note: this will reflect if you redefined the colour maps in SamSrf_defaults.json. 
%   
%   vartype :     'Polar_Lh', 'Polar_Rh', 'Eccen', 'Sigma', 'R^2', or 'Polar_Eeg'
%   scale :       Optional, determines the size of the image (defaults = 100%)
%   colmap :      Optional, instead of SamSrf_defaults uses this colour map.
%
% The output of this function is a bitmap which can be displayed or saved.
% 
% 31/10/2023 - Added option for M/EEG colour wheel (DSS)
%

if nargin < 2
    scale = 100;    % width of image
end

% load defaults
if nargin < 3
    samsrf_disp(['Using defaults in: ' which('SamSrf_defaults.json')]);
    SamSrfDefs = LoadSamSrfDefaults;
else
    SamSrfDefs.defs_cmap_angle = colmap;
    SamSrfDefs.defs_cmap_eccen = colmap;
    SamSrfDefs.defs_cmap_sigma = colmap;
    SamSrfDefs.defs_cmap_other = colmap;
end

% image pixels
[x y] = meshgrid(-scale:scale, -scale:scale);
[t r] = cart2pol(x,y);
t = t / pi * 180;   % convert to degrees
nr = ceil(r/scale*360); % normalized rho in degrees
nr(nr==0) = 1;

% image for output
imgR = zeros(2*scale+1, 2*scale+1);
imgG = zeros(2*scale+1, 2*scale+1);
imgB = zeros(2*scale+1, 2*scale+1);

if strcmpi(vartype, 'Polar_Eeg')
    t = t + 180;
    t = mod(ceil(t),360) + 1;   % ensure between 1-360

    % Colourmap
    Cmap = samsrf_cmap(SamSrfDefs.defs_cmap_angle, 360);    
    imgR(r<scale) = Cmap(t(r<scale),1);
    imgG(r<scale) = Cmap(t(r<scale),2);
    imgB(r<scale) = Cmap(t(r<scale),3);
elseif strcmpi(vartype, 'Polar_Lh')
    t = t + 270;  % calibrate angles for 0 degrees to be at 3 o'clock
    t = mod(ceil(t),360) + 1;   % ensure between 1-360

    % Colourmap
    Cmap = samsrf_cmap(SamSrfDefs.defs_cmap_angle, 360);        
    imgR(r<scale) = Cmap(t(r<scale),1);
    imgG(r<scale) = Cmap(t(r<scale),2);
    imgB(r<scale) = Cmap(t(r<scale),3);
elseif strcmpi(vartype, 'Polar_Rh')
    t = -t + 90;  % calibrate angles for 0 degrees to be at 9 o'clock
    t = mod(ceil(t),360) + 1;   % ensure between 1-360

    % Colourmap
    Cmap = samsrf_cmap(SamSrfDefs.defs_cmap_angle, 360);        
    imgR(r<scale) = Cmap(t(r<scale),1);
    imgG(r<scale) = Cmap(t(r<scale),2);
    imgB(r<scale) = Cmap(t(r<scale),3);
elseif strcmpi(vartype, 'Eccen')
    % Colourmap
    Cmap = samsrf_cmap(SamSrfDefs.defs_cmap_eccen, 360);        
    imgR(r<scale) = Cmap(nr(r<scale),1);
    imgG(r<scale) = Cmap(nr(r<scale),2);
    imgB(r<scale) = Cmap(nr(r<scale),3);
elseif strcmpi(vartype, 'Sigma')
    % Colourmap
    Cmap = samsrf_cmap(SamSrfDefs.defs_cmap_sigma, 400);        
    for c = 1:size(Cmap,1)
        imgR(c,abs(x(1,:))<50) = Cmap(c,1);
        imgG(c,abs(x(1,:))<50) = Cmap(c,2);
        imgB(c,abs(x(1,:))<50) = Cmap(c,3);
    end
elseif strcmpi(vartype, 'R^2') || strcmpi(vartype, 'nR^2')
    % Colourmap
    Cmap = samsrf_cmap(SamSrfDefs.defs_cmap_other, 800);        
    for c = 1:size(Cmap,1)
        imgR(c,abs(x(1,:))<50) = Cmap(c,1);
        imgG(c,abs(x(1,:))<50) = Cmap(c,2);
        imgB(c,abs(x(1,:))<50) = Cmap(c,3);
    end
else
    % Colourmap
    Cmap = samsrf_cmap(SamSrfDefs.defs_cmap_other, 400);        
    for c = 1:size(Cmap,1)
        imgR(c,abs(x(1,:))<50) = Cmap(c,1);
        imgG(c,abs(x(1,:))<50) = Cmap(c,2);
        imgB(c,abs(x(1,:))<50) = Cmap(c,3);
    end
end

img(:,:,1) = flipud(imgR);
img(:,:,2) = flipud(imgG);
img(:,:,3) = flipud(imgB);

if ~strcmpi(vartype, 'Polar_Lh') && ~strcmpi(vartype, 'Polar_Rh') && ~strcmpi(vartype, 'Eccen') && ~strcmpi(vartype, 'Polar_Eeg')
    img = img(:,abs(x(1,:))<scale/4,:);
end
imshow(img);
