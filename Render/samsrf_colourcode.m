function img = samsrf_colourcode(vartype, width, colmap)
%%
% img = samsrf_colourcode(vartype, [width, colmap])
%
% Displays the pseudo-colour codes for samsrf_surf. 
% Note: this will reflect if you redefined the colour maps in SamSrf_defaults.mat. 
%   
%   vartype :     'Polar_Lh', 'Polar_Rh', 'Eccen', 'Sigma', or 'R^2'
%   width :       Optional, width of the image
%   colmap :      Optional, instead of SamSrf_defaults uses this colour map.
%
% The output of this function is a bitmap which can be displayed or saved.
% 
% 15/09/2018 - Fixed errors with eccentricity & made polar consistent with 
%               samsrf_surf but numerically identical to previous version (DSS)
% 27/11/2018 - Modified to make more versatile for some colour schemes and 
%               added option to use colour maps without changing SamSrf_defaults (DSS)
%

if nargin < 2
    width = 100;    % width of image
end

% load defaults
if nargin < 3
    load('SamSrf_defaults.mat');
else
    def_cmap_angle = colmap;
    def_cmap_eccen = colmap;
    def_cmap_sigma = colmap;
    def_cmap_other = colmap;
end

% image pixels
[x y] = meshgrid(-width:width, -width:width);
[t r] = cart2pol(x,y);
t = t / pi * 180;   % convert to degrees
nr = ceil(r/width*360); % normalized rho in degrees
nr(nr==0) = 1;

% image for output
imgR = zeros(2*width+1, 2*width+1);
imgG = zeros(2*width+1, 2*width+1);
imgB = zeros(2*width+1, 2*width+1);

if strcmpi(vartype, 'Polar_Lh')
    t = t + 270;  % calibrate angles for 0 degrees to be at 3 o'clock
    t = mod(ceil(t),360) + 1;   % ensure between 1-360

    % Colourmap
    cstr = ['colormap(' def_cmap_angle '(360));'];
    Cmap = eval(cstr);        
    
    imgR(r<width) = Cmap(t(r<width),1);
    imgG(r<width) = Cmap(t(r<width),2);
    imgB(r<width) = Cmap(t(r<width),3);
elseif strcmpi(vartype, 'Polar_Rh')
    t = -t + 90;  % calibrate angles for 0 degrees to be at 9 o'clock
    t = mod(ceil(t),360) + 1;   % ensure between 1-360

    % Colourmap
    cstr = ['colormap(' def_cmap_angle '(360));'];
    Cmap = eval(cstr);        
    
    imgR(r<width) = Cmap(t(r<width),1);
    imgG(r<width) = Cmap(t(r<width),2);
    imgB(r<width) = Cmap(t(r<width),3);
elseif strcmpi(vartype, 'Eccen')
    % Colourmap
    cstr = ['colormap(' def_cmap_eccen '(360));'];
    Cmap = eval(cstr);        

    imgR(r<width) = Cmap(nr(r<width),1);
    imgG(r<width) = Cmap(nr(r<width),2);
    imgB(r<width) = Cmap(nr(r<width),3);
elseif strcmpi(vartype, 'Sigma')
    % Colourmap
    cstr = ['colormap(' def_cmap_sigma '(401));'];
    Cmap = eval(cstr);        
    
    for c = 1:size(Cmap,1)
        imgR(c,abs(x(1,:))<50) = Cmap(c,1);
        imgG(c,abs(x(1,:))<50) = Cmap(c,2);
        imgB(c,abs(x(1,:))<50) = Cmap(c,3);
    end
elseif strcmpi(vartype, 'R^2')
    % Colourmap
    cstr = ['colormap(' def_cmap_other '(801));'];
    Cmap = eval(cstr);        
    Cmap = Cmap(400:end,:);
    
    for c = 1:size(Cmap,1)
        imgR(c,abs(x(1,:))<50) = Cmap(c,1);
        imgG(c,abs(x(1,:))<50) = Cmap(c,2);
        imgB(c,abs(x(1,:))<50) = Cmap(c,3);
    end
else
    % Colourmap
    cstr = ['colormap(' def_cmap_other '(401));'];
    Cmap = eval(cstr);        

    for c = 1:size(Cmap,1)
        imgR(c,abs(x(1,:))<50) = Cmap(c,1);
        imgG(c,abs(x(1,:))<50) = Cmap(c,2);
        imgB(c,abs(x(1,:))<50) = Cmap(c,3);
    end
end

img(:,:,1) = flipud(imgR);
img(:,:,2) = flipud(imgG);
img(:,:,3) = flipud(imgB);

if ~strcmpi(vartype, 'Polar_Lh') && ~strcmpi(vartype, 'Polar_Rh') && ~strcmpi(vartype, 'Eccen')
    img = img(:,abs(x(1,:))<width/4,:);
end
imshow(img);
