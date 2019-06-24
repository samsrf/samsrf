function samsrf_showprf(Srf, Idx, Mode)
%
% samsrf_showprf(Srf, Idx, [Mode='C'])
%
% Produces a contour or surface plot of a receptive field profile.
% Use this for displaying reverse correlation pRFs, mean pRFs, or just to
%  visualise pRF profiles generated with the prf models.
%
%   Srf:    A Srf structure with reverse correlation data.
%           For plotting a pRF profile matrix, use a scalar with the mapping eccentricity.
%
%   Idx:    This defines what is to be plotted.
%               Scalar = Reverse correlation profile for vertex index Idx
%               Matrix = pRF profile matrix (e.g. from prf_gaussian_rf or samsrf_meanprf)
%
%   Mode:   Whether to plot a contour plot ('C') or a surface plot ('S').
%
% 09/08/2018 - SamSrf 6 version (DSS)
%

if nargin < 3
    Mode = 'C';
end
Mode = upper(Mode);

if ~isstruct(Srf) 
    % pRF profile matrix
    Mat = Idx; 
    % Half width of matrix
    hw = size(Idx,1) / 2;
    % Coordinate grid
    [Xc,Yc] = meshgrid([-hw:-1 1:hw]/hw, [-hw:-1 1:hw]/hw);
    Yc = flipud(Yc); % Because Matlab matrix is stupid
    % Convert to visual angle
    Xc = Xc * Srf;
    Yc = Yc * Srf; 
else
    % Expand Srf if necessary
    Srf = samsrf_expand_srf(Srf);
    % Reverse correlation plots
    Xc = Srf.X_coords;
    Yc = Srf.Y_coords;
    % Reshape vector into matrix
    Mat = reshape(Srf.Rmaps(:,Idx), size(Yc,1), size(Xc,2));
end

%% Plot data
if Mode == 'C'
    contourf(Xc, Yc, Mat);
    axis square
elseif Mode == 'S'
    surf(Xc, Yc, Mat, 'EdgeColor', 'none', 'FaceColor', 'interp');
else
    error('Unknown plotting mode specified.');
end
colormap hotcold
colorbar 

