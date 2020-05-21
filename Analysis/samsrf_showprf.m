function samsrf_showprf(SrfEcc, IdxMat, Mode)
%
% samsrf_showprf(SrfEcc, IdxMat, [Mode='C'])
%
% Produces a contour or surface plot of a receptive field profile.
% Use this for displaying reverse correlation pRFs, mean pRFs, or just to
%  visualise pRF profiles generated with the prf models.
%
%   SrfEcc:    A Srf structure with reverse correlation data or
%           for directly plotting a pRF profile matrix, a scalar with the mapping eccentricity.
%
%   IdxMat:    This defines what is to be plotted.
%               Scalar = Reverse correlation profile for vertex index Idx
%               Matrix = pRF profile matrix (e.g. from prf_gaussian_rf or samsrf_meanprf)
%
%   Mode:   Whether to plot a contour plot ('C'), a surface plot ('S'), or a 3D scatter plot ('D')
%
% 09/08/2018 - SamSrf 6 version (DSS)
% 13/05/2020 - Now gives clearer error when using pRF matrix (DSS)
% 18/05/2020 - Renamed input variables & edited help section (DSS)
%              Changes how colour scheme is scaled (DSS)
%              Added support for 3D scatter plot (DSS)
%

if nargin < 3
    Mode = 'C';
end
Mode = upper(Mode);

if ~isstruct(SrfEcc) 
    if ~isscalar(SrfEcc)
        error('If plotting pRF profile directly, first input must be scalar with mapping eccentricity!');
    end
    % Half width of matrix
    hw = size(IdxMat,1) / 2;
    % Coordinate grid
    [Xc,Yc] = meshgrid([-hw:-1 1:hw]/hw, [-hw:-1 1:hw]/hw);
    Yc = flipud(Yc); % Because Matlab matrix is stupid
    % Convert to visual angle
    Xc = Xc * SrfEcc;
    Yc = Yc * SrfEcc; 
else
    % Expand Srf if necessary
    SrfEcc = samsrf_expand_srf(SrfEcc);
    % Reverse correlation plots
    Xc = SrfEcc.X_coords;
    Yc = SrfEcc.Y_coords;
    % Reshape vector into matrix
    IdxMat = reshape(SrfEcc.Rmaps(:,IdxMat), size(Yc,1), size(Xc,2));
end

%% Plot data
if Mode == 'C'
    contourf(Xc, Yc, IdxMat, 100, 'EdgeColor', 'none');
    axis square
elseif Mode == 'S'
    surf(Xc, Yc, IdxMat, 'EdgeColor', 'none', 'FaceColor', 'interp');
elseif Mode == 'D'
    scatter3(Xc(:), Yc(:), IdxMat(:), 100, IdxMat(:), 'filled');
else
    error('Unknown plotting mode specified.');
end
Scale = max(abs([nanmin(IdxMat(:)) nanmax(IdxMat(:))]));
set(gca, 'Clim', [-1 +1]*Scale);
colormap hotcold
colorbar 
xlabel('Horizontal position (deg)');
ylabel('Vertical position (deg)');

