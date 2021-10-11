function samsrf_showprf(SrfEcc, IdxMat, Model, PlotType)
%
% samsrf_showprf(SrfEcc, IdxMat, Model, [PlotType='C'])
%
% Produces a contour or surface plot of a receptive field profile.
% Use this for displaying fitted pRF models, reverse correlation pRFs, 
%  mean pRFs, or just to visualise pRF profiles generated with the models.
%
%   SrfEcc:    A Srf structure or for directly plotting a pRF profile matrix, 
%               a scalar with the mapping eccentricity.
%
%   IdxMat:    This defines what is to be plotted.
%               Scalar = Reverse correlation profile for vertex index Idx
%               Matrix = pRF profile matrix (e.g. from prf_gaussian_rf or samsrf_meanprf). 
%
%   Model:     Model structure containing the pRF function etc if that's what your plotting.
%               If empty this isn't used.
%
%   PlotType:  Whether to plot a contour plot ('C'), a surface plot ('S'), or a 3D scatter plot ('D')
%
%
% 19/07/2020 - SamSrf 7 version (DSS)
% 06/08/2021 - Increase granularity of colour scheme (DSS)
% 11/10/2021 - Added option to plot pRF models from parameters (DSS)
%              Changed colour scheme to berlin (DSS) 
%

if nargin < 3
    Model = [];
end
if nargin < 4
    PlotType = 'C';
end
PlotType = upper(PlotType);

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
    if isempty(Model)
        % Reverse correlation plots
        Xc = SrfEcc.X_coords;
        Yc = SrfEcc.Y_coords;
        % Reshape vector into matrix
        IdxMat = reshape(SrfEcc.Rmaps(:,IdxMat), size(Yc,1), size(Xc,2));
    else
        % Use model parameters
        if ~isstruct(Model)
            error('If plotting pRF from model fit, Model structure must be provided!');
        end
        P = SrfEcc.Data(2:length(Model.Param_Names)+1, IdxMat); % Fit parameters
        P(Model.Scaled_Param==1) = P(Model.Scaled_Param==1) / Model.Scaling_Factor; % Rescale parameters if needed
        IdxMat = Model.Prf_Function(P,200); % Generate pRF profile from parameters
        % Half width of matrix
        hw = size(IdxMat,1) / 2;
        % Coordinate grid
        [Xc,Yc] = meshgrid([-hw:-1 1:hw]/hw, [-hw:-1 1:hw]/hw);
        Yc = flipud(Yc); % Because Matlab matrix is stupid
        % Convert to visual angle
        Xc = Xc * Model.Scaling_Factor*2;
        Yc = Yc * Model.Scaling_Factor*2; 
    end
end

%% Plot data
if PlotType == 'C'
    contourf(Xc, Yc, IdxMat, 500, 'EdgeColor', 'none');
    axis square
elseif PlotType == 'S'
    surf(Xc, Yc, IdxMat, 'EdgeColor', 'none', 'FaceColor', 'interp');
elseif PlotType == 'D'
    scatter3(Xc(:), Yc(:), IdxMat(:), 100, IdxMat(:), 'filled');
else
    error('Unknown plotting mode specified.');
end
if isstruct(Model)
    % If plotting model parameters, we must remove deadspace 
    axis([-1 +1 -1 +1] * Model.Scaling_Factor);
end
Scale = max(abs([nanmin(IdxMat(:)) nanmax(IdxMat(:))]));
set(gca, 'Clim', [-1 +1]*Scale);
colormap berlin
colorbar 
xlabel('Horizontal position (deg)');
ylabel('Vertical position (deg)');

