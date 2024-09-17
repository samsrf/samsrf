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
%              Alternatively, this can be a boolean (logical). 
%               If true, this incorporates the sign of the pRF (only works if a Srf is provided).
%
%   PlotType:  Whether to plot a contour plot ('C'), a surface plot ('S'), or a 3D scatter plot ('D')
%
%
% 14/03/2022 - Can now compute pRF profile for Srfs stripped of Rmaps (DSS)
% 20/04/2022 - SamSrf 8 version (DSS)
% 08/09/2023 - Now possible to invert sign if needed (DSS)
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
        samsrf_error('If plotting pRF profile directly, first input must be scalar with mapping eccentricity!');
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
    % Use Model structur?
    if isempty(Model) || islogical(Model)
        % Reverse correlation plots
        Xc = SrfEcc.X_coords;
        Yc = SrfEcc.Y_coords;
        % Stripped Srf?
        if isnan(SrfEcc.Rmaps)
            % Compute pRF profile 
            X = SrfEcc.Regs;  % Design matrix (regressors per pixel)
            if isfield(SrfEcc, 'Y_')
                Y = SrfEcc.Y_(:,IdxMat); % Observed time series
            else
                Y = SrfEcc.Y(:,IdxMat); % Observed time series
            end
            % Sign of pRF amplitude?
            sAmp = sign(SrfEcc.Data(5,IdxMat));
            % Regression on raw stimulus design
            warning off
            IdxMat = [Y ones(size(Y,1),1)] \ X; % Linear regression
            warning on
            IdxMat = IdxMat(1,:); % Remove intercept beta     
            if islogical(Model) && Model 
                IdxMat = IdxMat * sAmp;
            end
            IdxMat = reshape(IdxMat, [1 1] * sqrt(length(IdxMat))); % Reshape vector into matrix
            IdxMat = imresize(IdxMat, [size(Yc,1) size(Xc,2)]); % Down-sample r-map
        else
            % Retrieve pRF profile from Srf
            IdxMat = SrfEcc.Rmaps(:,IdxMat); % Data for this vertex
            IdxMat = reshape(IdxMat, size(Yc,1), size(Xc,2)); % Reshape vector into matrix
        end
    else
        % Use model parameters
        if ~isstruct(Model)
            samsrf_error('If plotting pRF from model fit, Model structure must be provided!');
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
% If empty matrix
if nansum(IdxMat(:)) == 0
    IdxMat(:) = 0;
    IdxMat(1,1) = .001; 
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
    samsrf_error('Unknown plotting mode specified.');
end
if isstruct(Model)
    % If plotting model parameters, we must remove deadspace 
    axis([-1 +1 -1 +1] * Model.Scaling_Factor);
end
Scale = max(abs([nanmin(IdxMat(:)) nanmax(IdxMat(:))]));
if Scale == 0
    Scale = 1;
end
set(gca, 'Clim', [-1 +1]*Scale);
colormap berlin
colorbar 
xlabel('Horizontal position (deg)');
ylabel('Vertical position (deg)');

