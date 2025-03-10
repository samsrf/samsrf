function [fParams, R2] = samsrf_fit2dprf(Rmap, PrfFcn, SeedParams, EccScaPars, ApFrm, FitParam)
%
% [fParams R2] = samsrf_fit2dprf(Rmap, PrfFcn, SeedParams, EccScaPars, ApFrm, FitParam)
%
% Fits a 2D pRF profile to the observed pRF profile. 
%
%   Rmap:           A reverse correlation pRF profile in square format (usually 50x50 pixels)
%                   You can use the prf_contour function to generate that.
%
%   PrfFcn:         Any 2D pRF function included in SamSrf or one you made yourself
%                   Use like in model-based pRF mapping so for example for standard 2D Gaussian:
%
%                       PrfFcn = @(P,ApWidth) prf_gaussian_rf(P(1), P(2), P(3), ApWidth)
%
%   SeedParams:     A row vector to seed the parameter optimisation procedure 
%                   Must have same number of parameters (ignoring betas which are always fit)
%
%                   Important: Must be in visual space & is internally divided 
%                   by your eccentricity or scaling factor!
%
%   EccScaPars:     A vector where the first entry defines the eccentricity/scaling factor
%                   and the following entries are booleans toggling whether a
%                   model parameter is scaled by the eccentricity/scaling factor
%
%   ApFrm:          Aperture frames (only used for masking the pRF profiles)
%
%   FitParam:       Defines the fitting algorithm & parameters for it:
%                     NaN:     Standard Nelder-Mead (fminsearch) fitting 
%                              If it is a vector then FitParam(2) defines parameter tolerance 
%                     Vector:  Initial step sizes per parameter for Hooke-Jeeves algorithm
%                               but the final two values define the number of iterations & shrinks!
%                               (It they are undefined this defaults to 1000 & 10 respectively)
%
% Returns the fitted parameters fParams and the goodness-of-fit R2.
%  The final two parameters are the betas for amplitude and intercept.
%
% If no output arguments are defined, the observed and modelled profiles are plotted.
%
% 15/02/2022 - Added option to fit to pRF coordinate data for CF fitting (DSS)
% 13/04/2022 - Added Hooke-Jeeves algorithm & adjustable Nelder-Mead tolerance (DSS)
% 14/04/2022 - Function defaults to using standard Nelder-Mead algorithm (DSS)
% 20/04/2022 - SamSrf 8 version (DSS)
%

% Use standard fit?
if nargin < 6
    FitParam = NaN;
end

% Check seed parameters are row vector 
if size(SeedParams,2) == 1
    SeedParams = SeedParams';
end
% Check scaling vector matches seed parameter vector
if length(SeedParams) ~= length(EccScaPars)-1
    samsrf_error('Length of seed parameter vector mismatches scaling vector!');
end

% Side length
dims = size(Rmap,1);

% Apertures provided?
if nargin < 5 || isempty(ApFrm)
    % No apertures defined so no masking applied
    Mask = 1;
else
    % Generate mask
    Scaled = NaN(dims, dims, size(ApFrm,3));
    for i = 1:size(ApFrm,3)
        % Rescale apertures to reverse correlation size
        Scaled(:,:,i) = imresize(ApFrm(:,:,i), [dims dims]);
    end
    Mask = sum(Scaled~=0,3)~=0;
end

% Which parameters must be scaled?
for p = 1:length(SeedParams)-2
    if EccScaPars(p+1)
        SeedParams(p) = SeedParams(p) / EccScaPars(1); % Scale this parameter to aperture space
    end
end

% Fit model
if isnan(FitParam(1))
    % Use Nelder-Mead algorithm
    OptimOpts = optimset('Display', 'off');
    if length(FitParam) > 1
        % Define parameter tolerance for Nelder-Mead alogrithm
        OptimOpts.TolX = FitParam(2);
    end
    [fParams, R] = fminsearch(@(P) errfun(PrfFcn,Mask,Rmap,P), [SeedParams 1 0], OptimOpts);
else
    % Use Hooke-Jeeves algorithm
    [fParams, R] = samsrf_hookejeeves(@(P) errfun(PrfFcn,Mask,Rmap,P), [SeedParams 1 0], [FitParam .1 .01], false(1,length(SeedParams)+2));
end
Resid = (Rmap - mean(Rmap(:))).^2; % Squared residuals
R2 = 1 - R/sum(Resid(:)); % Convert to R^2
Fmap = PrfFcn(fParams(1:end-2),dims).*Mask * fParams(end-1) + fParams(end); % Fitted profiles

% Which parameters must be scaled?
for p = 1:length(fParams)-2
    if EccScaPars(p+1)
        fParams(p) = fParams(p) * EccScaPars(1) / 2; % Scale this parameter to visual space
    end
end

% Plot fits?
if nargout == 0
    % Format figure
    set(gcf, 'Name', ['R^2 = ' num2str(R2)]);
    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1)-pos(3)*0.4 pos(2) pos(3)*1.6 pos(4)]);
    
    % Observed pRF profile
    subplot(1,2,1);
    m = nanmax(abs([nanmin(Rmap(:)) nanmax(Rmap(:))]));
    Rmap(1,1) = -m;
    Rmap(end,end) = +m;
    samsrf_showprf(EccScaPars(1), Rmap);
    title('Observed pRF profile');
    
    % Modelled pRF profile
    subplot(1,2,2);
    m = nanmax(abs([nanmin(Fmap(:)) nanmax(Fmap(:))]));
    Fmap(1,1) = -m;
    Fmap(end,end) = +m;
    samsrf_showprf(EccScaPars(1), Fmap);
    title('Modelled pRF profile');
end


%% Error function
function R = errfun(PrfFcn, Mask, ObsRfp, P)

ModelRfp = PrfFcn(P(1:end-2), size(Mask,1)) * P(end-1) + P(end); % Receptive field profile for these parameters
Ds = ObsRfp - ModelRfp.*Mask; % Difference between actual & modelled receptive field profile
R = sum(Ds(:).^2); % Sum of squared residuals