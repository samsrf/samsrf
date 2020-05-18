function [fParams, R2] = samsrf_fit2dprf(Rmap, PrfFcn, SeedParams, EccScaPars, ApFrm)
%
% [fParams R2] = samsrf_fit2dprf(Rmap, PrfFcn, SeedParams])
%
% Fits a 2D pRF profile to the observed pRF profile. 
%
%   Rmap:       A reverse correlation pRF profile in square format (usually 50x50 pixels)
%
%   PrfFcn:     Any 2D pRF function included in SamSrf or one you made yourself
%               Use like in model-based pRF mapping so for example for standard 2D Gaussian:

%                   PrfFcn = @(P,ApWidth) prf_gaussian_rf(P(1), P(2), P(3), ApWidth)
%
%   SeedParams: A row vector to seed the parameter optimisation procedure 
%               Must have same number of parameters (ignoring betas which are always fit)
%
%               Important: This is in aperture space so you need divide any location 
%                          or width parameters by your eccentricity/scaling factor!
%
%   EccScaPars: A vector where the first entry defines the eccentricity/scaling factor
%               and the following entries are booleans toggling whether a
%               model parameter is scaled by the eccentricity/scaling factor
%
%   ApFrm:      Aperture frames (only used for masking the pRF profiles)
%
% Returns the fitted parameters fParams and the goodness-of-fit R2.
%  The final two parameters are the betas for amplitude and intercept.
%
% If no output arguments are defined, the observed and modelled profiles are plotted.
%
% 15/05/2020 - Written (DSS)

% Check scaling vector matches seed parameter vector
if length(SeedParams) ~= length(EccScaPars)-1
    error('Mismatch betwen seed parameter vector & scaling vector!');
end

% Side length
dims = size(Rmap,1);

% Generate mask
Scaled = NaN(dims, dims, size(ApFrm,3));
for i = 1:size(ApFrm,3)
    % Rescale apertures to reverse correlation size
    Scaled(:,:,i) = imresize(ApFrm(:,:,i), [dims dims]);
end
Mask = sum(Scaled~=0,3)~=0;

% Which parameters must be scaled?
for p = 1:length(SeedParams)-2
    if EccScaPars(p+1)
        SeedParams(p) = SeedParams(p) / EccScaPars(1); % Scale this parameter to aperture space
    end
end

% Fit model
[fParams, R] = fminsearch(@(P) errfun(PrfFcn,Mask,Rmap,P), [SeedParams 1 0]);
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