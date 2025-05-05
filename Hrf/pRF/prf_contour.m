function R = prf_contour(Srf, v, Model)
%
% R = prf_contour(Srf, v, [Model])
%
% Returns the pRF profile for vertex v in Srf as a contour matrix. 
%
% If Model is undefined, the function requires reverse correlation data.
%
% If the optional input Model is defined it uses the model specification 
%  to create the pRF profile using the parameters in Srf.Data(:,v).
%  In this case, the matrix is truncated to remove the padding space that
%  the pRF model functions automatically add.
%
% 14/03/2022 - Can now compute pRF profile for Srfs stripped of Rmaps (DSS)
% 17/03/2022 - Fixed bug a nasty leprechaun added to break reverse-correlation! (DSS)
% 20/04/2022 - SamSrf 8 version (DSS)
%

if nargin > 2
    % Use model parameters
    if ~isstruct(Model)
        samsrf_error('Third input is not a Model structure!');
    end
    P = Srf.Data(2:length(Model.Param_Names)+1, v); % Fit parameters
    P(Model.Scaled_Param==1) = P(Model.Scaled_Param==1) / Model.Scaling_Factor; % Rescale parameters if needed
    Mat = Model.Prf_Function(P,200); % Generate pRF profile from parameters
    ql = size(Mat,1) / 4; % Quarter of the side length
    R = Mat(ql+1:ql*3, ql+1:ql*3); % Truncate padding from matrix
else
    % Use reverse correlation profile
    if ~isfield(Srf, 'Rmaps')
        samsrf_error('No reverse correlation data in this Srf!');
    end
    % Stripped Srf?
    if isnan(Srf.Rmaps)
        % Compute reverse correlation profile
		X = Srf.Regs;  % Design matrix (regressors per pixel)
        if isfield(Srf, 'Y_')
            Y = Srf.Y_(:,v); % Observed time series
        else
            Y = Srf.Y(:,v); % Observed time series
        end
        warning off
        Rmap = [Y ones(size(Y,1),1)] \ X; % Linear regression
        warning on
        Rmap = Rmap(1,:); % Remove intercept beta     	        
    else
        % Reverse correlation profile vector
        Rmap = Srf.Rmaps(:,v);
    end
    % Square side width
    Width = sqrt(length(Rmap));
    % Reshape vector into a map
    R = reshape(Rmap, Width, Width);
end
