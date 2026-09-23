function R = prf_contour(Srf, v, Model)
%
% R = prf_contour(Srf, v, [Model=NaN])
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
% Instead, Model can also simply be the vector or matrix of your HRF. 
%  This is to convolve the regressors with the HRF as is needed for
%  micro-probing or using vertex-wise HRFs.
%
% 14/03/2022 - Can now compute pRF profile for Srfs stripped of Rmaps (DSS)
% 17/03/2022 - Fixed bug a nasty leprechaun added to break reverse-correlation! (DSS)
% 20/04/2022 - SamSrf 8 version (DSS)
% 04/09/2026 - Now supports micro-probing data (DSS)
%              Also supports vertex-wise HRF data now (DSS)
%

if nargin < 3
    Model = NaN;
end
    
if isstruct(Model)
    % Use pRF model parameters
    P = Srf.Data(2:length(Model.Param_Names)+1, v); % Fit parameters
    P(Model.Scaled_Param==1) = P(Model.Scaled_Param==1) / Model.Scaling_Factor; % Rescale parameters if needed
    Mat = Model.Prf_Function(P,200); % Generate pRF profile from parameters
    ql = size(Mat,1) / 4; % Quarter of the side length
    R = Mat(ql+1:ql*3, ql+1:ql*3); % Truncate padding from matrix

else
    % Use reverse-correlation profile
    if ~isfield(Srf, 'Rmaps')
        samsrf_error('No reverse-correlation data in this Srf!');
    end

    % Stripped Srf?
    if isnan(Srf.Rmaps)
        % Compute reverse correlation profile
		X = Srf.Regs;  % Design matrix (regressors per pixel)
        % Convolve with HRF?
        if ~isnan(Model(1))
            % Vertex-wise HRF?
            if size(Model,2) > 1
                % Convolve with vertex HRF
                for i = 1:size(X,2) 
                    X(:,i) = prf_convolve_hrf(X(:,i), Model(:,v), 1); % Note that dowwnsampling has not been implemented!
                end
            else
                % Convolve with HRF
                for i = 1:size(X,2) 
                    X(:,i) = prf_convolve_hrf(X(:,i), Model, 1); % Note that dowwnsampling has not been implemented!
                end
            end
        end
        % Now generate profile
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
