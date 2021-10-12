function R = prf_contour(Srf, v, Model)
%
% R = prf_contour(Srf, v, [Model])
%
% Returns the pRF profile for vertex v in Srf as a contour matrix. 
% The matrix is flipped vertically to ensure it is correct in visual space.
% (See also samsrf_showprf for a very similar function for plotting)
%
% If Model is undefined, the function requires reverse correlation data.
%
% If the optional input Model is defined it uses the model specification 
%  to create the pRF profile using the parameters in Srf.Data(:,v).
%  In this case, the matrix is truncated to remove the padding space that
%  the pRF model functions automatically add.
%
% 19/07/2020 - SamSrf 7 version (DSS)
% 12/10/2021 - New option for plotting pRFs from model parameters.
%

if nargin > 2
    % Use model parameters
    if ~isstruct(Model)
        error('Third input is not a Model structure!');
    end
    P = Srf.Data(2:length(Model.Param_Names)+1, v); % Fit parameters
    P(Model.Scaled_Param==1) = P(Model.Scaled_Param==1) / Model.Scaling_Factor; % Rescale parameters if needed
    Mat = Model.Prf_Function(P,200); % Generate pRF profile from parameters
    ql = size(Mat,1) / 4; % Quarter of the side length
    R = Mat(ql+1:ql*3, ql+1:ql*3); % Truncate padding from matrix
else
    % Use reverse correlation profile
    if ~isfield(Srf, 'Rmaps')
        error('No reverse correlation data in this Srf!');
    end
    % Reverse correlation profile vector
    Rmap = Srf.Rmaps(:,v);
    % Square side width
    Width = sqrt(size(Rmap,1));
    % Reshape vector into a map
    R = reshape(Rmap, Width, Width);
end

% Flip vertically
R = flipud(R);