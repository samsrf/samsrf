function [Srf, vx] = samsrf_expand_srf(Srf)
%
% [Srf, vx] = samsrf_expand_srf(Srf)
%
% If the surface data Srf has not been compressed this function does nothing. 
% However, if only the data from a particular region of interest have been 
%  saved, this creates a full size surface structure and fills the appropriate 
%  vertices with the data. Such compressed data files are identified because 
%  the structure has a field called Srf.Roi that contains the vertices in the ROI.
%
% Moreover, if the anatomical surface meshes have been stored separately,
%  this function will automatically load them in. Such Srf structures are
%  identified by containing the field Srf.Meshes. This only works with Srfs
%  created in SamSrf 6 and higher.
%
% Due to various changes in previous SamSrf versions, this function is not
%  meant to be used for files created in versions prior to SamSrf 5. 
%
% The second output argument contains the list of vertices in the ROI.
%
% As of SamSrf 6, all functions that take Srfs as an input will check if
%  the Srf has been expanded and produce an error if it hasn't.
%
% 10/08/2018 - SamSrf 6 version (DSS)
% 28/11/2018 - Added support for noise ceiling (DSS)
% 30/11/2018 - Fixed show-stopping bug with noise ceiling (DSS)
% 01/12/2018 - Fixed platform-dependency with path when loading anatomy (DSS)
%

%% In case no values defined
if ~isfield(Srf, 'Values')
    Srf.Values = {};
    for r = 1:size(Srf.Data,1)
        Srf.Values{r,1} = ['Volume #' num2str(r)];
    end
end

%% Deal with anatomical meshes
if isfield(Srf, 'Meshes') && ~isfield(Srf, 'Vertices')
    disp('Loading anatomical meshes...');
    % Ensure path uses correct file separation
    if filesep == '/'
        sl = strfind(Srf.Meshes, '\');
    elseif filesep == '\'
        sl = strfind(Srf.Meshes, '/');
    else
        error('Unknown path separator in this OS!');
    end
    Srf.Meshes(sl) = filesep;
    
    % Load anatomical meshes
    load(Srf.Meshes);
    Srf.Version = Anat.Version;
    Srf.Structural = Anat.Structural;
    Srf.Hemisphere = Anat.Hemisphere;
    Srf.Vertices = Anat.Vertices;
    Srf.Faces = Anat.Faces;
    Srf.Normals = Anat.Normals;
    Srf.Pial = Anat.Pial;
    Srf.Inflated = Anat.Inflated;
    Srf.Sphere = Anat.Sphere;
    Srf.Curvature = Anat.Curvature;
    Srf.Area = Anat.Area;
    Srf.Thickness = Anat.Thickness;
end

%% Deal with region of interest 
vx = [];% Initialize output in case there is no ROI
if isfield(Srf, 'Roi')
    disp('Expanding surface data file...');
    
    % ROI vertices
    vx = Srf.Roi;
    nv = size(Srf.Vertices,1);
    
    % Create temporary data, expand, and store vertices
    D = Srf.Data;
    Srf.Data = zeros(size(D,1),nv);
    Srf.Data(:,vx) = D;
    if isfield(Srf, 'Y')
        Y = Srf.Y;
        Srf.Y = zeros(size(Y,1),nv);
        Srf.Y(:,vx) = Y;
    end
    if isfield(Srf, 'X')
        X = Srf.X;
        Srf.X = zeros(size(X,1),nv);
        Srf.X(:,vx) = X;
    end
    if isfield(Srf, 'Raw_Data')
        Rd = Srf.Raw_Data;
        Srf.Raw_Data = zeros(size(Rd,1),nv);
        Srf.Raw_Data(:,vx) = Rd;
    end
    if isfield(Srf, 'Rmaps')
        Rm = Srf.Rmaps;
        Srf.Rmaps = zeros(size(Rm,1),nv);
        Srf.Rmaps(:,vx) = Rm;
    end
    if isfield(Srf, 'ConFlds')
        Cf = Srf.ConFlds;
        Srf.ConFlds = zeros(size(Cf,1),nv);
        Srf.ConFlds(:,vx) = Cf;
    end
    if isfield(Srf, 'Noise_Ceiling')
        Nc = Srf.Noise_Ceiling;
        Srf.Noise_Ceiling = zeros(size(Nc,1),nv);
        Srf.Noise_Ceiling(:,vx) = Nc;
    end
    
    % Remove ROI field
    Srf = rmfield(Srf, 'Roi');
end

