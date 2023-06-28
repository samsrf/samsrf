function [Srf, vx] = samsrf_expand_srf(Srf, Use64bit)
%
% [Srf, vx] = samsrf_expand_srf(Srf, Use64bit)
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
% The optional second input toggles whether to force the function to use 
%  64 bit (double) data types, even if the default setting is to use 32 bit.
%  You can use this if you want to read in old data saved in 64 bit, but the
%  data will be numerically identcal. However, the data type could have some
%  downstream consequences on any analyses you are using.
%
% Due to various changes in previous SamSrf versions, this function is not
%  meant to be used for files created in versions prior to SamSrf 5. 
%
% The second output argument contains the list of vertices in the ROI.
%
% Note: this function only works for single-subject Srfs. There currently
%   is no support for expanding or compressing multi-subbject Srfs.
%
% 14/03/2022 - pRF reverse correlation profiles are not saved in data file by default (DSS)
% 20/04/2022 - SamSrf 8 version (DSS)
% 15/05/2022 - Now works with volumetric data (DSS)
% 16/10/2022 - Fixed possible backwards Matlab compatibility issue (DSS)
% 29/06/2023 - Data now converted into 32 bit unless fixed in SamSrf_defaults (DSS)
%              Also addded option to force using 64 bit data though (DSS)
%

%% In case no values defined
if ~isfield(Srf, 'Values')
    Srf.Values = {};
    for r = 1:size(Srf.Data,1)
        Srf.Values{r,1} = ['Volume #' num2str(r)];
    end
end

%% Do nothing if volumetric data! 
if ~strcmpi(Srf.Hemisphere, 'vol')
    
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
        Srf.Structural = Anat.Structural;
        Srf.Hemisphere = Anat.Hemisphere;
        Srf.Vertices = Anat.Vertices;
        Srf.Faces = Anat.Faces;
        if isfield(Anat, 'Normals')
            % Don't exist for anonymised Srfs
            Srf.Normals = Anat.Normals;
            Srf.Pial = Anat.Pial;
        end
        Srf.Inflated = Anat.Inflated;
        Srf.Sphere = Anat.Sphere;
        Srf.Curvature = Anat.Curvature;
        if isfield(Anat, 'Area')
            % Don't exist for anonymised Srfs
            Srf.Area = Anat.Area;
            Srf.Thickness = Anat.Thickness;
        end
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
            if ~isnan(Srf.Rmaps)
                Rm = Srf.Rmaps;
                Srf.Rmaps = zeros(size(Rm,1),nv);
                Srf.Rmaps(:,vx) = Rm;
            end
        end
        if isfield(Srf, 'ConFlds')
            if ~isnan(Srf.ConFlds)
                Cf = Srf.ConFlds;
                Srf.ConFlds = zeros(size(Cf,1),nv);
                Srf.ConFlds(:,vx) = Cf;
            end
        end
        if isfield(Srf, 'Noise_Ceiling')
            Nc = Srf.Noise_Ceiling;
            Srf.Noise_Ceiling = zeros(size(Nc,1),nv);
            Srf.Noise_Ceiling(:,vx) = Nc;
        end

        % Remove ROI field
        Srf = rmfield(Srf, 'Roi');
    end

else
    % For volumetric data return empty vertex field
    vx = [];
end

% Convert to 32 bit?
if nargin < 2 || ~Use64bit
    Srf = samsrf_32bit_srf(Srf);
else
    disp('Forcing any 64 bit data to remain!');
end