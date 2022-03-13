function samsrf_convert_old_srf(srfname, surfdir)
% 
% samsrf_convert_old_srf(srfname, surfdir)
%
% Converts an Srf file from an older version of SamSrf into one that is
% compatible with SamSrf version 6 by loading the structural surface data
% and then splitting it off into a separate mesh file in ../anatomy/.
% It then saves the Srf under the same name again and adds Srf.ConvertedFrom 
% to indicate the version from which it has been converted.
%
%   srfname:    Name of the SamSrf file (without .mat extension!)
%   surfdir:    Pathname to the 'surf' folder in the subject's Freesurfer folder
%
% 19/07/2020 - SamSrf 7 version (DSS)
% 14/10/2020 - Fixed typo in help section (DSS)
%              Will now load ASC files if binary files don't exist (DSS)
% 13/03/2022 - Now ensures file isn't just loaded from path (DSS)
%

%% Load old file
Mat = load(EnsurePath(srfname));

%% Load required surface data
P = fs_read_surf([surfdir filesep Mat.Srf.Hemisphere '.pial']); % Pial surface
I = fs_read_surf([surfdir filesep Mat.Srf.Hemisphere '.inflated']); % Inflated surface
S = fs_read_surf([surfdir filesep Mat.Srf.Hemisphere '.sphere']); % Spherical surface
try
    C = fs_read_curv([surfdir filesep Mat.Srf.Hemisphere '.curv']); % Cortical curvature 
    A = fs_read_curv([surfdir filesep Mat.Srf.Hemisphere '.area']); % Cortical surface area
    T = fs_read_curv([surfdir filesep Mat.Srf.Hemisphere '.thickness']); % Cortical thickness
catch
    % If no binary files available, hope there are ASC files, and use those
    C = Read_FreeSurfer([surfdir filesep Mat.Srf.Hemisphere '.curv.asc']); % Cortical curvature 
    A = Read_FreeSurfer([surfdir filesep Mat.Srf.Hemisphere '.area.asc']); % Cortical surface area
    T = Read_FreeSurfer([surfdir filesep Mat.Srf.Hemisphere '.thickness.asc']); % Cortical thickness
    C = C(:,5);
    A = A(:,5);
    T = T(:,5);
end

%% Add version number
if ~isfield(Mat.Srf, 'Version')
    Mat.Srf.Version = NaN;
end
Mat.Srf.ConvertedFrom = Mat.Srf.Version;
Mat.Srf.Version = samsrf_version;

%% Add fields for meshes
if ~isfield(Mat.Srf, 'Pial')
    Mat.Srf.Pial = P;
end
if ~isfield(Mat.Srf, 'Inflated')
    Mat.Srf.Inflated = I;
end
if ~isfield(Mat.Srf, 'Sphere')
    Mat.Srf.Sphere = S;
end

%% Add fields for surface stats
if ~isfield(Mat.Srf, 'Curvature')
    Mat.Srf.Curvature = C';
end
if ~isfield(Mat.Srf, 'Area')
    Mat.Srf.Area = A';
end
if ~isfield(Mat.Srf, 'Thickness')
    Mat.Srf.Thickness = T';
end

%% Remove Srf.Voxels if it exists
if isfield(Mat.Srf, 'Voxels')
    Mat.Srf = rmfield(Mat.Srf, 'Voxels');
end

%% Save file again
save(srfname, '-struct', 'Mat', '-v7.3');

%% Split off anatomical meshes
samsrf_anatomy_srf(srfname);