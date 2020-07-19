function samsrf_anatomy_srf(SrfName, AnatPath)
% 
% samsrf_anatomy_srf(SrfName, [AnatPath='../anatomy/'])
%
% Splits the anatomical meshes in SrfName and saves them as a separate 
%  file in the structure Anat. The file name is based on Srf.Structural and
%  the path defaults to '../anat/'. If you want the path to be something
%  else you can change this in the optional input AnatPath.
%
% If the anatomical mesh file already exists, the function shows a warning
%  and will not save it again. You would normally only need one such file 
%  for a given recon-all in FreeSurfer so you don't want to save a new one.
%  BUT IT IS YOUR RESPONSIBILITY TO CHECK THAT THIS IS THE CORRECT FILE!
%
% Srf is then saved without then anatomical meshes under the original name.
%
% Henceforth, Srf.Meshes contains the path name to the anatomical meshes.
%  Whenever you call samsrf_expand_srf it will automatically load them.
%  Conversely, when you call samsrf_compress_srf it will remove them.
%
% Note that Srf must be a SamSrf 6 Srf which contains all of the anatomical
%  data. If the Srf is from an older version, you will first need to
%  convert it using samsrf_convert_old_srf. Any new Srfs generated using
%  samsrf_vol2srf in SamSrf 6 will have automatically been separated though.
%
% 19/07/2020 - SamSrf 7 version (DSS)
%

%% Default path?
if nargin < 2
    AnatPath = ['..' filesep 'anatomy' filesep];
end

%% Load the full Srf 
F = load(SrfName);
F.Srf = samsrf_expand_srf(F.Srf);

%% Only continue if anatomy present
if isfield(F.Srf, 'Vertices')
    %% Create anatomical Srf
    Anat = struct;
    Anat.Version = F.Srf.Version;
    Anat.Structural = F.Srf.Structural;
    Anat.Hemisphere = F.Srf.Hemisphere;
    Anat.Vertices = F.Srf.Vertices;
    Anat.Faces = F.Srf.Faces;
    Anat.Normals = F.Srf.Normals;
    Anat.Pial = F.Srf.Pial;
    Anat.Inflated = F.Srf.Inflated;
    Anat.Sphere = F.Srf.Sphere;
    Anat.Curvature = F.Srf.Curvature;
    Anat.Area = F.Srf.Area;
    Anat.Thickness = F.Srf.Thickness;

    %% Remove fields from Srf
    F.Srf = rmfield(F.Srf, 'Vertices');
    F.Srf = rmfield(F.Srf, 'Faces');
    F.Srf = rmfield(F.Srf, 'Normals');
    F.Srf = rmfield(F.Srf, 'Pial');
    F.Srf = rmfield(F.Srf, 'Inflated');
    F.Srf = rmfield(F.Srf, 'Sphere');
    F.Srf = rmfield(F.Srf, 'Curvature');
    F.Srf = rmfield(F.Srf, 'Area');
    F.Srf = rmfield(F.Srf, 'Thickness');

    %% Anatomical mesh name
    AnatName = Anat.Structural;
    AnatName([strfind(AnatName,'/') strfind(AnatName,'\')]) = filesep;
    [~,AnatName,~] = fileparts(AnatName);
    % Add link to Srf
    F.Srf.Meshes = [AnatPath Anat.Hemisphere '_' AnatName];

    %% Save files
    save(SrfName, '-struct', 'F', '-v7.3'); % Save Srf without anatomy
    if ~exist(AnatPath, 'dir')
        mkdir(AnatPath); % Make anatomy folder
    end
    if exist([AnatPath Anat.Hemisphere '_' AnatName '.mat'], 'file')
        warning('Anatomical mesh file already exists so won''t save another one...');
    else
        save([AnatPath Anat.Hemisphere '_' AnatName], 'Anat', '-v7.3'); % Save anatomy
        disp(['Saved anatomical meshes in ' AnatPath Anat.Hemisphere '_' AnatName]);
    end
end