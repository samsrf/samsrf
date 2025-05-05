function samsrf_anatomy_srf(SrfName, AnatPath)
% 
% samsrf_anatomy_srf(SrfName, [AnatPath='../anatomy/'])
%
% Splits the anatomical meshes in SrfName and saves them as a separate 
%  file in the structure Anat. The file name is based on Srf.Structural and
%  the path defaults to '../anatomy/'. If you want the path to be something
%  else you can change this in the optional input AnatPath. If this input 
%  is empty, the function does nothing.
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
% 14/03/2022 - Now supports anonymised Srfs (DSS)
% 20/04/2022 - SamSrf 8 version (DSS)
% 19/12/2023 - Added option not to split off anatomy (DSS)
%

%% Default path?
if nargin < 2
    AnatPath = ['..' filesep 'anatomy' filesep];
end

%% Are we splitting anatomy?
if ~isempty(AnatPath)
    %% Load the full Srf
    F = load(EnsurePath(SrfName));
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
        if isfield(F.Srf, 'Normals')
            % Don't exist for anonymised Srfs
            Anat.Normals = F.Srf.Normals;
            Anat.Pial = F.Srf.Pial;
        end
        Anat.Inflated = F.Srf.Inflated;
        Anat.Sphere = F.Srf.Sphere;
        Anat.Curvature = F.Srf.Curvature;
        if isfield(F.Srf, 'Area')
            % Don't exist for anonymised Srfs
            Anat.Area = F.Srf.Area;
            Anat.Thickness = F.Srf.Thickness;
        end    

        %% Remove fields from Srf
        F.Srf = rmfield(F.Srf, 'Vertices');
        F.Srf = rmfield(F.Srf, 'Faces');
        if isfield(F.Srf, 'Normals')
            % Don't exist for anonymised Srfs
            F.Srf = rmfield(F.Srf, 'Normals');
            F.Srf = rmfield(F.Srf, 'Pial');
        end
        F.Srf = rmfield(F.Srf, 'Inflated');
        F.Srf = rmfield(F.Srf, 'Sphere');
        F.Srf = rmfield(F.Srf, 'Curvature');
        if isfield(F.Srf, 'Area')
            % Don't exist for anonymised Srfs
            F.Srf = rmfield(F.Srf, 'Area');
            F.Srf = rmfield(F.Srf, 'Thickness');
        end

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
            samsrf_disp('WARNING: Anatomical mesh file already exists so won''t save another one...');
        else
            save([AnatPath Anat.Hemisphere '_' AnatName], 'Anat', '-v7.3'); % Save anatomy
            samsrf_disp(['Saved anatomical meshes in ' AnatPath Anat.Hemisphere '_' AnatName]);
        end
    end
end