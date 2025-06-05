function Template2NativeMap(NatSrf, MeshFolder, TmpFolder)
%
% Template2NativeMap(NatSrf, MeshFolder, TmpFolder)
% 
% Warps the group average pRF map in the template lh/rh_pRF_fsaverage into 
% the native space belonging to NatSrf. The second input MeshFolder is the 
% path pointing to subject's surf folder which must contain the sphere.reg 
% linking the subject's native space with the fsaverage template. Finally,
% the third input TmpFolder points to the folder where the *h_pRF_fsaverage
% surface data files with the template map are located.
%
% If MeshFolder points to a GII file it loads this instead of the *h.sphere.reg.
% If using bilateral Srf files, the file name must contain wildcard *
% instead of the actual hemisphere letter.
%
% This procedure is a simple nearest neighbour transformation. It achieves
% very similar but non-identical results to FreeSurfer's tool.
%
% Saves a new Srf called lh/rh/bi_temp_map with the same number of vertices 
% as NatSrf but containing the map data from the group average template.
% Further, it saves the visual region labels from the average map in a
% separate folder called ROIs_temp_map.
%
% 19/09/2024 - Now accepts bilateral Srf input (DSS)
% 28/02/2025 - Fixed bug with too large files (DSS)
% 05/06/2025 - Can now also load GII registrations from HCP (DSS) 
%

% Load native map
NatSrfName = NatSrf;
warning off
load(EnsurePath(NatSrf));
warning on
NatSrf = samsrf_expand_srf(Srf);

% Determine folder 
CurDir = pwd;
[NatSrfDir, NatSrfName] = fileparts(NatSrfName);
if isempty(NatSrfDir)
    NatSrfDir = '.';
end
cd(NatSrfDir);

% Is bilateral Srf?
if strcmpi(NatSrf.Hemisphere, 'bi')
    samsrf_disp('Bilateral surface file: Running each hemisphere separately...');
    samsrf_newline;
    Hemis = {'lh' 'rh'};
    [Lsrf, Rsrf] = samsrf_hemi_srfs(NatSrf);
else
    Hemis = {Srf.Hemisphere};
end

% Loop thru hemispheres (if needed)
for h = 1:length(Hemis)
    % If bilateral Srf
    if length(Hemis) == 2
        if h == 1
            samsrf_disp('Left hemisphere:');
            NatSrf = Lsrf;
        else
            samsrf_disp('Right hemisphere:');
            NatSrf = Rsrf;
        end
    end

    % Load template
    load([TmpFolder filesep NatSrf.Hemisphere '_pRF_fsaverage']);
    TmpSrf = samsrf_expand_srf(Srf);
    TmpSrf.Data = [TmpSrf.Data; TmpSrf.Curvature]; % Add template curvature
    TmpSrf.Values{end+1} = 'TmpCurv';
    
    % Is there a noise ceiling in NatSrf?
    if isfield(NatSrf, 'Noise_Ceiling')
        NatSrf = rmfield(NatSrf, 'Noise_Ceiling'); % Remove this field
    end
    
    % Only if sphere data present
    if ~isfield(TmpSrf, 'Sphere')
        samsrf_error('TmpSrf contains no Sphere data!');
    end
    if ~isfield(NatSrf, 'Sphere')
        samsrf_error('NatSrf contains no Sphere data!');
    end
    
    % Vertices
    NsrcVx = size(TmpSrf.Sphere,1); % Number of source vertices
    NtgtVx = size(NatSrf.Sphere,1); % Number of target vertices
    TmpVx = TmpSrf.Sphere; % Source vertices
    NatVx = NatSrf.Sphere; % Target vertices
    [~,~,e] = fileparts(MeshFolder); % Extension of MeshFolder
    if strcmpi(e, '.gii')
        wc = strfind(MeshFolder, '*');
        if ~isempty(wc)
            MeshFolder(wc) = upper(NatSrf.Hemisphere(1)); % Replace wildcard
        end
        % Load GII registration
        RegVx = gifti(MeshFolder); 
        RegVx = RegVx.vertices; 
        MeshFolder(wc) = '*';
    else
        % Load FreeSurfer registration
        RegVx = fs_read_surf([MeshFolder filesep NatSrf.Hemisphere '.sphere.reg']); 
    end
    if size(NatVx,1) == size(RegVx,1)
        NatVx = RegVx;
    else
        samsrf_error('Number of registration vertices does not match native surface mesh!');
    end
    
    % Output Srf
    Srf = NatSrf;
    Srf.Functional = TmpSrf.Functional;
    Srf.Data = [];
    Srf.Values = TmpSrf.Values;
    
    %% Vertex assignments
    t0 = tic;
    samsrf_disp('Warping template into native surface...');
    samsrf_disp(' Please stand by...');
    Sva = NaN(1,NtgtVx); % Source vertex assignment
    parfor v = 1:NtgtVx
        % Vector from current target vertex to all source vertices
        xyz = TmpVx - repmat(NatVx(v,:), NsrcVx, 1);
        % Euclidian distances 
        ed = sqrt(xyz(:,1).^2 + xyz(:,2).^2 + xyz(:,3).^2);
        % Minimal distance
        mv = find(ed==min(ed),1);
        Sva(v) = mv;
    end
    
    %% Reassign vertices
    Srf.Source_Vertices = Sva;
    Srf.Data = TmpSrf.Data(:,Sva);
    
    %% Save Srf & ROIs
    save([Srf.Hemisphere '_temp_map'], 'Srf', '-v7.3');
    if ~isfield(TmpSrf, 'Roi_Names')
        samsrf_disp('WARNING: No ROI names defined in template map - assuming default names...');
        Rois = {'V1' 'V2' 'V3' 'V3A' 'V3B' 'V4' 'TO1' 'TO2'};
    else
        Rois = TmpSrf.Roi_Names; 
    end
    if h == 1 
        % Only need to do once
        mkdir('ROIs_temp_map');
    end
    for r = 1:length(Rois)
        rvx = find(Srf.Data(5,:)==r);
        if ~isempty(rvx)
            samsrf_srf2label(Srf, ['ROIs_temp_map' filesep Srf.Hemisphere '_' Rois{r}], 1, rvx);
        end
    end
    samsrf_disp(['Warping completed in ' num2str(toc(t0)/60) ' minutes.']);
    samsrf_newline;
end

% Combine hemisphere maps?
if length(Hemis) == 2
    samsrf_disp('Combining hemispheres...');
    L = load('lh_temp_map');
    R = load('rh_temp_map');
    Srf = samsrf_bilat_srf(L.Srf, R.Srf);
    save('bi_temp_map', 'Srf', '-v7.3');
    delete('lh_temp_map.mat');
    delete('rh_temp_map.mat');
    cd ROIs_temp_map
    for r = 1:length(Rois)
        samsrf_bilat_label(Srf, Rois{r});
    end
    delete('lh_*.label');
    delete('rh_*.label');
end
samsrf_newline;

cd(CurDir);