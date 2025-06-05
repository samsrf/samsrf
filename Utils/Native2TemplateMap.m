function Native2TemplateMap(NatSrf, MeshFolder, TmpFolder)
%
% Native2TemplateMap(NatSrf, MeshFolder, TmpFolder)
% 
% Warps whatever is in Srf.Data in the file NatSrf to the fsaverage template 
% using the lh/rh.sphere.reg files in the subject's surf folder (defined in MeshFolder).
% TmpFolder is the surf folder of the fsaverage template in the FreeSurfer subjects directory.
%
% If MeshFolder points to a GII file it loads this instead of the *h.sphere.reg.
% If using bilateral Srf files, the file name must contain wildcard *
% instead of the actual hemisphere letter.
%
% This procedure is a simple nearest neighbour transformation. It achieves
% very similar - but non-identical - results to FreeSurfer's own tool.
%
% Saves the spatially normalised Srf called the same as NatSrf with the suffix _sn.
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
    Hemis = {NatSrf.Hemisphere};
end

% Loop thru hemispheres (if needed)
for h = 1:length(Hemis)
    % If bilateral Srf
    if length(Hemis) == 2
        if h == 1
            samsrf_disp('Left hemisphere:');
            NatSrf = Lsrf; % Left hemisphere
        else
            samsrf_disp('Right hemisphere:');
            NatSrf = Rsrf; % Right hemisphere
        end
        NatSrfName = [Hemis{h} NatSrfName(3:end)]; % Update name
    end
    
    % Load surface vertices
    [V0 F] = fs_read_surf([TmpFolder filesep NatSrf.Hemisphere '.white']); % Grey-white surface
    P = fs_read_surf([TmpFolder filesep NatSrf.Hemisphere '.pial']); % Pial surface
    I = fs_read_surf([TmpFolder filesep NatSrf.Hemisphere '.inflated']); % Inflated surface
    S = fs_read_surf([TmpFolder filesep NatSrf.Hemisphere '.sphere']); % Spherical surface
    C = fs_read_curv([TmpFolder filesep NatSrf.Hemisphere '.curv']); % Cortical curvature 
    A = fs_read_curv([TmpFolder filesep NatSrf.Hemisphere '.area']); % Cortical surface area
    T = fs_read_curv([TmpFolder filesep NatSrf.Hemisphere '.thickness']); % Cortical thickness
    N = P - V0; % Cortical vectors for each vertex 
    
    % Create surface structure
    Srf = struct;
    Srf.Version = samsrf_version;
    Srf.Structural = NatSrf.Structural;
    Srf.Functional = NatSrf.Functional;
    Srf.Hemisphere = NatSrf.Hemisphere;
    Srf.Cortex_Steps = NatSrf.Cortex_Steps;
    Srf.Vertices = V0;
    Srf.Pial = P;
    Srf.Inflated = I;
    Srf.Sphere = S;
    Srf.Faces = F;
    Srf.Normals = N;
    Srf.Curvature = C';
    Srf.Area = A';
    Srf.Thickness = T';
    Srf.Data = [];
    Srf.Rule = 'Warped onto fsaverage';
    Srf.Values = NatSrf.Values;
        
    % Only if sphere data present
    if ~isfield(NatSrf, 'Sphere')
        samsrf_error('NatSrf contains no Sphere data!');
    end
    
    % Vertices
    NsrcVx = size(NatSrf.Sphere,1); % Number of source vertices
    NtgtVx = size(Srf.Sphere,1); % Number of target vertices
    NatVx = NatSrf.Sphere; % Source vertices
    TmpVx = Srf.Sphere; % Target vertices
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
    
    %% Vertex assignments
    t0 = tic;
    samsrf_disp('Warping native into template surface...');
    samsrf_disp(' Please stand by...');
    Sva = NaN(1,NtgtVx); % Source vertex assignment
    parfor v = 1:NtgtVx
        % Vector from current target vertex to all source vertices
        xyz = NatVx - repmat(TmpVx(v,:), NsrcVx, 1);
        % Euclidian distances 
        ed = sqrt(xyz(:,1).^2 + xyz(:,2).^2 + xyz(:,3).^2);
        % Minimal distance
        mv = find(ed==min(ed),1);
        Sva(v) = mv;
    end
    
    %% Reassign vertices
    Srf.Source_Vertices = Sva;
    Srf.Data = NatSrf.Data(:,Sva);
    
    %% Save Srf & ROIs
    save([NatSrfName '_sn'], 'Srf', '-v7.3');
    samsrf_disp(['Warping completed in ' num2str(toc(t0)/60) ' minutes.']);
    samsrf_newline;
end

% Combine hemisphere maps?
if length(Hemis) == 2
    samsrf_disp('Combining hemispheres...');
    L = load(['lh' NatSrfName(3:end) '_sn']);
    R = load(['rh' NatSrfName(3:end) '_sn']);
    Srf = samsrf_bilat_srf(L.Srf, R.Srf);
    if isfield(L, 'Model')
        % pRF or CF map
        Model = L.Model;
        save(['bi' NatSrfName(3:end) '_sn'], 'Srf', 'Model', '-v7.3');
    else
        % Other surface data
        save(['bi' NatSrfName(3:end) '_sn'], 'Srf', '-v7.3');
    end
    delete(['lh' NatSrfName(3:end) '_sn.mat']);
    delete(['rh' NatSrfName(3:end) '_sn.mat']);
end

cd(CurDir);