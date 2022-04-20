function Native2TemplateMap(NatSrf, MeshFolder, TmpFolder)
%
% Native2TemplateMap(NatSrf, MeshFolder, TmpFolder)
% 
% Warps whatever is in Srf.Data in the file NatSrf to the fsaverage template 
% using the lh/rh.sphere.reg files in the subject's surf folder (defined in MeshFolder).
% TmpFolder is the surf folder of the fsaverage template in the FreeSurfer subjects directory.
%
% This procedure is a simple nearest neighbour transformation. It achieves
% very similar - but non-identical - results to FreeSurfer's own tool.
%
% Saves the spatially normalised Srf called the same as NatSrf with the suffix _sn.
%
% 20/04/2022 - SamSrf 8 version (DSS)
%

% Load native map
NatSrfName = NatSrf;
load(EnsurePath(NatSrf));
NatSrf = samsrf_expand_srf(Srf);

% Create template Srf
curdir = cd;
cd(TmpFolder);

% Load surface vertices
[V0 F] = fs_read_surf([NatSrf.Hemisphere '.white']); % Grey-white surface
P = fs_read_surf([NatSrf.Hemisphere '.pial']); % Pial surface
I = fs_read_surf([NatSrf.Hemisphere '.inflated']); % Inflated surface
S = fs_read_surf([NatSrf.Hemisphere '.sphere']); % Spherical surface
C = fs_read_curv([NatSrf.Hemisphere '.curv']); % Cortical curvature 
A = fs_read_curv([NatSrf.Hemisphere '.area']); % Cortical surface area
T = fs_read_curv([NatSrf.Hemisphere '.thickness']); % Cortical thickness
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

cd ..
cd(curdir);

% Only if sphere data present
if ~isfield(NatSrf, 'Sphere')
    error('NatSrf contains no Sphere data!');
end

% Vertices
NsrcVx = size(NatSrf.Sphere,1); % Number of source vertices
NtgtVx = size(Srf.Sphere,1); % Number of target vertices
NatVx = NatSrf.Sphere; % Source vertices
TmpVx = Srf.Sphere; % Target vertices
RegVx = fs_read_surf([MeshFolder filesep NatSrf.Hemisphere '.sphere.reg']);
if size(NatVx,1) == size(RegVx,1)
    NatVx = RegVx;
else
    error('Number of registration vertices does not match native surface mesh!');
end

%% Vertex assignments
t0 = tic;
disp('Warping native into template surface...');
disp(' Please stand by...');
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
save([NatSrfName '_sn'], 'Srf');
disp(['Warping completed in ' num2str(toc(t0)/60) ' minutes.']);
new_line;
