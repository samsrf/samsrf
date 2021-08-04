function Template2NativeMap(NatSrf, MeshFolder)
%
% Template2NativeMap(NatSrf, MeshFolder)
% 
% Warps the group average pRF map in the template lh/rh_pRF_fsaverage into 
% the native space belonging to NatSrf. The second input MeshFolder is the 
% path pointing to subject's surf folder which must contain the sphere.reg 
% linking the subject's native space with the fsaverage template.
%
% This procedure is a simple nearest neighbour transformation. It achieves
% very similar but non-identical results to FreeSurfer's tool.
%
% Saves a new Srf called lh/rh_temp_map with the same number of vertices as 
% NatSrf but containing the map data from the group average template.
% Further, it saves the visual region labels from the average map in a
% separate folder called ROIs_temp_map.
%
% Note: You will need to download the group average map templates from the
%       SamSrf OSF website or create your own template and this must be on
%       the Matlab path, using the name lh/rh_pRF_fsaverage. 
%
% 07/08/2020 - Written (DSS)
% 28/03/2021 - Now automatically removes noise ceiling from NatSrf (DSS)
% 14/04/2021 - Exports anatomical surfaces automatically now (DSS)
% 10/06/2021 - Now only saves ROI labels if they exist (DSS) 
% 12/07/2021 - Added stand-by message since parallel progress reports are a pain (DSS)
%

% Load native map
load(NatSrf);
NatSrf = samsrf_expand_srf(Srf);
% Load template (on Matlab path)
load([NatSrf.Hemisphere '_pRF_fsaverage']);
TmpSrf = samsrf_expand_srf(Srf);

% Is there a noise ceiling in NatSrf?
if isfield(NatSrf, 'Noise_Ceiling')
    NatSrf = rmfield(NatSrf, 'Noise_Ceiling'); % Remove this field
end

% Only if sphere data present
if ~isfield(TmpSrf, 'Sphere')
    error('TmpSrf contains no Sphere data!');
end
if ~isfield(NatSrf, 'Sphere')
    error('NatSrf contains no Sphere data!');
end

% Vertices
NsrcVx = size(TmpSrf.Sphere,1); % Number of source vertices
NtgtVx = size(NatSrf.Sphere,1); % Number of target vertices
TmpVx = TmpSrf.Sphere; % Source vertices
NatVx = NatSrf.Sphere; % Target vertices
RegVx = fs_read_surf([MeshFolder filesep NatSrf.Hemisphere '.sphere.reg']);
if size(NatVx,1) == size(RegVx,1)
    NatVx = RegVx;
else
    error('Number of registration vertices does not match native surface mesh!');
end

% Output Srf
Srf = NatSrf;
Srf.Functional = TmpSrf.Functional;
Srf.Data = [];
Srf.Values = TmpSrf.Values;

%% Vertex assignments
t0 = tic;
disp('Warping template into native surface...');
disp(' Please stand by...');
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
save([Srf.Hemisphere '_temp_map'], 'Srf');
samsrf_anatomy_srf([Srf.Hemisphere '_temp_map']);
Rois = {'V1' 'V2' 'V3' 'V3A' 'V3B' 'V4' 'TO1' 'TO2'};
mkdir('ROIs_temp_map');
for r = 1:length(Rois)
    rvx = find(Srf.Data(5,:)==r);
    if ~isempty(rvx)
        samsrf_srf2label(Srf, ['ROIs_temp_map' filesep Srf.Hemisphere '_' Rois{r}], 1, rvx);
    end
end
disp(['Warping completed in ' num2str(toc(t0)/60) ' minutes.']);
new_line;
