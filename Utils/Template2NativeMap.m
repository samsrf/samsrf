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
% This procedure is a simple nearest neighbour transformation. It achieves
% very similar but non-identical results to FreeSurfer's tool.
%
% Saves a new Srf called lh/rh_temp_map with the same number of vertices as 
% NatSrf but containing the map data from the group average template.
% Further, it saves the visual region labels from the average map in a
% separate folder called ROIs_temp_map.
%
% 14/03/2022 - Now requires the pathname for the template map (DSS)
%              Warped map data now contains row with template curvatures (DSS)
%              Ensures now that random files aren't loaded from path (DSS)
% 25/07/2022 - Fixed incorrect help section (DSS)
%

% Load native map
load(EnsurePath(NatSrf));
NatSrf = samsrf_expand_srf(Srf);
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
save([Srf.Hemisphere '_temp_map'], 'Srf');
samsrf_anatomy_srf([Srf.Hemisphere '_temp_map']);
if ~isfield(TmpSrf, 'Roi_Names')
    warning('No ROI names defined in template map - assuming default names...');
    Rois = {'V1' 'V2' 'V3' 'V3A' 'V3B' 'V4' 'TO1' 'TO2'};
else
    Rois = TmpSrf.Roi_Names; 
end
mkdir('ROIs_temp_map');
for r = 1:length(Rois)
    rvx = find(Srf.Data(5,:)==r);
    if ~isempty(rvx)
        samsrf_srf2label(Srf, ['ROIs_temp_map' filesep Srf.Hemisphere '_' Rois{r}], 1, rvx);
    end
end
samsrf_disp(['Warping completed in ' num2str(toc(t0)/60) ' minutes.']);
samsrf_newline;
