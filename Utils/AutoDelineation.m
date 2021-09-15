function AutoDelineation(SrfName, NatMesh, TmpMesh, Thresh, Niter)
%
% AutoDelineation(Srf, NatMesh, TmpMesh, [Thresh, Niter])
%
% Fits the borders of visual regions V1-V3A/V4. It shows a movie of the 
% process so you can see how well it's going. Then it saves a delineation 
% file del_SrfName that you can use in DelineationTool for labelling ROIs.
% 
% THIS IS NOT A MAGIC WANT! You -must- inspect the maps in DealineationTool 
% & probably manually correct some smaller errors. Also, this algorithm 
% expects gentle retinotopic map gradients, so you probably need to heavily 
% smooth (spherical FWHM of 5-10mm?) your maps (see also samsrf_surfcalcs).
%
% The automatic delineation is based on the visual field meridians of a 
% group average map normalised to the fsaverage template. These initial
% borders are stored in lh/rh_WayPoints.mat in SamSrf/Utils. This contains
% waypoints for the meridians, ordered from the foveal confluence outwards
% into the periphery & by ascending visual area. 
%
% The algorithm works as follows:
%   1.  The waypoints comprising the normalised meridians are warped back
%       into native brain space (so you need the fsaverage template).
%   2.  Your retinotopic map is restricted to a particular range, using a
%       goodness-of-fit threshold & minimum & maximum eccentricity.
%   3.  The algorithm then searches the geodesic neighbourhood of each
%       waypoint for the maximum/minimum polar angle (depends on which 
%       meridian it is fitting). This search is constrained by being
%       limited to the eccentriity band of its original position (the width
%       of this band is adjustable). It is also constrained by the fact
%       that the algorithm expects a minimum Euclidean distance in sphere 
%       space from each waypoints to all others. 
%   4.  The search process is repeated for several iterations (adjustable).
%       The size of the search region is gradually reduced from 4 geodesic 
%       steps down to 1 geodesic step over the course of the iterations.
%   5.  Upon completion of the search process, the waypoints in each
%       meridian are connected into complete paths. 
%   6.  A path surrounding the whole thresholded region is also created.
%   7.  The end points of the meridian lines are then connected to this
%       surrounding path to complete the regional boundaries.
%   8.  Finally, the paths are expanded to be compatible with region
%       filling in the DelineationTool & the delineation file is saved.
%
% The default thresholds & number of iterations was chosen based on my test
% data. You may need to tweak this for your own needs. But this makes your
% initial auto-delineation reproducible. The parameters used are saved in
% the delineation file for posterity.
% 
%   SrfName:    Name of input surface data file with retinotopic map 
%   NatMesh:    Subject's surf folder which must contain lh/rh.sphere.reg
%   TmpMesh:    Template's surf folder which must contain lh/rh.sphere
%   Thresh:     Optional: 1x4 vector determining thresholds:
%                   (1) R^2 threshold (default = 0.1)
%                   (2) Central eccentricity limit (default = 0.5)
%                   (3) Peripheral eccentricity limit (default = 12)
%                   (4) With of eccentricity band (default = 0.5)
%   Niter:      Number of search iterations (default = 10)
%
% 14/09/2021 - Written first version (DSS)
% 15/09/2021 - Added distance constraint between neighbouring waypoints (DSS)
%

%% ROIs in delineation file (Change at your own leisure/peril!) %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% 
RoiList = {'V1' 'V2v' 'V3v' 'V4' 'V2d' 'V3d' 'V3A' 'V3B' 'LO1' 'LO2' 'VO1' 'VO2' 'TO1' 'TO2' 'V6' 'IPS0' 'IPS1' 'IPS2'}'; % ROI list
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% 

%% Parameters
RoiSeeds = NaN(length(RoiList),1); % Initialise list of ROI seed vertices
DefThresh = [0.1 0.5 12 0.5]; % Default thresholds (R^2, MinEcc, MaxEcc, EccBandwidth)
if nargin < 4 
    Thresh = DefThresh;
end
if length(Thresh) < 4
    Thresh = [Thresh DefThresh(length(Thresh)+1:end)];
end
if nargin < 5
    Niter = 20; % Number of search iterations
end
new_line;
disp('SamSrf - Automatic delineation tool');
disp([' R^2 threshold: ' num2str(Thresh(1))]);
disp([' Minimum eccentricity:   ' num2str(Thresh(2))]);
disp([' Maximum eccentricity:   ' num2str(Thresh(3))]);
disp([' Eccentricity bandwidth: ' num2str(Thresh(4))]);
disp([' Number of iterations:   ' num2str(Niter)]);

%% Load map
new_line;
disp(['Loading map ' SrfName '.mat...']);
load(SrfName, 'Srf');
Srf = samsrf_expand_srf(Srf);
if ~isfield(Srf, 'Values')
    error('This file does not contain a pRF map!');
end

%% Viewing region
load('SamSrf_defaults.mat', 'def_disproi');
if ~exist('def_disproi')
    def_disproi = NaN; 
end
if ~isempty(def_disproi) && (def_disproi(1) == '<' || def_disproi(1) == '>')
    % If ROI defined by coordinates
    if length(def_disproi) == 1
        error('You must define inflated mesh coordinate in def_disproi!');
    end
    switch def_disproi(2)
        case 'X'
            wv = Srf.Inflated(:,1);
        case 'Y'
            wv = Srf.Inflated(:,2);
        case 'Z'
            wv = Srf.Inflated(:,3);
        otherwise
            error('Invalid inflated mesh coordinate specified in def_disproi!');
    end
    if length(def_disproi) < 3
        error('You must define inflated mesh cut-off coordinate in def_disproi!');
    end
    wc = str2double(def_disproi(3:end));
    if def_disproi(1) == '<'
        ViewRoi = find(wv < wc);
    elseif def_disproi(1) == '>'
        ViewRoi = find(wv > wc);
    end
    disp(['Only displaying inflated mesh vertices with ' def_disproi]);
elseif isnan(def_disproi)
    % Use ROI from Srf
    disp('Using ROI from Srf if it exists');
else
    % Load region of interest
    Roi = [SrfName(1:2) '_' def_disproi];
    ViewRoi = samsrf_loadlabel(Roi);
    disp(['Only displaying ROI ' Roi]);
end
if ~isnan(ViewRoi)
    AllVx = true(size(Srf.Vertices,1),1);
    AllVx(ViewRoi) = false;
    Srf.Vertices(AllVx,:) = NaN;
end
new_line;

%% Sphere vertices in native & template 
NatVx = Srf.Sphere; % Target vertices
RegVx = fs_read_surf([NatMesh filesep Srf.Hemisphere '.sphere.reg']); % Registration vertices
TmpVx = fs_read_surf([TmpMesh filesep Srf.Hemisphere '.sphere']); % Template vertices
NnatVx = size(NatVx,1); % Number of target vertices
if size(NatVx,1) == size(RegVx,1)
    NatVx = RegVx;
else
    error('Number of registration vertices does not match native surface mesh!');
end

%% Load waypoints
Wpts = load([Srf.Hemisphere '_WayPoints']);
NatWpts = Wpts; % Waypoints after warping
PathTypes = fieldnames(Wpts);

%% Warp waypoints into native space
disp('Warping template waypoints into native space...');
% Loop thru path types
for t = 1:length(PathTypes)
    disp([' ' PathTypes{t}]);
    % Loop thru paths of this type
    for m = 1:length(Wpts.(PathTypes{t}))
        % Loop thru waypoints in each path
        for w = 1:length(Wpts.(PathTypes{t}){m})
            % Current waypoint vertex
            v = Wpts.(PathTypes{t}){m}(w);
            % Vector from current target vertex to all source vertices
            xyz = NatVx - repmat(TmpVx(v,:), NnatVx, 1);
            % Euclidian distances 
            ed = sqrt(xyz(:,1).^2 + xyz(:,2).^2 + xyz(:,3).^2);
            % Minimal distance
            mv = find(ed==min(ed),1);
            % Store warped vertex
            NatWpts.(PathTypes{t}){m}(w) = mv;
        end
        disp(['  ' num2str(m) '/' num2str(length(Wpts.(PathTypes{t}))) ': ' num2str(length(Wpts.(PathTypes{t}){m})) ' waypoints']);
    end
end

%% List of original waypoints
% Loop thru path types
Dots = []; % Path vertices 
for t = 1:length(PathTypes)
    % Loop thru paths of this type
    for m = 1:length(NatWpts.(PathTypes{t}))
        % Loop thru waypoints in each path
        for w = 1:length(NatWpts.(PathTypes{t}){m})
            % Current waypoint 
            v = NatWpts.(PathTypes{t}){m}(w); % Vertex index
            % Surround each waypoint for visibility
            Dots = [Dots; samsrf_borderpath(Srf, v)]; 
        end
    end
end

%% Process retinotopic map
Theta = atan2(Srf.Data(3,:), abs(Srf.Data(2,:))) / pi*180; % Polar angle mirrored into right hemifield
Rho = sqrt(Srf.Data(2,:).^2 + Srf.Data(3,:).^2); % Eccentricity
Gof = Srf.Data(1,:); % Goodness of fit
Mask = Gof <= Thresh(1) | Rho < Thresh(2) | Rho > Thresh(3); % Restrict to range
Theta(Mask) = NaN; % Mask out of range vertices
EccWpts = NatWpts; % To store original eccentricity
% Loop thru path types
for t = 1:length(PathTypes)
    % Loop thru paths of this type
    for m = 1:length(NatWpts.(PathTypes{t}))
        % Loop thru waypoints in each path
        for w = 1:length(NatWpts.(PathTypes{t}){m})
            v = NatWpts.(PathTypes{t}){m}(w); % Current waypoint 
            EccWpts.(PathTypes{t}){m}(w) = Rho(v); % Eccentricity of original position
        end
    end
end
% Display map
h = samsrf_surf(Srf, 'Sphere', [Thresh(1) 0 0 Thresh(2:3) -.1], {Dots; [NaN 1 1 1]}, '', 'Polar');
set(gcf, 'name', 'Initial waypoints');
pause(.5);

%% Fit each waypoint to map
disp('Fitting waypoints to empirical map...');
% Loop thru iterations
FitWpts = NatWpts; % Fitted waypoints
for i = 1:Niter
    % Adjust search radius?
    Radius = floor((Niter-i) / (Niter/4)) + 1;
    
    % Loop thru path types
    Dots = []; % Path vertices 
    cm = 1; % Count paths
    for t = 1:length(PathTypes)
        
        % Loop thru paths of this type
        for m = 1:length(NatWpts.(PathTypes{t}))
            
            % Loop thru waypoints in each path
            for w = 1:length(NatWpts.(PathTypes{t}){m})
                % Current waypoint 
                v = NatWpts.(PathTypes{t}){m}(w); % Vertex index
                % Only good points
                if ~isnan(v)
                    ecc = EccWpts.(PathTypes{t}){m}(w); % Original eccentricity 
                    roi = samsrf_georoi(v, Radius, Srf.Vertices, Srf.Faces); % Search region around waypoint
                    roi = roi(Rho(roi) > ecc-Thresh(4)/2 & Rho(roi) < ecc+Thresh(4)/2); % Limit to eccentricity range

                    % Constrain distance to other waypoints
                    roiVx = NatVx(roi,:); % Sphere coordinates of search region vertices
                    aw = Dots; % All waypoints
                    aw(isnan(aw)) = []; % Remove bad points
                    aw(aw == v) = []; % Remove current point 
                    % Are there other points?
                    if ~isempty(aw)
                        nbVx = NatVx(aw,:); % Sphere coordinates of other waypoints
                        goodroi = true(size(roi)); % Label search region vertices that are good
                        for j = 1:length(aw)
                            ds = sqrt(sum((roiVx - repmat(nbVx(j,:), length(roi), 1)).^2,2)); % Distances to this waypoint
                            goodroi(ds < 0.5) = false; % Untag search region vertices too close to this waypoint
                        end
                        roi = roi(goodroi);
                    end
                    
                    % Senseless to search a single vertex
                    if length(roi) > 1
                        % Current polar angles
                        T = Theta(roi); 
                        % Locate nearest fitting vertex
                        if strcmpi(PathTypes{t}, 'UpperVerticalMeridian') 
                            mv = roi(find(T == nanmax(T),1)); % Maximal polar angle          
                        elseif strcmpi(PathTypes{t}, 'LowerVerticalMeridian') 
                            mv = roi(find(T == nanmin(T),1)); % Minimal polar angle          
                        elseif strcmpi(PathTypes{t}, 'HorizontalMeridian') 
                            T = abs(T);
                            mv = roi(find(T == nanmin(T),1)); % Minimal absolute polar angle            
                        end

                        % Store fitted vertex
                        if ~isempty(mv)
                            FitWpts.(PathTypes{t}){m}(w) = mv;
                        else
                            FitWpts.(PathTypes{t}){m}(w) = NaN; % Nothing was found
                        end
                    end
                end
            end
            
            % Surround each waypoint for visibility
            Dots = [Dots; FitWpts.(PathTypes{t}){m}];  
            % Increase counter
            cm = cm + 1;  
        end
    end
    % Replace old waypoints
    NatWpts = FitWpts;
    
    % Show new waypoints
    samsrf_surf(Srf, 'Sphere', [Thresh(1) 0 0 Thresh(2:3) -.1], {samsrf_borderpath(Srf, Dots); [NaN 1 1 1]}, '', 'Polar', h);
    set(gcf, 'name', ['Iteration #' num2str(i)]);
    pause(.5);
end
new_line;

%% Create paths
Paths = {};

% Connecting waypoints
disp('Creating paths from waypoints...');
% Loop thru path types
cw = 1; % Count sub-paths
for t = 1:length(PathTypes)
    % Loop thru paths of this type
    for m = 1:length(NatWpts.(PathTypes{t}))
        % Loop thru waypoints in each path
        for w = 1:length(NatWpts.(PathTypes{t}){m})-1
            CurPts = NatWpts.(PathTypes{t}){m}(w:w+1); % Waypoints of current path
            CurPts(isnan(CurPts)) = []; % Remove bad points
            if length(CurPts) == 2 % Only if two points to connect
                Paths{cw} = ConnectWaypoints(CurPts, Srf.Sphere); % Connect path
                cw = cw + 1; % Increase counter
            end
        end
    end
end

% Add surrounding path
disp('Drawing path surrounding masked region');
sp = samsrf_borderpath(Srf, find(~Mask)); % Path surrounding masked region
% Paths{end+1} = sp; % Add surrounding path

% Connect end points to surrounding path
disp('Connecting end points with surrounding path');
Foveal = [];
Peripheral = [];
% Loop thru path types
for t = 1:length(PathTypes)
    % Loop thru end points
    for j = 1:2
        % Loop thru paths of this type
        for m = 1:3 
            if ~(m == 3 && strcmpi(PathTypes{t}, 'HorizontalMeridian')) % Only 2 horizontal meridians
                CurPts = NatWpts.(PathTypes{t}){m}; % Waypoints in this path
                CurPts(isnan(CurPts)) = []; % Remove bad points
                if ~isempty(CurPts)
                    if j == 1 
                        e = CurPts(1); % Foveal end point
                    else 
                        e = CurPts(end); % Peripheral end point
                    end
                    cvx = NatVx(e,:); % Current end point coordinates
                    svx = NatVx(sp,:); % Surrounding path coordinates
                    ds = sqrt((svx(:,1)-cvx(1)).^2 + (svx(:,2)-cvx(2)).^2 + (svx(:,3)-cvx(3)).^2); % Euclidian distance from end point to surrounding path vertices 
                    v = sp(find(ds == min(ds),1)); % Nearest vertex on the path
                    Paths{end+1} = ConnectWaypoints([e v], Srf.Sphere); % Connect path
                    % Store end points for later
                    if j == 1
                        Foveal = [Foveal v]; % Foveal end point
                    else
                        Peripheral = [Peripheral v]; % Foveal end point
                    end
                end
            end
        end 
    end
end

% Connect new foveal end points
disp('Connecting foveal end points...');
[~,x] = sort(NatVx(Foveal,3)); % Sort points based on Z-coordinate
Foveal = Foveal(x); % Reorder vertex indeces
for i = 1:length(Foveal)-1
    Paths{end+1} = ConnectWaypoints(Foveal(i:i+1), Srf.Sphere); % Connect path
end
% Connect new peripheral end points
disp('Connecting peripheral end points...');
[~,x] = sort(NatVx(Peripheral,3)); % Sort points based on Z-coordinate
Peripheral = Peripheral(x); % Reorder vertex indeces
for i = 1:length(Peripheral)-1
    Paths{end+1} = ConnectWaypoints(Peripheral(i:i+1), Srf.Sphere); % Connect path
end

% Paths for display
DrawPaths = Paths;
DrawPaths{end+1} = [NaN 1 1 1];
% Expand all paths
for p = 1:length(Paths)
    Paths{p} = [Paths{p}; samsrf_neighbours(Paths{p}, Srf.Faces)]; 
end

%% Save delineation
AutoDelinParams = struct;
AutoDelinParams.R2_Threshold = Thresh(1);
AutoDelinParams.Min_Eccentricity = Thresh(2);
AutoDelinParams.Max_Eccentricity = Thresh(3);
AutoDelinParams.Eccentricity_Bandwidth = Thresh(4);
AutoDelinParams.Number_of_Iterations = Niter;
% Vector with all path vertices 
Vs = [];
for i = 1:length(Paths)
    Vs = [Vs; Paths{i}];
end
% Saves everything to disc
new_line;
save(['del_' SrfName], 'SrfName', 'Vs', 'Paths', 'RoiList', 'RoiSeeds', 'AutoDelinParams');
disp(['Saved auto-delineation del_' SrfName '.mat']);

%% Display final delineation
samsrf_surf(Srf, 'Sphere', [Thresh(1) 0 0 Thresh(2:3) -.1], DrawPaths, '', 'Polar', h);
set(gcf, 'name', 'Inspect auto delineation');
rotate3d;
disp('Don''t forget to inspect & correct in DelineationTool before labelling the ROIs!');
new_line;



%% Draw connecting path
function Vs_path = ConnectWaypoints(Points, Vertices)
% Provide vector Points with two waypoints & Srf.Vertices
    Vs_path = [];
    v1 = Points(1); % 1st vertex
    v2 = Points(2); % 2nd vertex
    p2c = Vertices(v2,:) - Vertices(v1,:); % Vector from previous to current vertex
    vsp = []; % Vertices on new path
    for s = 0:.01:1
        cs = Vertices(v1,:) + p2c*s; % One small step along vector
        d = sqrt((Vertices(:,1)-cs(1)).^2 + (Vertices(:,2)-cs(2)).^2 + (Vertices(:,3)-cs(3)).^2); % Euclidian distance of all vertices from current step
        nv = find(d == min(d), 1); % New vertex on path
        vsp = [vsp; nv];
    end
    Vs_path = [Vs_path; vsp];
end

end
