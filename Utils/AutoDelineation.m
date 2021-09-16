function AutoDelineation(SrfName, NatMesh, TmpMesh, R2Thresh, MinEcc, MaxEcc, EccBw, Niter, InitRad)
%
% AutoDelineation(Srf, NatMesh, TmpMesh, [R2Thresh, MinEcc, MaxEcc, EccBw, Niter, InitRad])
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
%       The size of the search region is gradually reduced from InitRad 
%       geodesic steps down to 1 step over the course of Niter iterations.
%   5.  Upon completion of the search process, the waypoints in each
%       meridian are connected into complete paths. 
%   6.  A path surrounding the whole thresholded region is also created.
%   7.  The end points of the meridian lines are then connected to this
%       surrounding path to complete the regional boundaries.
%   8.  Finally, the paths are expanded to be compatible with region
%       filling in the DelineationTool & the delineation file is saved.
%
% The default thresholds & number of iterations was chosen based on my test
% data. You may need to tweak this for your own needs. It should be robust
% to relatively small differences in designs but it obviously cannot do a
% great job if you used a much larger field of view for instance.
%
% The parameters used are saved in the delineation file for posterity.
% Even though you may need to correct & adjust the automatic delineation,
% this makes your delineation much more reproducible. The auto-delineation
% is saved inside the delineation file for posterity so you can compare
% this to the final version.
% 
%   SrfName:    Name of input surface data file with retinotopic map 
%   NatMesh:    Subject's surf folder which must contain lh/rh.sphere.reg
%   TmpMesh:    Template's surf folder which must contain lh/rh.sphere
%
% Optional additional input arguments:
%   R2Thresh:   R^2 threshold (default = 0.1)
%   MinEcc:     Central eccentricity limit (default = 0.5)
%   MaxEcc:     Peripheral eccentricity limit (default = 9)
%   EccBw:      Width of eccentricity band (default = 0.5)
%   Niter:      Number of search iterations (default = 20)
%   InitRad:    Initial search radius (default = 2)
%
% 16/09/2021 - Completed (DSS)
%

%% ROIs in delineation file (Change at your own leisure/peril!) %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% 
RoiList = {'V1' 'V2v' 'V3v' 'V4' 'V2d' 'V3d' 'V3A' 'V3B' 'LO1' 'LO2' 'VO1' 'VO2' 'TO1' 'TO2' 'V6' 'IPS0' 'IPS1' 'IPS2'}'; % ROI list
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% 

%% Parameters
RoiSeeds = NaN(length(RoiList),1); % Initialise list of ROI seed vertices
if nargin < 4 
    R2Thresh = 0.1; % Goodness-of-fit threshold
end
if nargin < 5 
    MinEcc = 0.5; % Minimal eccentricity
end
if nargin < 6 
    MaxEcc = 9; % Maximal eccentricity
end
if nargin < 7 
    EccBw = 0.5; % Eccentricity band width
end
if nargin < 8
    Niter = 20; % Number of search iterations
end
if nargin < 9
    InitRad = 2; % Initial search radius
end
new_line;
disp('SamSrf - Automatic delineation tool');
disp([' R^2 threshold:          ' num2str(R2Thresh)]);
disp([' Minimum eccentricity:   ' num2str(MinEcc)]);
disp([' Maximum eccentricity:   ' num2str(MaxEcc)]);
disp([' Eccentricity bandwidth: ' num2str(EccBw)]);
disp([' Number of iterations:   ' num2str(Niter)]);
disp([' Initial search radius:  ' num2str(InitRad)]);
% mkdir('Png');

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
Mask = Gof <= R2Thresh | Rho < MinEcc | Rho > MaxEcc; % Restrict to range
Theta(Mask) = NaN; % Mask out of range vertices
EccWpts = NatWpts; % To store original eccentricity
% Loop thru path types
for t = 1:length(PathTypes)
    % Loop thru paths of this type
    for m = 1:length(NatWpts.(PathTypes{t}))
        % Loop thru waypoints in each path
        for w = 1:length(NatWpts.(PathTypes{t}){m})
            v = NatWpts.(PathTypes{t}){m}(w); % Current waypoint
            if strcmpi(PathTypes{t}, 'FovealBorder')
                EccWpts.(PathTypes{t}){m}(w) = MinEcc; % Minimal eccentricity
            elseif strcmpi(PathTypes{t}, 'PeripheralBorder')
                EccWpts.(PathTypes{t}){m}(w) = MaxEcc; % Maximal eccentricity
            else % Meridian paths
                EccWpts.(PathTypes{t}){m}(w) = Rho(v); % Eccentricity of original position
            end
        end
    end
end
% Display map
h = samsrf_surf(Srf, 'Sphere', [R2Thresh 0 0 MinEcc MaxEcc -.1], {Dots; [NaN 1 1 1]}, '', 'Polar');
set(gcf, 'name', 'Initial waypoints');
pause(.5);
new_line; 

%% Adjust peripheral border
disp('Adjusting peripheral border...');
sp = samsrf_borderpath(Srf, find(~Mask)); % Path surrounding masked region
t = find(strcmpi(PathTypes, 'PeripheralBorder')); % Select peripheral borders
% Loop thru peripheral border paths 
for m = 1:length(NatWpts.(PathTypes{t})) 
    % Loop thru peripheral border waypoints
    CurPts = NatWpts.(PathTypes{t}){m}; % Current peripheral border
    for w = 1:length(CurPts) 
        cvx = NatVx(CurPts(w),:); % Current end point coordinates
        svx = NatVx(sp,:); % Surrounding path coordinates
        ds = sqrt((svx(:,1)-cvx(1)).^2 + (svx(:,2)-cvx(2)).^2 + (svx(:,3)-cvx(3)).^2); % Euclidian distance from end point to surrounding path vertices 
        v = sp(find(ds == min(ds),1)); % Nearest vertex on the surrounding path
        CurPts(w) = v; % Store new vertex
    end
    NatWpts.(PathTypes{t}){m} = CurPts; % Store new peripheral border 
end

%% Fit each waypoint to map
disp('Fitting waypoints to empirical map...');
% Loop thru iterations
FitWpts = NatWpts; % Fitted waypoints
for i = 1:Niter    
    % Loop thru path types
    Dots = []; % Path vertices 
    cm = 1; % Count paths
    for t = 1:length(PathTypes)        
        % Adjust search radius
        Radius = floor((Niter-i) / (Niter/InitRad)) + 1;

        % Loop thru paths of this type
        for m = 1:length(NatWpts.(PathTypes{t}))            
            % Loop thru waypoints in each path
            for w = 1:length(NatWpts.(PathTypes{t}){m})
                % Current waypoint 
                v = NatWpts.(PathTypes{t}){m}(w); % Vertex index
                % Only good points
                if ~isnan(v)
                    % Create search region
                    ecc = EccWpts.(PathTypes{t}){m}(w); % Original eccentricity 
                    roi = samsrf_georoi(v, Radius, Srf.Vertices, Srf.Faces); % Search region around waypoint
                    roi = roi(Rho(roi) > ecc-EccBw/2 & Rho(roi) < ecc+EccBw/2); % Limit to eccentricity range

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
                        % Current pRF positions
                        T = Theta(roi); % Polar angle
                        R = Rho(roi); % Eccentricity
                        % Locate nearest fitting vertex
                        if strcmpi(PathTypes{t}, 'UpperVerticalMeridian') 
                            mv = roi(find(T == nanmax(T),1)); % Maximal polar angle          
                        elseif strcmpi(PathTypes{t}, 'LowerVerticalMeridian') 
                            mv = roi(find(T == nanmin(T),1)); % Minimal polar angle          
                        elseif strcmpi(PathTypes{t}, 'HorizontalMeridian') 
                            T = abs(T); % Ignoring upper vs lower hemifield
                            mv = roi(find(T == nanmin(T),1)); % Minimal absolute polar angle 
                        else % Eccentricity based borders
                            R = abs(R-ecc); % Distance from desired eccentricity
                            mv = roi(find(R == nanmin(R),1)); % Minimal absolute polar angle 
                        end

                        % Store fitted vertex
                        if ~isempty(mv)
                            FitWpts.(PathTypes{t}){m}(w) = mv; % Point was updated
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
    samsrf_surf(Srf, 'Sphere', [R2Thresh 0 0 MinEcc MaxEcc -.1], {samsrf_borderpath(Srf, Dots); [NaN 1 1 1]}, '', 'Polar', h);
    set(gcf, 'name', ['Iteration #' num2str(i)]);
    pause(.5);
end
new_line;

%% Create paths
Paths = {}; % List of sub-paths (for delineation file)
EccBorder = []; % Contiguous eccentricity border paths

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
                nm = ConnectWaypoints(CurPts, Srf.Sphere); % Newly connected path
                Paths{cw} = nm; % Add to path list
                % Add to contiguous eccentricity border?
                if strcmpi(PathTypes{t}, 'FovealBorder') || strcmpi(PathTypes{t}, 'PeripheralBorder')
                    EccBorder = [EccBorder; nm]; 
                end
                cw = cw + 1; % Increase counter
            end
        end
    end
end

% Close foveal borders
disp('Closing foveal borders...');
t = find(strcmpi(PathTypes, 'FovealBorder')); % Select foveal borders
% Loop thru peripheral border paths 
for m = 1:length(NatWpts.(PathTypes{t})) 
    % Loop thru peripheral border waypoints
    CurPts = NatWpts.(PathTypes{t}){m}; % Current foveal border
    % Can only do this if circle
    if length(CurPts) > 2
        nm = ConnectWaypoints(CurPts([1 end]), Srf.Sphere); % Newly connected path
        Paths{cw} = nm; % Add to path list
        EccBorder = [EccBorder; nm]; % Add to contiguous path
    end
end

% Connect end points to eccentricity borders
disp('Connecting end points with eccentricity borders...');
% Loop thru path types
for t = 1:length(PathTypes)
    % Loop thru paths of this type
    for m = 1:length(NatWpts.(PathTypes{t}))
        CurPts = NatWpts.(PathTypes{t}){m}; % Waypoints in this path
        CurPts(isnan(CurPts)) = []; % Remove bad points
        if ~isempty(CurPts)
            % Loop thru end points
            for j = 1:2
                if j == 1 
                    e = CurPts(1); % Posterior end point
                else 
                    e = CurPts(end); % Anterior end point
                    % Don't connect anterior end points beyond V3
                    if m > 2 
                        e = NaN;
                    end
                end
                if ~isnan(e)
                    cvx = NatVx(e,:); % Current end point coordinates
                    svx = NatVx(EccBorder,:); % Eccentricity path coordinates 
                    % Nearest path vertex to end point
                    ds = sqrt((svx(:,1)-cvx(1)).^2 + (svx(:,2)-cvx(2)).^2 + (svx(:,3)-cvx(3)).^2); % Euclidian distance 
                    v = EccBorder(find(ds == min(ds),1)); % Nearest vertex on the path
                    Paths{end+1} = ConnectWaypoints([e v], Srf.Sphere); % Connect path
                end
            end
        end
    end 
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
AutoDelinParams.R2_Threshold = R2Thresh;
AutoDelinParams.Min_Eccentricity = MinEcc;
AutoDelinParams.Max_Eccentricity = MaxEcc;
AutoDelinParams.Eccentricity_Bandwidth = EccBw;
AutoDelinParams.Number_of_Iterations = Niter;
AutoDelinParams.Initial_Search_Radius = InitRad;
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
samsrf_surf(Srf, 'Sphere', [R2Thresh 0 0 MinEcc MaxEcc -.1], DrawPaths, '', 'Polar', h);
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
