function AutoDelFile = AutoDelineation(SrfName, NatMesh, TmpMesh, Atlas, R2Thresh, MinEcc, MaxEcc, EccBw, Niter, InitRad)
%
% AutoDelFile = AutoDelineation(Srf, NatMesh, TmpMesh, [Atlas, R2Thresh, MinEcc, MaxEcc, EccBw, Niter, InitRad])
%
% Simple algorithm to fits the borders of visual regions. It shows a movie 
% of the search process so you can see how well it's going. Then it saves a 
% delineation file autodel_SrfName that you can load in DelineationTool. 
% 
% THIS IS NOT A MAGIC WANT! You -must- inspect the maps in DealineationTool 
% & probably manually correct some smaller errors. You still need to label
% the regions manually (there may be another tool for that though...).
% Also, this algorithm expects gentle retinotopic map gradients, so you 
% may need to smooth your maps beforehand (possibly quite heavily) & might
% want to apply other filtering (see also samsrf_surfcalcs for example). 
% However, too much smoothing can also go awry. TLDR: Check your results!
%
% The automatic delineation is based on the visual field atlas normalised 
% to the fsaverage template. These initial borders are stored in an atlas 
% file called AutoDelinAtlas_[Atlas].mat in SamSrf/Utils. These files 
% contain waypoints for cardinal meridians & some iso-eccentricity lines. 
%
% The algorithm works as follows:
%   1.  The waypoints comprising the normalised borders are warped back
%       into native brain space (so you need the fsaverage template).
%   2.  Your retinotopic map is restricted to a particular range, using a
%       goodness-of-fit threshold & minimum & maximum eccentricity.
%   3.  A path surrounding the whole thresholded region is also created to
%       estimate the initial location of the peripheral extent of the map.
%   4.  The algorithm then searches the geodesic neighbourhood of each
%       waypoint for the best pRF (depends on which border it is fitting). 
%       This search is constrained by being limited to the eccentriity band 
%       of its original position (the width of this band is adjustable). 
%       It is also constrained by expecting 0.5 mm minimum Euclidean 
%       distance in sphere space between all waypoints. 
%   5.  The search process is repeated for several iterations (adjustable).
%       The size of the search region is gradually reduced from InitRad 
%       geodesic steps down to 1 step over the course of Niter iterations.
%   6.  Upon completion of the search process, the waypoints in each path 
%       are connected into complete paths. 
%   7.  The end points of the meridian lines are then connected to the
%       foveal & peripheral borders as appropriate.
%   8.  Finally, the paths are expanded to be compatible with region
%       filling in the DelineationTool & the delineation file is saved.
%
% The default thresholds & number of iterations was defined by the atlas
% you are using. You may need to tweak this for your own needs. It should 
% be somewhat robust to relatively small differences in designs but it 
% obviously cannot do a great job if you used a much larger field of view
% than the atlas allows for.
%
% The parameters used are saved in the delineation file & the function
% returns the filename. Even though you may need to correct & adjust the 
% automatic delineation, this makes your delineation much more reproducible. 
% The auto-delineation is saved inside the delineation file for posterity 
% so you can compare this to the final version.
% 
%   SrfName:    Srf struct with a pRF map or
%                name of Srf file with retinotopic map 
%   NatMesh:    Subject's surf folder which must contain lh/rh.sphere.reg
%   TmpMesh:    Template's surf folder which must contain lh/rh.sphere
%
% Optional additional input arguments:
%   Atlas:      Atlas waypoints & parameters (default = 'InfernoSerpents')
%   R2Thresh:   R^2 threshold (default = 0.1)
%   MinEcc:     Central eccentricity limit (default = 0.5)
%   MaxEcc:     Peripheral eccentricity limit (default = 9)
%   EccBw:      Width of eccentricity band (default = 0.5)
%   Niter:      Number of search iterations (default = 20)
%   InitRad:    Initial search radius (default = 2)
%
% 19/09/2024 - Now also takes Srf struct as input map (DSS)
%              Returns file name of autodelineation (DSS)
% 07/10/2025 - Fixed bugs that shouldn't be happening anyway (DSS)
%

global SamSrfXPath

% ROI list for delineation file 
SamSrfDefs = LoadSamSrfDefaults;
if ~exist('SamSrfDefs.def_roilist')
    % ROI list if undefined in SamSrf_defaults
    RoiList = {'V1' 'V2v' 'V3v' 'V4' 'V2d' 'V3d' 'V3A' 'V3B' 'LO1' 'LO2' 'VO1' 'VO2' 'TO1' 'TO2' 'V6' 'IPS0' 'IPS1' 'IPS2'}'; % For backwards compatability
else
    % Pre-defined default ROI list
    RoiList = SamSrfDefs.def_roilist;
end

%% Parameters
if nargin < 4
    Atlas = 'InfernoSerpents';
end
if isempty(SamSrfXPath)
    % Matlab command window
    AtlasFile = ['AutoDelinAtlas_' Atlas '.mat'];
else
    % SamSrf X GUI
    AtlasFile = [SamSrfXPath filesep 'Utils' filesep 'AutoDelinAtlas_' Atlas '.mat'];
end
if ~exist(AtlasFile, 'file')
    samsrf_error(['Auto-delineation atlas ' Atlas ' does not exist!']);
end
AtlasData = load(AtlasFile);
% Assign default parameters from atlas file
RoiSeeds = NaN(length(RoiList),1); % Initialise list of ROI seed vertices
if nargin < 5 
    R2Thresh = AtlasData.R2Thresh; % Goodness-of-fit threshold
end
if nargin < 6 
    MinEcc = AtlasData.MinEcc; % Minimal eccentricity
end
if nargin < 7 
    MaxEcc = AtlasData.MaxEcc; % Maximal eccentricity
end
if nargin < 8 
    EccBw = AtlasData.EccBw; % Eccentricity band width
end
if nargin < 9
    Niter = AtlasData.Niter; % Number of search iterations
end
if nargin < 10
    InitRad = AtlasData.InitRad; % Initial search radius
end
samsrf_newline;

%% Welcome message
[vn,vd] = samsrf_version;
samsrf_clrscr; 
samsrf_disp('****************************************************************************');
samsrf_disp('   Welcome to the Seriously Annoying MatLab Surfer Auto-Delineation Tool!');
samsrf_disp('    by D. S. Schwarzkopf from the University of Auckland, New Zealand');
samsrf_newline;
samsrf_disp(['                 Version ' num2str(vn) ', Released on ' vd]);
samsrf_disp('      (see SamSrf/ReadMe.md for what is new in this version)');
samsrf_disp('****************************************************************************');
samsrf_newline;

%% Load map
samsrf_disp('Loading map...');
if isstruct(SrfName)
    % Map provided as Srf
    Srf = SrfName;
    SrfName = 'SamSrfX';
else
    % Map file defined
    load(SrfName, 'Srf'); % Load Srf
    SrfName = SrfName(4:end); % Remove prefix
end
Srf = samsrf_expand_srf(Srf); % Expand Srf
if ~isfield(Srf, 'Values') % Are there values?
    samsrf_error('This file does not contain a pRF map!');
end
% Which hemisphere(s)? 
if upper(Srf.Hemisphere(1)) == 'B' 
    % Bilateral Srf 
    Hemis = {'l' 'r'};
elseif upper(Srf.Hemisphere(1)) == 'L' || upper(Srf.Hemisphere(1)) == 'R'
    % Left or Right hemisphere Srf
    Hemis = {lower(Srf.Hemisphere(1))};
else
    % Not a valid brain surface file
    samsrf_error('This tool only works for valid surface data files!');
end
% If bilateral Srf
if length(Hemis) > 1
    samsrf_disp(' Splitting bilateral Srf into hemispheres...');
    [lh_Srf, rh_Srf] = samsrf_hemi_srfs(Srf); % Split into hemispheres
    samsrf_newline;
end

%% Report parameters
Rois = '';
for i = 1:length(AtlasData.Rois)
    Rois = [Rois ' ' AtlasData.Rois{i}];
end
samsrf_disp(['Atlas:                   ' Atlas ':' Rois]);
samsrf_newline;
samsrf_disp('Parameters:')
samsrf_disp([' R^2 threshold:          ' num2str(R2Thresh)]);
samsrf_disp([' Minimum eccentricity:   ' num2str(MinEcc)]);
samsrf_disp([' Maximum eccentricity:   ' num2str(MaxEcc)]);
samsrf_disp([' Eccentricity bandwidth: ' num2str(EccBw)]);
samsrf_disp([' Number of iterations:   ' num2str(Niter)]);
samsrf_disp([' Initial search radius:  ' num2str(InitRad)]);
samsrf_newline;

%% Loop thru all hemispheres to analyse
for h = 1:length(Hemis)
    figure;

    %% Which hemisphere & map?
    samsrf_disp(['Delineating map: ' Hemis{h} 'h_' SrfName]);
    % If bilateral Srf
    if length(Hemis) > 1
        if h == 1
            Srf = lh_Srf; % Delineate left hemisphere
        elseif h == 2
            Srf = rh_Srf; % Delineate right hemisphere
        end
    end

    %% Viewing region
    if ~exist('SamSrfDefs.def_disproi')
        SamSrfDefs.def_disproi = NaN; 
    end
    if ~isempty(SamSrfDefs.def_disproi) && (SamSrfDefs.def_disproi(1) == '<' || SamSrfDefs.def_disproi(1) == '>')
        % If ROI defined by coordinates
        if length(SamSrfDefs.def_disproi) == 1
            samsrf_error('You must define inflated mesh coordinate in SamSrfDefs.def_disproi!');
        end
        switch SamSrfDefs.def_disproi(2)
            case 'X'
                wv = Srf.Inflated(:,1);
            case 'Y'
                wv = Srf.Inflated(:,2);
            case 'Z'
                wv = Srf.Inflated(:,3);
            otherwise
                samsrf_error('Invalid inflated mesh coordinate specified in SamSrfDefs.def_disproi!');
        end
        if length(SamSrfDefs.def_disproi) < 3
            samsrf_error('You must define inflated mesh cut-off coordinate in SamSrfDefs.def_disproi!');
        end
        wc = str2double(SamSrfDefs.def_disproi(3:end));
        if SamSrfDefs.def_disproi(1) == '<'
            ViewRoi = find(wv < wc);
        elseif SamSrfDefs.def_disproi(1) == '>'
            ViewRoi = find(wv > wc);
        end
        samsrf_disp(['Only displaying inflated mesh vertices with ' SamSrfDefs.def_disproi]);
    elseif isnan(SamSrfDefs.def_disproi)
        % Use ROI from Srf
        samsrf_disp('Using ROI from Srf if it exists');
        if isfield(Srf, 'Roi')
            ViewRoi = Srf.Roi;
        else
            ViewRoi = NaN;
        end
    else
        % Load region of interest
        Roi = [Hemis{h} 'h_' SamSrfDefs.def_disproi];
        ViewRoi = samsrf_loadlabel(Roi);
        samsrf_disp(['Only displaying ROI ' Roi]);
    end
    if ~isnan(ViewRoi)
        AllVx = true(size(Srf.Vertices,1),1);
        AllVx(ViewRoi) = false;
        Srf.Vertices(AllVx,:) = NaN;
    end
    samsrf_newline;

    %% Sphere vertices in native & template 
    NatVx = Srf.Sphere; % Target vertices
    RegVx = fs_read_surf([NatMesh filesep Hemis{h} 'h.sphere.reg']); % Registration vertices
    TmpVx = fs_read_surf([TmpMesh filesep Hemis{h} 'h.sphere']); % Template vertices
    NnatVx = size(NatVx,1); % Number of target vertices
    if size(NatVx,1) == size(RegVx,1)
        NatVx = RegVx;
    else
        samsrf_error('Number of registration vertices does not match native surface mesh!');
    end

    %% Load waypoints
    if upper(Hemis{h}) == 'L'
        Wpts = AtlasData.Lhem; % Left hemisphere waypoints
    elseif upper(Hemis{h}) == 'R'
        Wpts = AtlasData.Rhem; % Right hemisphere waypoints
    else
        samsrf_error('Auto-delineation tool currently only works for left & right surface meshes!');
    end
    % Does atlas contain flexible peripheral border?
    if isfield(Wpts, 'IsoEccentricityLines')
        iel = abs(AtlasData.Eccens-MaxEcc); % Difference between iso-eccentricity lines & desired maximal eccentricity
        Wpts.PeripheralBorder = {Wpts.IsoEccentricityLines{iel == min(iel)}}; % Set desired iso-eccentricity line as peripheral border
        Wpts = rmfield(Wpts, 'IsoEccentricityLines'); % Remove these lines from process
    end
    NatWpts = Wpts; % Waypoints after warping
    PathTypes = fieldnames(Wpts);

    %% Warp waypoints into native space
    samsrf_disp('Warping template waypoints into native space...');
    % Loop thru path types
    for t = 1:length(PathTypes)
        samsrf_disp([' ' PathTypes{t}]);
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
            samsrf_disp(['  ' num2str(m) '/' num2str(length(Wpts.(PathTypes{t}))) ': ' num2str(length(Wpts.(PathTypes{t}){m})) ' waypoints']);
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
    ph = samsrf_surf(Srf, 'Sphere', [R2Thresh 0 0 MinEcc MaxEcc -.1], {Dots; [NaN 1 1 1]}, '', 'Polar');
    set(gcf, 'name', 'Initial waypoints');
    pause(.5);
    samsrf_newline; 

    %% Adjust peripheral border
    samsrf_disp('Adjusting peripheral border...');
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
    samsrf_disp('Fitting waypoints to empirical map...');
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
        samsrf_surf(Srf, 'Sphere', [R2Thresh 0 0 MinEcc MaxEcc -.1], {samsrf_borderpath(Srf, Dots); [NaN 1 1 1]}, '', 'Polar', ph);
        set(gcf, 'name', ['Iteration #' num2str(i)]);
        pause(.5);
    end
    samsrf_newline;

    %% Create paths
    Paths = {}; % List of sub-paths (for delineation file)
    EccBorder = []; % Contiguous eccentricity border paths

    % Connecting waypoints
    samsrf_disp('Creating paths from waypoints...');
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
    samsrf_disp('Closing foveal borders...');
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
    samsrf_disp('Connecting end points with eccentricity borders...');
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
    AutoDelin = struct;
    AutoDelin.Atlas = Atlas; % Atlas used for auto-delineation
    AutoDelin.Version = samsrf_version; % SamSrf version of this atlas
    AutoDelin.R2_Threshold = R2Thresh; % Goodness-of-fit threshold
    AutoDelin.Min_Eccentricity = MinEcc; % Minimal eccentricity
    AutoDelin.Max_Eccentricity = MaxEcc; % Maximal eccentricity
    AutoDelin.Eccentricity_Bandwidth = EccBw; % Eccentricity band width
    AutoDelin.Number_of_Iterations = Niter; % Number of iterations
    AutoDelin.Initial_Search_Radius = InitRad; % Initial search radius
    AutoDelin.Paths = DrawPaths(1:end-1); % Save autodelineated paths without expansion
    % Vector with all path vertices 
    Vs = [];
    for i = 1:length(Paths)
        Vs = [Vs; Paths{i}];
    end
    % Saves everything to disc
    AutoDelFile = ['autodel_' Hemis{h} 'h_' SrfName '.mat'];
    samsrf_newline;
    save(AutoDelFile, 'SrfName', 'Vs', 'Paths', 'RoiList', 'RoiSeeds', 'AutoDelin');
    samsrf_disp(['Saved auto-delineation ' AutoDelFile]);

    %% Display final delineation
    samsrf_surf(Srf, 'Sphere', [R2Thresh 0 0 MinEcc MaxEcc -.1], DrawPaths, '', 'Polar', ph);
    set(gcf, 'Name', 'Inspect auto delineation', 'CloseRequestFcn', 'closereq');
    clear global Type Data CurvGrey fh fv pv % Clearing all but Vertices to avoid failure in DelineationTool
    rotate3d;
    samsrf_disp('Don''t forget to inspect & correct in DelineationTool before labelling the ROIs!');
    samsrf_newline;
end

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
