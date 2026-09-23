function FastDelineationTool(SrfName, TmpName)
%
% FastDelineationTool(SrfName, TmpName)
%
% Fast delineation tool for retinotopic maps 
%
% It opens up the map in SrfName which must contain a gradient angle map 
% & probably also been smoothed & denoised etc. (so e.g. you have run it 
% through samsrf_denoisemap & samsrf_surfcalcs or samsrf_fieldsign).  
%
% Using the template maps ROI annotations & the empirical gradient angles
% the algorithm determines likely locations of all the ROIs in the template.
% You then have the option to reassign any ROIs that look wrong to new
% vertices in the map, using various maps to guide you (e.g. Polar angle,
% Eccentricity, Gradient Angle, Field Sign). Once you're happy with this
% assignment, you run a flood fill that labels ROIs based on fieldsigns.
%
% ROI labels are then saved in the subfolder ROIs_FastDelin.
%
%   SrfName:    Map file containing the gradient angle map            
%   TmpName:    Template map file warped into native space 
%
% 23/09/2026 - Completed functional version (DSS)

%% Initialise
global Vs 
Vs = [];
Bvx = {};

MaxSize = 300; % Maximum ROI size to be filled
FillTol = 15; % Filling gap tolerance to keep filling
CritGrd = 0.1; % Criterion gradient angle for filling

% Load Srfs
Srf = load(SrfName);
Srf = samsrf_expand_srf(Srf.Srf);
Srf.Data(end+1,:) = 0; % Add data row for ROI labels
Srf.Values{end+1} = 'ROI';
% Initialise seeds from template map or load saved delineation file
% Load template map
Tmp = load(TmpName);
Tmp = samsrf_expand_srf(Tmp.Srf);
if ~isfield(Tmp, 'Roi_Names')
    samsrf_error('This template map does not contain ROI names!');
end
RoiVal = find(strcmpi(Tmp.Values, 'ROI'));
if sum(Tmp.Data(RoiVal,:) < 0) == 0
    samsrf_disp('WARNING: No negative field sign - this may not be a suitable template map...');
    samsrf_newline;
end
Rois = Tmp.Roi_Names; % ROI names in template

%% Bilateral Srf?
if strcmpi(Srf.Hemisphere, 'bi')
    AdjRoiN = length(Rois); % Adjust ROI numbers for right hemisphere by this 
    bi_Rois = Rois;
    Rois = {};
    for r = 1:length(bi_Rois)*2
        if r <= length(bi_Rois)
            Hemi = 'lh';
        else
            Hemi = 'rh';
        end
        Rois{r} = [Hemi '_' bi_Rois{mod(r-1,length(bi_Rois))+1}];
    end

    % Ensure labels placed properly
    LabSph = Srf.Sphere;
    % Shift to origin
    LabSph(1:Srf.Nvert_Lhem,1) = LabSph(1:Srf.Nvert_Lhem,1) + 105; 
    LabSph(Srf.Nvert_Lhem+1:end,1) = LabSph(Srf.Nvert_Lhem+1:end,1) - 105; 
    % ROI labels outside sphere radius
    LabSph = LabSph * 1.2;
    % Shift back 
    LabSph(1:Srf.Nvert_Lhem,1) = LabSph(1:Srf.Nvert_Lhem,1) - 105; 
    LabSph(Srf.Nvert_Lhem+1:end,1) = LabSph(Srf.Nvert_Lhem+1:end,1) + 105; 
else
    AdjRoiN = 0; % Not relevant for single hemisphere
    Srf.Nvert_Lhem = size(Srf.Vertices,1);
    % Label ROIs with hemisphere 
    for r = 1:length(Rois)
        Rois{r} = [Srf.Hemisphere '_' Rois{r}];
    end

    % ROI labels outside sphere radius
    LabSph = Srf.Sphere * 1.2; 
end
FsVal = find(strcmpi(Srf.Values, 'Gradient Angle'),1); % Which row is field sign?
if isempty(FsVal)
    samsrf_error('This pRF map contains no gradient angle data!');
end

%% Extract ROIs from template
uR = unique(Tmp.Data(RoiVal,:)); % Unique ROI numbers (including sign)
[~,sa] = sort(abs(uR)); % Absolute ROI numbers in ascending sequence
uR = uR(sa); % Sort ROI numbers based on absolute value
% Loop thru hemispheres
for h = 1:2    
    % Loop thru ROIs
    for r = uR
        if r ~= 0
            % Template ROI label
            RoiVx = find(Tmp.Data(RoiVal,:) == r); % Vertices in this ROI
            % Which hemisphere?
            if h == 1
                RoiVx = RoiVx(RoiVx <= Srf.Nvert_Lhem); % Only left hemisphere vertices
            else
                RoiVx = RoiVx(RoiVx > Srf.Nvert_Lhem); % Only right hemisphere vertices
            end

            % Determine seed vertex for ROI
            RoiFs = sign(r) * Srf.Data(FsVal, RoiVx); % Field signs in this ROI
            sv = samsrf_cluster_midpt(RoiVx(RoiFs > 0), Srf.Sphere); % Midpoint of suprathreshold values in ROI
            
            % Store seeds for later
            Vs = [Vs; sv Srf.Data(FsVal, sv)]; % Seed vertex & field sign for this ROI
        end
    end
end

%% Show suggestions 
f = figure('WindowState', 'maximized', 'MenuBar', 'none', 'ToolBar', 'none');

% File menu
mmf = uimenu(f, 'Text', 'File');
mfl = uimenu(mmf, 'Text', 'Load ROI seeds');
mfl.MenuSelectedFcn = @LoadSeedsFcn;
mfs = uimenu(mmf, 'Text', 'Save ROI seeds');
mfs.MenuSelectedFcn = @SaveSeedsFcn;
% Map Menu
mm = uimenu(f, 'Text', 'Maps');
mmp = uimenu(mm, 'Text', 'Polar Angle');
mmp.MenuSelectedFcn = @ShowPolFcn;
mme = uimenu(mm, 'Text', 'Eccentricity');
mme.MenuSelectedFcn = @ShowEccFcn;
mms = uimenu(mm, 'Text', 'pRF size');
mms.MenuSelectedFcn = @ShowSigFcn;
mmg = uimenu(mm, 'Text', 'Gradient Angle');
mmg.MenuSelectedFcn = @ShowGrdFcn;
mmf = uimenu(mm, 'Text', 'Field Sign');
mmf.MenuSelectedFcn = @ShowFsFcn;
mmr = uimenu(mm, 'Text', 'ROI Labels');
mmr.MenuSelectedFcn = @ShowRoisFcn;
% mmr.Enable = 'off';
% ROIs menu
mr = uimenu(f, 'Text', 'ROIs');
mrf = uimenu(mr, 'Text', 'Fix ROI Position');
mrf.MenuSelectedFcn = @FixRoiFcn;
mra = uimenu(mr, 'Text', 'Add New Dummy');
mra.MenuSelectedFcn = @AddDummyFcn;
mrr = uimenu(mr, 'Text', 'Remove Last Dummy');
mrr.MenuSelectedFcn = @RemDummyFcn;
mrr.Enable = 'off';
% Rotation menu
mn = uimenu(f, 'Text', 'Rotation');
mn.MenuSelectedFcn = @RotModFcn;
% Delineation menu
md = uimenu(f, 'Text', 'Delineation');
mdm = uimenu(md, 'Text', ['Max ROI Size:       ' num2str(MaxSize)]);
mdm.MenuSelectedFcn = @SetMaxSizeFcn;
mdt = uimenu(md, 'Text', ['Gap Tolerance:      ' num2str(FillTol)]);
mdt.MenuSelectedFcn = @SetFillTolFcn;
mdc = uimenu(md, 'Text', ['Min Gradient Angle: ' num2str(CritGrd)]);
mdc.MenuSelectedFcn = @SetCritGrdFcn;
mdr = uimenu(md, 'Text', 'Fill ROIs');
mdr.MenuSelectedFcn = @FillRoisFcn;
mds = uimenu(md, 'Text', 'Smooth ROIs');
mds.MenuSelectedFcn = @SmoothRois;
mdl = uimenu(md, 'Text', 'Save ROI labels');
mdl.MenuSelectedFcn = @SaveRoisFcn;

% Set up ROI selection
datacursormode off
rotate3d

% Show polar map
p = samsrf_surf(Srf, 'Sphere', 0, '', '', 'Polar'); 
hold on
set(f, 'Name', 'Polar map');
% Label seeds
for i = 1:size(Vs,1)
    % Place label
    cv = Vs(i,1); % Current vertex
    quiver3(LabSph(cv,1), LabSph(cv,2), LabSph(cv,3), Srf.Sphere(cv,1)-LabSph(cv,1), Srf.Sphere(cv,2)-LabSph(cv,2), Srf.Sphere(cv,3)-LabSph(cv,3), 'off', 'Color', 'w', 'LineWidth', 2);
    text(LabSph(cv,1), LabSph(cv,2), LabSph(cv,3), Rois{i}, 'Interpreter', 'none', 'Backgroundcolor', 'w', 'EdgeColor', 'k');
    % Save seed vertex
    Vs(i,:) = [cv Srf.Data(FsVal,cv)];
end
% Stop if figure was closed
if ~isvalid(p)
    samsrf_disp('FastDelineationTool was closed by user.');
    return
end
 

%% % % % % % % % % % % % % %   Utility Functions   % % % % % % % % % % % % % % % % %%

    %% Load ROI seeds 
    function LoadSeedsFcn(~,~)
        delete(findobj(f, 'Type', 'DataTip')); % Remove data tips added
        title(''); % Remove any current title
        [Fn, Pn] = uigetfile('*.json', 'Load previous ROI seeds');
        if Fn 
            % Load parameters
            Pars = readstruct([Pn '/' Fn]);
            if ~isfield(Pars, 'ROI_Names')
                errordlg('Invalid Seed parameter file!', 'Error');
                return
            end
            if sum(~contains(Pars.ROI_Names, 'Dummy #')) ~= length(Rois) 
                errordlg('Seed parameters do not match number of ROIs in the template map!', 'Error');
                return
            end
            Vs = [Pars.ROI_Seeds' Pars.Gradient']; % ROI seed vertices assigned with gradient angles
            Rois = Pars.ROI_Names; % Load new ROI names 
            % Add dummy labels?
            for i = find(contains(Rois, 'Dummy #'))
                cv = Vs(i,1); % Current vertex
                quiver3(LabSph(cv,1), LabSph(cv,2), LabSph(cv,3), Srf.Sphere(cv,1)-LabSph(cv,1), Srf.Sphere(cv,2)-LabSph(cv,2), Srf.Sphere(cv,3)-LabSph(cv,3), 'off', 'Color', 'w', 'LineWidth', 2);
                text(LabSph(cv,1), LabSph(cv,2), LabSph(cv,3), Rois{i}, 'Interpreter', 'none', 'Backgroundcolor', 'w', 'EdgeColor', 'k', 'Color', 'r');
            end
            MaxSize = Pars.MaxSize; % Maximum ROI size
            FillTol = Pars.FillTol; % Filling gap tolerance 
            CritGrd = Pars.CritGrd; % Filling gap tolerance 
            % Update menu
            mdm.Text = ['Max ROI Size:       ' num2str(MaxSize)];
            mdt.Text = ['Gap Tolerance:      ' num2str(FillTol)];
            mdc.Text = ['Min Gradient Angle: ' num2str(CritGrd)];
            % Reassign labels
            t = findall(f, 'Type', 'Text');
            t = flipud(t(1:end-2)); % Invert order & remove two superfluous texts
            q = findall(f, 'Type', 'Quiver');
            q = flipud(q); % Invert order 
            for i = 1:length(q)
                v = Vs(i,1); % Current seed vertex
                t(i).Position = LabSph(v,:);
                q(i).XData = LabSph(v,1);
                q(i).YData = LabSph(v,2);
                q(i).ZData = LabSph(v,3);
                q(i).UData = Srf.Sphere(v,1) - LabSph(v,1);
                q(i).VData = Srf.Sphere(v,2) - LabSph(v,2);
                q(i).WData = Srf.Sphere(v,3) - LabSph(v,3);
                % Hide label if no field sign
                if Vs(i,2) == 0
                    t(i).Visible = 'off';
                    q(i).Visible = 'off';
                else
                    t(i).Visible = 'on';
                    q(i).Visible = 'on';
                end
            end
            samsrf_disp(['Loaded ROI seeds in: ' Pn '/' Fn]);
        end
    end

    %% Save ROI seeds 
    function Fn = SaveSeedsFcn(~,~)
        delete(findobj(f, 'Type', 'DataTip')); % Remove data tips added
        [~,Fn] = fileparts(SrfName);
        [Fn, Pn] = uiputfile([Fn '_FastDel.json'], 'Save ROI seeds');
        if Fn 
            % Save parameters
            Pars = struct;
            Pars.ROI_Names = Rois; % ROI names to compare against template
            Pars.ROI_Seeds = Vs(:,1); % ROI seed vertices assigned
            Pars.Gradient = Vs(:,2); % Gradient angles at seeds
            Pars.MaxSize = MaxSize; % Maximum ROI size
            Pars.FillTol = FillTol; % Filling gap tolerance 
            Pars.CritGrd = CritGrd; % Filling gap tolerance 
            writestruct(Pars, [Pn '/' Fn]);
            samsrf_disp(['Saved ROI seeds in: ' Pn '/' Fn]);
        end
    end

    %% Show polar map 
    function ShowPolFcn(~,~)
        samsrf_surf(Srf, 'Sphere', 0, Bvx, '', 'Polar', p);
        delete(findobj(f, 'Type', 'DataTip')); % Remove data tips added
    end

    %% Show eccentricity map 
    function ShowEccFcn(~,~)
        Ecc = sqrt(sum(Srf.Data(2:3,:).^2));
        samsrf_surf(Srf, 'Sphere', [0 0 max(Ecc)], Bvx, '', 'Eccentricity', p);
        delete(findobj(f, 'Type', 'DataTip')); % Remove data tips added
    end

    %% Show pRF size map 
    function ShowSigFcn(~,~)
        vs = find(strcmpi(Srf.Values, 'Sigma') | strcmpi(Srf.Values, 'Centre') | strcmpi(Srf.Values, 'Sigma1'));
        samsrf_surf(Srf, 'Sphere', [0 0 max(Srf.Data(vs,:))], Bvx, '', 'Eccentricity', p);
        delete(findobj(f, 'Type', 'DataTip')); % Remove data tips added
    end

    %% Show gradient angle map 
    function ShowGrdFcn(~,~)
        samsrf_surf(Srf, 'Sphere', [0 0 2], Bvx, '', 'Gradient Angle', p);
        delete(findobj(f, 'Type', 'DataTip')); % Remove data tips added
    end

    %% Show field sign map 
    function ShowFsFcn(~,~)
        samsrf_surf(Srf, 'Sphere', [0 0 2], Bvx, '', sign(Srf.Data(FsVal,:)), p);
        set(f, 'Name', 'Field sign');
        delete(findobj(f, 'Type', 'DataTip')); % Remove data tips added
    end

    %% Fix ROI position 
    function FixRoiFcn(~,~)
        global sr
        sr = listdlg('PromptString', 'Select ROI or Cancel to Finish', 'ListString', Rois, 'SelectionMode', 'single');
        if isempty(sr) 
            title('');
            datacursormode off
            delete(findobj(f, 'Type', 'DataTip')); % Remove data tips added
            rotate3d % Rotation back on
        else
            % Select new ROI seed
            title({'Select:'; Rois{sr}}, 'Interpreter', 'none');
            dcm_obj = datacursormode(f); % Data cursor mode object
            set(dcm_obj, 'UpdateFcn', @gvf); % Change selection function
            datacursormode on
        end
    end

    %% Add dummy ROI
    function AddDummyFcn(~,~)
        global sr
        Ndum = sum(contains(Rois, 'Dummy #'));
        Rois{end+1} = ['Dummy #' num2str(Ndum+1)];
        sr = length(Rois);
        % Select new ROI seed
        title({'Select:'; Rois{sr}}, 'Interpreter', 'none');
        dcm_obj = datacursormode(f); % Data cursor mode object
        set(dcm_obj, 'UpdateFcn', @gvf); % Change selection function
        datacursormode on
        mrr.Enable = 'on'; % Enable menu for removing dummies
    end
    
    %% Remove dummy ROI
    function RemDummyFcn(~,~)
        if contains(Rois{end}, 'Dummy #')
            Rois = Rois(1:end-1); % Remove dummy from ROI list
            Vs = Vs(1:end-1,:); % Remove from vertex list
            % Remove label
            t = findall(f, 'Type', 'Text');
            t = flipud(t(1:end-2)); % Invert order & remove two superfluous texts
            delete(t(end));
            q = findall(f, 'Type', 'Quiver');
            q = flipud(q); % Invert order 
            delete(q(end));
            % No more dummies left?
            if ~contains(Rois{end}, 'Dummy #')
                mrr.Enable = 'off'; % Disable menu for removing dummies
            end
        end                
    end

    %% Vertex selection 
    function txt = gvf(~,evnt)
        % Import vertex coordinates in mesh space 
        global Vertices sr     
        
        % Select new ROI seed
        pos = get(evnt, 'Position'); % Selected position in mesh space coordinates          
        v = find(Vertices(:,1) == pos(1) & Vertices(:,2) == pos(2) & Vertices(:,3) == pos(3), 1); % Selected vertex number
        Vs(sr,:) = [v Srf.Data(FsVal,v)]; % Collect seed vertices 
        txt = ''; % No data tip info
        
        % Report & label now ROI seed
        if contains(Rois{sr}, 'Dummy #')
            q = findall(f, 'Type', 'Quiver');
            if length(q) < length(Rois)
                quiver3(LabSph(v,1), LabSph(v,2), LabSph(v,3), Srf.Sphere(v,1)-LabSph(v,1), Srf.Sphere(v,2)-LabSph(v,2), Srf.Sphere(v,3)-LabSph(v,3), 'off', 'Color', 'w', 'LineWidth', 2);
                text(LabSph(v,1), LabSph(v,2), LabSph(v,3), Rois{sr}, 'Interpreter', 'none', 'Backgroundcolor', 'w', 'EdgeColor', 'k');
                samsrf_disp(['Assigned dummy ROI to #' num2str(v) ': Gradient Angle = ' num2str(Srf.Data(FsVal,v))]);
            else
                samsrf_disp(['Reassigned dummy ROI to #' num2str(v) ': Gradient Angle = ' num2str(Srf.Data(FsVal,v))]);
            end
        else
            samsrf_disp(['Reassigned ' Rois{sr} ' to #' num2str(v) ': Gradient Angle = ' num2str(Srf.Data(FsVal,v))]);
        end
        t = findall(f, 'Type', 'Text');
        t = flipud(t(1:end-2)); % Invert order & remove two superfluous texts
        q = findall(f, 'Type', 'Quiver');
        q = flipud(q); % Invert order 
        t(sr).Position = LabSph(v,:);
        q(sr).XData = LabSph(v,1);
        q(sr).YData = LabSph(v,2);
        q(sr).ZData = LabSph(v,3);
        q(sr).UData = Srf.Sphere(v,1) - LabSph(v,1);
        q(sr).VData = Srf.Sphere(v,2) - LabSph(v,2);
        q(sr).WData = Srf.Sphere(v,3) - LabSph(v,3);
        if contains(Rois{sr}, 'Dummy #')
            t(sr).Color = [1 0 0]; % Dummy ROIs are red
        end

        % Hide label if no field sign
        if Srf.Data(FsVal,v) == 0
            t(sr).Visible = 'off';
            q(sr).Visible = 'off';
        else
            t(sr).Visible = 'on';
            q(sr).Visible = 'on';
        end
    end

    %% Rotation mode 
    function RotModFcn(~,~)
        title('');
        datacursormode off
        delete(findobj(f, 'Type', 'DataTip')); % Remove data tips added
        rotate3d % Rotation back on
    end

    %% Set maximum ROI size 
    function SetMaxSizeFcn(~,~)
        Q = inputdlg('Maximum ROI size?', 'How many vertices?'); 
        Q = cell2mat(Q);
        if ~isempty(Q) 
            Q = str2double(Q);
            if ~isnan(Q)
                MaxSize = Q;
                mdm.Text = ['Max ROI Size:  ' num2str(MaxSize)];
            end
        end
    end

    %% Set fill gap tolerance 
    function SetFillTolFcn(~,~)
        Q = inputdlg('Minimum gap tolerance?', 'How many vertices?'); 
        Q = cell2mat(Q);
        if ~isempty(Q) 
            Q = str2double(Q);
            if ~isnan(Q)
                FillTol = Q;
                mdt.Text = ['Gap Tolerance:  ' num2str(FillTol)];
            end
        end
    end

    %% Set criterion gradient angle
    function SetCritGrdFcn(~,~)
        Q = inputdlg('Minimum gradient angle?', 'Filling criterion?'); 
        Q = cell2mat(Q);
        if ~isempty(Q) 
            Q = str2double(Q);
            if ~isnan(Q)
                CritGrd = Q;
                mdc.Text = ['Min Gradient Angle:  ' num2str(CritGrd)];
            end
        end
    end

    %% Fill ROIs 
    function FillRoisFcn(~,~)
        delete(findobj(f, 'Type', 'DataTip')); % Remove data tips added

        Rvx = {}; % Vertices per ROI
        LabMap = zeros(1,size(Srf.Vertices,1)); % Labelled vertices
        for r = 1:size(Vs,1)
            if Vs(r,1) > 0
                Rvx{r} = Vs(r,1); % Seed vertex for this region
                if contains(Rois{r}, 'Dummy #')
                    LabMap(Rvx{r}) = Inf; % Label as dummy
                else
                    LabMap(Rvx{r}) = sign(Vs(r,2)) * (r - (Vs(r,1) > Srf.Nvert_Lhem)*AdjRoiN);
                end
            else
                Rvx{r} = NaN;
            end    
        end
        
        % Keep growing regions 
        Bvx = {}; % Borders of ROIs
        for i = 1:30  
            % Stop if figure was closed
            if ~isvalid(p)
                samsrf_disp('FastDelineationTool was closed by user.');
                return
            end
            % Grow around labelled regions
            for r = 1:size(Vs,1)
                cB = samsrf_neighbours(Rvx{r}, Srf.Faces); % Neighbours around current ROI
                cB(sign(Vs(r,2))*Srf.Data(FsVal,cB) <= 0.1) = []; % Remove vertices below threshold 
                cB(ismember(cB, find(LabMap))) = []; % Remove already labelled
                % Only keep filling if new border sufficiently large
                if length(Rvx{r}) <= MaxSize || length(cB) > FillTol
                    Rvx{r} = [Rvx{r}; cB]; % Add vertices to ROI
                    Bvx{r} = samsrf_borderpath(Srf, Rvx{r}); % New border 
                    if contains(Rois{r}, 'Dummy #')
                        LabMap(cB) = Inf; % Label as dummy
                    else
                        LabMap(cB) = sign(Vs(r,2)) * (r - (Vs(r,1) > Srf.Nvert_Lhem)*AdjRoiN); % Label remaining border
                    end
                end
            end
            samsrf_surf(Srf, 'Sphere', 0, Bvx, '', 'Polar', p); % Show polar map 
            set(f, 'Name', 'Filling labels');
            pause(.01); % Short break for visualisation
        end
        % Remove dummy labels       
        LabMap(isinf(LabMap)) = 0;
        % Add ROIs to Srf
        Srf.Data(end,:) = LabMap;
        % Now show ROIs
        ShowRoisFcn; 
    end

    %% Show ROI labels
    function ShowRoisFcn(~,~)
        delete(findobj(f, 'Type', 'DataTip')); % Remove data tips added
        dSrf = Srf;
        dSrf.Data(1,Srf.Data(end,:)== 0) = 0;
        samsrf_surf(dSrf, 'Sphere', [.001 0 max(abs(Srf.Data(end,:)))], '', '', 'ROI', p); % Show ROI map
    end

    %% Smooth ROI map
    function SmoothRois(~,~)
        delete(findobj(f, 'Type', 'DataTip')); % Remove data tips added
        title('');
        if var(Srf.Data(end,:)) > 0
            mver = find(Srf.Data(1,:) > 0);
            NewLabMap = zeros(1,length(mver));
            title('Smoothing ROIs...');
            wb = waitbar(0, 'Smoothing ROIs...');
            for i = 1:length(mver)
                v = mver(i);
                sl = samsrf_sphere_roi(v, 5, Srf.Sphere);
                NewLabMap(i) = mode(Srf.Data(end,sl));
                waitbar(i/length(mver), wb);
            end
            Srf.Data(end,mver) = NewLabMap;
            title('');
            close(wb);
            % Now show ROIs 
            ShowRoisFcn; 
        else
            errordlg('No ROI have been labelled yet!', 'Error');
        end
    end

    %% Save ROI labels
    function SaveRoisFcn(~,~)
        delete(findobj(f, 'Type', 'DataTip')); % Remove data tips added
        Fn = SaveSeedsFcn; % Save seeds first
        [~,Fn] = fileparts(Fn); % Saved seeds filename
        mkdir(['ROIs_' Fn]); % Make folder for ROIs
        cd(['ROIs_' Fn]); % Make folder for ROIs
        title('Saving ROI labels...');
        for r = 1:length(Rois)
            if ~contains(Rois{r}, 'Dummy #') && Vs(r,2) ~= 0 % Only save ROIs not empty or dummies
                samsrf_srf2label(Srf, Rois{r}, 1, find(abs(Srf.Data(end,:))==r));
            end
        end
        % Bilateral Srf?
        if strcmpi(Srf.Hemisphere, 'bi')
            for r = 1:length(Rois)
                if contains(Rois{r}, 'lh_')
                    samsrf_bilat_label(Srf, Rois{r}(4:end)); % Combine hemispheres
                end
            end
            % Now delete hemisphere labels
            delete('lh_*.label');
            delete('rh_*.label');
        end
        title('');
        cd ..
    end

end
