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
% Eccentricity, Sigma, Gradient Angle, Field Sign). Once you're happy with 
% this assignment, you run a flood fill that labels ROIs based on fieldsigns.
%
% ROI labels are then saved in the subfolder ROIs_FastDelin.
%
%   SrfName:    Map file containing the gradient angle map            
%   TmpName:    Template map file warped into native space 
%
% 23/09/2026 - Completed functional version (DSS)
%              Bugfix for when no clusters can be found in a template ROI (DSS)
% 24/09/2026 - Opens window after successfully loading maps now (DSS)
%              ROI lables pointing at 0 field sign are now greyed out (DSS)
%              Removed waitbar from smoothing & report progess in window (DSS)
%              Fixed bug where bilateral labels were saved twice (DSS)
%              Now checks you only select labels in correct hemisphere (DSS) 
% 25/09/2026 - Checks now if field sign correct when selecting ROI (DSS)
%              ROI labels now contain field sign as well (DSS) 
% 26/09/2026 - ROI suggestions without good fieldsign now assigned zero (DSS)
% 28/09/2026 - If no inputs defined, opens GUI to select files (DSS)
% 29/09/2026 - Fixed bug with omitting ROIs in unilateral Srf (DSS)
%              Ensure info is visible in unilateral Srf (DSS)
% 12/10/2026 - Inconsequential cosmetic changes (DSS)
% 14/10/2026 - Added option to create duplicate ROIs (DSS)
%              You can now fill in specific ROIs to plug holes (DSS)
% 15/10/2026 - Plugging gaps now runs for one round only (DSS)
%              Now saves log of post-processing steps alongside ROIs (DSS)
% 16/10/2026 - You can now add new ROIs (DSS)
% 19/10/2026 - Seems to be working as intended - Famous last words... (DSS)

%% Initialise
global Vs 
Vs = [];
Bvx = {};

% Initial constants
MaxSize = 300; % Maximum ROI size to be filled
FillTol = 15; % Filling gap tolerance to keep filling
CritGrd = 0.1; % Criterion gradient angle for filling
PostProcLog = {}; % Log of post-processing steps like smoothing & hole plugging

%% Load pRF map 
if nargin == 0
    [fn, pn] = uigetfile('*.mat', 'Select empirical map');
    if fn ~= 0
        SrfName = [pn '/' fn];
    end
end
try
    Srf = load(SrfName);
catch
    samsrf_error(['Could not load map file ' SrfName]);
end
Srf = samsrf_expand_srf(Srf.Srf);
Srf.Data(end+1,:) = 0; % Add data row for ROI labels
Srf.Values{end+1} = 'ROI';

%% Load template map
if nargin <= 1
    [fn, pn] = uigetfile('*.mat', 'Select template map');
    if fn ~= 0
        TmpName = [pn '/' fn];
    end
end
try
    Tmp = load(TmpName);
catch
    samsrf_error(['Could not load template map file ' TmpName]);
end
if ~strcmpi(Srf.Hemisphere, Tmp.Srf.Hemisphere)
    samsrf_error('Empirical & template map not the same surface!');
end
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

% Open window
samsrf_newline;
f = figure('WindowState', 'maximized', 'MenuBar', 'none', 'ToolBar', 'none');

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

%% Retrieve gradient angle map
FsVal = find(strcmpi(Srf.Values, 'Gradient Angle'),1); % Which row is field sign?
if isempty(FsVal)
    samsrf_error('This pRF map contains no gradient angle data!');
end

%% Extract ROIs from template
uR = unique(Tmp.Data(RoiVal,:)); % Unique ROI numbers (including sign)
[~,sa] = sort(abs(uR)); % Absolute ROI numbers in ascending sequence
uR = uR(sa); % Sort ROI numbers based on absolute value
% Loop thru hemispheres
for h = 1:1+strcmpi(Srf.Hemisphere,'bi') % Only 2 if bilateral Srf  
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
            
            % Add field sign to ROI name for labelling
            Rois{(h-1)*(length(uR)-1)+abs(r)} = [Rois{(h-1)*(length(uR)-1)+abs(r)} ' ' signstr(r)];

            % Determine seed vertex for ROI
            RoiFs = sign(r) * Srf.Data(FsVal, RoiVx); % Field signs in this ROI
            sv = samsrf_cluster_midpt(RoiVx(RoiFs > 0), Srf.Sphere); % Midpoint of suprathreshold values in ROI
            % If cannot find anything use middle of template ROI
            if isempty(sv)
                sv = samsrf_cluster_midpt(RoiVx, Srf.Sphere); % Midpoint of template ROI
				Vs = [Vs; sv 0]; % Bad seed vertex & field sign set to zero
			else
				Vs = [Vs; sv Srf.Data(FsVal, sv)]; % Seed vertex & field sign for this ROI
            end
        end
    end
end

%% Show ROI suggestions 
% File menu
mmf = uimenu(f, 'Text', 'File');
mfl = uimenu(mmf, 'Text', 'Load ROI Seeds');
mfl.MenuSelectedFcn = @LoadSeedsFcn;
mfs = uimenu(mmf, 'Text', 'Save ROI Seeds');
mfs.MenuSelectedFcn = @SaveSeedsFcn;
% Map Menu
mm = uimenu(f, 'Text', 'Maps');
mmp = uimenu(mm, 'Text', 'Polar Angle');
mmp.MenuSelectedFcn = @ShowPolFcn;
mme = uimenu(mm, 'Text', 'Eccentricity');
mme.MenuSelectedFcn = @ShowEccFcn;
mms = uimenu(mm, 'Text', 'pRF Size');
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
mra = uimenu(mr, 'Text', 'Add New ROI');
mra.MenuSelectedFcn = @AddNewFcn;
mra = uimenu(mr, 'Text', 'Add New Dummy');
mra.MenuSelectedFcn = @AddDummyFcn;
mrd = uimenu(mr, 'Text', 'Duplicate ROI');
mrd.MenuSelectedFcn = @DuplicateFcn;
mrr = uimenu(mr, 'Text', 'Remove Last Extra');
mrr.MenuSelectedFcn = @RemLastRoiFcn;
mrr.Enable = 'off';
mro = uimenu(mr, 'Text', 'Omit ROIs');
mro.MenuSelectedFcn = @OmitRoiFcn;
% Rotation menu
mn = uimenu(f, 'Text', 'Rotation');
mn.MenuSelectedFcn = @RotModFcn;
% Delineation menu
md = uimenu(f, 'Text', 'Delineation');
mdm = uimenu(md, 'Text', ['Max Cluster Size:   ' num2str(MaxSize)]);
mdm.MenuSelectedFcn = @SetMaxSizeFcn;
mdt = uimenu(md, 'Text', ['Gap Tolerance:      ' num2str(FillTol)]);
mdt.MenuSelectedFcn = @SetFillTolFcn;
mdc = uimenu(md, 'Text', ['Min Gradient Angle: ' num2str(CritGrd)]);
mdc.MenuSelectedFcn = @SetCritGrdFcn;
mdr = uimenu(md, 'Text', 'Fill ROIs');
mdr.MenuSelectedFcn = @FillRoisFcn;
mdp = uimenu(md, 'Text', 'Plug Small Gaps');
mdp.MenuSelectedFcn = @PlugGapFcn;
mds = uimenu(md, 'Text', 'Smooth ROIs');
mds.MenuSelectedFcn = @SmoothRoisFcn;
mdl = uimenu(md, 'Text', 'Save ROI Labels');
mdl.MenuSelectedFcn = @SaveRoisFcn;

% Set up ROI selection
datacursormode off
rotate3d

% Show polar map
p = samsrf_surf(Srf, 'Sphere', 0, '', '', 'Polar'); 
hold on
title('');
if ~strcmpi(Srf.Hemisphere, 'bi')
    set(gca, 'TitleHorizontalAlignment', 'left');
end
set(f, 'Name', 'Polar map');
% Label seeds
for i = 1:size(Vs,1)
    % Place label
    cv = Vs(i,1); % Current vertex
    cq = quiver3(LabSph(cv,1), LabSph(cv,2), LabSph(cv,3), Srf.Sphere(cv,1)-LabSph(cv,1), Srf.Sphere(cv,2)-LabSph(cv,2), Srf.Sphere(cv,3)-LabSph(cv,3), 'off', 'Color', 'w', 'LineWidth', 2);
    ct = text(LabSph(cv,1), LabSph(cv,2), LabSph(cv,3), Rois{i}, 'Interpreter', 'none', 'Backgroundcolor', 'w', 'EdgeColor', 'k');
    % Save seed vertex
    Vs(i,:) = [cv Srf.Data(FsVal,cv)];
    % Grey out seed with empty field sign
    if Vs(i,2) == 0
        ct.Color = [.5 .5 .5];
        cq.Color = [.5 .5 .5];
    end
end
% Stop if figure was closed
if ~isvalid(p)
    samsrf_disp('FastDelineationTool was closed by user.');
    return
end
 

%% % % % % % % % % % % % % %   Utility Functions   % % % % % % % % % % % % % % % % %%

    %% Sign string
    function s = signstr(n)
        if sign(n) == -1
            s = '-';
        elseif sign(n) == 1
            s = '+';
        else
            s = '';
        end
    end

    %% Sign from char
    function n = strsign(s)
        if s == '-'
            n = -1;
        elseif s == '+'
            n = 1;
        else
            n = NaN;
        end
    end

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
            Vs = [Pars.ROI_Seeds' Pars.Gradient']; % ROI seed vertices assigned with gradient angles
            Rois = Pars.ROI_Names; % Load new ROI names 
            % Add dummy labels?
            if sum(contains(Rois, 'Dummy #') | contains(Rois, '~') | contains(Rois, '^')) > 0
                mrr.Enable = 'on'; % Enable menu for removing dummies
                for i = find(contains(Rois, 'Dummy #') | contains(Rois, '~') | contains(Rois, '^'))
                    cv = Vs(i,1); % Current vertex
                    quiver3(LabSph(cv,1), LabSph(cv,2), LabSph(cv,3), Srf.Sphere(cv,1)-LabSph(cv,1), Srf.Sphere(cv,2)-LabSph(cv,2), Srf.Sphere(cv,3)-LabSph(cv,3), 'off', 'LineWidth', 2);
                    text(LabSph(cv,1), LabSph(cv,2), LabSph(cv,3), Rois{i}, 'Interpreter', 'none', 'Backgroundcolor', 'w', 'EdgeColor', 'k');
                end
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
                % Grey out label if no field sign
                if Vs(i,2) == 0
                    t(i).FontAngle = 'italic';
                    t(i).Color = [.5 .5 .5];
                    q(i).Color = [.5 .5 .5];
                else
                    t(i).FontAngle = 'normal';
                    if contains(Rois{i}, 'Dummy #')
                        t(i).Color = [1 0 0];
                    else
                        t(i).Color = [0 0 0];
                    end
                    q(i).Color = [1 1 1];
                end
            end
            title({'ROI seeds'; 'loaded'});
            samsrf_disp(['Loaded ROI seeds in: ' Pn '/' Fn]);
        end
        datacursormode off
        rotate3d
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
            title({'ROI seeds'; 'saved'});
            samsrf_disp(['Saved ROI seeds in: ' Pn '/' Fn]);
        end
        datacursormode off
        rotate3d
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

    %% Add new ROI
    function AddNewFcn(~,~)
        global sr
        NewRoi = inputdlg('New ROI Name?');
        NewRoi = char(NewRoi);
        if isempty(NewRoi)  
            return
        end
        if length(NewRoi) < 4
            errordlg('ROI name must contain hemisphere prefix (lh_ or rh_)!', 'Error');
            datacursormode off
            rotate3d
            return
        end
        NewRoi = [NewRoi '^ ?']; % Placeholder field sign suffix
        if strcmpi(NewRoi(1:3), 'lh_') || strcmpi(NewRoi(1:3), 'rh_')
            Rois{end+1} = NewRoi;
            sr = length(Rois);
            % Select new ROI seed
            title({'Select:'; Rois{sr}}, 'Interpreter', 'none');
            dcm_obj = datacursormode(f); % Data cursor mode object
            set(dcm_obj, 'UpdateFcn', @gvf); % Change selection function
            datacursormode on
            mrr.Enable = 'on'; % Enable menu for removing dummies
        else
            errordlg('ROI name must contain hemisphere prefix (lh_ or rh_)!', 'Error');
            datacursormode off
            rotate3d
        end
    end

    %% Duplicate ROI
    function DuplicateFcn(~,~)
        global sr
        sr = listdlg('PromptString', 'Select ROI or Cancel to Finish', 'ListString', Rois(~contains(Rois, 'Dummy #') & ~contains(Rois, '~')), 'SelectionMode', 'single');
        if isempty(sr) 
            title('');
            datacursormode off
            delete(findobj(f, 'Type', 'DataTip')); % Remove data tips added
            rotate3d % Rotation back on
            return
        end
        
        Rois{end+1} = [Rois{sr}(1:end-2) '~' Rois{sr}(end-1:end)];
        sr = length(Rois);
        % Select new ROI seed
        title({'Select duplicate:'; Rois{sr}}, 'Interpreter', 'none');
        dcm_obj = datacursormode(f); % Data cursor mode object
        set(dcm_obj, 'UpdateFcn', @gvf); % Change selection function
        datacursormode on
        mrr.Enable = 'on'; % Enable menu for removing dummies
    end
    
    %% Remove last ROI
    function RemLastRoiFcn(~,~)
        % Only possible for dummies or duplicates
        if contains(Rois{end}, 'Dummy #') || contains(Rois{end}, '~') || contains(Rois{end}, '^')
            Rois = Rois(1:end-1); % Remove dummy from ROI list
            Vs = Vs(1:end-1,:); % Remove from vertex list
            % Remove label
            t = findall(f, 'Type', 'Text');
            t = flipud(t(1:end-2)); % Invert order & remove two superfluous texts
            delete(t(end));
            q = findall(f, 'Type', 'Quiver');
            q = flipud(q); % Invert order 
            delete(q(end));
            % No more dummies or duplicates left?
            if ~contains(Rois{end}, 'Dummy #') && ~contains(Rois{end}, '~') && contains(Rois{end}, '^')
                mrr.Enable = 'off'; % Disable menu for removing dummies
            end
        end                
        datacursormode off
        rotate3d
    end

    %% Omit ROIs
    function OmitRoiFcn(~,~)
        srs = listdlg('PromptString', 'Select ROI or Cancel to Finish', 'ListString', Rois);
        if ~isempty(srs) 
            for sr = srs
                fv = find(Srf.Sphere(:,2)' > 90 & Srf.Data(FsVal,:) == 0); % Frontopolar vertices with zero field sign
                if strcmpi(Srf.Hemisphere, 'bi')
                    if contains(Rois{sr}, 'lh_')
                        fv = fv(fv <= Srf.Nvert_Lhem); % Keep omitted ROI in left hemisphere
                    elseif contains(Rois{sr}, 'rh_')
                        fv = fv(fv > Srf.Nvert_Lhem); % Keep omitted ROI in right hemisphere
                    end
                end
                fv = fv(randperm(length(fv))); % Shuffled vertex list
                fv = fv(1); % Random vertex
                % Reassign ROI to empty vertex
                Vs(sr,1) = fv;
                Vs(sr,2) = 0;
                % Grey out ROI label
                t = findall(f, 'Type', 'Text');
                t = flipud(t(1:end-2)); % Invert order & remove two superfluous texts
                q = findall(f, 'Type', 'Quiver');
                q = flipud(q); % Invert order 
                t(sr).Position = LabSph(fv,:);
                q(sr).XData = LabSph(fv,1);
                q(sr).YData = LabSph(fv,2);
                q(sr).ZData = LabSph(fv,3);
                q(sr).UData = Srf.Sphere(fv,1) - LabSph(fv,1);
                q(sr).VData = Srf.Sphere(fv,2) - LabSph(fv,2);
                q(sr).WData = Srf.Sphere(fv,3) - LabSph(fv,3);
                t(sr).FontAngle = 'italic';
                t(sr).Color = [.5 .5 .5];
                q(sr).Color = [.5 .5 .5];
                title({['Omitting ' Rois{sr}]; '(Moved to front)'}, 'Interpreter', 'none');
            end
        end
        datacursormode off
        rotate3d
    end

    %% Vertex selection 
    function txt = gvf(~,evnt)
        % Import vertex coordinates in mesh space 
        global Vertices sr     
        
        % Select new ROI seed
        pos = get(evnt, 'Position'); % Selected position in mesh space coordinates          
        v = find(Vertices(:,1) == pos(1) & Vertices(:,2) == pos(2) & Vertices(:,3) == pos(3), 1); % Selected vertex number
        
        % In bilateral Srf only allow selecting in correct hemisphere
        if strcmpi(Srf.Hemisphere, 'bi')
            if strcmpi(Rois{sr}(1:3), 'lh_') && v > Srf.Nvert_Lhem
                txt = 'Wrong hemisphere!';
                return % Wrongly in right hemisphere!
            elseif strcmpi(Rois{sr}(1:3), 'rh_') && v <= Srf.Nvert_Lhem
                txt = 'Wrong hemisphere!';
                return % Wrongly in left hemisphere!
            end
        end

        % Check if correct field sign or zero for omitting
        if ~contains(Rois{sr}, 'Dummy #') && Srf.Data(FsVal,v) ~= 0 && Rois{sr}(end) ~= '?' && sign(Srf.Data(FsVal,v)) ~= strsign(Rois{sr}(end))
            txt = 'Wrong field sign!';
            return % Not correct field sign!
        end

        % Assign ROI seed
        Vs(sr,:) = [v Srf.Data(FsVal,v)]; % Collect seed vertices 
        txt = ''; % No data tip info

        % Adding ROI?
        if Rois{sr}(end) == '?'
            if Vs(sr,2) < 0
                Rois{sr}(end) = '-';
            elseif Vs(sr,2) > 0
                Rois{sr}(end) = '+';
            end
        end       

        % Report & label now ROI seed
        if contains(Rois{sr}, 'Dummy #')
            q = findall(f, 'Type', 'Quiver');
            if length(q) < length(Rois)
                quiver3(LabSph(v,1), LabSph(v,2), LabSph(v,3), Srf.Sphere(v,1)-LabSph(v,1), Srf.Sphere(v,2)-LabSph(v,2), Srf.Sphere(v,3)-LabSph(v,3), 'off', 'LineWidth', 2);
                text(LabSph(v,1), LabSph(v,2), LabSph(v,3), Rois{sr}, 'Interpreter', 'none', 'Backgroundcolor', 'w', 'EdgeColor', 'k');
                samsrf_disp(['Assigned dummy ROI to #' num2str(v) ': Gradient Angle = ' num2str(Srf.Data(FsVal,v))]);
            else
                samsrf_disp(['Reassigned dummy ROI to #' num2str(v) ': Gradient Angle = ' num2str(Srf.Data(FsVal,v))]);
            end
        else
            % If adding or duplicating ROI & no label yet
            if sr == length(Rois) && length(findall(f, 'Type', 'Quiver')) < sr
                quiver3(LabSph(v,1), LabSph(v,2), LabSph(v,3), Srf.Sphere(v,1)-LabSph(v,1), Srf.Sphere(v,2)-LabSph(v,2), Srf.Sphere(v,3)-LabSph(v,3), 'off', 'LineWidth', 2);
                text(LabSph(v,1), LabSph(v,2), LabSph(v,3), Rois{sr}, 'Interpreter', 'none', 'Backgroundcolor', 'w', 'EdgeColor', 'k');
                if contains(Rois{sr}, '+') || contains(Rois{sr}, '-')
                    % Duplicating ROI
                    samsrf_disp(['Duplicating ' Rois{sr} ' at #' num2str(v) ': Gradient Angle = ' num2str(Srf.Data(FsVal,v))]);
                else
                    % Adding ROI
                    samsrf_disp(['Adding ' Rois{sr} ' at #' num2str(v) ': Gradient Angle = ' num2str(Srf.Data(FsVal,v))]);
                end
            else
                samsrf_disp(['Reassigned ' Rois{sr} ' to #' num2str(v) ': Gradient Angle = ' num2str(Srf.Data(FsVal,v))]);
            end
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

        % Grey out label if no field sign
        if Srf.Data(FsVal,v) == 0
            t(sr).FontAngle = 'italic';
            t(sr).Color = [.5 .5 .5];
            q(sr).Color = [.5 .5 .5];
        else
            t(sr).FontAngle = 'normal';
            if contains(Rois{sr}, 'Dummy #')
                t(sr).Color = [1 0 0];
            else
                t(sr).Color = [0 0 0];
            end
            q(sr).Color = [1 1 1];            
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
                mdm.Text = ['Max Cluster Size:   ' num2str(MaxSize)];
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
                mdt.Text = ['Gap Tolerance:      ' num2str(FillTol)];
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
                mdc.Text = ['Min Gradient Angle: ' num2str(CritGrd)];
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
                    % Determine ROI number
                    d = strfind(Rois{r}, '~'); % Any duplicates?
                    if isempty(d)
                        lr = r; % No duplicates
                    else
                        lr = find(contains(Rois, Rois{r}(1:d-1)), 1); % Original ROI
                    end
                    % Label seed vertex
                    LabMap(Rvx{r}) = sign(Vs(r,2)) * (lr - (Vs(r,1) > Srf.Nvert_Lhem)*AdjRoiN);
                end
            else
                Rvx{r} = NaN;
            end    
        end
        
        % Keep growing regions 
        title({'Filling'; 'ROIs'});
        Bvx = {}; % Borders of ROIs
        for i = 1:30  
            % Stop if figure was closed
            if ~isvalid(p)
                samsrf_disp('FastDelineationTool was closed by user.');
                return
            end
            % Grow around labelled regions
            for r = 1:size(Vs,1)
                % Now start filling
                cB = samsrf_neighbours(Rvx{r}, Srf.Faces); % Neighbours around current ROI
                cB(sign(Vs(r,2))*Srf.Data(FsVal,cB) <= CritGrd) = []; % Remove vertices below threshold 
                cB(ismember(cB, find(LabMap))) = []; % Remove already labelled
                % Only keep filling if new border sufficiently large
                if length(Rvx{r}) <= MaxSize || length(cB) > FillTol
                    Rvx{r} = [Rvx{r}; cB]; % Add vertices to ROI
                    Bvx{r} = samsrf_borderpath(Srf, Rvx{r}); % New border 
                    if contains(Rois{r}, 'Dummy #')
                        LabMap(cB) = Inf; % Label as dummy
                    else
                        % Determine ROI number
                        d = strfind(Rois{r}, '~'); % Any duplicates?
                        if isempty(d)
                            lr = r; % No duplicates
                        else
                            lr = find(contains(Rois, Rois{r}(1:d-1)), 1); % Original ROI
                        end
                        % Label remaining border
                        LabMap(cB) = sign(Vs(r,2)) * (lr - (Vs(r,1) > Srf.Nvert_Lhem)*AdjRoiN); 
                    end
                end
            end
            samsrf_surf(Srf, 'Sphere', 0, Bvx, '', 'Polar', p); % Show polar map 
            pause(.01); % Short break for visualisation
        end
        title('');
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
        LabMap = dSrf.Data(end,:);
        % Renumber ROI labels by rank 
        ur = unique(dSrf.Data(end,:));
        [~,sx] = sort(abs(ur));
        ur = ur(sx);
        ur(ur == 0) = [];
        r = 1;
        for u = ur
            LabMap(dSrf.Data(end,:)==u) = r * sign(u); 
            r = r + 1;
        end        
        samsrf_surf(dSrf, 'Sphere', [.001 0 max(abs(LabMap))], '', '', LabMap, p); % Show ROI map
        rotate3d
    end

    %% Smooth ROI map
    function SmoothRoisFcn(~,~)
        delete(findobj(f, 'Type', 'DataTip')); % Remove data tips added
        title('');
        if var(Srf.Data(end,:)) > 0
            title({'Smoothing ROIs'; 'Standy by...'});
            pause(.01);
            mver = find(Srf.Data(1,:) > 0);
            NewLabMap = zeros(1,length(mver));
            samsrf_newline;
            samsrf_disp('Smoothing ROIs...')
            samsrf_progbar(0);
            for i = 1:length(mver)
                v = mver(i);
                sl = samsrf_sphere_roi(v, 5, Srf.Sphere);
                NewLabMap(i) = mode(Srf.Data(end,sl));
                samsrf_progbar(i/length(mver));
            end
            Srf.Data(end,mver) = NewLabMap;
            PostProcLog{end+1} = [char(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')) ': Smoothed ROI labels'];
            title({'Finished'; 'smoothing'});
            % Now show ROIs 
            ShowRoisFcn; 
        else
            errordlg('No ROIs have been labelled yet!', 'Error');
        end
        datacursormode off
        rotate3d
    end

    %% Plug specific gaps
    function PlugGapFcn(~,~)
        if var(Srf.Data(end,:)) == 0
            errordlg('No ROIs have been labelled yet!', 'Error');
            return
        end
        
        global sr
        sr = listdlg('PromptString', 'Select ROI or Cancel to Finish', 'ListString', Rois(~contains(Rois, 'Dummy #') & ~contains(Rois, '~')), 'SelectionMode', 'single');
        if isempty(sr) 
            title('');
            datacursormode off
            delete(findobj(f, 'Type', 'DataTip')); % Remove data tips added
            rotate3d % Rotation back on
            return
        end
        
        LabMap = Srf.Data(end,:);
        lr = sign(Vs(sr,2)) * (sr - (Vs(sr,1) > Srf.Nvert_Lhem)*AdjRoiN);
        cB = samsrf_neighbours(find(LabMap == lr), Srf.Faces); % Neighbours around current ROI
        % If bilateral only plug one hemisphere
        if strcmpi(Srf.Hemisphere,'bi')
            if Vs(sr,1) > Srf.Nvert_Lhem
                % Right hemisphere ROI
                cB(cB <= Srf.Nvert_Lhem) = []; % Remove left hemisphere vertices
            else
                % Left hemisphere ROI
                cB(cB > Srf.Nvert_Lhem) = []; % Remove right hemisphere vertices
            end
        end
        cB(Srf.Data(FsVal,cB) == 0) = []; % Remove empty vertices
        cB(LabMap(cB) ~= 0) = []; % Remove vertices from other ROIs
        LabMap(cB) = lr;                             
        Srf.Data(end,:) = LabMap; % Update map
        PostProcLog{end+1} = [char(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')) ': Plugged gaps in ' Rois{sr}];
        title({Rois{sr}; 'filled out'}, 'Interpreter', 'none');
        % Now show ROIs 
        ShowRoisFcn;             
    end

    %% Save ROI labels
    function SaveRoisFcn(~,~)
        if var(Srf.Data(end,:)) == 0
            errordlg('No ROIs have been labelled yet!', 'Error');
            return
        end
        
        delete(findobj(f, 'Type', 'DataTip')); % Remove data tips added
        Fn = SaveSeedsFcn; % Save seeds first
        [~,Fn] = fileparts(Fn); % Saved seeds filename
        mkdir(['ROIs_' Fn]); % Make folder for ROIs
        cd(['ROIs_' Fn]); % Make folder for ROIs
        title({'Saving'; 'ROIs'});
        samsrf_newline;
        % Bilateral Srf?
        if strcmpi(Srf.Hemisphere, 'bi')
            [Lsrf, Rsrf] = samsrf_hemi_srfs(Srf); 
        else
            cSrf = Srf;
        end
        % Loop thru all seeds
        for r = 1:size(Vs,1) 
             % Only save ROIs not dummies, duplicates, or empty
            if ~contains(Rois{r}, 'Dummy #') && ~contains(Rois{r}, '~') && Vs(r,2) ~= 0
                % Find ROI vertices
                lr = Srf.Data(end, Vs(r,1)); % Label number of current ROI
                rv = find(Srf.Data(end,:) == lr); % Find all ROI vertices
                
                % Bilateral Srf?
                if strcmpi(Srf.Hemisphere, 'bi')
                    if strcmpi(Rois{r}(1:3), 'lh_')
                        cSrf = Lsrf;
                        rv(rv > Srf.Nvert_Lhem) = []; % Remove right hemisphere indeces
                    elseif strcmpi(Rois{r}(1:3), 'rh_')
                        cSrf = Rsrf;
                        rv(rv <= Srf.Nvert_Lhem) = []; % Remove right hemisphere indeces
                        rv = rv - Srf.Nvert_Lhem; % Adjust vertex indeces for right hemisphere
                    end
                end

                % Include hemisphere prefix but without field sign & remove new ROI symbol ^
                samsrf_srf2label(cSrf, strrep(Rois{r}(1:end-2),'^',''), 1, rv);
            end
        end

        % Combine hemisphere labels for bilateral Srf
        if strcmpi(Srf.Hemisphere, 'bi')        
            % Find all ROIs
            uRois = dir('?h_*.label');
            uRois = {uRois.name}';
            % Strip prefix & extension
            for r = 1:length(uRois)
                uRois{r} = uRois{r}(4:end-6);
            end
            uRois = unique(uRois); % Unique ROIs
            % Create bilateral ROIs
            for r = 1:length(uRois)
                samsrf_bilat_label(Srf, uRois{r});
            end
            % Remove hemisphere ROIs
            if isunix
                !rm lh_*.label
                !rm rh_*.label
            else
                !del lh_*.label
                !del rh_*.label
            end
        end

        PostProcLog{end+1} = [char(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z')) ': ROI labels saved'];
        writecell(PostProcLog', 'Delineation.log', 'FileType', 'text');
        title('ROIs saved');
        datacursormode off
        rotate3d
        cd ..
    end

end
