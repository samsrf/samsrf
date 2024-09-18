function SamSrfAnalysis
% GUI for running SamSrf anlyses. Used for standalone app.

%% Version info
[vn,vd] = samsrf_version;

%% Load Model 
load([fileparts(mfilename('fullpath')) '/../Models/Standard_2D_Gaussian_pRF'], 'Algorithm', 'Model', 'xP');
Model = samsrf_model_defaults(Algorithm, Model); % Populate empty fields with defaults 

%% Create GUI
GuiFig = uifigure('Name', ['SamSrfAnalysis GUI v' num2str(vn)], 'Units', 'normalized', 'Position', [0.2 0.15 0.7 0.75]); % Open large GUI 
crfcn = @SamSrfCloseReq;
set(GuiFig, 'CloseRequestFcn', crfcn);
global GuiInfo wb % For use with SamSrf algorithms

% Model 
GuiModelMenu = uimenu(GuiFig,'Text','Model'); 
GuiModelNew = uimenu(GuiModelMenu,'Text','New Model'); % New model item
GuiModelNew.MenuSelectedFcn = @FcnModelNew;
GuiModelLoad = uimenu(GuiModelMenu,'Text','Load Model'); % Load model item
GuiModelLoad.MenuSelectedFcn = @FcnModelLoad;
GuiModelSave = uimenu(GuiModelMenu,'Text','Save Model'); % Save model item
GuiModelSave.MenuSelectedFcn = @FcnModelSave;

% Apertures
GuiApsMenu = uimenu(GuiFig,'Text','Apertures'); 
GuiApsLoad = uimenu(GuiApsMenu,'Text','Load Apertures'); % Select apertures
GuiApsLoad.MenuSelectedFcn = @FcnApsLoad;
GuiApsView = uimenu(GuiApsMenu,'Text','View Apertures'); % View apertures 
GuiApsView.MenuSelectedFcn = @FcnApsView;

% HRF
GuiHrfMenu = uimenu(GuiFig,'Text','HRF'); 
GuiHrfLoad = uimenu(GuiHrfMenu,'Text','Load HRF'); % Select apertures
GuiHrfLoad.MenuSelectedFcn = @FcnHrfLoad;

% Data 
GuiDataMenu = uimenu(GuiFig,'Text','Data'); 
GuiDataHemis = uimenu(GuiDataMenu,'Text','Select Hemisphere Data'); % Select GII files from one hemisphere 
GuiDataHemis.MenuSelectedFcn = @FcnDataHemis;
GuiDataBilat = uimenu(GuiDataMenu,'Text','Select Bilateral Data'); % Select GII files for both hemispheres 
GuiDataBilat.MenuSelectedFcn = @FcnDataBilat;
GuiDataVols = uimenu(GuiDataMenu,'Text','Select Volume Data'); % Select NII volume files
GuiDataVols.MenuSelectedFcn = @FcnDataVols;

% Surf folder 
GuiSurfMenu = uimenu(GuiFig,'Text','Surf'); 
GuiSurfSelect = uimenu(GuiSurfMenu,'Text','Select surf folder'); 
GuiSurfSelect.MenuSelectedFcn = @FcnSurfSelect;

% ROI
GuiRoiMenu = uimenu(GuiFig,'Text','ROI'); 
GuiRoiSelect = uimenu(GuiRoiMenu,'Text','Select ROI Label'); % Select ROI
GuiRoiSelect.MenuSelectedFcn = @FcnRoiSelect; 
GuiOccRoi = uimenu(GuiRoiMenu,'Text','Create occipital ROIs'); % Create occipital ROIs
GuiOccRoi.MenuSelectedFcn = @FcnOccRoi;

% Seed info for CF 
GuiCfMenu = uimenu(GuiFig,'Text','CF'); 
GuiCfMenu.Enable = 'off'; % Starts off inactive as undefined for standard 2D Gaussian
GuiSeedSelect = uimenu(GuiCfMenu,'Text','Select Seed ROI Label'); % Select Seed ROi
GuiSeedSelect.MenuSelectedFcn = @FcnSeedSelect; 
GuiTempSelect = uimenu(GuiCfMenu,'Text','Select Template Map'); % Select Template Map
GuiTempSelect.MenuSelectedFcn = @FcnTempSelect; 

% Analyse
GuiAnalyseMenu = uimenu(GuiFig,'Text','Analyse'); 
GuiAnalyseMenu.MenuSelectedFcn = @CompletenessCheck;
GuiAnalyseRun = uimenu(GuiAnalyseMenu,'Text','Run Analysis'); % Run Analysis
GuiAnalyseRun.MenuSelectedFcn = @FcnAnalyseRun; 
GuiAnalyseRun.Enable = 'off'; % Inactive until all necessary field filled in!

% Create grid layout
GuiGrid = uigridlayout(GuiFig, [6 3]); 
GuiGrid.RowHeight = {50 '1x' 50 '2x' 50 60}; 
GuiGrid.ColumnWidth = {400 '1x' '1x'};

% Algorithm info
GuiAlgo = uitextarea(GuiGrid);
GuiAlgo.Layout.Row = 1;
GuiAlgo.Layout.Column = 1;
GuiAlgo.Editable = 'off';
GuiAlgo.Value = {'Standard_2D_Gaussian_pRF'; ''; ['Algorithm: ' Algorithm]};
GuiAlgo.FontWeight = 'bold';

% Parameter & Model tables 
[GuiPars, GuiModel] = UpdateTables(Model);

% pRF profile
GuiPrf = uiaxes(GuiGrid);
GuiPrf.Layout.Row = [1 2];
GuiPrf.Layout.Column = 2;
axis(GuiPrf, 'off', 'square');
if isfield(Model, 'Prf_Function')
    title(GuiPrf, 'pRF profile');
    Rfp = flipud(Model.Prf_Function(xP, 200));
    contourf(GuiPrf, Rfp, 100, 'EdgeColor', 'none');
else
    title(GuiPrf, 'No visualisation of CF');
end
colormap(GuiPrf, 'berlin');
GuiPrf.CLim = [-1 1] * max(abs(Rfp(:)));
GuiPrf.Toolbar = [];
disableDefaultInteractivity(GuiPrf)

% HRF plot
GuiHrf = uiaxes(GuiGrid);
GuiHrf.YTick = [];
GuiHrf.Layout.Row = [1 2];
GuiHrf.Layout.Column = 3;
xlim(GuiHrf, [0 32]);
title(GuiHrf, 'HRF');
xlabel(GuiHrf, 'Time (s)');
ylabel(GuiHrf, 'Response (a.u.)');
GuiHrf.Toolbar = [];
disableDefaultInteractivity(GuiHrf)

% ROI label selection
GuiRoi = uitextarea(GuiGrid);
GuiRoi.Layout.Row = 3;
GuiRoi.Layout.Column = 3;
GuiRoi.Editable = 'off';
GuiRoi.Value = {'Region of Interest:'; ''; '< None selected >'}; 

% Model help text
GuiInfo = uitextarea(GuiGrid);
GuiInfo.Layout.Row = [3 6];
GuiInfo.Layout.Column = 2;
GuiInfo.Editable = 'off';
GuiModel.SelectionChangedFcn = @UpdateInfo;
GuiModel.CellEditCallback = @UpdateInfo;
GuiPars.SelectionChangedFcn = @UpdateInfo;

% Data selection
GuiFiles = uilistbox(GuiGrid);
GuiFiles.Layout.Row = 4;
GuiFiles.Layout.Column = 3;
GuiFiles.Items = {'< No files selected >'};
GuiFiles.Value = '< No files selected >';
GuiFiles.ClickedFcn = @HelpOnFiles;

% Surface folder
GuiSurf = uitextarea(GuiGrid);
GuiSurf.Layout.Row = 5;
GuiSurf.Layout.Column = 3;
GuiSurf.Editable = 'off';
GuiSurf.Value = {'Subject surf folder:'; ''; '< None selected >'};

% Data preprocessing
ProcTab = table;
ProcTab.Average = {true};
ProcTab.Normalise = {true};
ProcTab.Export = {false};
GuiProc = uitable(GuiGrid, 'Data', ProcTab, 'ColumnEditable', [true true true]);
clear ProcTab
GuiProc.Layout.Row = 6;
GuiProc.Layout.Column = 3;
GuiProc.SelectionChangedFcn = @UpdateInfo;

%% Initial update 
eventdata.EventName = 'SelectionChanged';
eventdata.DisplaySelection = [];
UpdateInfo(GuiModel, eventdata);

%% Welcome message
samsrf_clrscr; 
samsrf_disp('****************************************************************************');
samsrf_disp('     Kia ora to the Seriously Annoying MatLab Surfer Analysis Tool!');
samsrf_disp('     by D. S. Schwarzkopf from the University of Auckland, New Zealand');
samsrf_newline;
samsrf_disp(['                 Version ' num2str(vn) ', Released on ' vd]);
samsrf_disp('      (see SamSrf/ReadMe.md for what is new in this version)');
samsrf_disp('****************************************************************************');


    %% Info samsrf_display function
    function UpdateInfo(source, eventdata)
        % Update model help
        if strcmpi(eventdata.EventName, 'CellEdit') 
            % Check if troublesome fields were updated
            switch GuiModel.Data.Field{eventdata.DisplayIndices(1)}
                case 'Scaling_Factor'
                    if strcmpi(Algorithm, 'samsrf_fit_prf') && isfield(Model, 'SeedRoi')
                        % If pRF-from-CF analysis don't allow updating scaling factor
                        GuiModel.Data.Field{eventdata.DisplayIndices(1)} = {Inf};
                    end
                case 'Aperture_File'
                    if strcmpi(Algorithm, 'samsrf_fit_prf') && isfield(Model, 'SeedRoi')
                        % If pRF-from-CF analysis don't allow updating aperture file
                        GuiModel.Data.Field{eventdata.DisplayIndices(1)} = {''};
                    end
                case 'Hrf'
                    if ischar(eventdata.EditData)
                        % Only allow value scalar or empty input for HRF
                        if strcmpi(eventdata.EditData, 'Inf')
                            GuiModel.Data.Value(eventdata.DisplayIndices(1)) = {Inf};
                        elseif strcmpi(eventdata.EditData, '0')
                            GuiModel.Data.Value(eventdata.DisplayIndices(1)) = {0};
                        elseif strcmpi(eventdata.EditData, '1')
                            GuiModel.Data.Value(eventdata.DisplayIndices(1)) = {1};
                        elseif isempty(eventdata.EditData)
                            GuiModel.Data.Value(eventdata.DisplayIndices(1)) = {[]};
                        else
                            GuiModel.Data.Value(eventdata.DisplayIndices(1)) = {eventdata.PreviousData};
                        end
                    end
                case {'Param1' 'Param2' 'Param3' 'Param4' 'Param5' 'Param6' 'Param7' 'Param8' 'Param9' 'Param10'}
                    % Don't allow updating search grid definitions
                    GuiModel.Data.Value(eventdata.DisplayIndices(1)) = {eventdata.PreviousData};
            end
        elseif strcmpi(eventdata.EventName, 'SelectionChanged') && isempty(eventdata.DisplaySelection)
            % Nothing selected (Unsure if this can still happen...)
            GuiInfo.Value = '';
        else
            if strcmpi(source.Data.Properties.VariableNames{1}, 'Parameter')
                % Show model help text for this column in pRF parameter table
                switch eventdata.DisplaySelection(2)
                    case 1
                        GuiInfo.Value = ['Param_Names'; ''; string(ModelHelpText(Algorithm, 'Param_Names'))];
                    case 2
                        GuiInfo.Value = ['Scaled_Param'; ''; string(ModelHelpText(Algorithm, 'Scaled_Param'))];
                    case 3
                        GuiInfo.Value = ['Only_Positive'; ''; string(ModelHelpText(Algorithm, 'Only_Positive'))];
                end
            elseif strcmpi(source.Data.Properties.VariableNames{1}, 'Field')
                % Show model help text for this field
                GuiInfo.Value = [source.Data.Field{eventdata.DisplaySelection(1)}; ''; string(ModelHelpText(Algorithm, source.Data.Field{eventdata.DisplaySelection(1)}))];
                % If apertures selected
                if eventdata.DisplaySelection(2) == 2 && strcmpi(source.Data.Field{eventdata.DisplaySelection(1)}, 'Aperture_File') 
                    if ~isfield(Model, 'SeedRoi')
                        % Can select apertures directly if not pRF-from-CF analysis
                        FcnMatLoad('aps_*.mat', eventdata); % Select aperture file directly
                    end
                elseif eventdata.DisplaySelection(2) == 2 && (strcmpi(source.Data.Field{eventdata.DisplaySelection(1)}, 'Seed_Fine_Fit') || strcmpi(source.Data.Field{eventdata.DisplaySelection(1)}, 'Template')) 
                    FcnMatLoad('*.mat', eventdata); % Select Srf map file directly
                elseif eventdata.DisplaySelection(2) == 2 && strcmpi(source.Data.Field{eventdata.DisplaySelection(1)}, 'SeedRoi') 
                    FcnMatLoad('*.label', eventdata); % Select label file directly
                end
            elseif strcmpi(source.Data.Properties.VariableNames{1}, 'Average')
                % Show model help text for this column in processing table
                switch eventdata.DisplaySelection(2)
                    case 1
                        GuiInfo.Value = {'Average'; 
                                         ''; 
                                         'When checked, individual time series of runs (in the GII or NII files from the list) will be averaged into one run. This only makes sense if the sequence of time points is identical.'
                                         '';
                                         'When unchecked, the individual time series of all runs will be concatenated into one long run.'};
                    case 2
                        GuiInfo.Value = {'Normalise'; 
                                         ''; 
                                         'When checked, linear deterending & z-normalisation is applied to the time series of each run (in the GII or NII files from the list).'};
                    case 3
                        GuiInfo.Value = {'Export'; 
                                         ''; 
                                         'When checked, the individual pRF/CF parameter maps are exported in the same folder as the data, using the same file format as your input data (that is, surface GIIs or volume NIIs). You can then load these maps in the imaging analysis software of your choice.';
                                         '';
                                         'Irrespective of your choice here, the analysis will always save the maps in SamSrf MAT format as well.'};
                end
            end
        end

        % Update HRF plot
        if contains(Algorithm, '_cf') 
            % CF algorithms don't use HRF
            plot(GuiHrf, NaN(1,32));
            title(GuiHrf, 'No HRF in CF models!');
            GuiHrf.YLim = [-1.1 1.1];
            axis(GuiHrf, 'off');

            % Update CF profile plot
            if GuiModel.Data.Value{strcmpi(GuiModel.Data.Field, 'Fit_pRF')} == 0
                % Convex hull estimate
                title(GuiPrf, 'Convex hull pRF estimate');
                [chx, chy] = pol2cart([14 34 78 127 167 198 213 275 312 349]/180*pi, [.5 .4 .6 .5 .4 .4 .3 .5 .4 .4]);
                warning off
                plot(GuiPrf, polyshape(chx, chy));
                warning on
                hold(GuiPrf, 'on');
                scatter(GuiPrf, .05, .1, 'k', 'filled');
                hold(GuiPrf, 'off');
                axis(GuiPrf, [-1 1 -1 1]);
            elseif GuiModel.Data.Value{strcmpi(GuiModel.Data.Field, 'Fit_pRF')} == -1
                % Summary stats only
                title(GuiPrf, 'Summary stats only');
                scatter(GuiPrf, [0 -.05], [0 .05], [6000 20], [.7 .7 .9; 0 0 0], 'filled', 'MarkerEdgeColor', 'k');
                axis(GuiPrf, [-1 1 -1 1]);
            elseif GuiModel.Data.Value{strcmpi(GuiModel.Data.Field, 'Fit_pRF')} == 1
                % Fit 2D Gaussian CF 
                title(GuiPrf, '2D Gaussian CF fit');
                Rfp = prf_gaussian_rf(.2, .3, .4);
                contourf(GuiPrf, Rfp, 100, 'EdgeColor', 'none');
                axis(GuiPrf,[0 200 0 200]);
                GuiPrf.CLim = [-1 1] * max(abs(Rfp(:)));
            else
                title(GuiPrf, 'Invalid input to Fit_pRF!!!');
                plot(GuiPrf, 0);
            end
        else
            % pRF algorithms use HRF
            axis(GuiHrf, 'on');
            Hrf = GuiModel.Data.Value{strcmpi(GuiModel.Data.Field, 'Hrf')};
            TR = GuiModel.Data.Value{strcmpi(GuiModel.Data.Field, 'TR')};
            if isempty(Hrf)
                % de Haas canonical HRF
                ExaHrf = samsrf_hrf(TR);
                ExaHrf = ExaHrf/max(ExaHrf);
                plot(GuiHrf, 0:TR:32, ExaHrf, 'k', 'LineWidth', 2);
                title(GuiHrf, 'HRF: de Haas canonical');
                GuiHrf.YLim = [min(ExaHrf)-.1 1.1];
    
            elseif Hrf == 0
                % SPM canonical HRF
                ExaHrf = samsrf_doublegamma(TR, [6 16 1 1 6 0 32]);
                ExaHrf = ExaHrf/max(ExaHrf);
                plot(GuiHrf, 0:TR:32, ExaHrf, 'k', 'LineWidth', 2);
                title(GuiHrf, 'HRF: SPM canonical');
                GuiHrf.YLim = [min(ExaHrf)-.1 1.1];
    
            elseif isinf(Hrf)
                % Fit HRF on the fly
                if contains(Algorithm, 'revcor')
                    plot(GuiHrf, NaN(1,32), 'k', 'LineWidth', 2);
                    title(GuiHrf, 'Cannot fit HRF in reverse-correlation analysis!');
                    GuiHrf.YLim = [-1.1 1.1];
                    GuiModel.Data.Value{strcmpi(GuiModel.Data.Field, 'Hrf')} = {[]};
                else
                    ExaHrf = [samsrf_hrf(TR) samsrf_doublegamma(TR, [6 16 1 1 6 0 32]) samsrf_doublegamma(TR, [8 13 1 1 2 0 32])];
                    ExaHrf = ExaHrf./repmat(max(ExaHrf),size(ExaHrf,1),1);
                    plot(GuiHrf, 0:TR:32, ExaHrf, 'LineWidth', 2);
                    title(GuiHrf, 'Fitting HRF on the fly');
                    GuiHrf.YLim = [min(ExaHrf(:))-.1 1.1];
                end

            elseif Hrf == 1
                % No HRF convolution
                plot(GuiHrf, ones(1,32), 'k', 'LineWidth', 2);
                title(GuiHrf, 'No HRF convolution');
                GuiHrf.YLim = [0 2];
    
            elseif ischar(Hrf)
                % Use separate HRF fit
                load(Hrf, 'fP');
                if exist('fP', 'var')
                    ExaHrf = samsrf_doublegamma(Model.TR, [fP(1:2) 1 1 fP(3) 0 32])' * fP(4);
                    ExaHrf = ExaHrf/max(ExaHrf);
                    plot(GuiHrf, 0:TR:32, ExaHrf, 'k', 'LineWidth', 2);
                    title(GuiHrf, Hrf);
                    GuiHrf.YLim = [min(ExaHrf(:))-.1 1.1];
                end

            elseif strcmpi(class(Hrf),'double')
                % HRF stored as vector in CSV file
                Hrf = Hrf/max(Hrf);
                plot(GuiHrf, 0:TR:32, Hrf, 'k', 'LineWidth', 2);
                title(GuiHrf, 'User provided HRF vector');
                GuiHrf.YLim = [min(Hrf(:))-.1 1.1];
            else 
                % Cannot provide numbers via GUI
                plot(GuiHrf, NaN(1,32));
                title(GuiHrf, 'HRF undefined!!!');
                GuiHrf.YLim = [-1.1 1.1];
            end
    
            % Update pRF profile plot
            if sum(strcmpi(GuiModel.Data.Field, 'Prf_Function')) == 1
                % pRF function field exists
                if ~strcmpi(class(Model.Prf_Function), 'function_handle')
                    if GuiModel.Data.Value{strcmpi(GuiModel.Data.Field, 'Prf_Function')} == 0
                        % Convex hull estimate
                        title(GuiPrf, 'Convex hull pRF estimate');
                        [chx, chy] = pol2cart([14 34 78 127 167 198 213 275 312 349]/180*pi, [.5 .4 .6 .5 .4 .4 .3 .5 .4 .4]);
                        warning off
                        plot(GuiPrf, polyshape(chx, chy));
                        warning on
                        hold(GuiPrf, 'on');
                        scatter(GuiPrf, .05, .1, 'k', 'filled');
                        hold(GuiPrf, 'off');
                        axis(GuiPrf, [-1 1 -1 1]);
                    elseif GuiModel.Data.Value{strcmpi(GuiModel.Data.Field, 'Prf_Function')} == -1
                        % Summary stats only
                        title(GuiPrf, 'Summary stats only');
                        scatter(GuiPrf, [0 -.05], [0 .05], [6000 20], [.7 .7 .9; 0 0 0], 'filled', 'MarkerEdgeColor', 'k');
                        axis(GuiPrf, [-1 1 -1 1]);
                    else
                        title(GuiPrf, 'Invalid pRF function!!!');
                        plot(GuiPrf, 0);
                    end
                else
                    % pRF function defined
                    tstr = func2str(Model.Prf_Function);
                    tstr = strrep(tstr, '@(P,ApWidth)', '');
                    tstr = strrep(tstr, ',ApWidth)', ')');
                    Rfp = flipud(Model.Prf_Function(xP, 200));
                    contourf(GuiPrf, Rfp, 100, 'EdgeColor', 'none');
                    hold(GuiPrf, 'on');
                    axis(GuiPrf,[0 200 0 200]);
                    if contains(tstr, '(0,P(1),P(2)')
                        % Vertical tuning function
                        plot(GuiPrf, [100 100], [0 200], 'w:');
                    elseif contains(tstr, 'P(1),0,P(2)')
                        % Horizontal tuning function
                        plot(GuiPrf, [0 200], [100 100], 'w:');
                    elseif contains(tstr, 'cosd(P(1))/2,sind(P(1))/2')
                        % Circular tuning function
                        [ctx, cty] = pol2cart((0:360)/180*pi, 25);
                        plot(GuiPrf, ctx+100, cty+100, 'w:');
                        axis(GuiPrf,[50 150 50 150]);
                    else
                        % 2D pRF function
                        plot(GuiPrf, [0 200], [100 100], 'w:');
                        plot(GuiPrf, [100 100], [0 200], 'w:');
                    end
                    hold(GuiPrf, 'off');
                    GuiPrf.CLim = [-1 1] * max(abs(Rfp(:)));
                    for i = 1:size(GuiPars.Data.Parameter)
                        tstr = strrep(tstr, ['P(' num2str(i) ')'], GuiPars.Data.Parameter{i});
                    end
                    title(GuiPrf, tstr, 'Interpreter', 'none');
                end
            else
                % No pRF function so something is wrong
                title(GuiPrf, 'No pRF function defined!!!');
                plot(GuiPrf, 0);
                axis(GuiPrf,[0 200 0 200]);
            end
        end

        % (De-)activate CF & Apertures menus?
        if isfield(Model, 'SeedRoi') && isfield(Model, 'Template') 
            GuiCfMenu.Enable = 'on'; % CF menu available
            GuiApsMenu.Enable = 'off'; % No apertures menu
        else
            GuiApsMenu.Enable = 'on'; % Apertures menu available
        end

        % (De-)activate Analyse menu?
        GuiAnalyseRun.Enable = 'on';
        if isnan(GuiModel.Data.Value{strcmpi(GuiModel.Data.Field, 'Scaling_Factor')})
            % Scaling factor not yet defined
            GuiAnalyseRun.Enable = 'off';
        end
        if strcmpi(Algorithm, 'samsrf_fit_prf') || strcmpi(Algorithm, 'samsrf_revcor_prf')                        
            % Apertures define
            if isempty(GuiModel.Data.Value{strcmpi(GuiModel.Data.Field, 'Aperture_File')})
                GuiAnalyseRun.Enable = 'off';
            end
        end
        % Data files defined?
        if strcmpi(GuiFiles.Items{1}, '< No files selected >')
            GuiAnalyseRun.Enable = 'off';
        end
        % Surf folder selected?
        if strcmpi(GuiFiles.Items{1}(end-3:end), '.nii')
            GuiSurf.Value{3} = '< Not required for NII >';
            GuiSurfMenu.Enable = 'off';
        else
            GuiSurfMenu.Enable = 'on';
            if strcmpi(GuiSurf.Value{3}, '< None selected >') 
                GuiAnalyseRun.Enable = 'off';
            end
        end
        % Is pRF-from-CF analysis?
        if isfield(Model, 'SeedRoi') && isfield(Model, 'Template') && isempty(GuiModel.Data.Value{strcmpi(GuiModel.Data.Field, 'SeedRoi')}) && isempty(GuiModel.Data.Value{strcmpi(GuiModel.Data.Field, 'Template')})
            GuiAnalyseRun.Enable = 'off';
        end            
    end

    %% Completeness check info
    function CompletenessCheck(srf,event)
        CompCheckInfo = {'COMPLETENESS CHECK:'; ''};
        cl = 3;
        
        % ROI undefined
        if strcmpi(GuiRoi.Value{3}, '< None selected >')
            CompCheckInfo{cl,1} = 'Warning: No Region of Interest selected';
            cl = cl + 1;
            CompCheckInfo{cl,1} = '';
            cl = cl + 1;
        end

        % Scaling factor undefined
        if isnan(GuiModel.Data.Value{strcmpi(GuiModel.Data.Field, 'Scaling_Factor')})
            CompCheckInfo{cl} = 'Error: Scaling Factor not yet defined!';
            cl = cl + 1;
            CompCheckInfo{cl,1} = '';
            cl = cl + 1;
        end
        
        % Apertures undefined
        if strcmpi(Algorithm, 'samsrf_fit_prf') || strcmpi(Algorithm, 'samsrf_revcor_prf')                        
            if isempty(GuiModel.Data.Value{strcmpi(GuiModel.Data.Field, 'Aperture_File')})
                CompCheckInfo{cl,1} = 'Error: Apertures not yet defined!';
                cl = cl + 1;
                CompCheckInfo{cl,1} = '';
                cl = cl + 1;
            end
        end

        % Data files defined?
        if strcmpi(GuiFiles.Items{1}, '< No files selected >')
            CompCheckInfo{cl,1} = 'Error: No data files selected!';
            cl = cl + 1;
            CompCheckInfo{cl,1} = '';
            cl = cl + 1;
        end
        
        % Surf folder selected?
        if strcmpi(GuiFiles.Items{1}(end-3:end), '.gii')
            if strcmpi(GuiSurf.Value{3}, '< None selected >') 
                CompCheckInfo{cl,1} = 'Error: Surf folder not yet selected!';
                cl = cl + 1;
                CompCheckInfo{cl,1} = '';
                cl = cl + 1;
            end
        end
        
        % Is pRF-from-CF analysis?
        if isfield(Model, 'SeedRoi') && isfield(Model, 'Template') 
            if isempty(GuiModel.Data.Value{strcmpi(GuiModel.Data.Field, 'SeedRoi')}) 
                CompCheckInfo{cl,1} = 'Error: Seed ROI not yet selected!';
                cl = cl + 1;
                CompCheckInfo{cl,1} = '';
                cl = cl + 1;
            end
            if isempty(GuiModel.Data.Value{strcmpi(GuiModel.Data.Field, 'Template')})
                CompCheckInfo{cl,1} = 'Error: Template map not yet selected!';
                cl = cl + 1;
                CompCheckInfo{cl,1} = '';
            end
        end      
        if cl == 3
            CompCheckInfo{cl,1} = 'Ready to analyse!';
        end
        
        % Now update info field
        GuiInfo.Value = CompCheckInfo;
    end

    %% Turn Model struct into tables  
    function [GuiPars, GuiModel] = UpdateTables(Model)        
        % Create table for parameters?
        if isfield(Model, 'Param_Names')
            % Ensure boolean fields exist
            if ~isfield(Model, 'Scaled_Param')
                Model.Scaled_Param = false(1,length(Model.Param_Names));
            end
            % Create table
            ParTab = table; 
            % Ensuring column vectors
            ParTab.Parameter = Model.Param_Names(:);
            ParTab.Scaled = logical(Model.Scaled_Param(:));
            if isfield(Model, 'Only_Positive')
                ParTab.Positive = logical(Model.Only_Positive(:));
            end
        else
            % Create dummy table
            ParTab = table; 
            ParTab.Parameter = '';
            ParTab.Scaled = '';
            ParTab.Positive = '';
        end
        
        % Create table for other fields
        ModTab = table;
        mfn = fieldnames(Model); % All field names
        m = 1;
        % Loop thru all fields
        warning off
        for i = 1:length(mfn)
            % Only if not a pRF parameter
            if ~strcmpi(mfn{i},'Param_Names') && ~strcmpi(mfn{i},'Scaled_Param') && ~strcmpi(mfn{i},'Only_Positive') 
                % If parameter is not unused
                if contains(mfn{i}, 'Param') && length(Model.(mfn{i})) == 1 && Model.(mfn{i}) == 0
                    keep = false;
                else
                    keep = true;
                end
                if keep
                    ModTab.Field{m} = mfn{i}; % Name of this field
                    ModTab.Value{m} = Model.(mfn{i}); % Value of this field
                    m = m + 1;
                end
            end
        end
        warning on
    
        % pRF parameter table GUI
        GuiPars = uitable(GuiGrid, 'Data', ParTab, 'ColumnEditable', [true true true]);
        GuiPars.Layout.Row = 2;
        GuiPars.Layout.Column = 1;
        if ~isfield(Model, 'Param_Names')
            GuiPars.Enable = 'off';
        end
        
        % Model field table GUI
        GuiModel = uitable(GuiGrid, 'Data', ModTab, 'ColumnEditable', [false true]);
        GuiModel.Layout.Row = [3 6];
        GuiModel.Layout.Column = 1;
    end

    %% New Model file
    function FcnModelNew(srf,event)       
        load([fileparts(mfilename('fullpath')) '/../Models/Standard_2D_Gaussian_pRF'], 'Algorithm', 'Model', 'xP');
        Model = samsrf_model_defaults(Algorithm, Model); % Populate empty fields with defaults 
        GuiAlgo.Value = {'Standard_2D_Gaussian_pRF'; ''; ['Algorithm: ' Algorithm]};
        [GuiPars, GuiModel] = UpdateTables(Model);
        GuiModel.SelectionChangedFcn = @UpdateInfo;
        GuiModel.CellEditCallback = @UpdateInfo;
        GuiPars.SelectionChangedFcn = @UpdateInfo;
        eventdata.EventName = 'SelectionChanged';
        eventdata.DisplaySelection = [1 1];
        UpdateInfo(GuiModel, eventdata);        
    end

    %% Load Model file
    function FcnModelLoad(src,event)
        [fn, pn] = uigetfile('*.mat');
        if ~isscalar(fn)
            vs = whos('-file', [pn fn]);
            vs = {vs.name}';
            if sum(strcmpi(vs, 'Algorithm')) == 0 || sum(strcmpi(vs, 'xP')) == 0
                samsrf_error([pn fn ' was not saved by SamSrfAnalysis or has been corrupted!']);
            else
                warning off
                load([pn fn], 'Algorithm', 'Model', 'xP');
                warning on
                Model = samsrf_model_defaults(Algorithm, Model); % Populate empty fields with defaults 
                % If this is pRF-from-CF analysis, automatically set fields 
                if strcmpi(Algorithm, 'samsrf_fit_prf') && isfield(Model, 'SeedRoi')
                    Model.Aperture_File = '[Set by analysis]';
                    Model.Scaling_Factor = Inf;
                end
                
                GuiAlgo.Value = {fn(1:end-4) ; ''; ['Algorithm: ' Algorithm]};
                [GuiPars, GuiModel] = UpdateTables(Model);
                GuiModel.SelectionChangedFcn = @UpdateInfo;
                GuiModel.CellEditCallback = @UpdateInfo;
                GuiPars.SelectionChangedFcn = @UpdateInfo;
                eventdata.EventName = 'SelectionChanged';
                eventdata.DisplaySelection = [1 1];
                UpdateInfo(GuiModel, eventdata);        
            end
        end
    end
    
    %% Update Model struct from tables 
    function Model = UpdateModel(Model)
        % Fill in pRF parameters
        if isfield(Model, 'Param_Names')
            Model.Param_Names = GuiPars.Data.Parameter';
            Model.Scaled_Param = GuiPars.Data.Scaled';
        end
        if isfield(Model, 'Only_Positive')
            Model.Only_Positive = GuiPars.Data.Positive';
        end
        % Fill in remaining fields
        fs = GuiModel.Data.Field;
        for i = 1:length(fs)
            Model.(fs{i}) = GuiModel.Data.Value{i};
        end
    end

    %% Save Model file
    function FcnModelSave(src,event)
        Model = UpdateModel(Model);
        [fn, pn] = uiputfile('*.mat');
        if ~isscalar(fn)
            save([pn fn], 'Algorithm', 'Model', 'xP');
        end
    end

    %% Select aperture file
    function FcnApsLoad(src,event)
        [rn, rp] = uigetfile('aps_*.mat');
        if ~isscalar(rn)
            GuiModel.Data.Value(strcmpi(GuiModel.Data.Field, 'Aperture_File')) = {[rp rn(1:end-4)]}; 
            UpdateInfo(GuiModel, eventdata);
        end
    end

    %% View apertures 
    function FcnApsView(src,event)
        ViewApertures(GuiModel.Data.Value{strcmpi(GuiModel.Data.Field, 'Aperture_File')}); 
    end

    %% Select HRF file
    function FcnHrfLoad(src,event)
        [hn, hp] = uigetfile({'hrf_*.mat', 'HRF fitting parameters'; 'hrf_*.csv', 'Comma-separated HRF vector'});
        if ~isscalar(hn)
            if strcmpi(hn(end-3:end), '.csv')
                VecHrf = readtable([hp hn]);
                GuiModel.Data.Value(strcmpi(GuiModel.Data.Field, 'Hrf')) = {VecHrf.Var1};
            else
                GuiModel.Data.Value(strcmpi(GuiModel.Data.Field, 'Hrf')) = {[hp hn(1:end-4)]}; 
            end
            eventdata.EventName = 'SelectionChanged';
            eventdata.DisplaySelection = [find(strcmpi(GuiModel.Data.Field, 'Hrf')) 1];
            UpdateInfo(GuiModel, eventdata);
        end
    end

    %% Select hemispheric surface files 
    function FcnDataHemis(src,event)
        [DataFiles, DataPath] = uigetfile({'*.gii', 'GIfTI surface files'}, ...
                                           'Select files to analyse', 'MultiSelect', 'on');
        if ~isscalar(DataFiles)
            if ischar(DataFiles)
                DataFiles = {DataFiles};
            end
            DataFiles = DataFiles';
            % Add path
            for i = 1:length(DataFiles)
                DataFiles{i} = [DataPath DataFiles{i}];
            end
            % Update GUI
            GuiFiles.Items = DataFiles;
            UpdateInfo(GuiModel, eventdata);
        end
    end

    %% Select bilateral surface files 
    function FcnDataBilat(src,event)
        [DataFiles, DataPath] = uigetfile({'lh*.gii', 'GIfTI surface files'}, ...
                                           'Select files to analyse (must have "lh" suffix for combining hemispheres)', 'MultiSelect', 'on');
        if ~isscalar(DataFiles)
            if ischar(DataFiles)
                DataFiles = {DataFiles};
            end
            DataFiles = DataFiles';
            % Add path
            for i = 1:length(DataFiles)
                DataFiles{i} = [DataPath 'bi' DataFiles{i}(3:end)]; % Replace 'lh' with 'bi' for joining later
            end
            % Update GUI
            GuiFiles.Items = DataFiles;
            UpdateInfo(GuiModel, eventdata);
        end
    end

    %% Select volumetric files 
    function FcnDataVols(src,event)
        [DataFiles, DataPath] = uigetfile({'*.nii', 'NIfTI volume files'}, ...
                                           'Select files to analyse', 'MultiSelect', 'on');
        if ~isscalar(DataFiles)
            if ischar(DataFiles)
                DataFiles = {DataFiles};
            end
            DataFiles = DataFiles';
            % Add path
            for i = 1:length(DataFiles)
                DataFiles{i} = [DataPath DataFiles{i}];
            end
            % Update GUI
            GuiFiles.Items = DataFiles;
            UpdateInfo(GuiModel, eventdata);
        end
    end


    %% Select surf folder
    function FcnSurfSelect(src,event)
        sp = uigetdir('../surf', 'Select surf folder');
        if ~isscalar(sp)
            GuiSurf.Value = {'Subject surf folder:'; ''; sp}; 
        end
    end

    %% Select ROI label
    function FcnRoiSelect(src,event)
        [rn, rp] = uigetfile({'*.label', 'FreeSurfer label'; '*.nii', 'Binary NIfTI mask'});
        if ~isscalar(rn)
            GuiRoi.Value = {'Region of Interest:'; ''; [rp rn]}; 
            UpdateInfo(GuiModel, eventdata);
        end
    end

    %% Create occipital ROI
    function FcnOccRoi(src,event)
        dp = fileparts(GuiFiles.Items{1});
        if isempty(dp)
            samsrf_clrscr;
            samsrf_disp('ERROR: You need to select data files first!');
        else
            if ~strcmpi(GuiSurf.Value{3}, '< None selected >') && ~strcmpi(GuiSurf.Value{3}, '< Not required for NII >')
                samsrf_clrscr;
                samsrf_disp('Creating occipital ROIs...');
                cf = pwd;
                cd(dp);
                warning off
                MakeOccRoi(GuiSurf.Value{3})
                warning on
                cd(cf);
                samsrf_disp(' Occipital ROIs created.');
            elseif strcmpi(GuiSurf.Value{3}, '< None selected >')
                samsrf_clrscr;
                samsrf_disp('ERROR: No surf folder selected!');
            end
        end
    end

    %% Select seed ROI label
    function FcnSeedSelect(src,event)
        [sn, sp] = uigetfile({'*.label', 'FreeSurfer label'});
        if ~isscalar(sn)
            [sp, sn] = fileparts([sp sn]); % Remove extension
            GuiModel.Data.Value{strcmpi(GuiModel.Data.Field, 'SeedRoi')} = [sp filesep sn]; 
            UpdateInfo(GuiModel, eventdata);
        end
    end

    %% Select template map 
    function FcnTempSelect(src,event)
        [tn, tp] = uigetfile('*.mat');
        if ~isscalar(tn)
            [tp, tn] = fileparts([tp tn]); % Remove extension
            GuiModel.Data.Value{strcmpi(GuiModel.Data.Field, 'Template')} = [tp tn]; 
            UpdateInfo(GuiModel, eventdata);
        end
    end

    %% Run analysis
    function FcnAnalyseRun(src,event)
        % Should be unnecessary but update everything 
        UpdateInfo(GuiModel, eventdata);
        samsrf_clrscr;
        samsrf_disp(['Running ' Algorithm ' algorithm...']);
        samsrf_newline;

        % Update model fields
        Model = UpdateModel(Model);
        
        % Extract ROI label
        Roi = GuiRoi.Value{3};
        if strcmpi(Roi, '< None selected >')
            Roi = '';
        end

        % Convert data to Srf
        Srf = [];
        DataFiles = GuiFiles.Items;
        [dp,dn,de] = fileparts(DataFiles{1});
        % What file format?
        if strcmpi(de, '.nii')
            % Volumetric data
            if ~isempty(Roi)
                Roi = Roi(1:end-4); % Remove NII extension
            end
            samsrf_disp('Converting volumetric NII files to Srf...');
            samsrf_disp(DataFiles');
            Srf = samsrf_vol2mat(DataFiles, Roi, GuiProc.Data.Normalise{1}, GuiProc.Data.Average{1}, true); % ROI is used for creating data file
            samsrf_newline;
            Srf.Functional = DataFiles';
            Roi = ''; % Clear now
        elseif strcmpi(de, '.gii') 
            if ~isempty(Roi)
                Roi = Roi(1:end-6); % Remove label extension
            end
            if strcmpi(dn(1:2), 'bi')
                % Binocular data 
                LDataFiles = {};
                RDataFiles = {};
                for i = 1:length(DataFiles)
                    [dp,dn,de] = fileparts(DataFiles{i});
                    LDataFiles{i} = [dp filesep 'lh' dn(3:end)];
                    RDataFiles{i} = [dp filesep 'rh' dn(3:end)];
                end
                % GIfTI format
                samsrf_disp('Converting left hemisphere GII files to Srf...');
                samsrf_disp(LDataFiles');
                Lsrf = samsrf_gii2srf(LDataFiles, [GuiSurf.Value{3} filesep 'lh'], GuiProc.Data.Normalise{1}, GuiProc.Data.Average{1}, true, '');
                samsrf_newline;
                samsrf_disp('Converting right hemisphere GII files to Srf...');
                samsrf_disp(RDataFiles');
                Rsrf = samsrf_gii2srf(RDataFiles, [GuiSurf.Value{3} filesep 'rh'], GuiProc.Data.Normalise{1}, GuiProc.Data.Average{1}, true, '');
                samsrf_newline;
                % Combine hemispheres
                samsrf_disp('Combining hemisphere Srfs...');
                Srf = samsrf_bilat_srf(Lsrf, Rsrf);
                Srf.Functional = DataFiles';
                samsrf_newline;
                clear Lsrf Rsrf
            else
                % Hemispheric data 
                for i = 1:length(DataFiles)
                    [dp,dn,de] = fileparts(DataFiles{i});
                    DataFiles{i} = [dp filesep dn];
                end
                % GIfTI format
                samsrf_disp('Converting GII files to Srf...');
                samsrf_disp(DataFiles');
                Srf = samsrf_gii2srf(DataFiles, [GuiSurf.Value{3} filesep dn(1:2)], GuiProc.Data.Normalise{1}, GuiProc.Data.Average{1}, true, '');
                samsrf_newline;
                Srf.Functional = DataFiles';
            end
        end
        
        % Run analysis!
        if isempty(Srf)
            % Should not really happen but just in case
            samsrf_clrscr;
            samsrf_disp('ERROR: Something went wrong with importing data!');
        else
            cf = pwd;
            cd(dp);
            switch Algorithm
                case 'samsrf_fit_prf'
                    % Forward-model pRF analysis
                    if Model.Polar_Search_Space
                        p = 2;
                    else
                        p = 1;
                    end
                    % Rescale search space definitions if necessary
                    if ~isfield(Model, 'SeedRoi') % Unless pRF-from-CF analysis
                        for i = p:length(Model.Param_Names)
                            if Model.Scaled_Param(i)
                                Model.(['Param' num2str(i)]) = Model.(['Param' num2str(i)]) * Model.Scaling_Factor;
                            end
                        end
                    end
                    % Now ready for analysis
                    OutFile = samsrf_fit_prf(Model, Srf, Roi);
                case 'samsrf_revcor_prf'
                    % Reverse-correlation pRF analysis
                    OutFile = samsrf_revcor_prf(Model, Srf, Roi);
                case 'samsrf_revcor_cf'
                    % Reverse-correlation CF analysis
                    OutFile = samsrf_revcor_cf(Model, Srf, Roi);
            end
            
            % Export GII/NII files?
            if GuiProc.Data.Export{1}
                load(OutFile); % Load maps  
                if strcmpi(de, '.nii')
                    % Export as NIIs 
                    samsrf_mat2vol(Srf, Model.Name);
                else
                    % Export as GIIs
                    if strcmpi(Srf.Hemisphere, 'bi')
                        % Explore both hemispheres
                        [Lsrf, Rsrf] = samsrf_hemi_srfs(Srf);
                        samsrf_export_giis(Lsrf, [Lsrf.Hemisphere '_' Model.Name]);
                        samsrf_export_giis(Rsrf, [Rsrf.Hemisphere '_' Model.Name]);
                    elseif strcmpi(Srf.Hemisphere, 'vol') 
                        % Export NII files
                    else
                        % Export single hemisphere
                        samsrf_export_giis(Srf, [Srf.Hemisphere '_' Model.Name]);
                    end
                end
            end
            cd(cf);
        end
    end
    
    %% Help on Files 
    function HelpOnFiles(source, eventdata)
        GuiInfo.Value = {'Help for data file list';
                         '';
                         'Use the Data menu to select data files to analyse.'; 
                         '';
                         'You can either select GIfTI surface files from one hemisphere, or you can choose to combine hemispheres. In that case, you must select the left hemisphere files (with the suffix "lh") and the tool will then combine each file with the corresponding right hemisphere file (with suffix "rh").';
                         ''
                         'Instead, you can also select volumetric NIfTI files. In this case, the analysis will be conducted in volume space. Any surface-based analysis options are obviously unavailable in this situation.'};
    end    

    %% Select file
    function FcnMatLoad(fwc, eventdata)
        [rn, rp] = uigetfile(fwc);
        if ~isscalar(rn)
            [rp, rn] = fileparts([rp rn]); % Remove extension
            GuiModel.Data.Value(eventdata.DisplaySelection(1)) = {[rp filesep rn]}; 
        end
    end

    %% Clear global variables when figure is closed
    function SamSrfCloseReq(src, evnt)
        samsrf_newline;
        samsrf_disp('SamSrf says:'); 
        samsrf_disp(' "I hope you''ll come back to be annoyed by me again soon...');
        samsrf_disp('  Mā te wā!"');
        samsrf_newline;
        pause(.5);

        clear global GuiInfo wb
        delete(src);
    end
end

