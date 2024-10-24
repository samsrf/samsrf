function SamSrfX
% GUI for running SamSrf anlyses. Used for standalone app.

%% Version info
[vn,vd] = samsrf_version;

%% New Model 
[Model, Algorithm, xP] = Blank2dGModel;
Model = samsrf_model_defaults(Algorithm, Model); % Populate empty fields with defaults 

%% Create GUI
GuiFig = uifigure('Name', ['SamSrf X Analysis GUI v' num2str(vn)], 'Units', 'normalized', 'Position', [0.2 0.15 0.7 0.75]); % Open large GUI 
crfcn = @SamSrfCloseReq;
set(GuiFig, 'CloseRequestFcn', crfcn);
global SamSrfXPath GuiInfo wb % wb is for samsrf_progbar
SamSrfXPath = pwd; % Path to where SamSrfX is located

% Working folder
GuiFolderMenu = uimenu(GuiFig,'Text','Working Folder ~'); % Set current working directory
GuiFolderMenu.MenuSelectedFcn = @FcnChangeFolder;

% Data 
GuiDataMenu = uimenu(GuiFig,'Text','Data Files ~'); 
GuiDataHemis = uimenu(GuiDataMenu,'Text','Select Hemisphere Data'); % Select GII files from one hemisphere 
GuiDataHemis.MenuSelectedFcn = @FcnDataHemis;
GuiDataBilat = uimenu(GuiDataMenu,'Text','Select Bilateral Data'); % Select GII files for both hemispheres 
GuiDataBilat.MenuSelectedFcn = @FcnDataBilat;
GuiDataVols = uimenu(GuiDataMenu,'Text','Select Volume Data'); % Select NII volume files
GuiDataVols.MenuSelectedFcn = @FcnDataVols;

% Surf folder 
GuiSurfMenu = uimenu(GuiFig,'Text','Surf Folder ~'); 
GuiSurfSelect = uimenu(GuiSurfMenu,'Text','Select surf folder'); 
GuiSurfSelect.MenuSelectedFcn = @FcnSurfSelect;

% ROI
GuiRoiMenu = uimenu(GuiFig,'Text','Region of Interest ~'); 
GuiRoiSelect = uimenu(GuiRoiMenu,'Text','Select ROI Label'); % Select ROI
GuiRoiSelect.MenuSelectedFcn = @FcnRoiSelect; 
GuiOccRoi = uimenu(GuiRoiMenu,'Text','Create Occipital ROIs'); % Create occipital ROIs
GuiOccRoi.MenuSelectedFcn = @FcnOccRoi;

% Model 
GuiModelMenu = uimenu(GuiFig,'Text','Model Specification ~'); 
GuiModelNew = uimenu(GuiModelMenu,'Text','New Model'); % New model item
GuiModelNew.MenuSelectedFcn = @FcnModelNew;
GuiModelLoad = uimenu(GuiModelMenu,'Text','Load Model'); % Load model item
GuiModelLoad.MenuSelectedFcn = @FcnModelLoad;
GuiModelSave = uimenu(GuiModelMenu,'Text','Save Model'); % Save model item
GuiModelSave.MenuSelectedFcn = @FcnModelSave;
GuiModelInfo = uimenu(GuiModelMenu,'Text','Model Information'); % Display model info
GuiModelInfo.MenuSelectedFcn = @FcnModelInfo;

% Apertures
GuiApsMenu = uimenu(GuiFig,'Text','Apertures ~'); 
GuiApsLoad = uimenu(GuiApsMenu,'Text','Load Apertures'); % Select apertures
GuiApsLoad.MenuSelectedFcn = @FcnApsLoad;
GuiApsView = uimenu(GuiApsMenu,'Text','View Apertures'); % View apertures 
GuiApsView.MenuSelectedFcn = @FcnApsView;

% HRF
GuiHrfMenu = uimenu(GuiFig,'Text','Hemodynamic Response ~'); 
GuiHrfLoad = uimenu(GuiHrfMenu,'Text','Load HRF'); % Load prefit HRF 
GuiHrfLoad.MenuSelectedFcn = @FcnHrfLoad;

% Seed info for CF 
GuiCfMenu = uimenu(GuiFig,'Text','Connective Fields ~'); 
GuiCfMenu.Enable = 'off'; % Starts off inactive as undefined for standard 2D Gaussian
GuiSeedSelect = uimenu(GuiCfMenu,'Text','Select Seed ROI Label'); % Select seed ROI
GuiSeedSelect.MenuSelectedFcn = @FcnSeedSelect; 
GuiTempSelect = uimenu(GuiCfMenu,'Text','Select Template Map'); % Select template map
GuiTempSelect.MenuSelectedFcn = @FcnTempSelect; 

% Analyse
GuiAnalyseMenu = uimenu(GuiFig,'Text','Model Fitting ~'); 
GuiAnalyseMenu.MenuSelectedFcn = @CompletenessCheck;
GuiAnalyseRun = uimenu(GuiAnalyseMenu,'Text','Fit Model'); % Run analysis
GuiAnalyseRun.MenuSelectedFcn = @FcnAnalyseRun; 
GuiAnalyseRun.Enable = 'off'; % Inactive until all necessary field filled in!
GuiRepStr = uimenu(GuiAnalyseMenu,'Text','Replace String'); % Replace string
GuiRepStr.MenuSelectedFcn = @FcnRepStr; 
GuiRunBatch = uimenu(GuiAnalyseMenu,'Text','Run Batch Analysis'); % Run batch analysis
GuiRunBatch.MenuSelectedFcn = @FcnRunBatch; 

% Miscellaneous
GuiMiscMenu = uimenu(GuiFig,'Text','Miscellaneous'); 
GuiBensonMaps = uimenu(GuiMiscMenu,'Text','Convert Benson maps'); % Project Benson maps
GuiBensonMaps.MenuSelectedFcn = @FcnBensonMaps;
GuiTmp2Nat = uimenu(GuiMiscMenu,'Text','Warp Template > Native'); % Template 2 Native conversion
GuiTmp2Nat.MenuSelectedFcn = @FcnTmp2Nat;
GuiNat2Tmp = uimenu(GuiMiscMenu,'Text','Warp Native > Template'); % Native 2 Template conversion
GuiNat2Tmp.MenuSelectedFcn = @FcnNat2Tmp;
GuiDispMaps = uimenu(GuiMiscMenu,'Text','Map Display tool'); % Map display tool
GuiDispMaps.MenuSelectedFcn = @DisplayMaps;
GuiDelinTool = uimenu(GuiMiscMenu,'Text','Map Delineation Tool'); % Map delineation tool
GuiDelinTool.MenuSelectedFcn = @DelineationTool;
GuiEccenPlot = uimenu(GuiMiscMenu,'Text','Eccentricity Plots'); % Plot data by eccentricity
GuiEccenPlot.MenuSelectedFcn = @FcnPlotByEccen;
GuiVisFldCov = uimenu(GuiMiscMenu,'Text','Visual Field Coverage'); % Plot visual field coverage
GuiVisFldCov.MenuSelectedFcn = @FcnVisFldCov;


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
colormap(GuiPrf, samsrf_cmap('berlin'));
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
samsrf_disp('******************************************************************');
samsrf_disp(' Kia ora - this is Seriously Annoying MatLab Surfer!');
samsrf_disp(' by Sam Schwarzkopf at the University of Auckland, NZ');
samsrf_newline;
samsrf_disp([' Version ' num2str(vn) ', Released on ' vd]);
samsrf_disp('   (see SamSrf/ReadMe.md for what is new)');
samsrf_disp('******************************************************************');


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
                        FcnFileLoad('aps_*.mat', eventdata); % Select aperture file directly
                    end
                elseif eventdata.DisplaySelection(2) == 2 && (strcmpi(source.Data.Field{eventdata.DisplaySelection(1)}, 'Seed_Fine_Fit') || strcmpi(source.Data.Field{eventdata.DisplaySelection(1)}, 'Template')) 
                    FcnFileLoad('*.mat', eventdata); % Select Srf map file directly
                elseif eventdata.DisplaySelection(2) == 2 && strcmpi(source.Data.Field{eventdata.DisplaySelection(1)}, 'SeedRoi') 
                    FcnFileLoad('*.label', eventdata); % Select label file directly
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
        if sum(strcmpi(GuiModel.Data.Field, 'Scaling_Factor')) > 0 && isnan(GuiModel.Data.Value{strcmpi(GuiModel.Data.Field, 'Scaling_Factor')})
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
        CompCheckInfo{cl,1} = ['Current working folder: ' pwd];
        cl = cl + 1;
        CompCheckInfo{cl,1} = '';
        cl = cl + 1;

        % ROI undefined
        if strcmpi(GuiRoi.Value{3}, '< None selected >')
            CompCheckInfo{cl,1} = 'Warning: No Region of Interest selected';
            cl = cl + 1;
            CompCheckInfo{cl,1} = '';
            cl = cl + 1;
        end

        % Scaling factor undefined
        if ~strcmpi(Algorithm, 'samsrf_revcor_cf') && isnan(GuiModel.Data.Value{strcmpi(GuiModel.Data.Field, 'Scaling_Factor')})
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
        if cl == 5
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


    %% Change working directory
    function FcnChangeFolder(src,event)
        wp = uigetdir('.', 'Select working folder');
        figure(GuiFig);
        if ~isscalar(wp)
            cd(wp);
            samsrf_newline;
            samsrf_disp(['Changed working folder to: ' wp]);
        end
    end

    %% Blank 2D Gaussian model
    function [Model, Algorithm, xP] = Blank2dGModel
        Model.Name = 'pRF_Gaussian'; 
        Model.Prf_Function = @(P,ApWidth) prf_gaussian_rf(P(1), P(2), P(3), ApWidth); 
        Model.Param_Names = {'x0' 'y0' 'Sigma'}; 
        Model.Scaled_Param = [1 1 1];  
        Model.Only_Positive = [0 0 1]; 
        Model.Scaling_Factor = NaN; 
        Model.TR = 1; 
        Model.Hrf = 0; 
        Model.Aperture_File = ''; 
        Model.Polar_Search_Space = true; 
        Model.Param1 = 0:10:350; 
        Model.Param2 = 2.^(-5:.2:.6); 
        Model.Param3 = 2.^(-5.6:.2:1); 
        Algorithm = 'samsrf_fit_prf'; 
        xP = [0.4 0.5 0.2];
    end

    %% New Model file
    function FcnModelNew(srf,event)
        [Model, Algorithm, xP] = Blank2dGModel;
        Model = samsrf_model_defaults(Algorithm, Model); % Populate empty fields with defaults 
        GuiAlgo.Value = {'Standard_2D_Gaussian_pRF'; ''; ['Algorithm: ' Algorithm]};
        [GuiPars, GuiModel] = UpdateTables(Model);
        GuiModel.SelectionChangedFcn = @UpdateInfo;
        GuiModel.CellEditCallback = @UpdateInfo;
        GuiPars.SelectionChangedFcn = @UpdateInfo;
        eventdata.EventName = 'SelectionChanged';
        eventdata.DisplaySelection = [1 1];
        UpdateInfo(GuiModel, eventdata);        
        figure(GuiFig);
    end

    %% Load Model file
    function FcnModelLoad(src,event)
        [fn, pn] = uigetfile('*.mod');
        if ~isscalar(fn)
            vs = whos('-file', [pn fn]);
            vs = {vs.name}';
            if sum(strcmpi(vs, 'Algorithm')) == 0 || sum(strcmpi(vs, 'xP')) == 0
                samsrf_error([pn fn ' was not saved by SamSrfX or has been corrupted!']);
            else
                warning off
                load([pn fn], '-mat', 'Algorithm', 'Model', 'xP', 'DataFiles', 'SurfFolder', 'Roi', 'Avg', 'Nrm', 'Exp');
                warning on
                Model = samsrf_model_defaults(Algorithm, Model); % Populate empty fields with defaults 
                % If this is pRF-from-CF analysis, automatically set fields 
                if strcmpi(Algorithm, 'samsrf_fit_prf') && isfield(Model, 'SeedRoi')
                    Model.Aperture_File = '[Set by analysis]';
                    Model.Scaling_Factor = Inf;
                end

                % If Model Spec contains file info
                samsrf_newline;
                if exist('DataFiles', 'var')
                    samsrf_disp('Saved Model contains data file & ROI info...');
                    GuiFiles.Items = DataFiles;
                    GuiSurf.Value{3} = SurfFolder;
                    GuiRoi.Value{3} = Roi;                    
                    cd(pn);
                    samsrf_newline;
                    samsrf_disp(['Changed working folder to: ' pn]);
                else
                    samsrf_disp('Predefined Model only contains parameters...');
                end
                
                % If Model Spec contained flags
                if exist('Avg', 'var')
                    GuiProc.Data.Average{1} = Avg;
                    GuiProc.Data.Normalise{1} = Nrm;
                    GuiProc.Data.Export{1} = Exp;
                end

                GuiAlgo.Value = {fn(1:end-4) ; ''; ['Algorithm: ' Algorithm]};
                [GuiPars, GuiModel] = UpdateTables(Model);
                GuiModel.SelectionChangedFcn = @UpdateInfo;
                GuiModel.CellEditCallback = @UpdateInfo;
                GuiPars.SelectionChangedFcn = @UpdateInfo;
                eventdata.EventName = 'SelectionChanged';
                eventdata.DisplaySelection = [1 1];
                UpdateInfo(GuiModel, eventdata);        
                figure(GuiFig);
            end
        end
        FcnModelInfo([],[]);
        figure(GuiFig);
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
        [fn, pn] = uiputfile('*.mod');
        if ~isscalar(fn)
            DataFiles = GuiFiles.Items;
            SurfFolder = GuiSurf.Value{3};
            Roi = GuiRoi.Value{3};
            Avg = GuiProc.Data.Average{1};
            Nrm = GuiProc.Data.Normalise{1};
            Exp = GuiProc.Data.Export{1};
            save([pn fn], '-mat', 'Algorithm', 'Model', 'xP', 'DataFiles', 'SurfFolder', 'Roi', 'Avg', 'Nrm', 'Exp');
            GuiAlgo.Value{1} = fn(1:end-4);
            figure(GuiFig);
        end
    end

    %% Select aperture file
    function FcnApsLoad(src,event)
        [rn, rp] = uigetfile('aps_*.mat');
        if ~isscalar(rn)
            GuiModel.Data.Value(strcmpi(GuiModel.Data.Field, 'Aperture_File')) = {[rp rn(1:end-4)]}; 
            UpdateInfo(GuiModel, eventdata);
        end
        figure(GuiFig);
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
        figure(GuiFig);
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
        figure(GuiFig);
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
        figure(GuiFig);
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
        figure(GuiFig);
    end

    %% Select surf folder
    function FcnSurfSelect(src,event)
        sp = uigetdir('../surf', 'Select surf folder');
        figure(GuiFig);
        if ~isscalar(sp)
            GuiSurf.Value = {'Subject surf folder:'; ''; sp}; 
        end
        figure(GuiFig);
    end

    %% Select ROI label
    function FcnRoiSelect(src,event)
        [rn, rp] = uigetfile({'*.label', 'FreeSurfer label'; '*.nii', 'Binary NIfTI mask'}, 'Select Region of Interest');
        if ~isscalar(rn)
            GuiRoi.Value = {'Region of Interest:'; ''; [rp rn]}; 
            UpdateInfo(GuiModel, eventdata);
        end
        figure(GuiFig);
    end

    %% Create occipital ROI
    function FcnOccRoi(src,event)
        if ~strcmpi(GuiSurf.Value{3}, '< None selected >') && ~strcmpi(GuiSurf.Value{3}, '< Not required for NII >')
            samsrf_clrscr;
            samsrf_disp('Creating occipital ROIs...');
            warning off
            MakeOccRoi(GuiSurf.Value{3})
            warning on
            samsrf_disp(' Occipital ROIs created.');
            samsrf_done;
        elseif strcmpi(GuiSurf.Value{3}, '< None selected >')
            samsrf_clrscr;
            samsrf_disp('ERROR: No surf folder selected!');
        end
        figure(GuiFig);
    end

    %% Select seed ROI label
    function FcnSeedSelect(src,event)
        [sn, sp] = uigetfile({'*.label', 'FreeSurfer label'}, 'Select seed ROI');
        if ~isscalar(sn)
            [sp, sn] = fileparts([sp sn]); % Remove extension
            GuiModel.Data.Value{strcmpi(GuiModel.Data.Field, 'SeedRoi')} = [sp filesep sn]; 
            UpdateInfo(GuiModel, eventdata);
        end
        figure(GuiFig);
    end

    %% Select template map 
    function FcnTempSelect(src,event)
        [tn, tp] = uigetfile('*.mat', 'Select template map');
        if ~isscalar(tn)
            [tp, tn] = fileparts([tp tn]); % Remove extension
            GuiModel.Data.Value{strcmpi(GuiModel.Data.Field, 'Template')} = [tp filesep tn]; 
            UpdateInfo(GuiModel, eventdata);
        end
        figure(GuiFig);
    end

    %% Replace string
    function FcnRepStr(src,event)
        % If strings provided?
        if ischar(src)
            Old = src;
            New = event;
        else
            % Which string to replace?
            Old = cell2mat(inputdlg('Find string:', 'Old string'));
        end

        if ~isempty(Old)
            wp = pwd;
            % Check if there is anything to replace
            if ~contains(wp, Old)
                samsrf_error([Old ' does not appear in working path!']);
            end
            % What to replace with?
            if ~ischar(New)
                New = cell2mat(inputdlg('Replace with:', 'New string'));
            end
            if ~isempty(New)    
                wp = strrep(wp, Old, New);
                if ~exist(wp, 'dir') 
                    samsrf_error([New ' does not exist!']);
                end
                cd(wp);
                % Replace all instance with wildcard
                for i = 1:length(GuiFiles.Items)
                    GuiFiles.Items{i} = strrep(GuiFiles.Items{i}, Old, New);
                end
                % Replace in surf folder
                GuiSurf.Value{3} = strrep(GuiSurf.Value{3}, Old, New);
                % Replace in ROI label
                GuiRoi.Value{3} = strrep(GuiRoi.Value{3}, Old, New);
                % Replace in all strings in Model
                Model = UpdateModel(Model);
                Fs = fieldnames(Model);
                for i = 1:length(Fs)
                    if ischar(Model.(Fs{i})) 
                        Model.(Fs{i}) = strrep(Model.(Fs{i}), Old, New);
                    end
                end
                [GuiPars, GuiModel] = UpdateTables(Model);
                UpdateInfo(GuiModel, eventdata);
                figure(GuiFig);
                samsrf_clrscr;
                samsrf_disp(['Replaced all instances of "' Old '" with "' New '".']);
                samsrf_newline;
                samsrf_disp(['Working folder is now: ' wp]);
            end
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
                samsrf_done;
            end
        end
    end
    
    %% Run analysis batch
    function FcnRunBatch(src,event)
        % Current working folder
        wp = pwd;
        % Select main folder
        bp = uigetdir('../..', 'Select main folder');
        if ~isscalar(bp)
            % Select subjects
            DataBatch = dir([bp filesep '*.']);
            DataBatch = {DataBatch.name}';
            x = listdlg('ListString', DataBatch, 'PromptString', 'Select subjects to analyse');
            if ~isempty(x)
                % Selected subjects but ignore . & ..
                DataBatch = DataBatch(x);
                DataBatch = DataBatch(~strcmpi(DataBatch, '.') & ~strcmpi(DataBatch, '..'));

                % Current subject ID
                CurSub = '';
                for i = 1:length(DataBatch)
                    if contains(wp, DataBatch{i})
                        CurSub = DataBatch{i};
                    end
                end

                % Ready to run batch?
                if isempty(CurSub)
                    samsrf_error('Cannot find this subject in your batch! Are you sure you''re in the correct working folder?');
                else
                    % Open batch window
                    Bf = uifigure('Name', 'SamSrfX Batch Analysis');
                    Bf.Position(3) = 300; % Narrower than normal
                    Bg = uigridlayout(Bf, [1 1]);
                    Bt = table;
                    Bt.Subject = DataBatch;     
                    Bt = uitable(Bg, 'Data', Bt, 'ColumnEditable', false);
                    Bt.Enable = 'inactive';
    
                    % Loop thru subjects
                    for i = 1:length(DataBatch)
                        % Replace string with current subject
                        FcnRepStr(CurSub, DataBatch{i});
                        Bt.Selection = [i 1];
                        CurSub = DataBatch{i};
                        figure(Bf);
                        pause(.1);
                        % Analyse this subject
                        FcnAnalyseRun([],[]);
                        % Record as completed
                        Bt.Data.Subject{i} = [Bt.Data.Subject{i} ' completed'];
                        figure(Bf);
                    end
                    samsrf_newline;
                    samsrf_disp(['Completed batch analysis of ' num2str(length(DataBatch)) ' subjects!']);
                end
            end
        end
    end

    %% Create Benson maps
    function FcnBensonMaps(src,event)
        bp = uigetdir('../benson', 'Select folder with Benson GII files');
        figure(GuiFig);
        if ~isscalar(bp) 
            if strcmpi(GuiSurf.Value{3}, '< None selected >')
                samsrf_error('Surf folder not defined!');
            else
                cf = pwd;
                cd(bp);
                samsrf_clrscr;

                % Left hemisphere 
                if exist('lh_benson.gii', 'file')
                    samsrf_benson2srf('lh_benson', GuiSurf.Value{3});
                else
                    samsrf_error('lh_benson.gii does not exist!');
                end
                samsrf_newline;
                % Right hemisphere 
                if exist('rh_benson.gii', 'file')
                    samsrf_benson2srf('rh_benson', GuiSurf.Value{3});
                else
                    samsrf_error('rh_benson.gii does not exist!');
                end
                % Combine hemispheres
                samsrf_newline;
                samsrf_disp('Combining hemispheres...');
                L = load('lh_benson.mat');
                R = load('rh_benson.mat');
                Srf = samsrf_bilat_srf(L.Srf, R.Srf);
                save('bi_benson', 'Srf');
                samsrf_disp('Saved bi_benson.mat');
                cd ROIs_Benson
                samsrf_newline;
                samsrf_bilat_label(Srf, 'V1');
                samsrf_bilat_label(Srf, 'V2');
                samsrf_bilat_label(Srf, 'V3');
                % Clear & go back
                clear L R Srf
                cd(cf);
                samsrf_done;
            end
        end
    end

    %% Template to Native conversion
    function FcnTmp2Nat(src,event)
        if strcmpi(GuiSurf.Value{3}, '< None selected >')
            samsrf_error('Surf folder not defined!');
        else
            [nn, np] = uigetfile('*.mat', 'Select native data');
            tp = uigetdir('*.mat', 'Select template map folder');            
            figure(GuiFig);
            if ~isscalar(nn) 
                samsrf_clrscr;
                cf = cd;
                cd(np);
                Template2NativeMap(nn, GuiSurf.Value{3}, tp);
                cd(cf);
                samsrf_done;
            end
        end
    end

    %% Native to Template conversion
    function FcnNat2Tmp(src,event)
        if strcmpi(GuiSurf.Value{3}, '< None selected >')
            samsrf_error('Surf folder not defined!');
        else
            [nn, np] = uigetfile('*.mat', 'Select native data');            
            figure(GuiFig);
            if ~isscalar(nn) 
                samsrf_clrscr;
                SamSrfDefs = LoadSamSrfDefaults;
                if isfield(SamSrfDefs, 'def_fsaverage') && ~isempty(SamSrfDefs.def_fsaverage)                
                    cf = cd;
                    cd(np);
                    Native2TemplateMap([np nn(1:end-4)], GuiSurf.Value{3}, [SamSrfDefs.def_fsaverage filesep 'surf']);
                    cd(cf);
                    samsrf_done;
                else
                    samsrf_error('No template brain defined!');
                end
            end
        end
    end

    %% Plot by eccentricity
    function FcnPlotByEccen(srf,event)
        [mn, mp] = uigetfile('*.mat', 'Select map file');
        if ~isscalar(mn)
            % Load map
            load([mp mn]);
            Srf = samsrf_expand_srf(Srf);
            
            % Select ROIs
            [Rois, rp] = uigetfile('*.label', 'Select ROI labels', 'MultiSelect', 'on');
            if ischar(Rois)
                Rois = {Rois(1:end-6)};
            elseif isscalar(Rois)
                Rois = {''};
            end
            if length(Rois) > 1
                Colours = lines(length(Rois));
                for r = 1:length(Rois)
                    Rois{r} = Rois{r}(1:end-6); % Remove label extension
                end
            else
                Colours = [0 0 0];
            end
            
            % Define R^2 threshold
            Thrsh = inputdlg('R^2 threshold: ', 'Select R^2 threshold');
            Thrsh = str2double(cell2mat(Thrsh));
            if isnan(Thrsh)
                samsrf_error('Invalid R^2 threshold entered!');
                Thrsh = 0;
            end
           
            % Select data field
            ValNum = listdlg('ListString', Srf.Values, 'SelectionMode', 'single', 'PromptString', 'Which subject?');
            if isempty(ValNum)
                samsrf_disp('You must pick a field! Selecting 1st one...');
                ValNum = 1;
            end
            
            % Define bins
            Bins = inputdlg('Eccentricity bins: ', 'Select eccentricity bins');
            try
                Bins = eval(['[' cell2mat(Bins) ']']);
            catch
                samsrf_error('Invalid eccentricity bins entered!');
                Bins = [0 Inf];
            end

            % Loop thru ROIs
            figure('Name', [Srf.Values{ValNum} ' by eccentricity']);
            hs = [];
            for r = 1:length(Rois)
                T = table;
                if length(Bins) > 2
                    % Binned medians
                    [Res,h] = samsrf_plot(Srf, Srf.Values{ValNum}, Srf, 'Eccentricity', Bins, [rp Rois{r}], Thrsh, 'Median', Colours(r,:));
                    T.Eccentricity_Bin = Res(:,1);
                    T.(['Median_' Srf.Values{ValNum}]) = Res(:,2);
                    T.Lower_Bound = Res(:,3);
                    T.Upper_Bound = Res(:,4);
                    T.Num_Data_Pts = Res(:,5);
                else
                    % Scatter plot
                    [Res,h] = samsrf_plot(Srf, Srf.Values{ValNum}, Srf, 'Eccentricity', Bins, [rp Rois{r}], Thrsh, 'Scatter', Colours(r,:));
                    T.Eccentricity = Res(:,1);
                    T.(Srf.Values{ValNum}) = Res(:,2);
                end
                hs = [hs h];
                % Save results
                writetable(T, [mn(1:end-4) '_' Srf.Values{ValNum} '-vs-Eccentricity_' Rois{r} '_R^2=' num2str(Thrsh*100) '%.csv']);
            end
            title([Srf.Values{ValNum} ' by eccentricity']);
            legend(hs, Rois);
        end        
    end

    %% Visual field coverage plots
    function FcnVisFldCov(src,event)
        [mn, mp] = uigetfile('*.mat', 'Select map file');
        if ~isscalar(mn)
            % Load map
            load([mp mn]);
            Srf = samsrf_expand_srf(Srf);
            
            % Select ROIs
            [Rois, rp] = uigetfile('*.label', 'Select ROI labels', 'MultiSelect', 'on');
            if isscalar(Rois)
                samsrf_error('You must select a ROI!');
            elseif ischar(Rois)
                Rois = {Rois(1:end-6)};
            end
            if length(Rois) > 1
                for r = 1:length(Rois)
                    Rois{r} = Rois{r}(1:end-6); % Remove label extension
                end
            end

            % Define R^2 threshold
            Thrsh = inputdlg('R^2 threshold: ', 'Select R^2 threshold');
            Thrsh = str2double(cell2mat(Thrsh));
            if isnan(Thrsh)
                samsrf_error('Invalid R^2 threshold entered!');
                Thrsh = 0;
            end
            
            % Loop thru ROIs
            for r = 1:length(Rois)
                figure('Name', ['Visual Field Coverage - ' Rois{r}]);
                samsrf_vfcoverage(Srf, Model.Scaling_Factor, [rp Rois{r}], Thrsh);
                axis([-1 1 -1 1] * Model.Scaling_Factor);
                colormap(samsrf_cmap('batlow'));
            end
        end
    end

    %% Model information
    function FcnModelInfo(src,event)
        if strcmpi(Algorithm, 'samsrf_fit_prf') && ~isfield(Model, 'SeedRoi')          
            GuiInfo.Value = {'Forward-modelling pRF analysis'; ''; 
                             'The standard approach for pRF analysis based on the original description by Dumoulin & Wandell, 2008. See the Introduction to pRF document for more detail.'; ''; 
                             'Essentially, the model works by predicting the neural response of an assumed pRF model to the stimuli as defined in your apertures file. This is then usually convolved with an HRF (which may be predefined or fit concurrently) and may also involve a compressive nonlinearity of the response.'; ''; 
                             'The model then attempts to find the pRF parameters that provide the best fit of this predicted time series with the actually measured response. This can be done via a coarse-to-fine fitting approach or by initialising the fit using a map, e.g. a coarse-fit you already ran before, a group-average atlas map, or the maps predicted by the DeepRetinotopy model.'; '';
                             'This approach probably yields the best estimates of pRF size, that is, how selective a given voxel is to stimulus position. However, this is based on a lot of assumptions (such as the shape and parameters of your pRF model) and can take a long time to fit.'};

        elseif strcmpi(Algorithm, 'samsrf_revcor_prf')
            GuiInfo.Value = {'Reverse correlation pRF analysis'; ''; 
                             'This approach uses reverse correlation to fit a pRF like in our study Urale et al., 2022, Human Brain Mapping.'; ''; 
                             'The analysis convolves the stimulus apertures with a predefined HRF and the correlates these convolved apertures with the actual measured response for each voxel. This yields a correlation profile of where in the stimulus space the voxel responded. By default, this only quantifies basic geometric properties of this profile but you can also fit a pRF model to this profile.'; '';
                             'This approach makes far fewer assumptions and is -much- faster than forward-modelling pRF analysis. In our blind spot mapping study we found that it is also much more accurate in delineating scotomas. This makes sense because it captures irreguar shapes of responses without forcing a pRF model on them (unless you deliberately fit a pRF model to the reverse correlation profiles of course).'; '';
                             'Estimates of pRF size are not directly comparable to those from forward-modelling analysis especially when you use a stimulus design with strong spatio-temporal correlations (i.e. most people).'};

        elseif strcmpi(Algorithm, 'samsrf_revcor_cf')
            GuiInfo.Value = {'Reverse correlation CF analysis'; ''; 
                             'This approach uses reverse correlation to fit a connective field like in our study Tangtartharakul et al., 2023, Human Brain Mapping.'; ''; 
                             'The analysis takes activity in a seed ROI (such as V1) and correlates this with the measured time series for a voxel. This yields a correlation profiles across the seed ROI. Using a template retinotopic map (e.g. the Benson templates, DeepRetinotopy prediction, or an empirical pRF map) these correlations are then projected into stimulus space to produce a correlation profile (e.g. a reflection of the CF in the visual field). By default, only basic geometric properties of this profile are estimated, although you can also fit a 2D Gaussian model to these.'; '';
                             'Estimates of the CF size are not directly comparable to pRF size from forward-modelling analysis. They reflect different things: while pRFs show how a voxel''s response is modulated as the stimulus traverse space, the CF reflects how neighbouring voxels coactivate. There is clearly a relationship between these measures but they are no equivalent.'};

        elseif strcmpi(Algorithm, 'samsrf_fit_prf') && isfield(Model, 'SeedRoi')          
            GuiInfo.Value = {'Deriving forward-modelling pRF estimates from CF profiles'; ''; 
                             'This approach uses standard forward-modelling pRF analysis to fit CF models directly in stimulus space.'; ''; 
                             'The analysis first uses a template retinotopic map (e.g. the Benson templates, DeepRetinotopy prediction, or an empirical pRF map) of a seed ROI (such as V1) to project brain activity during your runs into the stimulus space. These backprojections are then used in the same way stimulus apertures are used in normal forward-modelling pRF analysis.'; '';
                             'This approach may yield better guestimates of pRF size from CF data than reverse correlation CF analysis but the caveats for CF analysis also apply here (anyway, watch this space - this analysis is still experimental!). Like forward-modelling pRF analysis it is based on a lot of assumptions (such as the shape and parameters of your pRF model) and due to the amount of data and noisy responses it typically takes a long time to fit.'};

        else
            GuiInfo.Value = {'Something is amiss! I do not know anything about this algorithm!'}; 
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
    function FcnFileLoad(fwc, eventdata)
        [rn, rp] = uigetfile(fwc);
        figure(GuiFig);
        if ~isscalar(rn)
            [rp, rn] = fileparts([rp rn]); % Remove extension
            GuiModel.Data.Value(eventdata.DisplaySelection(1)) = {[rp filesep rn]}; 
        end
    end

    %% Clear global variables when figure is closed
    function SamSrfCloseReq(src, evnt)
        samsrf_clrscr;
        samsrf_newline;
        samsrf_disp('SamSrf says:'); 
        samsrf_disp(' "I hope you''ll come back to be annoyed by me again soon...');
        samsrf_disp('  Mā te wā!"');
        samsrf_newline;
        pause(1);

        clear global SamSrfXPath GuiInfo wb
        delete(src);
    end
end

