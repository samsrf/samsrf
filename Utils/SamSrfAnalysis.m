function SamSrfAnalysis(Model)
% GUI for running SamSrf anlyses. Used for standalone app.

%% Welcome message
[vn,vd] = samsrf_version;
clc; 
disp('****************************************************************************');
disp('     Welcome to the Seriously Annoying MatLab Surfer Analysis Tool!');
disp('     by D. S. Schwarzkopf from the University of Auckland, New Zealand');
new_line;
disp(['                 Version ' num2str(vn) ', Released on ' vd]);
disp('      (see SamSrf/ReadMe.md for what is new in this version)');
disp('****************************************************************************');

Algorithm = 'samsrf_fit_prf';

%% Populate model defaults 
Model = samsrf_model_defaults(Algorithm, Model);

%% Create table for parameters?
if isfield(Model, 'Param_Names')
    % Ensure boolean fields exist
    if ~isfield(Model, 'Scaled_Param')
        Model.Scaled_Param = false(1,length(Model.Param_Names));
    end
    if ~isfield(Model, 'Only_Positive')
        Model.Only_Positive = false(1,length(Model.Only_Positive));
    end
    % Create table
    ParTab = table; 
    % Ensuring column vectors
    ParTab.Parameter = Model.Param_Names(:);
    ParTab.Scaled = logical(Model.Scaled_Param(:));
    ParTab.Positive = logical(Model.Only_Positive(:));
end

%% Create table for other parameters
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

%% Create GUI
GuiFig = uifigure('Name', ['SamSrfAnalysis GUI v' num2str(vn)], 'Units', 'normalized', 'Position', [0.1 0.15 0.8 0.75]); % Open large GUI 
GuiModelMenu = uimenu(GuiFig,'Text','Model'); % Create menu
GuiModelMenuLoad = uimenu(GuiModelMenu,'Text','Load'); % Load model
GuiModelMenuSave = uimenu(GuiModelMenu,'Text','Save'); % Save model
GuiDataMenu = uimenu(GuiFig,'Text','Data'); % Create menu
% Create grid layout
GuiGrid = uigridlayout(GuiFig, [3 4]); 
GuiGrid.RowHeight = {30 '1x' '2x'}; 
GuiGrid.ColumnWidth = {400 '1x' '1x' '1x'};
% Algorithm info
GuiAlgo = uitextarea(GuiGrid);
GuiAlgo.Layout.Row = 1;
GuiAlgo.Layout.Column = 1;
GuiAlgo.Editable = 'off';
% pRF parameter table
GuiPars = uitable(GuiGrid, 'Data', ParTab, 'ColumnEditable', [true true true]);
GuiPars.Layout.Row = 2;
GuiPars.Layout.Column = 1;
% Model field table
GuiModel = uitable(GuiGrid, 'Data', ModTab, 'ColumnEditable', [false true]);
GuiModel.Layout.Row = 3;
GuiModel.Layout.Column = 1;
% pRF profile
GuiPrf = uiaxes(GuiGrid);
GuiPrf.Layout.Row = [1 2];
GuiPrf.Layout.Column = 2;
axis(GuiPrf, 'off');
title(GuiPrf, 'pRF profile');
% HRF plot
GuiHrf = uiaxes(GuiGrid);
GuiHrf.YTick = [];
GuiHrf.Layout.Row = [1 2];
GuiHrf.Layout.Column = 3;
xlim(GuiHrf, [0 32]);
title(GuiHrf, 'HRF');
xlabel(GuiHrf, 'Time (s)');
ylabel(GuiHrf, 'Response (a.u.)');
% Model help text
GuiInfo = uitextarea(GuiGrid);
GuiInfo.Layout.Column = [2 3];
GuiInfo.Editable = 'off';
GuiModel.SelectionChangedFcn = @UpdateInfo;
% Data selection
GuiFiles = uilistbox(GuiGrid);
GuiFiles.Layout.Row = [1 3];
GuiFiles.Layout.Column = 4;
GuiFiles.Items = {'Select files to analyse'};
GuiFiles.Value = 'Select files to analyse';

clear ParTab ModTab


    %% Info display function
    function UpdateInfo(source, eventdata)
        % Update model help
        if isempty(eventdata.DisplaySelection)
            % Nothing selected 
        else
            % Display model help text for this field
            GuiInfo.Value = [source.Data.Field{eventdata.DisplaySelection(1)}; ''; string(ModelHelpText(Algorithm, source.Data.Field{eventdata.DisplaySelection(1)}))];
        end

        % Update HRF plot
        Hrf = GuiModel.Data.Value{strcmpi(GuiModel.Data.Field, 'Hrf')};
        TR = GuiModel.Data.Value{strcmpi(GuiModel.Data.Field, 'TR')};
        if isempty(Hrf)
            ExaHrf = samsrf_hrf(TR);
            ExaHrf = ExaHrf/max(ExaHrf);
            plot(GuiHrf, 0:TR:32, ExaHrf, 'k', 'LineWidth', 2);
            title(GuiHrf, 'HRF: Canonical de Haas');
            GuiHrf.YLim = [min(ExaHrf)-.1 1.1];
        elseif isnan(Hrf)
            ExaHrf = samsrf_doublegamma(TR, [6 16 1 1 6 0 32]);
            ExaHrf = ExaHrf/max(ExaHrf);
            plot(GuiHrf, 0:TR:32, ExaHrf, 'k', 'LineWidth', 2);
            title(GuiHrf, 'HRF: Canonical SPM');
            GuiHrf.YLim = [min(ExaHrf)-.1 1.1];
        elseif isinf(Hrf)
            ExaHrf = [samsrf_hrf(TR) samsrf_doublegamma(TR, [6 16 1 1 6 0 32]) samsrf_doublegamma(TR, [8 13 1 1 2 0 32])];
            ExaHrf = ExaHrf./repmat(max(ExaHrf),size(ExaHrf,1),1);
            plot(GuiHrf, 0:TR:32, ExaHrf, 'LineWidth', 2);
            title(GuiHrf, 'Fitting HRF on the fly');
            GuiHrf.YLim = [min(ExaHrf(:))-.1 1.1];
        elseif Hrf == 1
            plot(GuiHrf, ones(1,32), 'k', 'LineWidth', 2);
            title(GuiHrf, 'No HRF convolution');
            GuiHrf.YLim = [0 2];
        else 
            plot(GuiHrf, NaN(1,32), 'k', 'LineWidth', 2);
            title(GuiHrf, 'HRF undefined!');
            GuiHrf.YLim = [-1.1 1.1];
        end
    end

end
