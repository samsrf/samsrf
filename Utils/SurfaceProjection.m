function SurfaceProjection
%
% ******* SamSrf surface projection tool ******* 
%
% Use this tool to project NII data to the cortical surfaces created by recon-all.
%
% 10/08/2018 - SamSrf 6 version (DSS)
% 12/09/2018 - Fixed bug with Rule switch (DSS)
% 19/09/2018 - Fixed bug with default cortical steps (DSS)
% 02/04/2020 - Added new line after default cortical sampling steps (DSS)
% 04/04/2020 - Added overlay for the functional data (DSS)
% 08/04/2020 - Now allows user to skip out at file selection & normalisation dialog (DSS)
%

%% Select paths 
CurrPath = pwd;
SubjPath = uigetdir('.', 'Select subject folder'); 
hemfolder = [SubjPath filesep 'surf' filesep]; % Surface folder

%% Check registration of white-matter surface
cd([SubjPath filesep 'label']);
strimg = samsrf_checkreg;
cd(CurrPath);

%% Welcome message
[vn vd] = samsrf_version; 
new_line;
disp('****************************************************************************');
disp('     Welcome to the Seriously Annoying MatLab Surfer Projection Tool!');
disp('    by D. S. Schwarzkopf from the University of Auckland, New Zealand');
new_line;
disp(['                 Version ' num2str(vn) ', Released on ' vd]);
disp('      (see SamSrf/ReadMe.txt for what is new in this version)');
disp('****************************************************************************');
new_line;

%% Select functional data files
[ef,ep] = uigetfile('*.nii', 'Select 4D-Nifty files', 'MultiSelect', 'on');
if isnumeric(ef) && ef == 0
    disp('No functional scans selected.');
    return
end
if ischar(ef)
    ef = {ef};
end
for i = 1:length(ef)
    % Add path & drop extension
    ef{i} = [ep ef{i}(1:end-4)];
end

%% Cortex sampling steps
ctxsteps = inputdlg('Cortex sampling steps?', '');
ctxsteps = eval(['[' cell2mat(ctxsteps) ']']); % Cortical sampling steps
if isempty(ctxsteps)
    disp('Choosing default cortex sampling of 0.5.');
    ctxsteps = 0.5;
    new_line;
end

%% Overlay for functional scans
% Get current figure handle
fig = gcf;
if ~verLessThan('matlab', '8.5')
    fig = fig.Number;
end
% Convert label into functional voxel mask
samsrf_label2nii([SubjPath filesep 'label' filesep 'lh.cortex'], ef{1}, strimg, [hemfolder 'lh'], ctxsteps);
samsrf_label2nii([SubjPath filesep 'label' filesep 'rh.cortex'], ef{1}, strimg, [hemfolder 'rh'], ctxsteps);
% Overlay functional voxels
spm_orthviews('RemoveBlobs', fig); % Remove previous overlay to avoid display error
spm_orthviews('AddColouredImage', fig, [SubjPath filesep 'label' filesep 'lh.cortex.nii'], [0 1 0]); % Left hemisphere
spm_orthviews('AddColouredImage', fig, [SubjPath filesep 'label' filesep 'rh.cortex.nii'], [0 1 0]); % Right hemisphere
spm_orthviews('Redraw'); % Redraw so the overlay appears

%% Normalise data?
nrmls = questdlg('Normalise time series?', '', 'Yes', 'No', 'Yes'); 
if isempty(nrmls)
    disp('Surface projection aborted by user.');
    return
end
nrmls = strcmpi(nrmls, 'Yes'); 

%% How to handle multiple runs?
if length(ef) > 1
    avrg = questdlg('How to handle multiple runs?', '', 'Average', 'Concatenate', 'Separate', 'Average'); 
else
    avrg = 'Separate';
end

%% Sampling rule to be used
RuleStrs = {'Mean' 'Median' 'Maximum' 'Minimum' 'Sum' 'Geomean' 'None'};
if length(ctxsteps) > 1
    rulenum = listdlg('ListString', RuleStrs, 'SelectionMode', 'single');
    if isempty(rulenum)
        error('No sampling rule specified!');
    end
else
    % If only one cortex sampling step, there is no point defining the rule
    rulenum = 1;
end
rule = RuleStrs{rulenum};
if strcmpi(rule, 'None')
    rule = '';
end

%% Loop thru hemispheres
hemsurf = {'lh' 'rh'};
for h = 1:2
    %% Convert NII to surface images
    if strcmpi(avrg, 'Average')
        % Average all the runs after projection
        samsrf_vol2srf(ef, strimg, [hemfolder hemsurf{h}], ctxsteps, rule, nrmls, true);
    elseif strcmpi(avrg, 'Concatenate')
        % Concatenate all the runs after projection
        samsrf_vol2srf(ef, strimg, [hemfolder hemsurf{h}], ctxsteps, rule, nrmls, false);
    else
        % Project each run separately
        for i = 1:length(ef)
            samsrf_vol2srf(ef{i}, strimg, [hemfolder hemsurf{h}], ctxsteps, rule, nrmls);
        end
    end
end

%% Return home
cd(CurrPath);
close all