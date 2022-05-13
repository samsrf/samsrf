function SurfaceProjection
%
% ******* SamSrf surface projection tool ******* 
%
% Use this tool to project NII data to the cortical surfaces created by recon-all.
% This creates SamSrf surface data files using samsrf_vol2srf.
%
% You can -also- use this tool to convert MGH or GII files created by other tools into SamSrf files.
% If doing the latter, you will need to run this twice for files for each cortical hemisphere
% unless you are converting MGH/GII files separately (no averaging or concatenation).
%
% If using MGH/GII files, you must follow a strict naming convention where your files 
% for the left & right hemispheres are prefixed 'lh' & 'rh', respectively.
% 
% 20/04/2022 - SamSrf 8 version (DSS)
% 14/05/2022 - Added option to read GIfTI files (DSS)
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
[ef,ep] = uigetfile({'*.nii'; '*.mgh'; '*.gii'}, 'Select 4D-NIfTI/MGH/GII files', 'MultiSelect', 'on');
if isnumeric(ef) && ef == 0
    disp('No functional scans selected.');
    return
end
if ischar(ef)
    ef = {ef};
end
% Work out file type
[~,~,ft] = fileparts(ef{1});
ft = ft(2:end);
for i = 1:length(ef)
    % Add path & drop extension
    ef{i} = [ep ef{i}(1:end-4)];
end

%% Only for volume data
if strcmpi(ft,'nii')
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
end

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

%% Only for volume data
if strcmpi(ft,'nii')
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
end

%% Convert NII volumes to Srf files
if strcmpi(ft,'nii')
    %% Loop thru hemispheres
    hemsurf = {'lh' 'rh'};
    for h = 1:2
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
%% Convert MGH/GII surfaces to Srf files    
else
    if strcmpi(avrg, 'Average')
        % Average all the runs after projection
        [~,fn] = fileparts(ef{1});
        hemsurf = fn(1:2); % Use hemisphere as indicated by -first- file name!
        if strcmpi(ft,'mgh')
            samsrf_mgh2srf(ef, [hemfolder hemsurf], nrmls, true); % MGH files
        elseif strcmpi(ft,'gii')
            samsrf_gii2srf(ef, [hemfolder hemsurf], nrmls, true); % GIfTI files
        end
    elseif strcmpi(avrg, 'Concatenate')
        % Concatenate all the runs after projection
        [~,fn] = fileparts(ef{1});
        hemsurf = fn(1:2); % Use hemisphere as indicated by -first- file name!
        if strcmpi(ft,'mgh')            
            samsrf_mgh2srf(ef, [hemfolder hemsurf], nrmls, false); % MGH files
        elseif strcmpi(ft,'gii')
            samsrf_gii2srf(ef, [hemfolder hemsurf], nrmls, false); % GIfTI files
        end
    else
        % Project each run separately
        for i = 1:length(ef)
            [~,fn] = fileparts(ef{i});
            hemsurf = fn(1:2); % Use hemisphere as indicated by -this- file name!
            if strcmpi(ft,'mgh')            
                samsrf_mgh2srf(ef{i}, [hemfolder hemsurf], nrmls); % MGH files
            elseif strcmpi(ft,'gii')
                samsrf_gii2srf(ef{i}, [hemfolder hemsurf], nrmls); % GIfTI files
            end
        end
    end
end

%% Return home
cd(CurrPath);
close all