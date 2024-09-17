function strf = samsrf_checkreg(labpath)
%
% strf = samsrf_checkreg(labpath)
%
% This tool allow you to check the coregistration between FreeSurfer's
% mysterious coordinate system and what SamSrf believes it should be doing.
%
% labpath must point to the 'label' folder in your FreeSurfer subject's folder.
% You can omit this if the present directory is the 'label' folder.
%
% The tool creates nifti masks of the ?h.cortex.label files in the 'label' 
% folder of your subject. It then uses the 'surf' folder and the T1 scan
% you used for recon-all. This must be present as a NII file either in the 
% 'mri/orig' or in the 'mri' folder. It checks in 'mri/orig' first and then
% in 'mri'. It fails if there is more than one NII file present.
%
% The function returns the name of the structural NII file.
%
% It then loads an SPM window of the T1 scan overlaid with the two cortex
% labels. You will have to wait for it to finish loading and the overlay 
% will only appear once you click on the SPM window. You can use this overlay
% to check if the coregistration and the FreeSurfer reconstruction make sense. 
% The red surface should be well aligned with the grey-white matter boundary 
% and it shouldn't miss large parts of the cortex.
%
% IMPORTANT: You must run this script from inside the 'label' folder 
%            of the FreeSurfer subject you want to analyse!
%
% 20/04/2022 - SamSrf 8 version (DSS)
%

% Check whether you are in a 'label' folder
if nargin == 0
    labpath = pwd;
end
wd = labpath(end-4:end);
if ~strcmp(wd,'label')
    error('Not in a ''label'' directory!')
end

% Find the structural scan
fs = dir([labpath filesep '..' filesep 'mri' filesep 'orig' filesep '*.nii']);
if isempty(fs)
    fs = dir([labpath filesep '..' filesep 'mri' filesep '*.nii']);
    if isempty(fs)
        error('Cannot locate NII for the reconstructed T1!');
    end
    strf = [labpath filesep '..' filesep 'mri' filesep fs(1).name];
else
    strf = [labpath filesep '..' filesep 'mri' filesep 'orig' filesep fs(1).name];
end
if length(fs) > 1
    error('NII file in mri folder is ambiguous!');
end
samsrf_newline;
samsrf_disp(['Using ' fs(1).name ' as template.'])
samsrf_newline;

% Call check registration with the structural
if ~exist('spm', 'file')
    samsrf_disp('Sorry but I need SPM to show you these brains :(');
    return
end
spm_check_registration(strf);
strf = strf(1:end-4); % Truncate extension

% Get current figure handle
fig = gcf;
if ~verLessThan('matlab', '8.5')
    fig = fig.Number;
end
    
% Create volume masks
hem = {'lh' 'rh'};
for h = 1:2
    % Label file for this hemisphere
    labf = [labpath filesep hem{h} '.cortex'];
    % Convert label into nifti mask
    samsrf_label2nii(labf, strf, strf, [labpath filesep '..' filesep 'surf' filesep hem{h}], 0);
    % Superimpose the label mask
    spm_orthviews('AddColouredImage', fig, [labf '.nii'], [1 0 0]);
    % Redraw so the overlay appears
    spm_orthviews('Redraw');
end