function V = samsrf_checkreg_temp(ctxsmp)
%
% Creates a nifti binary mask of a cortical label. All files are selected 
% through dialogs. Use this to check whether your coregistration between
% the FreeSurfer files and SamSurfer is correct. This script directly calls 
% spm_check_registration and superimposes the binary mask of the label.
% The optional input ctxsmp defines cortical sampling steps (default=0.5).
%
% The files you want to select are as follows:
%
%   1. The label file you want to check (you could just pick the
%       lh.cortex.label or rh.cortex.label from the label folder).
%
%   2. The structural you used for reconstruction in FreeSurfer.
%
%   3. The template you want to use for creating the binary mask. 
%       This determines the space and dimensions of the mask, so if you
%       want to check that the general alignment of your structural is
%       correct, you want to chose the same structural as in 2 here. 
%       If you want to check that the functional data that you believe to 
%       be coregistered to your structural scan are truly in register, you 
%       want to select one of your functional scans here.
%
%   4. The surface folder of your subject. It's generally called 'surf' but
%       if you moved your files about it might be something else.
%   
%   5. Finally you need to decide which cortical hemisphere this is from.
%       For most things the toolbox can work it out on its own but in this
%       case it isn't that clever.
%
% In the SPM window you can then check if the overlaid mask actually falls 
% into the grey matter of your structural (see samsrf_checkreg.m for more 
% details). This is most easily seen when using the structural as a template 
% because this has the best resolution.
%
% 20/04/2022 - SamSrf 8 version (DSS)
%

if nargin == 0
    ctxsmp = 0.5;
end

[labf, labp] = uigetfile('*.label', 'Select label');
[strf, strp] = uigetfile('*.nii', 'Select structural');
[tmpf, tmpp] = uigetfile('*.nii', 'Select template');
labf = labf(1:end-6);
strf = strf(1:end-4);
tmpf = tmpf(1:end-4);
srfp = uigetdir('.', 'Select surface folder');
hem = questdlg('Hemisphere?','Question','Left','Right','Left');
if strcmpi(hem, 'Left')
    hem = 'lh';
elseif strcmpi(hem, 'Right')
    hem = 'rh';
end

% Convert label into nifty mask
samsrf_label2nii([labp labf], [tmpp tmpf], [strp strf], [srfp filesep hem], ctxsmp);
% Call check registration with the structural
if ~exist('spm', 'file')
    samsrf_error('Sorry but I need SPM to show you these brains :(');
end
spm_check_registration([strp strf '.nii']);

% Get current figure handle
fig = gcf;
if ~verLessThan('matlab', '8.5')
    fig = fig.Number;
end

% Superimpose the label mask
spm_orthviews('AddColouredImage', fig, [labp labf '.nii'], [1 0 0]);
% Redraw so the overlay appears
spm_orthviews('Redraw');
