function samsrf_write_nii(NiiFile, hdr, img)
%
% samsrf_write_nii(NiiFile)
%
% Writes a NIfTI file NiiFile (without extension) using the reader defined
% by def_nii_io in SamSrf_defaults which is a case-invariant char:
%
%   'spm':      Uses SPM's NII reading function (default)
%   'cbi':      Jonas Larsson's cbiNifti toolbox
%   'matlab':   MATLAB's native NII reading function (version >=9.3)
%
% 20/09/2019 - Wrote this function (DSS)

%% Load default parameters
load('SamSrf_defaults.mat');
if ~exist('def_nii_io', 'var')
    def_nii_io = 'spm';
end

%% Add extension
NiiFile = [NiiFile '.nii'];

%% Load NII file
if strcmpi(def_nii_io, 'spm')
    hdr.fname = NiiFile;
    spm_write_vol(hdr, img);
end
