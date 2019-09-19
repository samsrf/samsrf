function [img, hdr] = samsrf_load_nii(NiiFile)
%
% [img, hdr] = samsrf_load_nii(NiiFile)
%
% Reads in NIfTI file NiiFile (without extension) using the reader defined
% by def_nii_io in SamSrf_defaults which is a case-invariant char:
%
%   'spm':      Uses SPM's NII reading function (default)
%   'cbi':      Jonas Larsson's cbiNifti toolbox
%   'matlab':   MATLAB's native NII reading function (version >=9.3)
%
% Returns the image data in img and the NII header in hdr.
%
% 19/09/2019 - Wrote this function (DSS)

%% Load default parameters
load('SamSrf_defaults.mat');
if ~exist('def_nii_io', 'var')
    def_nii_io = 'spm';
end

%% Add extension
NiiFile = [NiiFile '.nii'];

%% Load NII file
if strcmpi(def_nii_io, 'spm')
    hdr = spm_vol(NiiFile);
    img = spm_read_vols(hdr);
end
