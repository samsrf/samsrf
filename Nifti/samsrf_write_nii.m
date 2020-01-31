function samsrf_write_nii(NiiFile, hdr, img)
%
% samsrf_write_nii(NiiFile)
%
% Writes a NIfTI file NiiFile (without extension) using the reader defined
% by def_nii_io in SamSrf_defaults which is a case-invariant char:
%
%   'spm':      Uses SPM's NIfTI writing function (default)
%   'cbi':      Jonas Larsson's cbiNifti toolbox
%   'matlab':   MATLAB's native NIfTI writing function (version >=9.3)
%
% 20/09/2019 - Started writing this function (DSS)
%

%% Load default parameters
load('SamSrf_defaults.mat');
if ~exist('def_nii_io', 'var')
    def_nii_io = 'spm';
end

%% Add extension
NiiFile = [NiiFile '.nii'];

%% Load NII file
if strcmpi(def_nii_io, 'spm')
    if ~exist('spm_write_vol.m', 'file')
      error('SPM is not installed on the path!');
    end
    hdr.fname = NiiFile;
    spm_write_vol(hdr, img);

elseif strcmpi(def_nii_io, 'cbi')
    if ~exist('cbiReadNifti.m', 'file')
      error('Larsson''s cbiNifti toolbox is not installed on the path!');
    end
    
elseif strcmpi(def_nii_io, 'matlab')    
    if verLessThan('matlab', '9.3')
      error('Native NIfTI reader not available in this MatLab version!');
    end

end
