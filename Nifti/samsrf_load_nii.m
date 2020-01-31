function [hdr, img] = samsrf_load_nii(NiiFile)
%
% [hdr, img] = samsrf_load_nii(NiiFile)
%
% Reads in NIfTI file NiiFile (without extension) using the reader defined
% by def_nii_io in SamSrf_defaults, which is a case-invariant char:
%
%   'spm':      Uses SPM's NIfTI reading function (default)
%   'cbi':      Jonas Larsson's cbiNifti toolbox
%   'matlab':   MatLab's native NIfTI reading function (version >=9.3)
%
% *** If you have other suggestions for NIfTI loading/writing code to be included,
%     please contact Sam at s.schwarzkopf@auckland.ac.nz (no promises though... ;)
%
% Returns the header in hdr and the image data in img. 
% If img is undefined, no image data is loaded.
%
% 19/09/2019 - Started writing this function (DSS)
% 29/01/2020 - Added functionality for Larsson & Matlab methods (DSS)

%% Load default parameters
load('SamSrf_defaults.mat');
if ~exist('def_nii_io', 'var')
    def_nii_io = 'spm';
end

%% Add extension
NiiFile = [NiiFile '.nii'];

%% Load NII file
if strcmpi(def_nii_io, 'spm')
    if ~exist('spm_vol.m', 'file')
      error('SPM is not installed on the path!');
    end
    % Load header
    hdr = spm_vol(NiiFile);
    % Load image data also?
    if nargout > 1
      img = spm_read_vols(hdr);
    end
    % Nifti origin
    hdr.nii_orig = hdr.mat(1:3,4);

elseif strcmpi(def_nii_io, 'cbi')
    if ~exist('cbiReadNifti.m', 'file')
      error('Larsson''s cbiNifti toolbox is not installed on the path!');
    end
    % Load file 
    warning off % Some MatLab inconsistencies in that toolbox...
    if nargout == 1
      % Only load header
      hdr = cbiReadNiftiHeader(NiiFile);
    else
      % Load image data also
      [img,hdr] = cbiReadNifti(NiiFile);     
    end
    warning on
    hdr.mat = hdr.qform44;
    hdr.dim = hdr.dim; %%% STILL NEEDS TO BE DONE!!!
    % Nifti origin
    hdr.nii_orig = hdr.mat(1:3,4);
    hdr.nii_orig(hdr.nii_orig > 0) = hdr.nii_orig(hdr.nii_orig > 0) + 1;
    hdr.nii_orig(hdr.nii_orig < 0) = hdr.nii_orig(hdr.nii_orig < 0) - 1;    
    
elseif strcmpi(def_nii_io, 'matlab')
    if exist('OCTAVE_VERSION', 'builtin') || verLessThan('matlab', '9.3') 
      error('Native MatLab NIfTI reader not available!');
    end
    % Load header
    hdr = niftiinfo(strimg);
    hdr.mat = hdr.Transform.T';
    m = sum(abs(hdr.mat(:,1:3)), 2);
    hdr.mat(:,4) = hdr.mat(:,4) + m .* sign(hdr.mat(:,4));     
    hdr.dim = hdr.ImageSize;
    % Load image data?
    if nargout > 1
        img = niftiread(NiiFile);
    end
    % Nifti origin
    hdr.nii_orig = hdr.mat(1:3,4);
    hdr.nii_orig(hdr.nii_orig > 0) = hdr.nii_orig(hdr.nii_orig > 0) + 1;
    hdr.nii_orig(hdr.nii_orig < 0) = hdr.nii_orig(hdr.nii_orig < 0) - 1;  
    
end
