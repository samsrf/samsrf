function ApsMat = MakeApsMat(ApsName)
%
% ApsMat = MakeApsMat(ApsName)
%
% Loads the apertures in ApsName & if they are not in MAT format, saves 
% them as a MAT file in SamSrf format. Returns the name of the MAT file.
% Note that this file is saved in the current working directory.
%
% Accepts the following file formats:
%
% Volumetric NII apertures must be formatted with x & y coordinates containing pixels,
% & the z-dimension containing individual frames.
%
% Animated GIF apertures should be formated as a movie of frames (self-explanatory?).
%
% 08/10/2025 - As it was written (DSS)

% File parts info
[~,n,e] = fileparts(ApsName);

% Convert NII?
if strcmpi(e,'.nii')
    if exist('spm', 'file') % Use SPM 
	    hdr = spm_vol(EnsurePath(ApsName));
        ApFrm = spm_read_vols(hdr);
        ApsMat = ['aps_' n];
        save(ApsMat, 'ApFrm');
        samsrf_disp('Volumetric NII apertures saved as:');
        samsrf_disp(['   ' ApsMat '.mat']);
    else 
        samsrf_error('Sorry but I need SPM or NIfTI-patch to load NII files :(');
    end
elseif strcmpi(e,'.gif')
    ApFrm = double(squeeze(imread(EnsurePath(ApsName), 'frames', 'all'))) / 255;
    ApsMat = ['aps_' n];
    save(ApsMat, 'ApFrm');
    samsrf_disp('Animated GIF apertures saved as:');
    samsrf_disp(['   ' ApsMat '.mat']);
else
    if strcmpi(e, '.mat') || isempty(e)
        % Already a MAT file 
        ApsMat = ApsName;
    else
        samsrf_error('Invalid aperture format!');
    end
end