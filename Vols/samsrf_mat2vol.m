function samsrf_mat2vol(SrfName)
% 
% samsrf_mat2vol(SrfName)
%
% Converts a volumetric SamSrf file back into a NII functional volume.
%
%   SrfName:    Name of SamSrf file (without extension).
%                 Saves each row in Srf.Data as SrfName_Value.nii after
%                 removing the vol_ prefix. If Srf.Values doesn't exist
%                 it simply saves each row as _tr#.
%
% 09/08/2018 - SamSrf 6 version (DSS)
%

% Load fake surface data
load(SrfName);

% Header information
hdr = Srf.NiiHeader;
% Volume image to be created
dim = hdr.dim;
% Remove header data info
hdr = rmfield(hdr, 'pinfo');

% Loop thru time points
for v = 1:size(Srf.Data,1)
    % Reshape volume image
    img = reshape(Srf.Data(v,:), dim);
    % What is this volume called?
    if isfield(Srf, 'Values')
        VolStr = Srf.Values{v};
    else
        VolStr = ['tr' num2str(v)];
    end
    % Save volume file
    hdr.fname = [SrfName(5:end) '_' VolStr '.nii'];
    spm_write_vol(hdr, img);
end

% Finished!
disp(['Saved ' SrfName ' as volumetric files.']);
