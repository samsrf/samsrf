function samsrf_mat2vol(SrfName, OutFile)
% 
% samsrf_mat2vol(SrfName, [OutFile])
%
% Converts a volumetric SamSrf file back into a NII functional volume.
%
%   SrfName:    Either the name of a SamSrf file (without extension) 
%               or a Srf structure. If the latter, you must define OutFile.
%   OutFile:    Name of the files to be saved. If undefined, uses SrfName,
%               unless a Srf was provided as input, then this is mandatory.
%
% Saves each row in Srf.Data as [FileName]_Value.nii after removing the vol_ prefix. 
% [Filename] is either SrfName or OutFile.
%
% If Srf.Values doesn't exist it simply saves each row as _tr#.
%
% 15/09/2024 - Modified to better align with samsrf_export_giis (DSS)
%

% Load fake surface data
if ischar(SrfName)
    load(EnsurePath(SrfName));
    if nargin > 1 
        SrfName = OutFile;
    else
        SrfName = SrfName(5:end);
    end
else
    Srf = SrfName;
    SrfName = OutFile;
end
Srf = samsrf_expand_srf(Srf);

% Header information
hdr = Srf.NiiHeader(1); % Use first volume if there are more in header
% Volume image to be created
dim = hdr.dim;
% Remove header data info
hdr = rmfield(hdr, 'pinfo');

% Expand Srf.Data to contain whole image
if isfield(Srf, 'Roi') && ~isempty(Srf.Roi)
    RoiData = Srf.Data;
    Srf.Data = zeros(size(RoiData,1), prod(dim));
    % Loop thru volumes/time points
    for v = 1:size(Srf.Data,1)
        Srf.Data(v, Srf.Roi) = RoiData(v,:);
    end
end

% Enforce float32 data type
hdr.dt = [16 0]; 

% Loop thru volumes/time points
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
    hdr.fname = [SrfName '_' VolStr '.nii'];
    
    if exist('spm', 'file')
        spm_write_vol(hdr, img);
        samsrf_disp(['Saved ' hdr.fname ' as volumetric data.']);
    else
        samsrf_error('Sorry but I need SPM or NIfTI-patch to load NII files :(');
    end
end


