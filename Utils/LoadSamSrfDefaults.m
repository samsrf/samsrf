function json = LoadSamSrfDefaults
%
% json = LoadSamSrfDefaults
% 
% Reads the JSON file containing the SamSrf default parameters as a struct.
%
% NOTE: This replaces the SamSrf_defaults.mat file used prior to SamSrf X!
%
% 19/09/2024 - Written (DSS)

if exist('SamSrf_defaults.mat', 'file')
    samsrf_disp('WARNING: SamSrf_defaults.mat found - this is obsolete!');
end

% Load JSON parameters
if exist('SamSrf_defaults.json', 'file')
    f = fopen('SamSrf_defaults.json');
    raw = fread(f,Inf);
    str = char(raw');
    json = jsondecode(str);
    fclose(f);
else
    samsrf_error('No SamSrf_default.json found!');
end