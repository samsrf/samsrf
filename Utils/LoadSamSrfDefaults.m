function json = LoadSamSrfDefaults
%
% json = LoadSamSrfDefaults
% 
% Reads the JSON file containing the SamSrf default parameters as a struct.
%
% NOTE: This replaces the SamSrf_defaults.mat file used prior to SamSrf X!
%
% 19/09/2024 - Written (DSS)
% 24/09/2024 - Adjustments to work both in deployed mode & command window (DSS)
% 07/10/2025 - Now also works for standalone functions outside SamSrfX GUI (DSS)

if exist('SamSrf_defaults.mat', 'file')
    samsrf_disp('WARNING: SamSrf_defaults.mat found - this is obsolete!');
end

% Default parameters filename
if isdeployed
    global SamSrfXPath
    if isempty(SamSrfXPath)
        SamSrfXPath = pwd;
    end
	SamSrfDefsName = [SamSrfXPath filesep 'SamSrf_defaults.json'];
else
    SamSrfDefsName = 'SamSrf_defaults.json';
end

% Load JSON parameters
if exist(SamSrfDefsName, 'file')
    f = fopen(SamSrfDefsName);
    raw = fread(f,Inf);
    str = char(raw');
    json = jsondecode(str);
    fclose(f);
else
    samsrf_error('No SamSrf_defaults.json found!');
end
