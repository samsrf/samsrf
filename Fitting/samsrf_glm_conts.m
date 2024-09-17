function samsrf_glm_conts(SrfFile, Conts, ContNames)
%
% samsrf_glm_conts(SrfFile, Conts, ContNames)
%
% Calculates GLM contrasts on the GLM surface data file defined by SrfFile.
% The contrasts are defined in the matrix Conts. Each row is a different
% contrast vector. ContNames is a cell string defining the names of each
% contrast. The contrasts are then saved in a file name called the same as
% SrfFile by appended with '_conts'.
%
% Srf.Data in the saved file contains the t-map for each contrast per row.
% Srf.Values contains the names of the contrasts. The Srf contains a field
% Crit_Ts that defines the critical t-statistics for p=0.05, p=0.001, and
% p=0.05 (Bonferroni corrected).
%
% 20/04/2022 - SamSrf 8 version (DSS)
%

% Convert ContNames into cell if it's a string
if ischar(ContNames)
    ContNames = {ContNames};
end

% Ensure contrasts are in rows
if size(Conts,1) > 1 && size(Conts,2) == 1
    Conts = Conts';
end

% Load data
load(EnsurePath(SrfFile));
[Srf,vx] = samsrf_expand_srf(Srf);
bSrf = Srf;
Srf.Data = [];
Srf.Values = {};

% Add version number
Srf.Version = samsrf_version;

% Check that contrasts are correct
if size(Conts,2) ~= size(bSrf.Data,1)-1
    samsrf_error('Contrast contains not the right number of columns!');
end

% Critical t-statistics
nver = size(bSrf.Data,2); % Number of vertices
if isfield(bSrf, 'Y')
    Srf.df = size(bSrf.Y,1) - size(bSrf.X_glm,2) - 1; % Degrees of freedom
    Srf.Crit_Ts = tinv(1 - [0.05 0.001 0.05/nver]/2, Srf.df); % Critical T's (doesn't work for z-scored data)
    Srf = rmfield(Srf, 'Y'); % Remove time courses
end
Srf = rmfield(Srf, 'X_glm'); % Remove regressors
if isfield(Srf, 'Raw_Data') 
    Srf = rmfield(Srf, 'Raw_Data'); % Remove unsmoothed data 
end

% Calculate contrasts
Resid = bSrf.Data(end,:); % Residuals
bSrf.Data = bSrf.Data(1:end-1,:);
for c = 1:size(Conts,1)
    curcon = Conts(c,:)'; % Current contrast
    % Calculate difference of betas
    D = bSrf.Data .* repmat(curcon, 1, nver);
    D = sum(D); 
    T = D ./ Resid;  % Calculate t-ratio of contrast over error
    Srf.Data = [Srf.Data; T];
    Srf.Values{c} = ContNames{c};
end
Srf.Values = Srf.Values';
Srf.Functional = 'GLM contrasts';

% Compress & save contrast file
Srf = samsrf_compress_srf(Srf,vx);
save([SrfFile '_conts'], 'Srf');
