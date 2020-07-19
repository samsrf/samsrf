function ApName = CreateRoiRespApertures(DataSrfs, MapSrfs, RoiNames, Eccen, ApWidth)
%
% ApName = CreateRoiRespApertures(DataSrfs, MapSrfs, RoiNames, Eccen, [ApWidth=100])
%
% Uses the retinotopic maps (could be Benson prediction maps) in MapSrfs
% to backproject brain activity in DataSrfs into visual space. Each of
% these can be a vector of Srfs that must be of equal length. For each Srf
% in the vector you must also define a ROI in the cell array RoiNames.
% Typically, you would have a Srf for each hemisphere.
%
% To make this work across both Benson and empirical maps, all Sigmas in
% the MapSrfs will be set to 1. The backprojection does not take pRF size
% into account at all because it uses samsrf_backproj_srclt.
%
% Eccen defines the maximum eccentricity up to which to include data.
% 
% The optional input ApWidth (default=100) defines the height/width of the
% apertures to be generated.
%
% Saves the backprojected movie of the time course as an aperture that will
% be named ApName in the format aps_[RoiNames].
%
% 19/07/2020 - SamSrf 7 version (DSS)
%

if nargin < 5
    ApWidth = 100;
end

% Aperture file name
ApName = 'aps_';

% Number of Srfs
N = length(DataSrfs);
if length(MapSrfs) ~= N
    error('Inconsistent number of map and data Srfs!');
end
if length(RoiNames) ~= N
    error('Not enough ROIs defined!');
end

%% Prepare data
Data = []; % Data matrix
Map = []; % Map matrix
% Loop thru Srfs
for i = 1:N
    % Add ROI name
    if i > 1
        ApName = [ApName '+']; % Add plus sign if not the first ROI
    end
    ApName = [ApName RoiNames{i}];
    
    % Expand all data 
    CurDataSrf = samsrf_expand_srf(DataSrfs(i));
    CurMapSrf = samsrf_expand_srf(MapSrfs(i));
    
    % Load ROI
    Roi = samsrf_loadlabel(RoiNames{i});   

    % Prepare data 
    Data = [Data CurDataSrf.Data(:,Roi)]; % Data matrix
    Map = [Map CurMapSrf.Data(:,Roi)]; % Map matrix
end
% Set sigmas to 1
Map(4,:) = 1; 

%% Backprojection
Tc = samsrf_backproj_srclt(Data, Map, Eccen, 0);

%% Rescale images
Tc = 0.5 * Tc / nanmax(abs(Tc(:))) + 0.5; % Rescale into 0-1 range where 0.5 is zero
ApFrm = NaN(ApWidth, ApWidth, size(Tc,3)); % Aperture movie
for i = 1:size(Tc,3)
    ApFrm(:,:,i) = imresize(Tc(:,:,i), [1 1]*ApWidth) - 0.5;
end

%% Save apertures
save(ApName, 'ApFrm');