function samsrf_anat2map
%
% samsrf_anat2map
%
% Use this function to create an anatomical map based on sphere coordinates
% relative to an origin (for example, the foveal centre). Such anatomical
% maps can be useful for various situations, for instance as a template map
% in connective field analysis.
%
% The function will ask you to select a SamSrf file with a pRF/pTC map. 
% It will then display the mesh contained within and show the polar map, 
% if possible. If no polar map can be generated, it will show an empty map.
%
% Use the data cursor mode to select the vertex you want for your origin.
% Enter the desired vertex number in the prompt in the command line.
% The function will then calculate the relative anatomical coordinates, 
% and display these as a polar map. The Srf is then saved as ?h_anatmap.
%
% 20/04/2022 - SamSrf 8 version (DSS)
%

%% Load Srf
[f,p] = uigetfile('*h_*.mat', 'Select Srf file');
load([p f]);
[Srf,vx] = samsrf_expand_srf(Srf);

%% Display mesh
if isfield(Srf, 'Values')
    if strcmpi(Srf.Values{2}, 'x0') && strcmpi(Srf.Values{3}, 'y0')
        samsrf_surf(Srf, 'Inflated', 0, '', [], 'Polar'); % Polar exists
    else
        samsrf_surf(Srf, 'Inflated', 0.05, '', [], 1); % Show 1st data row
    end
else
    samsrf_error('Selected Srf must be a pRF or pTC map!');
end
set(gcf, 'name', 'Select origin vertex');

%% Input origin vertex
Origin = input('Enter origin vertex number: ');
close gcf

%% Calculate map
Srf.Functional = ['Anatomical maps (origin = ' num2str(Origin) ')'];
Srf.Data = [];
Srf.Values = {'R^2'; 'x0'; 'y0'};
% Label vertices inside the ROI as R^2 = 1
Srf.Data(1,:) = zeros(1,size(Srf.Vertices,1));
Srf.Data(1,vx) = 1;
% Calculate anatomical coordinates
Sph = Srf.Sphere(:,[1 3]);
% Coordinates relative to origin
Sph(:,1) = Sph(:,1) - Sph(Origin,1); % Left-Right coordinate   
Sph(:,2) = Sph(:,2) - Sph(Origin,2); % Inferior-Superior coordinate
Sph(:,2) = -Sph(:,2); % Invert to match visual field maps
Srf.Data = [Srf.Data; Sph']; % Add coordinates to Srf

%% Remove unnecessary fields
if isfield(Srf, 'Raw_Data')
    Srf = rmfield(Srf, 'Raw_Data');
end
if isfield(Srf, 'X')
    Srf = rmfield(Srf, 'Y');
end
if isfield(Srf, 'Y')
    Srf = rmfield(Srf, 'Y');
end
if isfield(Srf, 'Rmaps')
    Srf = rmfield(Srf, 'Rmaps');
end
if isfield(Srf, 'ConFlds')
    Srf = rmfield(Srf, 'ConFlds');
end

%% Display anatomical maps
samsrf_surf(Srf, 'Inflated', 0, '', [], 'Polar'); % Polar exists
set(gcf, 'name', Srf.Functional);

%% Compress Srf & save
Srf = samsrf_compress_srf(Srf,vx);
save([Srf.Hemisphere '_anatmap'], 'Srf', '-v7.3');

  