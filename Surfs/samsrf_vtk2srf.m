function samsrf_vtk2srf(Vtks)
%
% Convert Leipzig 7T VTK files in the cell array Vtks into a SamSrf file.
% The first VTK must be the segmentation surface containing the curvature.
% The following (optional) ones can contain betas, individual TRs, etc.
% Obviously, these should all be of the same surface model.
%
% The saved Srf file will be named after Vtks{1}. This file name must
% contain either 'lcr' or 'rcr' to identify the hemisphere.
%
% 20/04/2022 - SamSrf 8 version (DSS)
% 29/06/2023 - Added conversion to 32 bit (single) data (DSS)
%

if ischar(Vtks)
    Vtks = {Vtks};
end

% Load surface & curvature
[V, F, D] = Read_vtk(Vtks{1});
V(:,1) = -V(:,1);
V(:,3) = -V(:,3);
samsrf_disp(['Surface & curvature from: ' Vtks{1}]);

% Create Srf
Srf = struct;
Srf.Structural = 'n/a';
Srf.Functional = 'Converted from Leipzig-VTK';
Srf.Hemisphere = '';
Srf.Cortex_Steps = NaN;
Srf.Vertices = V;
Srf.Curvature = D;
Srf.Faces = F;
Srf.Normals = NaN;
Srf.Data = [];
Srf.Rule = 'n/a';
Srf.Values = Vtks(2:end)';
if strfind(Vtks{1}, 'lcr')
    Srf.Hemisphere = 'lh';
elseif strfind(Vtks{1}, 'rcr')
    Srf.Hemisphere = 'rh';
else
    samsrf_error('Can''t figure out which hemisphere this is!');    
end

% Loop through VTK files
for i = 2:length(Vtks)
    [~,~,D] = Read_vtk(Vtks{i});
    Srf.Data = [Srf.Data; D'];
    samsrf_disp(['Read in #' num2str(i-1) ': ' Vtks{i}]);
end

% Convert to 32 bit?
Srf = samsrf_32bit_srf(Srf);

% Save Srf
save([Srf.Hemisphere '_' Vtks{1}(1:end-4)], 'Srf');


