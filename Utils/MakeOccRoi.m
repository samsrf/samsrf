function MakeOccRoi(MeshFolder, Y)
%
% MakeOccRoi(MeshFolder, [Y = -35])
%
% Creates an occipital ROI for pRF analysis for each hemisphere.
% The first input argument is the MeshFolder (e.g. '..\surf'). 
% The optional second input Y defines the border of the ROI in 
%  anterior-posterior coordinates (default = -35).
%
% Note that the vertex coordinates saved in the label file are 
%  from the spherical mesh, -not- the white matter mesh as normally.
%
% 15/09/2024 - Now also automatically creates binocular ROI (DSS)
%

if nargin < 2
    Y = -35;
end

Hemis = {'lh' 'rh'};
% Loop thru hemispheres
for h = 1:2
    [Vs Fs] = fs_read_surf([MeshFolder filesep Hemis{h} '.inflated']);
    Srf = struct;
    Srf.Hemisphere = Hemis{h};
    Srf.Vertices = Vs;
    Srf.Faces = Fs;
    Srf.Data = zeros(1,size(Srf.Vertices,1));
    vx = find(Vs(:,2) <= Y);
    samsrf_srf2label(Srf, [Hemis{h} '_occ'], 1, vx);
    if h == 1
        Lsrf = Srf;
    else
        Rsrf = Srf;
    end
end

% Bilateral Srf
Srf = samsrf_bilat_srf(Lsrf, Rsrf);
% Bilateral ROI
samsrf_bilat_label(Srf, 'occ');
