function MakeOccRoi(MeshFolder, Y)
%
% MakeOccRoi(MeshFolder, [Y = -35])
%
% Creates an occipital ROI for pRF analysis for each hemisphere.
% The first input argument is the MeshFolder (e.g. '..\surf'). 
% The optional second input Y defines the border of the ROI in 
%  anterior-posterior coordinates (default = -35).
%

if nargin < 2
    Y = -35;
end

Hemis = {'lh' 'rh'};
% Loop thru hemispheres
for h = 1:2
    Vs = fs_read_surf([MeshFolder filesep Hemis{h} '.inflated']);
    Srf = struct;
    Srf.Vertices = Vs;
    Srf.Data = zeros(1,size(Srf.Vertices,1));
    vx = find(Vs(:,2) <= Y);
    samsrf_srf2label(Srf, [Hemis{h} '_occ'], 1, vx);
end

