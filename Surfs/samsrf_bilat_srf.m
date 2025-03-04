function Srf = samsrf_bilat_srf(SrfL, SrfR)
%
% Srf = samsrf_bilat_srf(SrfL, SrfR)
%
% Returns a Srf structure that combines the data from both hemispheres. 
%  This is useful when you want to combine hemispheres but an analysis function 
%  requires a Srf as input. You can also display these combined hemispheres with 
%  the rendering functions and in theory you can even delineate them together,
%  but this is probably not recommended and currently you cannot actually save 
%  combined labels with the DelineationTool.
%
% The left hemisphere is always first in the resulting data fields. Therefore,
%  that the vertex indices for the right hemisphere are increased by the number 
%  of vertices in the left hemisphere. So right hemisphere vertex N is now:
%       N_combined = N + size(SrfL.Vertices,1) 
%  The combined Srf contains a filed called Srf.Nvert_Lhem with that number. 
%
% By convention, surface data files for combined hemispheres are prefixed 'bi_'
% If the original surface structures contained a meshes field, this is removed.
% You can split the anatomical surfaces using samsrf_anatomy_srf as usual.
%
% You can create combined ROI labels using the function samsrf_bilat_label.
%  This will save a label without the hemisphere prefix (so e.g. V1.label) and 
%  the vertex indices are now the ones based on the combined Srf.
%
% SrfL/SrfR: Srf structures for left & right hemispheres, respectively.
%             Automatically expanded but you may want to denoise them first.
%
%
% Warning: May fail with very large data files due to lack of memory.
%
% 20/04/2022 - SamSrf 8 version (DSS)
% 29/11/2022 - Now supports Srfs with missing surface fields (DSS)
%

% Expand Srfs
SrfL = samsrf_expand_srf(SrfL);
SrfR = samsrf_expand_srf(SrfR);

% Shift hemispheres
if isfield(SrfL, 'Pial')
    SrfL.Pial(:,1) = SrfL.Pial(:,1)-3;
    SrfR.Pial(:,1) = SrfR.Pial(:,1)+3;
else
    samsrf_disp('Warning: Pial surfaces not in Srf...');
end

if isfield(SrfL, 'Inflated')
    SrfL.Inflated(:,1) = SrfL.Inflated(:,1)-50;
    SrfR.Inflated(:,1) = SrfR.Inflated(:,1)+50;
else
    samsrf_disp('Warning: Inflated surfaces not in Srf...');
end

if isfield(SrfL, 'Sphere')
    SrfL.Sphere(:,1) = SrfL.Sphere(:,1)-105;
    SrfR.Sphere(:,1) = SrfR.Sphere(:,1)+105;
else
    samsrf_disp('Warning: Sphere surfaces not in Srf...');
end

% Combine hemispheres
Srf = SrfL;
Srf.Hemisphere = 'bi';
Srf.Nvert_Lhem = size(SrfL.Vertices,1);
Srf.Faces = [SrfL.Faces; SrfR.Faces + Srf.Nvert_Lhem];

Srf.Vertices = [SrfL.Vertices; SrfR.Vertices];

if isfield(SrfL, 'Normals')
    Srf.Normals = [SrfL.Normals; SrfR.Normals];
else
    samsrf_disp('Warning: Normals data not in Srf...');
end

if isfield(SrfL, 'Pial')
    Srf.Pial = [SrfL.Pial; SrfR.Pial];
end

if isfield(SrfL, 'Inflated')
    Srf.Inflated = [SrfL.Inflated; SrfR.Inflated];
end

if isfield(SrfL, 'Sphere')
    Srf.Sphere = [SrfL.Sphere; SrfR.Sphere];
end

if isfield(SrfL, 'Curvature')
    Srf.Curvature = [SrfL.Curvature SrfR.Curvature];
else
    samsrf_disp('Warning: Curvature data not in Srf...');
end
if isfield(SrfL, 'Area')
    Srf.Area = [SrfL.Area SrfR.Area];
    Srf.Thickness = [SrfL.Thickness SrfR.Thickness];
else
    samsrf_disp('Warning: Cortical dimensions not in Srf...');
end
Srf.Data = [SrfL.Data SrfR.Data];

if isfield(SrfL, 'Raw_Data')
    Srf.Raw_Data = [SrfL.Raw_Data SrfR.Raw_Data];
end
if isfield(SrfL, 'X')
    Srf.X = [SrfL.X SrfR.X];
end
if isfield(SrfL, 'Y')
    Srf.Y = [SrfL.Y SrfR.Y];
end
if isfield(SrfL, 'Rmaps')
    if ~isnan(SrfL.Rmaps)
        Srf.Rmaps = [SrfL.Rmaps SrfR.Rmaps];
    end
end
if isfield(SrfL, 'ConFlds')
    if ~isnan(SrfL.ConFlds) && ~isnan(SrfR.ConFlds)
        Srf.ConFlds = [SrfL.ConFlds SrfR.ConFlds];
    else
        Srf.ConFlds = NaN;
    end
end
if isfield(SrfL, 'SeedVx')
    Srf.SeedVx = [SrfL.SeedVx; SrfR.SeedVx + Srf.Nvert_Lhem];
end
if isfield(SrfL, 'Noise_Ceiling')
    Srf.Noise_Ceiling = [SrfL.Noise_Ceiling SrfR.Noise_Ceiling];
end

% Remove mesh field if it exists
if isfield(Srf, 'Meshes')
    Srf = rmfield(Srf, 'Meshes');
end