function Srf = samsrf_combine_srf(SrfL, SrfR)
%
% Srf = samsrf_combine_srf(SrfL, SrfR)
%
% Returns a Srf structure that combines the data from both hemispheres. 
% This is useful when you want to combine hemispheres but the function 
% requires a Srf as input. The left hemisphere will always be first in the
% resulting data fields. 
%
% SrfL/SrfR: Srf structures for left & right hemispheres, respectively.
%             Automatically expanded but you may want to denoise them first.
%
% Roi:       ROI label name, e.g. 'V1'
%
% 28/10/2020 - Written (DSS)
%

% Expand Srfs
SrfL = samsrf_expand_srf(SrfL);
SrfR = samsrf_expand_srf(SrfR);

% Shift hemispheres
if isfield(SrfL, 'Sphere')
    SrfL.Pial(:,1) = SrfL.Pial(:,1)-3;
    SrfR.Pial(:,1) = SrfR.Pial(:,1)+3;
    SrfL.Inflated(:,1) = SrfL.Inflated(:,1)-50;
    SrfR.Inflated(:,1) = SrfR.Inflated(:,1)+50;
    SrfL.Sphere(:,1) = SrfL.Sphere(:,1)-105;
    SrfR.Sphere(:,1) = SrfR.Sphere(:,1)+105;
else
    warning('Anatomical surfaces not in Srf...');
end

% Combine hemispheres
Srf = SrfL;
Srf.Hemisphere = 'both';
Srf.Num_Left_Hemi = size(SrfL.Vertices,1);
Srf.Faces = [SrfL.Faces; SrfR.Faces + Srf.Num_Left_Hemi];

Srf.Vertices = [SrfL.Vertices; SrfR.Vertices];
Srf.Normals = [SrfL.Normals; SrfR.Normals];

if isfield(SrfL, 'Pial')
    Srf.Pial = [SrfL.Pial; SrfR.Pial];
    Srf.Inflated = [SrfL.Inflated; SrfR.Inflated];
    Srf.Sphere = [SrfL.Sphere; SrfR.Sphere];
end

if isfield(SrfL, 'Curvature')
    Srf.Curvature = [SrfL.Curvature SrfR.Curvature];
    Srf.Area = [SrfL.Area SrfR.Area];
    Srf.Thickness = [SrfL.Thickness SrfR.Thickness];
else
    warning('Cortical folding data not in Srf...');
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
    Srf.Rmaps = [SrfL.Rmaps SrfR.Rmaps];
end
if isfield(SrfL, 'ConFlds')
    Srf.ConFlds = [SrfL.ConFlds SrfR.ConFlds];
end
if isfield(SrfL, 'Noise_Ceiling')
    Srf.Noise_Ceiling = [SrfL.Noise_Ceiling SrfR.Noise_Ceiling];
end