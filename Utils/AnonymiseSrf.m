function AnonymiseSrf(SrfFile, RoiVx, Curv)
%
% AnonymiseSrf(SrfFile, [RoiVx=[], Curv=[]])
%
% Removes potentially sensitive anatomical data from the Srf structure in
% SrfFile and replaces the Srf.Vertices field with the Srf.Inflated field.
%
% This function will redact the structural T1 file name. 
% HOWEVER, the names of the functional files are NOT redacted! 
% Please ensure that these do not contain identifying information!
%
% Please note that the function saves whatever else is in the original file
% so if this contains any identifying info (e.g. in the Model structure)
% then this is still present afterwards... 
%
% The second input RoiVx contains the vertex indeces of the ROI to compress 
% which you can use to reduce disc space (see samsrf_compress_srf).
%
% If the third input Curv is defined this contains curvature statistics.
% This must be a vector of the same length as the number of vertices in the
% surface model with the curvature. The idea is to use the nativised stats
% from the fsaverage template instead of the native curvatures, so that the
% unique anatomical information is removed but this generic curvature can
% be used to navigate in the deidentified surface model.
%
% Saves a new Srf file with the suffix _anon.
%
% 14/03/2022 - Written (DSS)
% 20/04/2022 - SamSrf 8 version (DSS)
%

if nargin < 2
    RoiVx = [];
end
if nargin < 3
    Curv = [];
end

% Load data file
load(EnsurePath(SrfFile));
% Expand Srf
Srf = samsrf_expand_srf(Srf);

% Remove meshes if it exists
if isfield(Srf, 'Meshes')
    Srf = rmfield(Srf, 'Meshes');
end

% Redact T1 file name 
Srf.Structural = 'Redacted-T1';

% Reorganise anatomy
Srf.Vertices = Srf.Inflated; % Replace white surface with inflated
% Remove unwanted fields 
Srf = rmfield(Srf, 'Pial');
Srf = rmfield(Srf, 'Normals');
Srf = rmfield(Srf, 'Area');
Srf = rmfield(Srf, 'Thickness');
% Replace or remove curvature?
if isempty(Curv)
    Srf.Curvature = zeros(1,size(Srf.Vertices,1));
    curvstr = ' without curvature.';
else
    if length(Curv) ~= size(Srf.Vertices,1)
        samsrf_error('Curvature statistics don''t match anatomy!');
    end
    Srf.Curvature = Curv;
    curvstr = ' with curvature.';
end

% Compress the bastard
Srf = samsrf_compress_srf(Srf, RoiVx);

% Save new file
samsrf_disp(['Anonymised surface data file ' SrfFile curvstr]);
clear RoiVx Curv curvstr
save([SrfFile '_anon'], '-v7.3');
% Separate anatomy
samsrf_anatomy_srf([SrfFile '_anon']);