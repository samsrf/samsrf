function Srf = samsrf_compress_srf(Srf,vx)
%
% Srf = samsrf_compress_srf(Srf,vx)
%
% Compresses the surface data Srf so that only the vertices indexed by vx
% (for example from a particular region of interest) are saved. This adds a
% field named Srf.Roi to the structure to mark this data structure as a
% compressed version. You can use samsrf_expand_srf to expand it again.
% Only if the length of vx is shorter than the number of vertices overall
% the compression takes place, unless vx is empty.
%
% Note: this function only works for single-subject Srfs. There currently
%   is no support for expanding or compressing multi-subbject Srfs.
%
% 16/07/2020 - SamSrf 7 version (DSS)
%

%% Remove data outside of region of interest
if ~isempty(vx) && length(vx) < size(Srf.Vertices,1)
    disp('Compressing surface data file...');
    
    % ROI vertices
    Srf.Roi = vx;
    
    % Remove unwanted vertices
    Srf.Data = Srf.Data(:,vx);
    if isfield(Srf, 'Y')
        Srf.Y = Srf.Y(:,vx);
    end
    if isfield(Srf, 'X')
        Srf.X = Srf.X(:,vx);
    end
    if isfield(Srf, 'Raw_Data')
        Srf.Raw_Data = Srf.Raw_Data(:,vx);
    end
    if isfield(Srf, 'Rmaps')
        Srf.Rmaps = Srf.Rmaps(:,vx);
    end
    if isfield(Srf, 'ConFlds')
        Srf.ConFlds = Srf.ConFlds(:,vx);
    end
    if isfield(Srf, 'Noise_Ceiling')
        Srf.Noise_Ceiling = Srf.Noise_Ceiling(:,vx);
    end
end

%% If anatomy is saves separately, remove it from Srf
if isfield(Srf, 'Meshes')
    disp('Removing anatomical meshes...');
    Srf = rmfield(Srf, 'Vertices');
    Srf = rmfield(Srf, 'Faces');
    Srf = rmfield(Srf, 'Normals');
    Srf = rmfield(Srf, 'Pial');
    Srf = rmfield(Srf, 'Inflated');
    Srf = rmfield(Srf, 'Sphere');
    Srf = rmfield(Srf, 'Curvature');
    Srf = rmfield(Srf, 'Area');
    Srf = rmfield(Srf, 'Thickness');
end
