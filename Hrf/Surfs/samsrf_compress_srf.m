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
% 14/02/2022 - CF correlation profiles are not saved in data file by default (DSS)
% 14/03/2022 - pRF reverse correlation profiles are not saved in data file by default (DSS)
%              Now supports anonymised Srfs (DSS)
% 20/04/2022 - SamSrf 8 version (DSS)
% 15/05/2022 - Now works with volumetric data (DSS)
% 29/06/2023 - Defaults to saving 32 bit data unless fixed in SamSrf_defaults (DSS)  
% 22/08/2023 - Now supports EEG/MEG data files (DSS)
%

%% Do nothing if volumetric data! 
if ~strcmpi(Srf.Hemisphere, 'vol') && ~strcmpi(Srf.Hemisphere, 'eeg')

    %% Remove data outside of region of interest
    if ~isempty(vx) && length(vx) < size(Srf.Vertices,1)
        samsrf_disp('Compressing surface data file...');

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
            if ~isnan(Srf.Rmaps)
                Srf.Rmaps = Srf.Rmaps(:,vx);
            end
        end
        if isfield(Srf, 'ConFlds')
            if ~isnan(Srf.ConFlds)
                Srf.ConFlds = Srf.ConFlds(:,vx);
            end
        end
        if isfield(Srf, 'Noise_Ceiling')
            Srf.Noise_Ceiling = Srf.Noise_Ceiling(:,vx);
        end
    end

    %% If anatomy is saved separately, remove it from Srf
    if isfield(Srf, 'Meshes')
        samsrf_disp('Removing anatomical meshes...');
        Srf = rmfield(Srf, 'Vertices');
        Srf = rmfield(Srf, 'Faces');
        if isfield(Srf, 'Normals')
            % Don't exist for anonymised Srfs
            Srf = rmfield(Srf, 'Normals');
            Srf = rmfield(Srf, 'Pial');
        end
        Srf = rmfield(Srf, 'Inflated');
        Srf = rmfield(Srf, 'Sphere');
        Srf = rmfield(Srf, 'Curvature');
        if isfield(Srf, 'Area')
            % Don't exist for anonymised Srfs
            Srf = rmfield(Srf, 'Area');
            Srf = rmfield(Srf, 'Thickness');
        end
    end
end

%% Convert to 32 bit?
Srf = samsrf_32bit_srf(Srf);
