function mSrf = samsrf_average_srf(SrfCell, mSrfFile)
%
% mSrf = samsrf_average_srf(SrfCell, mSrfFile)
%
% Averages the data in the surface data files listed in the cell SrfCell.
% This makes most sense in the context of averaging time courses but there
% may be other uses for it.
%
% The function returns the averaged Srf in mSrf.
% If mSrfFile is defined, the averaged Srf is saved in this file 
% (name is prepended by the hemisphere as with all Srf files).
%
% 20/04/2022 - SamSrf 8 version (DSS)
%

if nargin == 1
    mSrfFile = '';
end

% Initialise averaged Srf
load(EnsurePath(SrfCell{1}));
Srf = samsrf_expand_srf(Srf);
mSrf = Srf;
mSrf.Data = NaN(size(Srf.Data,1), size(Srf.Data,2), length(SrfCell));
if isfield(mSrf, 'Raw_Data')
    mSrf.Raw_Data = NaN(size(Srf.Raw_Data,1), size(Srf.Raw_Data,2), length(SrfCell));
end

% Loop through surface files
for i = 1:length(SrfCell)
    % Load current Srf
    load(EnsurePath(SrfCell{i}));
    [Srf,vx] = samsrf_expand_srf(Srf);
    mSrf.Data(:,:,i) = Srf.Data;
    if isfield(mSrf, 'Raw_Data')
        mSrf.Raw_Data(:,:,i) = Srf.Raw_Data;
    end
end

% Average data
mSrf.Data = nanmean(mSrf.Data,3);
if isfield(mSrf, 'Raw_Data')
    mSrf.Raw_Data = nanmean(mSrf.Raw_Data,3);
end
mSrf.Functional = SrfCell;

% Compress Srf again
mSrf = samsrf_compress_srf(mSrf,vx);

% Save averaged file
if ~isempty(mSrfFile)
    savdir = fileparts(mSrfFile);
    if isempty(savdir)
        savdir = '.';
    end
    mSrfFile = [savdir filesep mSrf.Hemisphere '_' mSrfFile];
    Srf = mSrf;
    clear mSrf;
    save(mSrfFile, 'Srf');
    disp(['Saved average Srf ' mSrfFile]);
end
