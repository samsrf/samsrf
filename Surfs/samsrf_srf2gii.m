function samsrf_srf2gii(SrfData, Value)
%
% samsrf_srf2asc(SrfData, [Value=1])
%
% Saves data in a Srf as a GIfTI surface file. Note that this requires
% hemispheric Srf data so if yours are bilateral they must be split first.
%
% NOTE: This function requires SPM12 for GIfTI functionality.
%
%   SrfData:    Either a Srf variable or a surface data file 
%
%   Value:      Which row in Srf.Data to export (defaults to 1)
%                This can either be the name (e.g. 'Sigma') or the number.
%                The name in Srf.Values for this row determines the GII name.
%
%               If Value is Inf, all rows in Srf.Data are saved in a matrix.
%                In that case the file is named after the SrfData file.
%
%               If you provided a Srf variable, the file is simply named lh/rh_Srf.gii
%                and you will probably want to rename it afterwards.
%
% 14/09/2022 - Thus it was writ (DSS)
%

% Load/prepare data
if ischar(SrfData)
    load(EnsurePath(SrfData));
else
    Srf = SrfData;
    SrfData = [Srf.Hemisphere '_Srf'];
end
Srf = samsrf_expand_srf(Srf);

% Which row?
if nargin < 2
    Value = 1;
end

% Name of file
if isscalar(Value)
    if isinf(Value)
        valrow = 1:size(Srf.Data,1);
        giiname = [SrfData '.gii'];
    else
        valrow = Value;
        giiname = [Srf.Hemisphere '_' Srf.Values{valrow} '.gii'];
    end
else
    valrow = find(strcmpi(Srf.Values, Value));
    giiname = [Srf.Hemisphere '_' Srf.Values{valrow} '.gii'];
end

% Save GII file
G = gifti;
G.cdata = Srf.Data(valrow,:)';
save(G, giiname);
disp(['Saved ' giiname '.']);
