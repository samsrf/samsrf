function samsrf_srf2gii(SrfData, Value)
%
% samsrf_srf2gii(SrfData, [Value=1])
%
% Saves data in a Srf as a GIfTI surface file. Note that this requires
% hemispheric Srf data so if yours are bilateral they should be split first.
% Technically, it can save the bilateral GII file of course - it is just
% that probably no other program be able to read them properly.
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
% 07/10/2025 - Can now also export time series fields X & Y (DSS)
% 08/10/2025 - Ensured that it can be read by FreeView (DSS)
% 09/10/2025 - Updated help section for clarification (DSS)
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
if isscalar(Value) && ~ischar(Value)
    if isinf(Value)
        valrow = 1:size(Srf.Data,1);
        giiname = [SrfData '.gii'];
    else
        valrow = Value;
        giiname = [Srf.Hemisphere '_' Srf.Values{valrow} '.gii'];
    end
else
    if strcmp(Value, 'X') || strcmp(Value, 'Y') 
        valrow = Inf;
        giiname = [Srf.Hemisphere '_' Value '.gii'];
    else
        valrow = find(strcmpi(Srf.Values, Value));
        giiname = [Srf.Hemisphere '_' Srf.Values{valrow} '.gii'];
    end
end

% Save GII file
G = gifti;
if isinf(valrow)
	G.cdata = Srf.(Value)';
else
	G.cdata = Srf.Data(valrow,:)';
end
SaveGiiFile(G, giiname); % Save GII in FreeView readable format
samsrf_disp(['Saved ' giiname '.']);


%% Save GII file
function SaveGiiFile(Sav, Fname)
% Saves GII & then changes data dimensions so can be loaded by Freeview

save(Sav, Fname, 'GZipBase64Binary'); % Save GII
f = fopen(Fname); % Open GII
S = char(fread(f)'); % Read GII
S = strrep(S, "ColumnMajorOrder", "RowMajorOrder"); % Fix dimensions
fclose(f); % Close file

f = fopen(Fname, 'w'); % Now open for writing
fwrite(f, S); % Write data
fclose(f); % Close file

