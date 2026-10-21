function samsrf_export_giis(Srf, OutFile, Raw, Maps)
%
% samsrf_export_giis(Srf, OutFile, [Raw=false, Maps=[]])
%
% Saves GIfTI files for all the value fields in Srf (it doesn't do anything 
% if there is no Srf.Values field). If Srf contains x0 and y0 values it can 
% also save labels containing polar angle and eccentricity maps. You also have 
% to define a file name to prefix to each label.
%
% The optional input Raw is a boolean that toggles whether to save the
% raw Srf.Raw_Data fields or the standard Srf.Data fields (default).
%
% The optional input Maps defines the list of maps to export. This can
% either be a cell array for the names in Srf.Values or a vector with the
% rows in Srf.Data. Defaults to all rows (including Polar & Eccentricity).
%
% You can use this function to load maps in Fernanda Ribeiro's annotation web app, 
% in Freeview, tksurfer, or any other tools that reads GII files.
%
% IMPORTANT: This function requires SPM12 for GIfTI functionality!
%
% This is a modern equivalent of samsrf_export_labels which we kept in the
% toolbox for backwards compatibility.
%
% Note this function is very similar to samsrf_srf2gii but that one is for
% a more general use as it can save time series information. 
%
% 30/07/2024 - Created (DSS)
% 10/06/2025 - Added check if SPM is on path for GII functionality (DSS)
% 14/09/2026 - New option to only export some of the maps (DSS)
%              Removed complex outputs & Mu conversion as obsolete (DSS)
%

if ~exist('spm', 'file') % Use SPM 
    samsrf_error('Requires SPM for writing GII files!');
end
if strcmpi(Srf.Hemisphere, 'bi')
    samsrf_disp('WARNING: This is a bilateral Srf so most tools won''t read this GII...');
end

% Default
if nargin < 3 
    Raw = 0;
end
if nargin < 4 || isempty(Maps)
    Maps = 1:length(Srf.Values);
end
if ~iscell(Maps)
    Maps = Srf.Values(Maps);
end

% Expand if necessary
Srf = samsrf_expand_srf(Srf);

% Pick raw values if requested
if Raw
    for i = 1:size(Srf.Raw_Data,1)
        Srf.Data(i,:) = Srf.Raw_Data(i,:);
    end
end

if isfield(Srf, 'Values')
    %% Initialise GIfTI
    Sav = gifti;
    
    %% If this is a pRF retinotopic map
    if cell2mat(strfind(Srf.Values,'x0'))
        
        % Extract coordinates
        X = Srf.Data(2,:);  % Horizontal coordinate
        Y = Srf.Data(3,:);  % Vertical coordinate

        % Polar angle
        P = atan2(Y,X);         
        % Polar angle for FreeSurfer
        P = P / pi * 180;
        
        % Eccentricity 
        E = sqrt(Srf.Data(2,:).^2 + Srf.Data(3,:).^2); 
        E(Srf.Data(1,:) == 0) = 0;
        
        %% Save surfaces?       
        % Polar angle
        if sum(strcmpi(Maps, 'Polar'))
            Fname = [OutFile '_pol.gii'];
            Sav.cdata = P';
            SaveGiiFile(Sav, Fname)
            samsrf_disp(['Saved ' Fname]);
        end     
        % Eccentricity
        if sum(strcmpi(Maps, 'Eccentricity'))
            Fname = [OutFile '_ecc.gii'];
            Sav.cdata = E';
            SaveGiiFile(Sav, Fname)
            samsrf_disp(['Saved ' Fname]);    
        end
    end
    
    %% Save surfaces for each field
    for i = 1:length(Srf.Values)
        if sum(strcmpi(Maps, Srf.Values{i}))
            Fname = [OutFile '_' lower(Srf.Values{i}) '.gii'];
            sp = strfind(Fname, ' ');
            Fname(sp) = '_';
            Sav.cdata = Srf.Data(i,:)';
            SaveGiiFile(Sav, Fname)
            samsrf_disp(['Saved ' Fname]);
        end
    end
end
samsrf_newline;


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

