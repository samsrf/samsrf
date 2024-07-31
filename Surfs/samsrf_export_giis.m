function samsrf_export_giis(Srf, OutFile, Raw)
%
% samsrf_export_giis(Srf, OutFile, [Raw=false])
%
% Saves GIfTI files for all the value fields in Srf (it doesn't do anything 
% if there is no Srf.Values field). If Srf contains x0 and y0 values it also
% saves labels containing polar angle and eccentricity maps. You also have 
% to define a file name to prefix to each label.
%
% The optional input Raw is a boolean that toggles whether to save the
% raw Srf.Raw_Data fields or the standard Srf.Data fields (default).
%
% You can use this function to load maps in Freeview, tksurfer, or other tools.
%
% NOTE: This function requires SPM12 for GIfTI functionality.
%
% This is a modern equivalent of samsrf_export_labels which we kept in the
% toolbox for backwards compatibility.
%
% 30/07/2024 - Created (DSS)
%

% Default
if nargin < 3 
    Raw = 0;
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
        [Pr Pi] = pol2cart(P, ones(size(P)) .* Srf.Data(1,:)); % Complex numbers
        Pr(Srf.Data(1,:) == 0) = 0;
        Pi(Srf.Data(1,:) == 0) = 0;
        P = P / pi * 180;
        if Srf.Hemisphere == 'rh'
            Pr = -Pr;
        end
        
        % Eccentricity 
        E = sqrt(Srf.Data(2,:).^2 + Srf.Data(3,:).^2); 
        % Eccentricity for FreeSurfer
        [Er Ei] = pol2cart(E * 2*pi, ones(size(E)) .* Srf.Data(1,:)); % Complex numbers
        Er(Srf.Data(1,:) == 0) = 0;
        Ei(Srf.Data(1,:) == 0) = 0;
        % In the 1D label, negative numbers will be below half maximum eccentricity.
        E = E - 0.5; 
        E(Srf.Data(1,:) == 0) = 0;
        
        %% Save surfaces       
        % Polar angle
        Fname = [OutFile '_pol.gii'];
        Sav.cdata = P';
        SaveGiiFile(Sav, Fname)
        disp(['Saved ' Fname]);
        
        % Eccentricity
        Fname = [OutFile '_ecc.gii'];
        Sav.cdata = E';
        SaveGiiFile(Sav, Fname)
        disp(['Saved ' Fname]);        
    end
    
    %% If this is a tuning curve map
    if cell2mat(strfind(Srf.Values, 'Mu'))
        
        % Extract values
        R2 = Srf.Data(find(strcmpi(Srf.Values, 'nR^2')), :);
        if isempty(R2)
            R2 = Srf.Data(find(strcmpi(Srf.Values, 'R^2')), :);
            disp('Using raw R^2 values as magnitude.');
        else
            disp('Using normalised R^2 values as magnitude.');
        end
        Mu = Srf.Data(find(strcmpi(Srf.Values, 'Mu')), :);
        
        % Convert Mu into complex map
        Blank = ones(size(Mu));
        P = atan2(Mu, Blank);
        [Pr, Pi] = pol2cart(P, ones(size(P)) .* R2);
        Pr(R2 == 0) = 0;
        Pi(R2 == 0) = 0;
        P = P ./ Pi * 180;
                        
        %% Save surface       
        % Mu
        Fname = [OutFile '_mu.gii'];
        Sav.cdata = P';
        SaveGiiFile(Sav, Fname)
        disp(['Saved ' Fname]);        
    end

    %% Save surfaces for each field
    for i = 1:length(Srf.Values)
        Fname = [OutFile '_' lower(Srf.Values{i}) '.gii'];
        sp = strfind(Fname, ' ');
        Fname(sp) = '_';
        Sav.cdata = Srf.Data(i,:)';
        SaveGiiFile(Sav, Fname)
        disp(['Saved ' Fname]);
    end
end
    
new_line;

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

