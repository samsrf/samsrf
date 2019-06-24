function samsrf_export_labels(Srf, Ecc, OutFile, Raw)
%
% samsrf_export_labels(Srf, Ecc, OutFile, Raw)
%
% Saves labels for all the value fields in Srf (it doesn't do anything if
% there is no Srf.Values field). If Srf contains x0 and y0 values it also
% saves labels containing polar angle and eccentricity maps. You must to 
% define the maximal eccentricity you used for your stimulus in Ecc as 
% otherwise the eccentricity map wouldn't know what the scale is. Finally,
% you also have to define a file name to prefix to each label.
%
% You can use this function to load maps in FreeSurfer.
%
% 09/08/2018 - SamSrf 6 version (DSS)
%

% Default
if nargin < 4 
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
  
    % Vertex indices
    mver = 1:size(Srf.Vertices,1);

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
        E = sqrt(Srf.Data(2,:).^2 + Srf.Data(3,:).^2) / Ecc; 
        % Eccentricity for FreeSurfer
        [Er Ei] = pol2cart(E * 2*pi, ones(size(E)) .* Srf.Data(1,:)); % Complex numbers
        Er(Srf.Data(1,:) == 0) = 0;
        Ei(Srf.Data(1,:) == 0) = 0;
        % In the 1D label, negative numbers will be below half maximum eccentricity.
        E = E - 0.5; 
        E(Srf.Data(1,:) == 0) = 0;
        
        %% Save labels
        Sav = Srf;
        % Polar angle
        Fname = [OutFile '_pol'];
        Sav.Data = P;
        samsrf_srf2label(Sav, Fname, 1, mver);
        % Polar angle - real component
        Fname = [OutFile '_pol_r'];
        Sav.Data = Pr;
        samsrf_srf2label(Sav, Fname, 1, mver);
        % Polar angle - imaginary component
        Fname = [OutFile '_pol_i'];
        Sav.Data = Pi;
        samsrf_srf2label(Sav, Fname, 1, mver);
        % Eccentricity
        Fname = [OutFile '_ecc'];
        Sav.Data = E;
        samsrf_srf2label(Sav, Fname, 1, mver);
        % Eccentricity - real component
        Fname = [OutFile '_ecc_r'];
        Sav.Data = Er;
        samsrf_srf2label(Sav, Fname, 1, mver);
        % Eccentricity - imaginary component
        Fname = [OutFile '_ecc_i'];
        Sav.Data = Ei;
        samsrf_srf2label(Sav, Fname, 1, mver);
    end
    
    %% If this is a tuning curve map
     if cell2mat(strfind(Srf.Values, 'Mu'))
        
        % Extract values
        R2 = Srf.Data(find(strcmpi(Srf.Values, 'R^2')), :);
        Mu = Srf.Data(find(strcmpi(Srf.Values, 'Mu')), :);
        
        % Convert Mu into complex map
        Blank = ones(size(Mu));
        P = atan2(Mu, Blank);
        [Pr, Pi] = pol2cart(P, ones(size(P)) .* R2);
        Pr(R2 == 0) = 0;
        Pi(R2 == 0) = 0;
        P = P ./ Pi * 180;
                
        %% Save labels
        Sav = Srf;
        
        % Mu
        Fname = [OutFile '_mu'];
        Sav.Data = P;
        samsrf_srf2label(Sav, Fname, 1, mver);
        
        % Mu - real component
        Fname = [OutFile '_mu_r'];
        Sav.Data = Pr;
        samsrf_srf2label(Sav, Fname, 1, mver);
        
        % Mu - imaginary component
        Fname = [OutFile '_mu_i'];
        Sav.Data = Pi;
        samsrf_srf2label(Sav, Fname, 1, mver);
    end

    %% Save labels for each field
    for i = 1:length(Srf.Values)
        Fname = [OutFile '_' lower(Srf.Values{i})];
        sp = strfind(Fname, ' ');
        Fname(sp) = '_';
        samsrf_srf2label(Srf, Fname, i, mver);
    end
end
    
new_line;

