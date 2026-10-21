function hs = samsrf_showlabels(Srf, Mesh, Rois)
% 
% hs = samsrf_showlabels(Srf, Mesh, Rois)
%
% Shows labels for ROIs in the cell array Rois on a surface using the model Mesh as char (e.g. 'Sphere').
% This assumes a surface Srf has been displayed using samsrf_surf in the current figure.
%
% Path names & extensions from label names are automatically removed to avoid clutter.
% If a label has been defined via a vector of vertices, it is labelled with a simple descriptor.
%
% If Rois is empty nothing happens.
%
% Returns a struct hs that contains the handles to all the labels & their quiver arrows.
%
% 25/09/2026 - Written (DSS)

if isempty(Rois)
    hs = [];
    return
end
figure(gcf); hold on

%% Mesh vertices
Srf = samsrf_expand_srf(Srf); % Expand just in case
Mesh = lower(Mesh);
Mesh(1) = upper(Mesh(1));
if isfield(Srf, Mesh)
    Vs = Srf.(Mesh);
elseif strcmpi(Mesh, 'Fake-Flat')
    samsrf_disp('Labelling looks odd with fake-flat so this is not an option');
    return
elseif strcmpi(Mesh, 'Inflated-50%')
    Vs = Srf.Vertices + (Srf.Inflated-Srf.Vertices)*0.5;
elseif strcmpi(Mesh, 'Pial-50%')
    Vs = Srf.Pial + (Srf.Inflated-Srf.Pial)*0.5;
else
    if strcmpi(Mesh, 'White')
        Vs = Srf.Vertices;
    else
        samsrf_disp(['Warning: Unknown mesh ' Mesh ' specified, so using Vertices...']);
        Vs = Srf.Vertices;
    end
end
% Shift bilateral Srf to origin
if strcmpi(Srf.Hemisphere, 'bi')
    nvl = Srf.Nvert_Lhem;
    Lvx = 1:nvl;
    Rvx = nvl + 1:size(Srf.Vertices,1);
    if contains(Mesh, 'Pial')
        Vs(Lvx,1) = Vs(Lvx,1)+3;
        Vs(Rvx,1) = Vs(Rvx,1)-3;
    elseif contains(Mesh, 'Inflated')
        Vs(Lvx,1) = Vs(Lvx,1)+50;
        Vs(Rvx,1) = Vs(Rvx,1)-50;
    elseif strcmp(Mesh, 'Sphere')
        Vs(Lvx,1) = Vs(Lvx,1)+105;
        Vs(Rvx,1) = Vs(Rvx,1)-105;
    end 
else
    nvl = size(Srf.Vertices,1); 
end
% Place label properly
LabVs = Vs * 1.2; % Label at 20% above mesh surface
% Shift bilateral Srf back where it belongs 
if strcmpi(Srf.Hemisphere, 'bi')
    if contains(Mesh, 'Pial')
        Vs(Lvx,1) = Vs(Lvx,1)-3;
        Vs(Rvx,1) = Vs(Rvx,1)+3;
        LabVs(Lvx,1) = LabVs(Lvx,1)-3;
        LabVs(Rvx,1) = LabVs(Rvx,1)+3;
    elseif contains(Mesh, 'Inflated')
        Vs(Lvx,1) = Vs(Lvx,1)-50;
        Vs(Rvx,1) = Vs(Rvx,1)+50;
        LabVs(Lvx,1) = LabVs(Lvx,1)-50;
        LabVs(Rvx,1) = LabVs(Rvx,1)+50;
    elseif strcmp(Mesh, 'Sphere')
        Vs(Lvx,1) = Vs(Lvx,1)-105;
        Vs(Rvx,1) = Vs(Rvx,1)+105;
        LabVs(Lvx,1) = LabVs(Lvx,1)-105;
        LabVs(Rvx,1) = LabVs(Rvx,1)+105;
    end
end

%% Ensure cell array
if ~iscell(Rois)
    Rois = {Rois};
end

%% Loop thru ROIs/labels
ln = 1;
for i = 1:length(Rois)
    % Only do this if not colour value
    if ~isnan(Rois{i}(1))
        %% Retrieve label vertices
        if ischar(Rois{i})
            if contains(Rois{i}, '.label')
                CurVx = samsrf_loadlabel(Rois{i}(1:end-6)); % Label file
            else
                CurVx = samsrf_loadpaths(Rois{i}); % Path file
            end
        else
            CurVx = Rois{i}; % Vector of vertices
        end
        if ischar(Rois{i})
            % Label ROIs or paths by name
            [~,LabTxt] = fileparts(Rois{i});
        else
            % Label vectors with count
            LabTxt = ['Vector #' num2str(ln)];
            ln = ln + 1;
        end

        %% 1st hemisphere (left if bilateral)
        % Label midpoint & neighbours
        cv = samsrf_cluster_midpt(CurVx(CurVx <= nvl), Vs); 
        % Place label
        quiver3(LabVs(cv,1), LabVs(cv,2), LabVs(cv,3), Vs(cv,1)-LabVs(cv,1), Vs(cv,2)-LabVs(cv,2), Vs(cv,3)-LabVs(cv,3), 'off', 'Color', 'w', 'LineWidth', 2);
        text(LabVs(cv,1), LabVs(cv,2), LabVs(cv,3), LabTxt, 'Interpreter', 'none', 'Backgroundcolor', 'w', 'EdgeColor', 'k');

        %% Right hemisphere if bilateral
        % Label midpoint & neighbours
        cv = samsrf_cluster_midpt(CurVx(CurVx > nvl), Vs); 
        if ~isempty(cv)
            % Place label
            quiver3(LabVs(cv,1), LabVs(cv,2), LabVs(cv,3), Vs(cv,1)-LabVs(cv,1), Vs(cv,2)-LabVs(cv,2), Vs(cv,3)-LabVs(cv,3), 'off', 'Color', 'w', 'LineWidth', 2);
            text(LabVs(cv,1), LabVs(cv,2), LabVs(cv,3), LabTxt, 'Interpreter', 'none', 'Backgroundcolor', 'w', 'EdgeColor', 'k');        
        end
    end
end

%% Return handles
hs = struct;
t = findall(gcf, 'Type', 'Text');
t = flipud(t(1:end-2)); % Invert order & remove two superfluous texts
q = findall(gcf, 'Type', 'Quiver');
q = flipud(q); % Invert order 
hs.Labels = t;
hs.Quivers = q;

hold off
