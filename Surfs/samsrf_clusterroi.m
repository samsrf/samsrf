function Vs = samsrf_clusterroi(v, Statistic, Threshold, Faces)
%
% Vs = samsrf_clusterroi(v, Statistic, Threshold, Faces)
%
% Selects a contiguous cluster surounding vertex v of vertices for which 
% Statistic is greater than Threshold. Statstic is a vector with the 
% relevant values per vertex. If you want to select using significance 
% (p-values) you need to invert the sign of Statistic.
% The function requires the list of faces in Faces.
%
% Returns a list of vertex indices in Vs.
%
% 08/04/2022 - Written (DSS)
%

% Is core vertex above threshold?
if Statistic(v) > Threshold
    % Core vertex is first selection
    Vs = v;
    Statistic(v) = -Inf; % Remove core vertex from future selections
else
    % Nothing left to do
    Vs = [];
    return
end

% Loop until nothing left
nb = v; % Initialise floodfill
while ~isempty(nb)
    nb = samsrf_neighbours(nb, Faces); % Vertices neighbouring previous neighbours
    if ~isempty(nb)
        cS = Statistic(nb); % Statistic of current neighbour vertices
        nb = nb(cS > Threshold); % Remove neighbours below threshold or already selected
        Statistic(nb) = -Inf; % Remove new vertices from future selection
        Vs = [Vs; nb]; % Add new vertices to output vector
    end
end