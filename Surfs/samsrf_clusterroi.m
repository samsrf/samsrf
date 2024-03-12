function Vs = samsrf_clusterroi(v, Statistic, Threshold, Faces)
%
% Vs = samsrf_clusterroi(v, Statistic, Threshold, Faces)
%
% Selects a contiguous cluster surounding vertex v of vertices 
% for which Statistic is greater than Threshold. The function also requires 
% the list of faces in Faces.
%
% Statstic is a vector with the relevant values per vertex. If you want to 
% select using significance (log p-values) you need to invert the sign. 
%
% If Threshold is a vector, the 2nd element contains the stopping criterion.
% This defines the minimum of surrounding suprathreshold vertices necessary 
% in each iteration needed for filling to continue. 
%   0 (Default) = Filling as long as -any- neighbours are above threshold
%   0.5 = Filling as long as half of neighbours are suprathreshold
% 	1 = Filling only continues as long as all neighbours are suprathreshold
%   <0 = Filling whie absolute number of neighbours are suprathreshold
%
% Returns a list of vertex indices in Vs.
%
% 20/04/2022 - SamSrf 8 version (DSS)
% 13/03/2024 - Added stopping parameter for filling (DSS)
%

% Tolerance?
if length(Threshold) == 1
    Threshold = [Threshold 0];
end
if Threshold(2) > 1
    error('Filling criterion must not be greater than 1!');
end

% Is core vertex above threshold?
if Statistic(v) > Threshold(1)
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
        % Any filling criterion?
        if Threshold(2) > 0
            % Use proportion of neighbours
            if mean(cS > Threshold(1)) <= Threshold(2)
                cS(cS > Threshold(1)) = -Inf; % Unlabel if doesn't meet criterion
            end
        elseif Threshold(2) < 0 
            % Use number of neighbours
            if length(nb) > -Threshold(2) 
                if sum(cS > Threshold(1)) <= -Threshold(2)
                    cS(cS > Threshold(1)) = -Inf; % Unlabel if doesn't meet criterion
                end
            else
                cS(nb) = -Inf; % Not yet enough neighbours so keep filling
            end
        end
        nb = nb(cS > Threshold(1)); % Remove neighbours below threshold or already selected
        Statistic(nb) = -Inf; % Remove new vertices from future selection
        Vs = [Vs; nb]; % Add new vertices to output vector
    end
end