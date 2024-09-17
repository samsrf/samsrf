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
% The 3rd element of Threshold defines the tolerance. This is the
% proportion of vertices -below- threshold in an iteration that will be 
% filled in. If this is zero (default), no subthreshold vertices are filled. 
% If it is >0 then some vertices will be ignored. A tolerance of 1 would
% continue indefinitely so this is unwise.
%
% Returns a list of vertex indices in Vs.
%
% 20/04/2022 - SamSrf 8 version (DSS)
% 13/03/2024 - Added stopping parameter for filling (DSS)
% 16/03/2024 - Added tolerance parameter for filling subthreshold (DSS)
%              Fixed incorrect comment (DSS)
%

% Criterion?
if length(Threshold) == 1
    Threshold = [Threshold 0];
end
% Tolerance?
if length(Threshold) == 2
    Threshold = [Threshold 0];
end
if Threshold(2) > 1
    samsrf_error('Criterion must not be greater than 1!');
end
if Threshold(3) < 0 || Threshold(3) >= 1
    samsrf_error('Tolerance must be greater than 0 & smaller than 1!');
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
    nb = nb(nb ~= v); % Remove core vertex
    
    if ~isempty(nb)
        cS = Statistic(nb); % Statistic of current neighbour vertices
        % Any filling criterion?
        if Threshold(2) > 0
            % Use proportion of neighbours
            if mean(cS > Threshold(1)) <= Threshold(2)
                cS(cS > Threshold(1)) = -Inf; % Label already selected
            end
        elseif Threshold(2) < 0 
            % Use number of neighbours
            if length(nb) > -Threshold(2) 
                if sum(cS > Threshold(1)) <= -Threshold(2)
                    cS(cS > Threshold(1)) = -Inf; % Label already selected
                end
            else
                cS(nb) = -Inf; % Not yet enough neighbours so keep filling
            end
        end
        
        % Any filling tolerance?
        if Threshold(3) > 0
            % Proportion unlabelled subthreshold neighbours within tolerance?
            if mean(cS <= Threshold(1) & isinf(cS)) <= Threshold(3)
                cS(cS <= Threshold(1)) = Threshold(1) + 1; % Reassign as suprathreshold
            end
        end
        
        % Fill in suprathreshold vertices
        nb = nb(cS > Threshold(1)); % Remove neighbours below threshold or already selected
        Statistic(nb) = -Inf; % Remove new vertices from future selection
        Vs = [Vs; nb]; % Add new vertices to output vector
    end
end