function Ds = samsrf_geomatrix(V, F, Vs, MaxDist)
%
% Ds = samsrf_geomatrix(V, F, Vs, [MaxDist=20])
%
% Calculates geodesic distances between the vertices defined by the indeces 
% in vector Vs using the vertex coordinates V and the faces F. The function
% determines the geodesic distance in terms of steps between neighbouring
% vertices. All distances greater than MaxDist are set to Inf.
%
% Returns a N x N matrix Ds containing the cortical distances between each 
% of the N vertices and all the others.
%
% 05/05/2020 - Written (DSS)
%

%% Default parameters
if nargin < 4
    MaxDist = 20;
end
% Check if waitbar is to be used
wb = samsrf_waitbarstatus;

%% Constants & output matrix
% Number of vertices
N = length(Vs);
% Output matrix
Ds = NaN(N,N);

%% Loop thru vertices
if wb h = waitbar(0, 'Calculating distance matrix...', 'Units', 'pixels', 'Position', [100 100 360 70]); end
disp('Calculating distance matrix...');
for v = 1:N
    [Nv,Nd] = samsrf_georoi(Vs(v), MaxDist, V, F); % Neighbourhood vertices & distances
    X = Inf(1,size(V,1)); % Dummy vector of vertices
    X(Nv) = Nd; % Distances on surface
    Ds(v,:) = X(Vs); % Store in matrix
    if wb waitbar(v/N, h); end
end
if wb close(h); end