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
% 20/04/2022 - SamSrf 8 version (DSS)
%

%% Default parameters
if nargin < 4
    MaxDist = 20;
end

%% Constants & output matrix
% Number of vertices
N = length(Vs);
% Output matrix
Ds = NaN(N,N);

%% Loop thru vertices
samsrf_disp('Calculating distance matrix...');
samsrf_disp(' Please stand by...');
parfor v = 1:N
    [Nv,Nd] = samsrf_georoi(Vs(v), MaxDist, V, F); % Neighbourhood vertices & distances
    X = Inf(1,size(V,1)); % Dummy vector of vertices
    X(Nv) = Nd; % Distances on surface
    Ds(v,:) = X(Vs); % Store in matrix
end