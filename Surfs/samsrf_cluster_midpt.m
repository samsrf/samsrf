function v = samsrf_cluster_midpt(Vx, Mesh)
%
% v = samsrf_cluster_midpt(Vx, Mesh)
% 
% Finds the mean position of the cluster defined by vertex indeces Vx
% assuming the 3d Mesh (3-column matrix). The idea is to use the sphere mesh 
% but you can use other meshes, but results may not be what you want.
%
% 20/04/2026 - Written (DSS)
%

D = sum(abs(zscore(Mesh(Vx,:))),2); % Z-score of each vertex in each dimension
v = Vx(find(D == min(D), 1)); % Vertex closest to centre   