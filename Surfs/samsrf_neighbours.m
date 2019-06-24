function Nb = samsrf_neighbours(v, F)
%
% Nb = samsrf_neighbours(v, F)
%
% Finds the neighbours of a particular vertex v (or a vector of vertices).
% You also need to provide the faces of the mesh in F.
%
% 09/08/2018 - SamSrf 6 version (DSS)
%
 
Nf = ismember(F,v);
Nb = F(sum(Nf,2)>0,:);
Nb = unique(Nb(:));
Nb(ismember(Nb,v)) = [];
