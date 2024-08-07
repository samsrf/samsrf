function [Vs, Ds] = samsrf_georoi(v, s, V, F)
%
% [Vs, Ds] = samsrf_georoi(v, s, V, F)
%
% Finds the vertices in the neighbourhood of v within s steps to create a
% geodesic ROI. Requires vertices V and faces F. Returns the vertices in 
% Vs and in Ds the distances in steps from the centre.
%
% 20/04/2022 - SamSrf 8 version (DSS)
%

n = size(V,1);
Ds = zeros(n,1);
Ds(v) = 0;
Vs = [];

Nb = v;
for i = 1:s+1
    Vs = [Vs; Nb];
    cn = samsrf_neighbours(Nb, F);
    Ds(Nb) = Ds(Nb) + 1;
    Ds(cn) = Ds(cn) + 1;
    Nb = cn;
end

Ds = s - Ds + 1;
Vs = unique(Vs);
Ds = Ds(Vs);