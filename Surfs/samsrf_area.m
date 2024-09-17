function A = samsrf_area(V, F, Vs)
%
% A = samsrf_area(V, F, [Vs])
%
% Calculates the area for vertices Vs (optional, default is all) based on 
% vertex coordinates V and faces F. V can contain the 3D coordinates from 
% FreeSurfer segmentation or the 2D coordinates from the pRF mapping. 
% This is only really needed for determining the area of visual space 
% associated with a set of vertices which you need for example to calculate 
% cortical magnification. FreeSurfer already calculates the cortical surface 
% area for each vertex (but see the note in samsrf_vertexarea.m about this!)
%
% 20/04/2022 - SamSrf 8 version (DSS)
%

if nargin < 3
    Vs = 1:size(V,1);
end

A = 0; 
samsrf_disp('  Calculating area...'); 
for i = 1:length(Vs)
    A = A + samsrf_vertexarea(Vs(i),V,F); 
end
