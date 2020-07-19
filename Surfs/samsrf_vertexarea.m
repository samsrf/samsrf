function A = samsrf_vertexarea(v, V, F)
%
% A = samsrf_vertexarea(v, V, F)
%
% Calculates the area for vertex v based on vertex coordinates V and faces F.
% V can contain the 3D coordinates from FreeSurfer segmentation or the 2D
% coordinates from the pRF mapping.
%
% We only really need to run this function for the area of visual space 
% associated with each vertex as the FreeSurfer area already provides the 
% area of the cortical representation.
%
% 19/07/2020 - SamSrf 7 version (DSS)
%

F = F(sum(F==v,2)~=0,:);
if size(V,2) == 2
    V(:,3) = 0;
end

A = 0;
for i = 1:size(F,1)
    a = V(F(i,2),:)-V(F(i,1),:);    % Vector from vertex 1 to vertex 2
    b = V(F(i,3),:)-V(F(i,1),:);    % Vector from vertex 1 to vertex 3
    c = norm(cross(a,b)) / 2;       % Half the magnitude of the cross product 
    A = A + c;  % Add this area to the sum
end

% Divide by 3 because of overlap
A = A / 3;