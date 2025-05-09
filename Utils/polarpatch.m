function h = polarpatch(X, Y, C)
%
% h = polarpatch(X, Y, C)
% 
% Plots a surface patch containing a polar plot. X and Y contain
% coordinates, and C contains the associated intensity values.
%
% 02/03/2022 - Written (DSS)
% 20/04/2022 - SamSrf 8 version (DSS)
% 12/02/2024 - Fixed bug when using 32-bit data format (DSS)
% 17/09/2024 - Removed warning of duplicate points (DSS)
%

% Ensure column vectors
if size(X,2) > 1
    X = X';
end
if size(Y,2) > 1
    Y = Y';
end
if size(C,2) > 1
    C = C';
end

% Ensure double format
X = double(X);
Y = double(Y);
C = double(C);

% Plot patch
mx = range(X)/100; % Margins for X axis
my = range(Y)/100; % Margins for Y axis
warning off
tri = delaunay(X,Y); % Delaunay triangulation
warning on
h = trisurf(tri, X, Y, zeros(size(C)), 'FaceVertexCData', C, 'EdgeColor', 'none', 'FaceColor', 'interp'); % Draw patch
axis square
mx = range(X)/100; % Margins for X axis
my = range(Y)/100; % Margins for Y axis
axis([min(X)-mx max(X)+mx min(Y)-my max(Y)+my]); % Restrict axes
view([0 90]); % Turn camera for flat view
