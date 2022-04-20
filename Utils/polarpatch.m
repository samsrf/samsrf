function h = polarpatch(X, Y, C)
%
% h = polarpatch(X, Y, C)
% 
% Plots a surface patch containing a polar plot. X and Y contain
% coordinates, and C contains the associated intensity values.
%
% 02/03/2022 - Written (DSS)
% 20/04/2022 - SamSrf 8 version (DSS)
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

% Plot patch
mx = range(X)/100; % Margins for X axis
my = range(Y)/100; % Margins for Y axis
tri = delaunay(X,Y); % Delaunay triangulation
h = trisurf(tri, X, Y, zeros(size(C)), 'FaceVertexCData', C, 'EdgeColor', 'none', 'FaceColor', 'interp'); % Draw patch
axis square
mx = range(X)/100; % Margins for X axis
my = range(Y)/100; % Margins for Y axis
axis([min(X)-mx max(X)+mx min(Y)-my max(Y)+my]); % Restrict axes
view([0 90]); % Turn camera for flat view
