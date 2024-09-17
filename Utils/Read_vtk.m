function [vertex,face,points] = Read_vtk(filename, verbose)
%
% td_read_vtk - read data from Leipzig 7T VTK file.
%
%   [vertex,face,points] = Read_vtk(filename, verbose);
%
%   'vertex' is a 'nb.vert x 3' array specifying the position of the vertices.
%   'face' is a 'nb.face x 3' array specifying the connectivity of the mesh.
%   'points' is a 'nb.vert x 1' vector specifying the point data of each vertex
%
%   Copyright (c) Mario Richtsfeld
%
% 20/04/2022 - SamSrf 8 version (DSS)
%

if nargin<2
    verbose = 1;
end

fid = fopen(filename,'r');
if( fid==-1 )
    samsrf_error('Can''t open the file.');
end

str = fgets(fid);   % -1 if eof
if ~strcmp(str(3:5), 'vtk')
    samsrf_error('The file is not a valid VTK one.');    
end

%%% read header %%%
str = fgets(fid);
str = fgets(fid);
str = fgets(fid);
str = fgets(fid);
nvert = sscanf(str,'%*s %d %*s', 1);

% read vertices
[A,cnt] = fscanf(fid,'%f %f %f', 3*nvert);
if cnt~=3*nvert
    samsrf_disp('WARNING: Problem in reading vertices.');
end
A = reshape(A, 3, cnt/3);
vertex = A;

% read polygons
str = fgets(fid);
str = fgets(fid);

info = sscanf(str,'%c %*s %*s', 1);

if((info ~= 'P') && (info ~= 'V'))
    str = fgets(fid);    
    info = sscanf(str,'%c %*s %*s', 1);
end

if(info == 'P')
    
    nface = sscanf(str,'%*s %d %*s', 1);

    [A,cnt] = fscanf(fid,'%d %d %d %d\n', 4*nface);
    if cnt~=4*nface
        samsrf_disp('WARNING: Problem in reading faces.');
    end

    A = reshape(A, 4, cnt/4);
    face = A(2:4,:)+1;
end

if(info ~= 'P')
    face = 0;
end

% read vertex ids
if(info == 'V')
    
    nv = sscanf(str,'%*s %d %*s', 1);

    [A,cnt] = fscanf(fid,'%d %d \n', 2*nv);
    if cnt~=2*nv
        samsrf_disp('WARNING: Problem in reading faces.');
    end

    A = reshape(A, 2, cnt/2);
    face = repmat(A(2,:)+1, 3, 1);
end

if((info ~= 'P') && (info ~= 'V'))
    face = 0;
end

% read POINT_Data
str = fgets(fid);

npoints = sscanf(str, '%*s %d', 1);

str = fgets(fid);
str = fgets(fid);

[A, cnt] = fscanf(fid, '%f', npoints);

points = A;

fclose(fid);

vertex = vertex';
face = face';

return

