function [vertices, faces, border] = fs_read_patch(fpatch, fsurf)
% [vertices, faces, border] = fs_read_patch(fpatch, fsurf)
%
% Inputs
%   fpatch      [string] name of binary patch file in FreeSurfer format
%   fsurf       [string] name of a matching surface file (e.g. lh.orig)            
%
% Output
%   vertices    [matrix] X,Y,Z coordinates of vertex locations
%   faces       [matrix] m x 3 matrix of triangular faces
%   border      [vector] 1 = border vertex
%                        0 = not a border vertex
%
% % Changelog %
%
% 04/04/2017    Written
%
% Ivan Alvarez
% FMRIB Centre, University of Oxford

% Read patch file (modified from read_patch)
fid = fopen(fpatch,'r');
if (fid == -1)
   error('could not open file %s', fpatch) ;
end
ver = fread(fid, 1, 'int', 0, 'b');
if (ver ~= -1)
   error('incorrect version # %d (not -1) found in file',ver) ;
end
patch.npts = fread(fid, 1, 'int', 0, 'b') ;
for i = 1:patch.npts
    ind = fread(fid, 1, 'int', 0, 'b') ;
    border = ind < 0;
    if (ind < 0)
       ind = -ind - 1 ;
    else
       ind = ind - 1 ;
    end
    patch.ind(i) = ind ;
    patch.x(i) = fread(fid, 1, 'float', 0, 'b') ;
    patch.y(i) = fread(fid, 1, 'float', 0, 'b') ;
    patch.z(i) = fread(fid, 1, 'float', 0, 'b') ;
    patch.border(i) = border ;
end
fclose(fid);

% Read surface file
[surf.vertex, surf.faces] = fs_read_surf(fsurf);

% Matlab counts from 1, not 0
patch.ind = patch.ind + 1;

% Vertices in the patch
vertices = nan(size(surf.vertex));
vertices(patch.ind, :) = [patch.x; patch.y; patch.z]';

% Faces are the same
faces = surf.faces;

% Border vertices
border = patch.border;

% Done
%