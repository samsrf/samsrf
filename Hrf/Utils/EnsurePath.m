function newname = EnsurePath(oldname)
% 
% Ensures that if no path is defined the pathname is to the current folder.
% Use this to avoid random files being loaded from the MatLab path in error.
% (This particularly happens on public machines when people edit the path...)
%
% 20/04/2022 - SamSrf 8 version (DSS)
%

[p, n, e] = fileparts(oldname); % Deconstruct file name
% Path undefined
if isempty(p)
    p = '.';
end
newname = [p filesep n e]; % Reconstruct file name
