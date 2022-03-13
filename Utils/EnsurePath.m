function newname = EnsurePath(oldname)
% 
% Ensures that if no path is defined the pathname is to the current folder.
% Use this to avoid random files being loaded from the MatLab path because
% some idiot to put their data files on the path on a public computer...
%

[p, n, e] = fileparts(oldname); % Deconstruct file name
% Path undefined
if isempty(p)
    p = '.';
end
newname = [p filesep n e]; % Reconstruct file name
