function fs_bin2asc(subj)
% Converts the surface binaries in subj/surf to ASC format so SamSrf can 
% read them. This is a wrapper for a shell script and should work from
% anywhere as long as SamSrf is on the Matlab path.
%
% The input subj can contain a list of subjects.
%

% Check if this is Windows
if ispc
    error('I''m sorry but this won''t work under Windows!');
end

% Determine TCL script folder  
tcldir = fileparts(mfilename('fullpath'));

% Run script
unix([tcldir filesep 'fs_bin2asc.sh ' subj]);