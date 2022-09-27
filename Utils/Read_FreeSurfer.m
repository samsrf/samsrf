function data = Read_FreeSurfer(fname)
%data = Read_FreeSurfer(fname)
%
% Reads the contents of FreeSurfer ASCII files.
% The matrix is organized into five colums:
%   1:      vertex     vertex index number
%   2-4:    x,y,z      coordinates of vertices in world space
%   5:      values     contains overlay values etc
%
% If the file name has the extension '.label' the program reads a label
% (in which the first 2 lines are header information and thus ignored).
% If no file extension is given, the program assumes no headers.
%
% 13/03/2022 - Now ensures that files aren't loaded from path (DSS)
%              Produces an error if file cannot be loaded now (DSS)
%              Reports full pathname of loaded file (DSS)
% 20/04/2022 - SamSrf 8 version (DSS)
% 27/09/2022 - Removed loading error as this breaks dependent functions (DSS)
%

fname = EnsurePath(fname); % Ensure correct path
[p, n, e] = fileparts(fname); % Deconstruct file name
fid = fopen(fname);
try 
    if strcmpi(e,'.label')
        c = textscan(fid, repmat('%n',1,6), 'headerlines', 2);
    else
        c = textscan(fid, repmat('%n',1,5));
    end
    disp(['Loaded ' fname]);
    fclose(fid);
    data = cell2mat(c);
catch
    data = NaN;
end

