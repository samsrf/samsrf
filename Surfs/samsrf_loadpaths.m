function [Vs, Paths] = samsrf_loadpaths(PathFileName)
%
% Loads the vertex indeces of FreeSurfer paths in file PathName into Vs.
% The second output argument Paths contains the Paths cell array if PathFileName
% is a delineation file from DelineationTool. Otherwise this is empty.
%
% 20/04/2022 - SamSrf 8 version (DSS)
%

if ~isempty(PathFileName)
    [~,~,ext] = fileparts(PathFileName);
    if strcmpi(ext, '.mat')
        load(EnsurePath(PathFileName));
    else
        fid = fopen(PathFileName);
        cp = textscan(fid, repmat('%n',1,4), 'headerlines',4);
        cp = cell2mat(cp);
        Vs = cp(:,4)+1;

        while size(cp,1) > 0
            cp = textscan(fid, repmat('%n',1,4), 'headerlines',3);
            cp = cell2mat(cp);
            Vs = [Vs; cp(:,4)+1];   
        end

        fclose(fid);
        Paths = {};
    end
else
    Vs = [];
    Paths = {};
end