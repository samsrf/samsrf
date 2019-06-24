function Vs = samsrf_loadpaths(Paths)
%
% Loads the vertex indeces of FreeSurfer paths in file Paths into Vs.
%
% 09/08/2018 - SamSrf 6 version (DSS)
%

[~,~,ext] = fileparts(Paths);
if strcmpi(ext, '.mat')
    load(Paths);
else
    fid = fopen(Paths);
    cp = textscan(fid, repmat('%n',1,4), 'headerlines',4);
    cp = cell2mat(cp);
    Vs = cp(:,4)+1;

    while size(cp,1) > 0
        cp = textscan(fid, repmat('%n',1,4), 'headerlines',3);
        cp = cell2mat(cp);
        Vs = [Vs; cp(:,4)+1];   
    end

    fclose(fid);
end