function Vs = samsrf_loadpaths(Paths)
%
% Loads the vertex indeces of FreeSurfer paths in file Paths into Vs.
%
% 19/07/2020 - SamSrf 7 version (DSS)
%

if ~isempty(Paths)
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
else
    Vs = [];
end