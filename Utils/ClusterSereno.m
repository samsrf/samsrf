function ClusterSereno(Hemis, RoiPath)
%
% After generating the Sereno Atlas ROIs (see atlas template on OSF) 
% you can use this function to cluster these ROIs into more manageable 
% larger ROI clusters based on anatomical criteria. 
% 
% Hemis defines the hemisphere ('L' or 'R'). Omit for bilateral labels.
% 
% RoiPath determines the path with the Sereno labels, which defaults to '.' 
%
% 25/08/2022 - Written (PWBU & DSS)
% 17/09/2022 - Fixed bug when directing to different path (DSS)
% 14/12/2023 - Bugfix when providing a path (DSS)
% 09/06/2025 - Added cluster for dlPFC (DSS)
%

if nargin == 0
    Hemis = '';
end
if nargin <= 1
    RoiPath = '.';
end
if upper(Hemis) == 'L'
    Hemis = 'lh_';
elseif upper(Hemis) == 'R'
    Hemis = 'rh_';
else
    Hemis = '';
end

homecd = pwd; 
cd(RoiPath);

%% Define complexes
cluster_struct = struct;

complexes = {'V1' ,'V2', 'V3', 'V3AB', 'VO', 'LO', 'MT', 'V6', 'ProSt', 'STVPIC', ...
    'V4', 'PeriFus', 'FEFs', 'IntrParSul', 'AntPre', 'EarlyIntP' 'dlPFC'};

% V1
cluster_struct.(complexes{1}) = {'V1_lower', 'V1_upper'};

% V2
cluster_struct.(complexes{2}) = {'V2_lower', 'V2_upper'};

% V3
cluster_struct.(complexes{3}) = {'V3' 'VP'};

% V3A+B
cluster_struct.(complexes{4}) = {'V3A', 'V3B' 'DI'};

% VO+
cluster_struct.(complexes{5}) = {'VO2', 'VO1' 'V8'};

% LO+
cluster_struct.(complexes{6}) = {'LO1', 'LO2', 'LO3' 'PGp' 'OPA'};

% MT+
cluster_struct.(complexes{7}) = {'MSTd', 'MSTv', 'MT_lower', 'MT_upper', 'MTc' 'FSTd'};

% V6+
cluster_struct.(complexes{8}) = {'V6', 'V6A' 'aPOS', 'POm'};

% Prostriata
cluster_struct.(complexes{9}) = {'ProS1', 'ProS2' };

% STV+PIC
cluster_struct.(complexes{10}) = {'STV1', 'STV2' '7b_PICv', '7b_PICvs'};

% V4
cluster_struct.(complexes{11}) = {'V4v' 'hV4'} ;

% Perifusiform
cluster_struct.(complexes{12}) = {'FFC' 'VVC' 'PH'};

% FEFs
cluster_struct.(complexes{13}) = {'FEF' 'dmFEF'} ;

% Intraparietal sulcus area
cluster_struct.(complexes{14}) = {'IPS4', 'IPS5'  'VIP1v' 'VIP2v' 'VIP1vs' 'VIP2vs'};

% Anterior pre-cuneus
cluster_struct.(complexes{15}) = {'aPCu1', 'aPCu2'};

% Early intraparietal
cluster_struct.(complexes{16}) = {'LIP0', 'LIP1' 'V7' 'cIPS' 'PEc'};

% Dorso-lateral prefrontal cortex
cluster_struct.(complexes{17}) = {'DLPFC', 'DLPFCA' 'DLPFCaud'};

%% Add hemisphere prefix if needed
for i = 1:length(complexes)
    for j = 1:length(cluster_struct.(complexes{i}))
        cluster_struct.(complexes{i}){j} = [Hemis cluster_struct.(complexes{i}){j}];
    end
end

%% Create cluster ROIs
mkdir('Sereno_clusters'); 

for i = 1:length(complexes)
    samsrf_combine_labels(cluster_struct.(complexes{i}),... 
        ['Sereno_clusters', filesep, Hemis, complexes{i}]);
end

cd(homecd);

end

