function SamSrfAnalysis(Algorithm, ModelJson, DataWildcard, varargin); %Roi, SurfFolder, Normalise, Average)
%
% Standalone command line tool for analysis. Run without inputs for help information.
%
% 08/10/2025 - At the beginning was the word (DSS)
% 09/10/2025 - Further updates & command line help section (DSS) 
% 10/10/2025 - Now allows flexible inputs & can save results in different folder than data (DSS)

%% Help section for command line
samsrf_clrscr;
if nargin == 0
    samsrf_disp(['Usage:' newline ...
                 newline ... 
                 '   SamSrfAnalysis Algorithm ModelJson DataPath -outpath= -roi= -surf= -norm=1 -aver=1' newline ...
                 newline ...
                 'Important Notes:' newline ...
                 '  The tool must be run from the folder where the executable lives!' newline ...
                 '  The first three arguments are mandatory & must be defined in this order!' newline ...
                 '  All paths must be defined explicitly or relative to where the -data- files are!' newline ...
                 newline ...
                 '      Algorithm:  prf / rc / cf (for forward-model pRF, reverse-correlation pRF, or reverse-correlation CF)' newline ...
                 '      ModelJson:  JSON file with the model description but do not include extension.' newline ...
                 '      DataPath:   Path & filename (allows wildcards) to the data files, which can be in NII or GII format' newline ...
                 newline ...
                 'The remaining input arguments are optional & can be defined in any order.' newline ...
                 'They must include the equation symbol = to define the input:' newline ...
                 newline ...
                 '      -outpath:   Path where you want to save the results files, if not the data path.' newline ... 
                 '      -roi:       Name of ROI definition to restrict analysis without file extension.' newline ... 
                 '                  Must be FreeSurfer LABEL for surface analysis or NII binary mask for volumetric analysis' newline ...
                 '                  (Can also be -roi= if you don''t want to provide a ROI)' newline ...
                 '      -surf: 		Path to the FreeSurfer surface folder (can be -surf= for no surf)' newline ...
                 '      -norm:  	1 (default) to detrend & z-score time series in each run or 0 for no normalisation' newline ...
                 '      -aver: 		1 (default) if you want to average time series across runs or 0 for concatenating runs' newline ...
                 newline ... 
                 'Examples:' newline ...
                 newline ... 
                 '   SamSrfAnalysis prf 2dGaussian /data/001/func/lh_*.gii -roi=lh_occ -surf=../surf' newline ...
                 '      Runs a 2D Gaussian forward-model pRF analysis specified by 2dGaussian.json of the data in /data/001/lh_*.gii,' newline ...
                 '      restricted to ROI lh_occ.label & using the surface meshes in ../surf' newline ...
                 newline ... 
                 '   SamSrfAnalysis prf 2dGaussian /data/001/func/lh_*.gii' newline ...
                 '      Runs the same analysis as above but simply using the data contained in the lh_*.gii files' newline ...
                 '      but without treating it as a surface & not restricted to any ROI.' newline ...
                 newline ... 
                 '   SamSrfAnalysis rc RevCor /data/001/func/lh_*.gii -roi= -surf=../surf/lh.inflated.gii' newline ...
                 '      Runs a reverse-correlation pRF analysis specified by RevCor.json on the same data' newline ...
                 '      but uses the inflated surface provided by the GII & not restricted to any ROI' newline ...
                 newline ... 
                 'Demo using Example Dataset:' newline ... 
                 '		Download the example dataset called X001 from the SamSrf OSF repository.' newline ... 
                 '		For convenience, copy the 3 example JSON files from SamSrf/Models/Json into the X001/prf folder.' newline ... 
                 '		You can then run the 3 standard analyses as follows.' newline ... 
                 newline ... 
                 ]);
    return
end

%% Determine input parameters
switch lower(Algorithm)
    case 'prf'
        samsrf_disp('Forward-modelling pRF analysis');
    case 'rc'
        samsrf_disp('Reverse-correlation pRF analysis');
    case 'cf'
        samsrf_disp('Reverse-correlation CF analysis');
end
disp([' Specified in ' ModelJson '.JSON']);
[DataPath, dn, de] = fileparts(DataWildcard);
if isempty(DataPath)
	DataPath = './';
end
disp([' - Data files are in: 	' DataPath])
[ResultsPath, Roi, SurfFolder, Normalise, Average] = InputParams(varargin);

%% Ensure path
SamSrfDefs = LoadSamSrfDefaults; % Dummy call to ensure paths
HomeFolder = pwd;
if ~isempty(DataPath)
	cd(DataPath);
end
DataWildcard = [dn de];

%% Convert data to Srf
Srf = [];
% What file format?
if contains(lower(DataWildcard), '.nii') 
    % Volumetric data
    samsrf_disp('Converting volumetric NII files to Srf...');
    Srf = samsrf_vol2mat(strrep(lower(DataWildcard), '.nii', ''), Roi, Normalise, Average, true); % ROI is used for creating data file
    samsrf_newline;
    Roi = ''; % Clear now
elseif contains(lower(DataWildcard), '.gii') 
    % GIfTI format
    if ~isempty(SurfFolder)
        if contains(lower(DataWildcard), 'bi')
            LSurfFld = [SurfFolder filesep 'lh'];
            RSurfFld = [SurfFolder filesep 'rh'];
        elseif contains(lower(DataWildcard), 'lh')
            SurfFolder = [SurfFolder filesep 'lh'];
        elseif contains(lower(DataWildcard), 'rh')
            SurfFolder = [SurfFolder filesep 'rh'];
        else
            if contains(lower(SurfFolder), '.gii')
                samsrf_disp(['Using surface mesh in ' SurfFolder]);
            else
                samsrf_error('Cannot work out what hemisphere to analyse!');
            end
        end
    end
    if contains(lower(DataWildcard), 'bi')
        % Binocular data 
	    samsrf_disp('Converting left hemisphere GII files to Srf...');
        Lsrf = samsrf_gii2srf(strrep(strrep(lower(DataWildcard), '.gii', ''), 'bi', 'lh'), LSurfFld, Normalise, Average, true, '');
        samsrf_newline;
        samsrf_disp('Converting right hemisphere GII files to Srf...');
        Rsrf = samsrf_gii2srf(strrep(strrep(lower(DataWildcard), '.gii', ''), 'bi', 'rh'), RSurfFld, Normalise, Average, true, '');
        samsrf_newline;
        % Combine hemispheres
        samsrf_disp('Combining hemisphere Srfs...');
        Srf = samsrf_bilat_srf(Lsrf, Rsrf);
        samsrf_newline;
        clear Lsrf Rsrf
    else
        % Hemispheric data 
        samsrf_disp('Converting GII files to Srf...');        
        Srf = samsrf_gii2srf(strrep(lower(DataWildcard), '.gii', ''), SurfFolder, Normalise, Average, true, '');
        samsrf_newline;
    end
else
    samsrf_error('Data format not supported!');
end

%% Run analysis!
if ~isempty(ResultsPath)
	cd(ResultsPath);
end
if isempty(Srf)
    % Should not really happen but just in case
    samsrf_clrscr;
    samsrf_disp('ERROR: Something went wrong with importing data!');
else
    switch lower(Algorithm)
        case 'prf'
            % Forward-modelling pRF analysis
            OutFile = samsrf_fit_prf(ModelJson, Srf, Roi);
        case 'rc'
            % Reverse-correlation pRF analysis
            OutFile = samsrf_revcor_prf(ModelJson, Srf, Roi);
        case 'cf'
            % Reverse-correlation CF analysis
            OutFile = samsrf_revcor_cf(ModelJson, Srf, Roi);
    end
    
    %% Export GII/NII files
    load(OutFile); % Load maps  
    if contains(lower(DataWildcard), '.nii')
        % Export as NIIs 
        samsrf_mat2vol(Srf, OutFile);
    else
        % Export as GIIs
	    usc = strfind(OutFile, '_');
        if strcmpi(Srf.Hemisphere, 'bi')
            % Explore both hemispheres
            [Lsrf, Rsrf] = samsrf_hemi_srfs(Srf);
            samsrf_export_giis(Lsrf, [Lsrf.Hemisphere '_' OutFile(usc(1)+1:end)]);
            samsrf_export_giis(Rsrf, [Rsrf.Hemisphere '_' OutFile(usc(1)+1:end)]);
        else
            % Export single hemisphere
            samsrf_export_giis(Srf, [Srf.Hemisphere '_' OutFile(usc(1)+1:end)]);
        end
    end
    samsrf_done;
end

cd(HomeFolder);


%% Determine input arguments
function [ResultsPath, Roi, SurfFolder, Normalise, Average] = InputParams(Args)

% Defaults
ResultsPath = '';
Roi = '';
SurfFolder = '';
Normalise = 1;
Average = 1;

% Loop thru inputs in case they are different
for i = 1:length(Args)
    % Must start with -
    if Args{i}(1) ~= '-'
    	samsrf_error('Optional input arguments must start with minus sign!');
    end
	% Must contain =
	e = strfind(Args{i}, '=');
    if isempty(e)
    	samsrf_error('Optional input arguments must contain equals sign!');
    end
    % If not empty
    if e < length(Args{i})
        Val = Args{i}(e+1:end);
	    if contains(lower(Args{i}), '-outpath')
            ResultsPath = Val;
        elseif contains(lower(Args{i}), '-roi')
        	Roi = Val;
        elseif contains(lower(Args{i}), '-surf')
        	SurfFolder = Val;
        elseif contains(lower(Args{i}), '-norm')
        	Normalise = str2double(Val);
        elseif contains(lower(Args{i}), '-aver')
        	Average = str2double(Val);
        end
    end
end

% Report back
if isempty(Roi)
	samsrf_disp(' - No region of Interest');
else
	samsrf_disp([' - Region of Interest:	' Roi]);
end
if isempty(SurfFolder)
	samsrf_disp(' - No surface folder');
else
	samsrf_disp([' - Surface folder: 	' SurfFolder]);
end
if Normalise
	samsrf_disp(' - De-trending & z-normalising time series');
else
	samsrf_disp(' - No temporal normalisation');
end
if Average
	samsrf_disp(' - Averaging individual runs, if more than one');
else
	samsrf_disp(' - Concatenating individual runs, if more than one');
end
if isempty(ResultsPath)
	samsrf_disp(' - Results are saved in data folder');
else
	samsrf_disp([' - Path for results: 	' ResultsPath]);
end
samsrf_newline;