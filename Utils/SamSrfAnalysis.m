function SamSrfAnalysis(Algorithm, ModelJson, DataWildcard, Roi, SurfFolder, Normalise, Average)
%
% Standalone command line tool for analysis.
%
%	Algorithm:	
% 		prf : 	Forward modelling pRF
% 		rc  :	Reverse correlation pRF
%		cf  :	Reverse correlation CF
%
% 	ModelJson: 	
% 		JSON file with model specification but don't include extension!
%       Pathname is relative to where the data are unless specified explicitly.
%
%	DataWildcard:
%		Which data files to analyse, e.g. 
% 		lh*.gii    bi_Run*.gii    Run*.nii
%
%       The file extension determines whether analysis is done on surface or in volume.
%
%       This wildcard should also include the pathname of the data files.
%       All analysis will be conducted in that folder. So all the other
%       files you specify in this script must be defined relative to there.
%
%       To analyse bilateral surface data the wildcard must contain 'bi',
%        so ensure that the file name does not inadvertently contain 'bi'.
%        In this case, it will search for the files using this wildcare but
%        replacing 'bi' with 'lh' & 'rh', respectively. You must ensure
%        that this does not find any files you don't want to include!
%
%	Roi: 
%		Region of interest to analyse as FreeSurfer .label for surface analysis or as NII for volumetric analysis.
%		In either case, you only provide the file name but not the extension!
%		Defaults to '-noroi' (In Matlab, you can also use '')
%       Pathname is relative to where the data are unless specified explicitly.
%
% 	SurfFolder:
%		Pathname of subject's surf folder. Defaults to '-nosurf' (In Matlab, you can also use '') 
% 		In volumetric analysis this is ignored because no surf folder is needed. 
%       Pathname is relative to where the data are unless specified explicitly.
%
%	Normalise:
%		Toggle 1 (default) or 0 to detrend & z-score
%
%	Average:
%		Toggle 1 (default) for averaging or 0 for concatenating runs.
%
% 08/10/2025 - At the beginning was the word (DSS)
% 09/10/2025 - Further updates & command line help section (DSS) 

%% Help section for command line
if nargin == 0
    samsrf_disp(['Usage:' newline ...
                 newline ... 
                 '   SamSrfAnalysis Algorithm ModelJson DataPath Roi=-noroi SurfFolder=-nosurf Normalise=1 Average=1' newline ...
                 '      Algorithm:  prf / rc / cf (for forward-model pRF, reverse-correlation pRF, or reverse-correlation CF)' newline ...
                 '      ModelJson:  JSON file with the model description. Path is relative to where the data files are ' newline ...
                 '      DataPath:   Path & filename (allows wildcards) to the data files, which can be in NII or GII format' newline ...
                 '      Roi:        Name of ROI definition to restrict analysis without file extension.' newline ... 
                 '                  Must be FreeSurfer LABEL for surface analysis or NII binary mask for volumetric analysis' newline ...
                 '                  (Can also be -noroi if you don''t want to provide a ROI)' newline ...
                 '      SurfFolder: Path to the FreeSurfer surface folder (can be -nosurf)' newline ...
                 '      Normalise:  1 (default) to detrend & z-score time series in each run or 0 for no normalisation' newline ...
                 '      Average:    1 (default) if you want to average time series across runs or 0 for concatenating runs' newline ...
                 newline ... 
                 'Examples:' newline ...
                 newline ... 
                 '   SamSrfAnalysis prf 2dGaussian /data/001/func/lh_*.gii lh_occ ../surf' newline ...
                 '      Runs a 2D Gaussian forward-model pRF analysis specified by 2dGaussian.json of the data in /data/001/lh_*.gii,' newline ...
                 '      restricted to ROI lh_occ.label & using the surface meshes in ../surf' newline ...
                 newline ... 
                 '   SamSrfAnalysis prf 2dGaussian /data/001/func/lh_*.gii' newline ...
                 '      Runs the same analysis as above but simply using the data contained in the lh_*.gii files' newline ...
                 '      but without treating it as a surface & not restricted to any ROI' newline ...
                 newline ... 
                 '   SamSrfAnalysis rc RevCor /data/001/func/lh_*.gii -noroi ../surf/lh.inflated.gii' newline ...
                 '      Runs a reverse-correlation pRF analysis specified by RevCor.json on the same data' newline ...
                 '      but uses the inflated surface provided by the GII & not restricted to any ROI' newline ...
                 newline ... 
                 ]);
    return
end

%% Ensure inputs are formatted correctly & defaults set
if nargin < 4
    Roi = '';
end
if nargin < 5
    SurfFolder = '';
end
if nargin < 6
    Normalise = true;
end
if nargin < 7 
    Average = true;
end
if strcmpi(Roi, '-noroi')
    Roi = '';
end
if strcmpi(SurfFolder, '-nosurf')
    SurfFolder = '';
end
if ischar(Normalise)
    Normalise = str2double(Normalise);
end
if ischar(Average)
    Average = str2double(Average);
end

%% Ensure path
SamSrfDefs = LoadSamSrfDefaults;
[dp dn de] = fileparts(DataWildcard);
cf = pwd;
cd(dp);
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
            SurfFolder = LSurfFld;
        elseif contains(lower(DataWildcard), 'rh')
            SurfFolder = RSurfFld;
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
    
    % Export GII/NII files
    load(OutFile); % Load maps  
    if contains(lower(DataWildcard), '.nii')
        % Export as NIIs 
        samsrf_mat2vol(Srf, OutFile);
    else
        % Export as GIIs
        if strcmpi(Srf.Hemisphere, 'bi')
            % Explore both hemispheres
            [Lsrf, Rsrf] = samsrf_hemi_srfs(Srf);
            samsrf_export_giis(Lsrf, [Lsrf.Hemisphere '_' OutFile(4:end)]);
            samsrf_export_giis(Rsrf, [Rsrf.Hemisphere '_' OutFile(4:end)]);
        else
            % Export single hemisphere
            samsrf_export_giis(Srf, [Srf.Hemisphere '_' OutFile(4:end)]);
        end
    end
    samsrf_done;
end

cd(cf);

