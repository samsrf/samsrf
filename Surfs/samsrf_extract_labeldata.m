function Data = samsrf_extract_labeldata(SrfName, Roi)
%
% Data = samsrf_extract_labeldata(SrfName, Roi)
%
% Returns the data in SrfName that is in label Roi. SrfName can either be a 
% Srf variable or a surface data file. If a file name is given, the function
% also expands the Srf if necessary. The function always extracts all rows 
% from Srf.Data, so if you want the unsmoothed data you must first replace
% Srf.Data with Srf.Raw_Data.
%
% 19/07/2020 - SamSrf 7 version (DSS)
%

% Load label vertices
Vs = samsrf_loadlabel(Roi);

% Load surface data if a file is given
if ischar(SrfName)
    load(EnsurePath(SrfName));
else   
    Srf = SrfName;
end

% Expand Srf if necessary
Srf = samsrf_expand_srf(Srf);
    
% Extract ROI data
Data = Srf.Data(:,Vs);