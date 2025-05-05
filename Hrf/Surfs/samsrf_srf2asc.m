function samsrf_srf2asc(srfdata, ascname, srfdir, valnum)
%
% samsrf_srf2asc(srfdata, ascname, srfdir, [valnum])
%
% Saves srfdata as a Freesurfer ASCII overlay file in ascname.
% The input can be either a filename or a Srf structure.
%
% The output name will automatically be prefixed by '<hemi>.'
%
% You must define the surf folder in srfdir because the overlay must be
% saved there in order to be useful to FreeSurfer.
%
% You can define with row of data to convert with valnum (default is 1).
%
% This function is essentially the same as samsrf_srf2label but it instead
% saves a FreeSurfer ASCII file that can be converted into binary format
% using mris_convert. It always saves the entire hemisphere.
%
% 20/04/2022 - SamSrf 8 version (DSS)
%

if ischar(srfdata)
    load(EnsurePath(srfdata));
else
    Srf = srfdata;
end
Srf = samsrf_expand_srf(Srf);

if nargin < 3
    valnum = 1;
end

% Number of vertices 
nver = size(Srf.Vertices,1);

% Change NaNs to 0
Srf.Data(isnan(Srf.Data)) = 0;

% Write ASCII data
fid = fopen([srfdir filesep Srf.Hemisphere '.' ascname '.asc'], 'w');
for v = 1:nver
    fprintf(fid, '%d %5.3f %5.3f %5.3f %f\n', v-1, Srf.Vertices(v,1), Srf.Vertices(v,2), Srf.Vertices(v,3), Srf.Data(valnum,v));
end
fclose(fid);
samsrf_disp(['Saved ' srfdir filesep Srf.Hemisphere '.' ascname '.asc.']);
