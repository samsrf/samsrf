function samsrf_srf2label(srfdata, labelname, valnum, Vx)
%
% samsrf_srf2label(srfdata, labelname, [valnum, Vx])
%
% Saves srfdata as a Freesurfer label in labelname.
% This can be either a filename or a Srf structure.
% You can define with row of data to convert with valnum (default is 1).
% You can also define the vertex indeces you want to include in your label
% to constrain its size (default is all vertices).
%
% 20/04/2022 - SamSrf 8 version (DSS)
% 24/10/2023 - Added support for M/EEG and volumetric data files (DSS)
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

% If vertex indeces undefined
if nargin < 4
    nver = size(Srf.Data,2);
    Vx = 1:nver;
else
    if size(Vx,2) == 1
        Vx = Vx';
    end
    nver = length(Vx);
end

% Change NaNs to 0
Srf.Data(isnan(Srf.Data)) = 0;

% If volumetric data
if size(Srf.Vertices,2) == 1
    Srf.Vertices = repmat(Srf.Vertices,[1 3]);
end
% If M/EEG data
if size(Srf.Vertices,2) == 1
    Srf.Vertices = [Srf.Vertices zeros(size(Srf.Vertices,1),1)];
end

fid = fopen([labelname '.label'], 'w');
fprintf(fid, '#! Converted from SamSurfer.\n');
fprintf(fid, '%d\n', nver);
for v = Vx
    fprintf(fid, '%d %5.3f %5.3f %5.3f %f\n', v-1, Srf.Vertices(v,1), Srf.Vertices(v,2), Srf.Vertices(v,3), Srf.Data(valnum,v));
end
fclose(fid);
disp(['Saved ' labelname '.label']);