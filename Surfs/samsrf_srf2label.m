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
% 19/07/2020 - SamSrf 7 version (DSS)
% 08/10/2021 - Turned command line report back on (DSS) 
%

if ischar(srfdata)
    load(srfdata);
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

fid = fopen([labelname '.label'], 'w');
fprintf(fid, '#! Converted from SamSurfer.\n');
fprintf(fid, '%d\n', nver);
for v = Vx
    fprintf(fid, '%d %5.3f %5.3f %5.3f %f\n', v-1, Srf.Vertices(v,1), Srf.Vertices(v,2), Srf.Vertices(v,3), Srf.Data(valnum,v));
end
fclose(fid);
disp(['Saved ' labelname '.label']);