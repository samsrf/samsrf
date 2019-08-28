function AvgPrf = samsrf_meanprf(SrfDv, SrfIv, Roi, Eccs, Thr)
%
% AvgPrf = samsrf_meanprf(SrfDv, SrfIv, Roi, Eccs, [Thr=0.01])
%
% Calculates the mean pRF profile in SrfDv within the eccentricity band Eccs 
%  of region of interest Roi, based on the eccentricities in SrfIv.
%
% Thr is optional and defines the R^2 threshold for vertices to include.
%
% This function requires Srf's analysed with reverse correlation!
%
% 09/08/2018 - SamSrf 6 version (DSS)
% 18/08/2019 - Corrected variable names in the help section (DSS) 
%

if nargin < 5
    Thr = 0.01;
end

%% Expand Srfs if necessary
SrfDv = samsrf_expand_srf(SrfDv);
SrfIv = samsrf_expand_srf(SrfIv);

%% Load ROI label
rv = samsrf_loadlabel(Roi);
if ~isnan(rv)
    SrfDv.Data = SrfDv.Data(:,rv);
    SrfDv.Rmaps = SrfDv.Rmaps(:,rv);
    SrfIv.Data = SrfIv.Data(:,rv);
else
    error(['ROI ' Roi ' not found!']);
end

%% Remove rubbish & calculate measures
sv = SrfIv.Data(1,:) >= Thr(1);
SrfDv.Data = SrfDv.Data(:,sv);
SrfDv.Rmaps = SrfDv.Rmaps(:,sv);
SrfIv.Data = SrfIv.Data(:,sv);
E = sqrt(SrfIv.Data(2,:).^2 + SrfIv.Data(3,:).^2);

%% Restrict eccentricity range
SrfDv.Data = SrfDv.Data(:, E >= Eccs(1) & E < Eccs(2)); 
SrfDv.Rmaps = SrfDv.Rmaps(:, E >= Eccs(1) & E < Eccs(2)); 

%% Calculate average profile
dims = sqrt(size(SrfDv.Rmaps,1));
Prf = NaN(dims, dims, size(SrfDv.Rmaps,2));
[x,y] = meshgrid(1:dims, 1:dims);
for i = 1:size(SrfDv.Rmaps,2)
    R = prf_contour(SrfDv, i, x, y); %% CURRENTLY DOESN'T WORK!!!
    if sum(double(R(:) ~= 0)) > 0
        cR = prf_centre_prf(R);
        Prf(:,:,i) = cR;
    end
end
Prf = Prf(:,:,~isnan(Prf(1,1,:)));
AvgPrf = mean(Prf,3) / size(Prf,3);

