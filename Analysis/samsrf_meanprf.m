function AvgPrf = samsrf_meanprf(SrfDv, SrfIv, Roi, Eccs, Thrsh, Nrmls)
%
% AvgPrf = samsrf_meanprf(SrfDv, SrfIv, Roi, Eccs, [Thrsh=0.01, Nrmls='None'])
%
% Calculates the mean pRF profile in SrfDv within the eccentricity band defined 
%  by the 1x2 vector Eccs in region Roi, based on eccentricities in SrfIv. 
%
% Thrsh is optional and defines the R^2 threshold for vertices to include.
%  If this is a 1x2 vector, the second entry determines the proportion
%  relative to the peak response above which to include in the average.
%  This defaults to 0, which means all pixels are included.
%
% Nrmls is optional and determines what spatial normalisation is applied to 
%  the pRF profile before averaging (only the first letter is required):
%   'None':     No normalisation
%   'Centre':   Centre pRFs so their peak is at origin
%   'Rotate':   Rotate pRFs to fall on 3 o'clock axis
%   'Both':     Rotate pRFs and then centre them 
%
% Note: This function requires Srf's analysed with reverse correlation!
%
% 09/08/2018 - SamSrf 6 version (DSS)
% 18/08/2019 - Corrected variable names in the help section (DSS) 
% 14/05/2020 - Added alignment to 3 o'clock (DSS)
%              Made centering of pRF profile optional (DSS)
%              Now uses nanmean instead of mean (DSS)
%

if nargin < 5
    Thrsh = 0.01;
end
if nargin < 6
  Nrmls = 'None';
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
sv = SrfIv.Data(1,:) >= Thrsh(1);
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
for i = 1:size(SrfDv.Rmaps,2)
    R = prf_contour(SrfDv, i); 
    % Thrshold intensities?
    if length(Thrsh) > 1
        R(R < nanmax(R(:)) * Thrsh(2)) = 0;
    end
    % Spatial normalisations?
    if sum(double(R(:) ~= 0)) > 0
        switch upper(Nrmls(1))
            case 'N'
                % Do nothing
            case 'C' % Centre pRFs only
                R = prf_centre_prf(R); 
            case 'R' % Rotate pRFs only
                R = prf_rotate_prf(R);
            case 'B' % Rotate & then centre pRFs
                R = prf_rotate_prf(R);
                R = prf_centre_prf(R);
            otherwise
                error('Invalid normalisation method specified!');
        end
        Prf(:,:,i) = R;
    end
end
AvgPrf = nanmean(Prf,3);

