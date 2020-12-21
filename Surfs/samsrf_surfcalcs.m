function Srf = samsrf_surfcalcs(InSrf, Roi, R2Thrsh, Eccens, Method, Fwhms)
%
% Srf = samsrf_surfcalcs(InSrf, [Roi='', R2Thrsh=0.01, Eccens=[0 Inf], Method='S', Fwhms=[10 3]])
% 
% Runs a range of surface calculations on a pRF map:
% First, it smoothes the map with a large kernel and calculates the field sign.
% Next, it smoothes the map again with a smaller kernel and calculates the CMF.
% Finally, it rectifies the now smoothed field sign map and returns the Srf.
%
%   InSrf:      Input Srf
%   Roi:        ROI label to restrict analysis
%   R2Thrsh:    R^2 threshold
%   Eccens:     Eccentricity range (1x2 vector for minimum and maximum)
%   Method:     Smoothing method (only first letter needed) from
%                'Sphere':      Distance determined on the spherical surface (SamSrf standard)
%                'Geodeosic':   Distance determined by geodesic steps (FreeSurfer standard)
%                'Dijkstra':    Distance determined by Dijkstra's geodesic distance (best method, but slow)
%   Fwhms:      Smoothing kernels (1x2 vector, first kernel is for field sign, second is for everything else)
%
% 05/07/2020 - SamSrf 7 version (DSS)
% 21/12/2020 - Inconsequential bug fix (DSS)
% 22/12/2020 - Fixed typo in help section (DSS)
%

if nargin < 2
    Roi = '';
end
if nargin < 3
    R2Thrsh = 0.01;
end
if nargin < 4
    Eccens = [0 Inf];
end
if nargin < 5
    Method = 'S';
end
if nargin < 6
    Fwhms = [10 3];
end

% Expand if needed & store raw data
[Srf, RoiVx] = samsrf_expand_srf(InSrf);
Srf = samsrf_smooth_sphere(Srf, 0); % Remove previous smoothing
Raw = Srf.Data; % Save raw data

% Filter as required
Ecc = sqrt(Srf.Data(2,:).^2 + Srf.Data(3,:).^2); % Eccentricities
Srf.Data(1,Ecc < Eccens(1) | Ecc > Eccens(2)) = 0; % Remove anything outside eccentricity range

% Smooth heavily & calculate field sign
switch upper(Method(1))
    case 'S'
        Srf = samsrf_smooth_sphere(Srf, Fwhms(1), Roi, R2Thrsh);
    case 'G'
        Srf = samsrf_smooth_georoi(Srf, Fwhms(1), Roi, R2Thrsh);
    case 'D'
        Srf = samsrf_smooth_dijkstra(Srf, Fwhms(1), Roi, R2Thrsh);
end
Srf = samsrf_fieldsign(Srf, 5, Roi, R2Thrsh);
LocFs = size(Srf.Data,1); % Which row is field sign?

% Smooth again, calculate CMF & rectify smoothed field sign
Srf.Data = [Raw; Srf.Data(end,:)]; % Add field sign to original raw data 
Srf = rmfield(Srf, 'Raw_Data'); % Clear heavily smoothed data as it won't be needed any more
switch upper(Method(1))
    case 'S'
        Srf = samsrf_smooth_sphere(Srf, Fwhms(2), Roi, R2Thrsh);
    case 'G'
        Srf = samsrf_smooth_georoi(Srf, Fwhms(2), Roi, R2Thrsh);
    case 'D'
        Srf = samsrf_smooth_dijkstra(Srf, Fwhms(2), Roi, R2Thrsh);
end
Srf = samsrf_cortmagn(Srf, Roi); % Calculate CMF
Srf.Data(LocFs,:) = sign(Srf.Data(LocFs,:)); % Rectify field sign

% Compress data if needed
Srf = samsrf_compress_srf(Srf, RoiVx);
