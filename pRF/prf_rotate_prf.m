function tR = prf_rotate_prf(R, alpha)
%
% tR = prf_rotate_prf(R, [alpha=Inf])
%
% Rotate the reverse correlation pRF profile in the r-map R by angle alpha. 
% It also sets all NaNs to zero. If alpha is set to Inf, it rotates the pRF 
% so that the peak is in the 3 o'clock position.
%
% 20/04/2022 - SamSrf 8 version (DSS)
%

if nargin < 2
    alpha = Inf;
end

% Transformation parameters
dims = size(R,1);  % Dimensions of r-map

% Rotate the matrix
if isinf(alpha)
    mR = max(R(:));  % Peak correlation
    [r c] = find(R == mR, 1);  % Peak matrix indeces
    dx = c - dims/2;  % Column indeces relative to centre
    dy = r - dims/2;  % Row indeces relative to centre
    theta = atan2(-dy,dx) / pi*180;  % Polar angle of peak
    tR = imrotate(R, -theta, 'bicubic', 'crop'); % Rotate backwards to 3 o'clock
else
    tR = imrotate(R, alpha, 'bicubic', 'crop'); % Rotate by alpha
end

% Remove NaNs
tR(isnan(tR)) = 0;
