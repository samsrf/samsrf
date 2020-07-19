function cR = prf_centre_prf(R)
%
% cR = prf_centre_prf(R)
%
% Centres the pRF in the r-map R on zero. It also sets all NaNs to zero. 
%
% Note: This function requires a Srf analysed with reverse correlation.
%
% 19/07/2020 - SamSrf 7 version (DSS)
%

% Transformation parameters
dims = size(R,1);  % Dimensions of r-map
mR = max(R(:));  % Peak correlation
[r c] = find(R == mR);  % Peak matrix indeces
[mr mc] = meshgrid(1:dims, 1:dims);  % Matrix indeces
dx = mc-c;  % Column indeces relative to peak
dy = mr-r;  % Row indeces relative to peak

% Centre the matrix on zero
cR = zeros(dims,dims);
for i = 1:numel(dx) 
    if abs(dx(i)) < 25 && abs(dy(i)) < 25 
        cR(dy(i)+25, dx(i)+25) = R(r+dy(i),c+dx(i)); 
    end
end

% Remove NaNs
cR(isnan(cR)) = 0;
