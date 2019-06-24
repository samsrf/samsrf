function Srf = samsrf_rotate(InSrf, theta)
%
% Srf = samsrf_rotate(InSrf, theta)
%
% Rotates the pRFs in Srf by angle theta.
%

%% Expand if necessary
InSrf = samsrf_expand_srf(InSrf);

%% Coordinates of pRFs
X = InSrf.Data(2,:);
Y = InSrf.Data(3,:);
% Unsmoothed data?
if isfield(InSrf, 'Raw_Data')
    rX = InSrf.Raw_Data(2,:);
    rY = InSrf.Raw_Data(3,:);
end

%% Rotate pRFs
nX = X .* cosd(theta) - Y .* sind(theta);
nY = Y .* cosd(theta) + X .* sind(theta);
% Unsmoothed data?
if isfield(InSrf, 'Raw_Data')
    nrX = rX .* cosd(theta) - rY .* sind(theta);
    nrY = rY .* cosd(theta) + rX .* sind(theta);
end

%% Add to output
Srf = InSrf;
Srf.Data(2,:) = nX;
Srf.Data(3,:) = nY;
% Unsmoothed data?
if isfield(InSrf, 'Raw_Data')
    Srf.Raw_Data(2,:) = nrX;
    Srf.Raw_Data(3,:) = nrY;
end
