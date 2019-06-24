function Srf = samsrf_lopass(InSrf, StdDev)
%
% Srf = samsrf_lopass(InSrf, StdDev)
%
% Temporally smoothes the time courses in InSrf.Data using a Gaussian kernel
% with StdDev to remove high-frequency noise. StdDev is defined in terms of
% the number of volumes (TRs). Filtered data is returned in Srf.Data but 
% raw data are not kept so save this separately if needed.
%
% 30/10/2018 - Created this function (DSS)
% 01/11/2018 - Removed raw data because this breaks other functions (DSS)
%

Srf = InSrf;

% Create smoothing kernel
x = -round(StdDev*5):round(StdDev*5);
Kernel = exp(-(x.^2/(2*StdDev.^2)));

% Add information
if iscellstr(Srf.Functional)
    for iStr = 1:length(Srf.Functional)
        Srf.Functional{iStr} = [Srf.Functional{iStr} ' (Low pass filtered with SD=' num2str(StdDev) ')'];
    end
else
    Srf.Functional = [Srf.Functional ' (Low pass filtered with SD=' num2str(StdDev) ')'];
end

% Convolution 
Srf.Data = NaN(size(InSrf.Data));
for v = 1:size(InSrf.Data,2)
    Srf.Data(:,v) = conv(InSrf.Data(:,v), Kernel, 'same');
end