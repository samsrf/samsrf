function cSrf = samsrf_bandpass(Srf, Hz, Tr)
%
% cSrf = samsrf_bandpass(Srf, Hz, Tr)
%
% Uses fast fourier transform to band-pass filter the time series data 
%   in Srf.Data. Returns the "cleaned" Srf. Hz is a two component vector
%   defining the high-pass & low-pass cutoff, respectively, in Hz. 
%   In order to calculate this, the TR must be defined in seconds. 
%
% 20/04/2022 - SamSrf 8 version (DSS)
%

%% Band-pass filter
cSrf = Srf; % Cleaned Srf
Y = Srf.Data; % Time courses
Nt = size(Y,1); % Number of time points
Cps = round(Hz .* (Nt/Tr)) + 1; % Cycles per time series for each cutoff
Cps(Cps > Nt) = Nt; % Ensure it stays within bounds

disp('Band-pass filtering time series...'); 
cF = fft(Y); % Fast fourier transform 
cF([1:Cps(1) Cps(2):end],:) = 0; % Remove out-of-range cycles
tY = real(ifft(cF)); % Transform back
cSrf.Data = tY;
