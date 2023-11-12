function Hrf = samsrf_plothrf(Srf, v)
%
% Hrf = samsrf_plothrf(Srf, v)
%
% If HRF parameters were estimated during pRF modelling of map Srf, this 
% function plots the HRF for vertex v. For comparison it also plots the
% canonical HRFs from SamSrf, SPM12, and mrVista. Returns the HRF.
%
% 13/11/2023 - Written (DSS)
%

% Which data rows are HRF parameters?
HrfParams = [find(strcmpi(Srf.Values, 'RLat')) ...
             find(strcmpi(Srf.Values, 'ULat')) ...
             find(strcmpi(Srf.Values, 'RDisp')) ...
             find(strcmpi(Srf.Values, 'UDisp')) ... 
             find(strcmpi(Srf.Values, 'R/U'))];

% Generate HRFs
Hrf = samsrf_doublegamma(.1, Srf.Data(HrfParams,v)); % Fit HRF
CanHrf = samsrf_hrf(.1); % SamSrf        
SpmHrf = samsrf_doublegamma(.1); % SPM12
VisHrf = samsrf_doublegamma(.1, [5.8034  12.7082  0.6733  0.5663  2.5227]); % mrVista

% Plot HRFs with normalised peaks        
figure; hold on
plot(0:.1:32, Hrf/max(Hrf)*max(CanHrf), 'k', 'linewidth', 2); % Fit HRF
plot(0:.1:32, CanHrf, 'r', 'linewidth', 2); % SamSrf 
plot(0:.1:32, SpmHrf/max(SpmHrf)*max(CanHrf), 'b--', 'linewidth', 2); % SPM12 
plot(0:.1:32, VisHrf/max(VisHrf)*max(CanHrf), 'm:', 'linewidth', 2); % mrVista
set(gca, 'fontsize', 12);
xlabel('Time (s)');
ylabel('Response');
title(['Vertex = ' num2str(v)]);
legend({'Fit HRF' 'SamSrf' 'SPM12' 'Vista'});
