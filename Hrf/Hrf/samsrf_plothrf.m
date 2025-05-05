function [Hrf, Params] = samsrf_plothrf(Srf, vx, R2Thr)
%
% [Hrf, Params] = samsrf_plothrf(Srf, [vx=[], R2Thr=0])
%
% If HRF parameters were estimated during pRF modelling of map Srf, this 
% function plots the HRF for vertices vx. The 3rd optional input R2Thr
% defines the R^2 threshold for vertices to include.
%
% For comparison it also plots canonical HRFs by de Haas (previous SamSrf default), 
% from SPM12, and mrVista. If vx is empty (default) all vertices are included.
%
% If an output argument is defined, the function returns the HRF but does 
% not plot anything. Useful for analysis. The second output contains the
% actual HRF parameters of the fit.
%
% 13/11/2023 - Written (DSS)
% 12/10/2024 - Can now analyse list of vertices (DSS)
%              Now has option not to plot figure (DSS) 
% 20/10/2024 - Added R^2 threshold option (DSS)
% 21/10/2024 - Outputs HRF parameters as well now (DSS)
% 02/11/2024 - HRF output is now determined by average parameters (DSS)
%

% Which data rows are HRF parameters?
HrfParams = [find(strcmpi(Srf.Values, 'RLat')) ...
             find(strcmpi(Srf.Values, 'ULat')) ...
             find(strcmpi(Srf.Values, 'RDisp')) ...
             find(strcmpi(Srf.Values, 'UDisp')) ... 
             find(strcmpi(Srf.Values, 'R/U'))];
if isempty(HrfParams)
    samsrf_error('This data file does not contain HRF fits!');
end

if nargin < 2
    vx = [];
end
if nargin < 3
    R2Thr = 0;
end
if isempty(vx)
    vx = 1:size(Srf.Data,1);
end

% Expand Srf
Srf = samsrf_expand_srf(Srf);

% Generate HRFs
Hrf = [];
Params = [];
for x = 1:length(vx)
    v = vx(x); % Current vertex
    if Srf.Data(1,v) > R2Thr
        Params = [Params Srf.Data(HrfParams,v)]; % HRF parameters for this vertex
    end
end
% Average parameters across vertices
Params = nanmean(Params,2); 
% HRF with these parameters
Hrf = samsrf_doublegamma(.1, Params); 

% Plot HRFs?
if nargout == 0
    % Canonical HRFs for comparison
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
    if length(vx) == 1
        title(['Vertex = ' num2str(vx)]);
    else
        title('Multiple vertices');
    end
    legend({'Fit HRF' 'de Haas' 'SPM12' 'Vista'});
end