function samsrf_simvsfit_ecc(Srf, Bins, Threshold)
%
% samsrf_simvsfit_ecc(Srf, Bins, [Threshold=-Inf])
%
% Plots a comparison of simulated ground truth pRFs & noisy model fits in Srf.
% Specifically, it works out the unique pRF sizes in the ground truths.
% Then in plots for each true sigma the modelled eccentricity against the
% ground truth eccentricity as in a classical binning analysis.
%
% Bins defines the eccentricity bins to calculate. That depends on how you
%   chose your ground truths so this is up to you to define. In order to
%   help you define that, if Bins is undefined the function reports the 
%   unique eccentricity values in the ground truth set instead.  
%
% Threshold defines the R^2 threshold of the model fits to include in 
%   the comparison. Defaults to -Inf (includes all).
%
% 20/04/2022 - SamSrf 8 version (DSS)
%

if nargin < 3
    Threshold = -Inf;
end

%% Unique ground truths
uE = unique(sqrt(Srf.Ground_Truth(1,:).^2 + Srf.Ground_Truth(2,:).^2));
uS = unique(Srf.Ground_Truth(3,:));
Cmap = lines(length(uS));
if nargin < 2
    samsrf_newline;
    samsrf_disp('Unique ground truth eccentricities');
    samsrf_disp(uE);
    return
end

%% Ground truth as new Srf
Trf= Srf;
Trf.Data = [Trf.Data(1,:); Trf.Ground_Truth];
Trf = rmfield(Trf, 'Ground_Truth');

%% Plot data
i = 1;
Lgnd = {};
Hs = [];
for s = uS
    Lgnd{i} = ['\sigma = ' num2str(s)];
    curM = Srf; curM.Data = curM.Data(:,Trf.Data(4,:) == s);
    curT = Trf; curT.Data = curT.Data(:,Trf.Data(4,:) == s);
    [~,h] = samsrf_plot(curM, 'Eccentricity', curT, 'Eccentricity', Bins, '', Threshold, 'Mean', Cmap(i,:));
    Hs = [Hs h];
    i = i + 1;
end
axis square
xlabel('Ground truth eccentricity');
ylabel('Modelled eccentricity');
title(['R^2 > ' num2str(Threshold)]);
line(xlim, xlim, 'color', 'k', 'linewidth', 2, 'linestyle', '--');
legend(Hs, Lgnd, 'location', 'northwest');