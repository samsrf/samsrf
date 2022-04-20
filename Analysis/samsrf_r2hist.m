function [p75, med, r2h] = samsrf_r2hist(Srf, Roi)
%
% [p75, med, r2h] = samsrf_r2hist(Srf, [Roi=''])
%
% Returns in p75 the 75th percentile and in med the median of the R^2 values 
% in Srf and plots the histogram of all reasonable vertices (R^2>0.01). 
% Can be restricted to region of interest Roi. It doesn't check if the
% values are in fact R^2 values but simply takes Srf.Data(1,:).
% The X-axis label states what the name of this field is so if this doesn't
% say R^2 or nR^2 you aren't plotting goodness-of-fit values...
%
% 20/04/2022 - SamSrf 8 version (DSS)
%

if nargin < 2
    Roi = '';
end

% Expand Srf if necessary
Srf = samsrf_expand_srf(Srf);

% Load ROI
if isempty(Roi)
    Vs = 1:size(Srf.Vertices,1);
    Roi = 'All vertices';
else
    Vs = samsrf_loadlabel(Roi);
end

% R^2 values
r2h = Srf.Data(1,Vs);

% Remove rubbish
r2h = r2h(r2h > 0.01 & ~isnan(r2h));

% 90th percentile & median
p75 = prctile(r2h, 75);
med = median(r2h);

% Plot histogram
hist(r2h, -.05:.05:1.1);
hold on
h1 = line([med med], ylim, 'color', 'r');
h2 = line([p75 p75], ylim, 'color', 'r', 'linestyle', '--');
xlim([-.05 1.1]);
legend([h1 h2], {'Median' '90^t^h %ile'});
title(Roi);
xlabel(Srf.Values{1});
ylabel('Number of data points');
