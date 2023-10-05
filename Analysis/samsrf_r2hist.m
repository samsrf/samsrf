function [p75, med, r2h] = samsrf_r2hist(Srf, Roi, DataRow)
%
% [p75, med, r2h] = samsrf_r2hist(Srf, [Roi=''], [DataRow=1])
%
% Returns in p75 the 75th percentile and in med the median of the values 
% in Srf and plots the histogram of all positive & non-NaN vertices. 
%
% Can be restricted to region of interest Roi. It doesn't check if the
% values are in fact R^2 values but simply takes Srf.Data(1,:) or whatever
% row you are indicating in input DataRow. 
%
% The X-axis label states what the name of this field is so if this doesn't
% say R^2 or nR^2 you aren't plotting goodness-of-fit values...
%
% 20/04/2022 - SamSrf 8 version (DSS)
% 05/10/2023 - Added option to plot any row in Srf.Data you want (DSS)
%              Adjusted maximum when plotting normalised R^2 (DSS)
%              Changed colour scheme of histogram (DSS)
%              Fixed error in legend (DSS)
%

if nargin < 2
    Roi = '';
end
if nargin < 3
    DataRow = 1;
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
r2h = Srf.Data(DataRow,Vs);

% Remove rubbish
r2h = r2h(r2h> 0 & ~isnan(r2h));

% 90th percentile & median
p75 = prctile(r2h, 75);
med = median(r2h);

% Plot histogram
if strcmpi(Srf.Values{DataRow}, 'nR^2')
    m = 1.6;
else
    m = 1.1;
end
hist(r2h, -.05:.05:m);
h = findobj(gca,'Type','patch');
h.FaceColor = [.8 .8 .8];
hold on
h1 = line([med med], ylim, 'color', 'r');
h2 = line([p75 p75], ylim, 'color', 'r', 'linestyle', '--');
xlim([-.05 m]);
legend([h1 h2], {'Median' '75^t^h %ile'});
title(Roi);
xlabel(Srf.Values{DataRow});
ylabel('Number of data points');
