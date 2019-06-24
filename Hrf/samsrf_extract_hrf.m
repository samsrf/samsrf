function samsrf_extract_hrf(fname, tr, roi, ntrials)
%
% samsrf_extract_hrf([fname, tr, roi, ntrials])
%
% Extracts the HRF for a data set. 
%
%   fname:      Surface data file name of the HRF measurement
%   tr:         TR of the scan (default = 2)
%   roi:        ROI label to restrict the measurement (default = '')
%   ntrials:    Number of trials in the scan (default = 10)
%
% Saves a matlab file containing the fitted points in Hrf, the function 
% parameters fP, and the raw data points in Raw_avg and Raw_err.
% This fit is done for all active vertices.
%
% 09/08/2018 - SamSrf 6 version (DSS)
%

if nargin < 1
    [fname, pname] = uigetfile('*.mat');
    fname = [pname fname];
    tr = 2;
    roi = '';
    ntrials = 10;
elseif nargin < 2
    fname = [fname '.mat'];
    tr = 2;
    roi = '';
    ntrials = 10;
elseif nargin < 3
    fname = [fname '.mat'];
    roi = '';
    ntrials = 10;
elseif nargin < 4
    fname = [fname '.mat'];
    ntrials = 10;
end

% Load data
load(fname);
Srf = samsrf_expand_srf(Srf);
img = Srf.Data;

% Volumes per trial
nvols = size(img,1);
vols_per_trial = nvols / ntrials;

% Remove outliers
ol = abs(nanmean(img,2)) > 1.5;
vm = mean(img(~ol,:));
for v = 1:length(ol)
    if ol(v)
        img(ol(v),:) = vm;
    end
end

% Average HRF 
Raw_hrf = NaN(vols_per_trial, size(img,2), ntrials);
t = 0;
for v = 1:vols_per_trial:nvols
    t = t + 1;
    Raw_hrf(:,:,t) = img(v:v+vols_per_trial-1,:);
end
Raw_avg = mean(Raw_hrf,3);
Raw_err = std(Raw_hrf,0,3) / sqrt(ntrials);

% Restrict to ROI?
if ~isempty(roi)
    mver = samsrf_loadlabel(roi);
else
    mver = 1:size(Raw_avg,2);
end

% Average across vertices
can = spm_hrf(tr);
vis_vol = 2:floor(vols_per_trial/2);  % Beginning of HRF
gv = mean(Raw_avg(vis_vol,mver))-mean(Raw_err(vis_vol,mver)) > 0; % Only visually responsive vertices
Raw_avg = nanmean(Raw_avg(:,mver(gv)),2);
Raw_err = nanmean(Raw_err(:,mver(gv)),2);

% Fit a function 
fP = fminsearch(@(P)hrf_errfun(tr, P, Raw_avg), [6 16 6 1]);
Hrf = spm_hrf(tr, [fP(1:2) 1 1 fP(3) 0 32]) * fP(4);
Hrf = max(can) / max(Hrf) * Hrf;
sHrf = spm_hrf(.1, [fP(1:2) 1 1 fP(3) 0 32]) * fP(4);
if mean(sHrf(1:floor(length(sHrf)/2))) < 0
    sHrf = -sHrf;
end
sHrf = max(Raw_avg) / max(sHrf) * sHrf;

% Plot the HRF
figure; hold off;
plot((0:vols_per_trial-1)*tr, Raw_avg, 'ko', 'linestyle', 'none', 'linewidth', 2, 'markerfacecolor', 'k');
hold on; 
t = 0:.1:32;
plot(t, sHrf, 'color', [.5 .5 .5], 'linewidth', 2);
yl = ylim;
line([0 tr], [1 1]*yl(1), 'color', 'r', 'linewidth', 3);
set(gca, 'fontsize', 15);
legend('Response', 'Fit', 'Stimulus');
xlabel('Time (s)');
ylabel('Response (z)');
xlim([-.1 vols_per_trial*tr+.1]);

% Save the HRF data
[p f e] = fileparts(fname);
if ~isempty(roi)
    roi = ['_' roi];
end
save(['hrf_' f roi], 'Hrf', 'sHrf', 'fP', 'Raw_avg', 'Raw_err', 'tr', '-v7.3');


% Error function to fit HRF
function err = hrf_errfun(TR,P,Y)
X = spm_hrf(TR, [P(1:2) 1 1 P(3) 0 32])' * P(4); 
X = X';
if length(X) > length(Y)
    X = X(1:length(Y));
elseif length(X) < length(Y)
    Y = Y(1:length(X));
end
err = sum((Y-X).^2);
