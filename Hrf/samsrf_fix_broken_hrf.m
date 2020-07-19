function samsrf_fix_broken_hrf(fname, BadVols)
%
% samsrf_fix_broken_hrf([fname, BadVols])
%
% Loads a fitted HRF data set for review and allows dropping of bad data points. 
%
%   fname:      Surface data file name of the HRF measurement
%
% Saves the fixed data file under the same name suffixed by '_fixed'.
%
% 19/07/2020 - SamSrf 7 version (DSS)
%

if nargin < 1
    [fname, pname] = uigetfile('*.mat');
    fname = [pname fname];
end
load(fname);
vols_per_trial = length(Raw_avg);

if nargin < 2
    % Review old data 
    figure; hold off;
    plot(1:vols_per_trial, Raw_avg, 'ko', 'linestyle', 'none', 'linewidth', 2, 'markerfacecolor', 'k');
    hold on; 
    plot(max(Raw_avg) / max(Hrf) * Hrf, 'color', 'b', 'linewidth', 2);
    yl = ylim;
    line([0 1], [1 1]*yl(1), 'color', 'r', 'linewidth', 3);
    set(gca, 'fontsize', 15);
    legend('Response', 'Fit', 'Stimulus');
    xlabel('Volumes (#)');
    ylabel('Response (z)');
    xlim([-.1 vols_per_trial+.1]);
    % Define bad volumes
    BadVols = inputdlg('Drop these volumes: ','',1);
    BadVols = eval(['[' cell2mat(BadVols) ']']);
end

% Drop bad volumens
Raw_avg(BadVols) = NaN;
close(gcf);

% Fit a function 
can = samsrf_doublegamma(tr);
fP = fminsearch(@(P)hrf_errfun(tr, P, Raw_avg), [6 16 6 1]);
Hrf = samsrf_doublegamma(tr, [fP(1:2) 1 1 fP(3) 0 32]) * fP(4);
Hrf = max(can) / max(Hrf) * Hrf;
sHrf = samsrf_doublegamma(.1, [fP(1:2) 1 1 fP(3) 0 32]) * fP(4);
if mean(sHrf(1:floor(length(sHrf)/2))) < 0
    sHrf = -sHrf;
end
sHrf = max(Raw_avg) / max(sHrf) * sHrf;

% Plot the new HRF
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
save([f '_fixed'], 'Hrf', 'sHrf', 'fP', 'Raw_avg', 'Raw_err', 'tr', '-v7.3');


% Error function to fit HRF
function err = hrf_errfun(TR,P,Y)
X = samsrf_doublegamma(TR, [P(1:2) 1 1 P(3) 0 32])' * P(4); 
X = X';
X = X(~isnan(Y));
Y = Y(~isnan(Y));
if length(X) > length(Y)
    X = X(1:length(Y));
elseif length(X) < length(Y)
    Y = Y(1:length(X));
end
err = sum((Y-X).^2);
