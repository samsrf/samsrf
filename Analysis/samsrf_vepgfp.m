function Gfp = samsrf_vepgfp(Srf)
%
% Gfp = samsrf_vepgfp(Srf)
%
% Plots the global field power across time averaged for all sensors.
% Each curve depicts the power for one stimulus condition.
% This requires a Srf file with EEG/MEG data. 
%
% 08/09/2023 - Written (DSS)
% 19/09/2023 - Now returns the global field power curve for each row (DSS)
%

%% Check if M/EEG data
if ~strcmpi(Srf.Hemisphere, 'eeg')
    error('This Srf does not contain MEG/EEG data!');
end

%% Average for each time point 
uT = unique(Srf.TimePts); % Unique time points
Gfp = []; % Average global field power per stimulus

% Loop thru time points
for t = uT
    cur = Srf.Data(:, Srf.TimePts==t).^2; % Squared powers for this time point (conditions in rows, sensors in columns)
    Gfp = [Gfp mean(cur,2)]; % Add average to output
end

%% Plot data
plot(uT, Gfp');
hold on
plot(uT, mean(Gfp), 'k', 'linewidth', 2);
set(gcf, 'Units', 'normalized');
set(gcf, 'Position', [.1 .1 .8 .4]);
set(gca, 'fontsize', 12);
xlabel('Time');
ylabel('Mean squared power');
legend(Srf.Values, 'location', 'EastOutside');
