
function Y = prf_predict_timecourse(Rfp, ApFrm, Znorm, DispTc, Hrf)
%
% Y = prf_predict_timecourse(Rfp, ApFrm, Znorm, [DispTc, Hrf])
%
% Predicts the time course resulting from receptive field profile Rfp and stimulus mask movie ApFrm. 
% If Znorm is true the time series is z-normalised. When simulating data, this must be turned off!
%
%   If DispTc is true the time course is plotted as a movie. 
%   IMPORTANT: If Hrf is defined the time series is then convolved with the HRF but this is for display purposes only -
%              for actual model fitting you will need to convolve Y with the HRF separately!
%
% 20/08/2018 - SamSrf 6 version (date added - no changes from v6) (DSS)
% 25/11/2019 - Fixed bug with missing widerthantall function for displaying (DSS)
% 27/05/2020 - Added option to remove z-score normalisation (DSS) 
%

if nargin < 4
    DispTc = false;
end

% Output time course vector
Y = NaN(size(ApFrm,3),1); 

if DispTc    
    % Colour map
    cmap = [linspace(1,0,10) linspace(0,0,10) linspace(0,0,10) linspace(0,1,10) linspace(1,1,10) linspace(1,1,10); ...    % Red
            linspace(0,0,10) linspace(0,1,10) linspace(1,0,10) linspace(0,0,10) linspace(0,1,10) linspace(1,1,10); ...    % Green
            linspace(1,1,10) linspace(1,0,10) linspace(0,0,10) linspace(0,0,10) linspace(0,0,10) linspace(0,1,10)]';      % Blue

    % Plot apertures    
    figure; 
    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1)-pos(3)*0.4 pos(2) pos(3)*1.6 pos(4)]);
    subplot(2,3,1); 
    imshow(Rfp/2+.5);
    colormap(gca, cmap);
    set(gca,'fontsize',10);
    title('Receptive field profile');    
end

% Predict response at each time point
for i = 1:size(ApFrm,3)
    [y o m] = prf_predict_response(ApFrm(:,:,i),Rfp); 
    Y(i) = y; 
    
    if DispTc
        subplot(2,3,2); 
        set(gca,'fontsize',16);
        imshow(m/2+.5); 
        colormap(gca, cmap);
        title('Stimulus mask');

        subplot(2,3,3); 
        set(gca,'fontsize',16);
        imshow(o/2+.5); 
        colormap(gca, cmap);
        title('Overlap with pRF');

        subplot(2,1,2); 
        set(gca,'fontsize',10);
        hold on 
        plot(Y,'k','linewidth', 2); 
        xlim([0 size(ApFrm,3)]); 
        ylim([-100 200]);
        pause(0.01);
        set(gca, 'yticklabel', '');
        xlabel('Volume #');
        ylabel('Predicted response')
        title('Predicted time series');
    end
end

% Normalize time series?
if Znorm
    Y = zscore(Y);
end

if DispTc
    subplot(2,1,2); 
    set(gca,'fontsize',15);
    hold off
    plot(Y,'k','linewidth',2); 
    xlim([0 size(ApFrm,3)]); 
    ylim([-4 4]);
    set(gca, 'yticklabel', '');
    xlabel('Volume #');
    ylabel('Predicted response')
    title('Predicted time series');
    
    % Are we convolving with HRF?
    if nargin > 3
        cY = conv(Y, Hrf);
        cY = cY(1:length(Y)); % Truncate back to original length
        hold on
        plot(cY,'r','linewidth',2);
        hold off
        legend('Neuronal only','Convolved with HRF');
    end
end

