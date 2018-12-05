function Y = prf_predict_timecourse(Rfp, ApFrm, DispTc, Hrf)
%
% Y = prf_predict_timecourse(Rfp, ApFrm, [DispTc, Hrf])
%
% Predicts the time course resulting from receptive field profile Rfp 
%   and stimulus mask movie ApFrm. If DispTc is true the time course 
%   is plotted as a movie. If Hrf is defined the time series is then
%   convolved with the HRF (but in model fitting this is done separately).
%
% 20/08/2018 - SamSrf 6 version (data added - no changes from v6) (DSS)
%

if nargin < 3
    DispTc = false;
end

Y = NaN(size(ApFrm,3),1); 

if DispTc    
    % Colour map
    cmap = [linspace(1,0,10) linspace(0,0,10) linspace(0,0,10) linspace(0,1,10) linspace(1,1,10) linspace(1,1,10); ...    % Red
            linspace(0,0,10) linspace(0,1,10) linspace(1,0,10) linspace(0,0,10) linspace(0,1,10) linspace(1,1,10); ...    % Green
            linspace(1,1,10) linspace(1,0,10) linspace(0,0,10) linspace(0,0,10) linspace(0,0,10) linspace(0,1,10)]';      % Blue

    % Plot apertures    
    figure; widerthantall; 
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

% Normalize time series
Y = zscore(Y);

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