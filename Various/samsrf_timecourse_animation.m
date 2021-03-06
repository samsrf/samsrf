
function samsrf_timecourse_animation(Rfp, ApFrm, Hrf)
%
% samsrf_timecourse_animation(Rfp, ApFrm, [Hrf])
%
% Predicts the time course resulting from receptive field profile Rfp and stimulus mask movie ApFrm. 
% The time course is plotted as a movie and at the end it is convolved with the HRF. 
% Hrf is a vector with the HRF by volume, e.g. samsrf_hrf(1)
%
% 02/06/2020 - SamSrf 7 version (DSS) 
%

% Output time course vector
Y = NaN(size(ApFrm,3),1); 

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

% Predict response at each time point
Y = prf_predict_timecourse(Rfp, ApFrm);
for i = 1:size(ApFrm,3)
    [y, o, m] = prf_predict_response(ApFrm(:,:,i),Rfp); 
    
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
    plot(Y(1:i),'k','linewidth', 2); 
    xlim([0 size(ApFrm,3)]); 
    pause(0.01);
    set(gca, 'yticklabel', '');
    xlabel('Volume #');
    ylabel('Predicted response')
    title('Predicted time series');
end

% Display final time course?
subplot(2,1,2); 
set(gca,'fontsize',15);
hold off
plot(Y,'k','linewidth',2); 
xlim([0 size(ApFrm,3)]); 
set(gca, 'yticklabel', '');
xlabel('Volume #');
ylabel('Predicted response')
title('Predicted time series');

% Are we convolving with HRF?
if nargin >= 3
    cY = prf_convolve_hrf(Y, Hrf);
    hold on
    plot(cY,'r','linewidth',2);
    hold off
    legend('Neuronal only','Convolved with HRF');
end

