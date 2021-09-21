
function samsrf_timecourse_animation(Rfp, ApFrm, Hrf, Downsampling)
%
% samsrf_timecourse_animation(Rfp, ApFrm, [Hrf, Downsampling])
%
% Predicts the time course resulting from receptive field profile Rfp and stimulus mask movie ApFrm. 
% The time course is plotted as a movie and at the end it is convolved with the HRF. 
% Optional input Hrf is a vector with the HRF by volume for convolution, e.g. samsrf_hrf(1). Default = 1
% Optional input Downsampling defines the factor by which the time series is downsampled. Default = 1
%
% 02/06/2020 - SamSrf 7 version (DSS) 
% 22/09/2021 - Added option to downsample predictions (DSS)
%

if nargin < 3
    Hrf = 1;
end
if nargin < 4
    Downsampling = 1;
end

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

% Convolving with HRF & downsample if desired
cY = prf_convolve_hrf(Y, Hrf, Downsampling);
hold on
plot(1:Downsampling:length(Y), cY, 'r', 'linewidth',2);
hold off
legend('Neuronal only','Convolved (& downsampled?)');

