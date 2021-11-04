function samsrf_clonecam(Z)
%
% samsrf_clonecam([Z])
%
% Clones the camera angle from the current figure in all other figures.
% You shouldn't have any figures open that aren't surface rendering
% figures. MatLab doesn't care if your plot is a 3D-mesh or a flat scatter
% plot. It will change the camera regardless...
%
% Also, this function isn't very clever. It only updates the camera view.
% The optional input argument Z defines the zoom to apply to each figure.
%
% 19/07/2020 - SamSrf 7 version (DSS)
% 22/10/2021 - Now also adjusts axis limits (DSS)
%

% Get current view
c = get(gca, 'CameraPosition');
a = axis;
% Loop through figures
h = findobj('Type','Figure');
for i = 1:length(h)
    figure(h(i));
    if ~strcmpi(get(gcf,'Name'), 'DisplayMaps')        
        set(gca, 'CameraPosition', c);
        axis(a);
        if nargin > 0
            zoom(Z);
        end
    end
end
