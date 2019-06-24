function samsrf_lighting(tog)
%
% samsrf_lighting(tog)
%
% Toggles the camera lighting in the surface renderer 'on' or 'off'.
%
% 09/08/2018 - SamSrf 6 version (DSS) 
%

% Loop through figures
h = findobj('Type','Figure');
for i = 1:length(h)
    figure(h(i));
    if strcmpi(tog, 'on')
        lighting gouraud
        material dull
        camlight headlight
    elseif strcmpi(tog, 'off')
        lighting none
    end
end
