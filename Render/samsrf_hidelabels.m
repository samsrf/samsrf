function samsrf_hidelabels
% 
% samsrf_hidelabels
%
% Removes all ROI labels in a surface plot of displayed in current figure by samsrf_surf.
%
% 25/09/2026 - Written (DSS)

% Find handles
figure(gcf);
t = findall(gcf, 'Type', 'Text');
t = flipud(t(1:end-2)); % Invert order & remove two superfluous texts
q = findall(gcf, 'Type', 'Quiver');
q = flipud(q); % Invert order 
delete(t);
delete(q);

