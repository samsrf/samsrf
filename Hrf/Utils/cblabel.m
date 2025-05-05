function h = cblabel(LabelStr)
% Adds a label to a colour bar in the current figure axes and returns the handle to it.

h = colorbar;
set(get(h,'label'), 'string', LabelStr);
