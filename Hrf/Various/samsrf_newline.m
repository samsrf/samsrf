function samsrf_newline(str)
% Adds a new line - if SamSrfAnalysis was run this is in the GUI
% otherwise it is in the Matlab command window.

global GuiInfo

% Command window or GUI?
if isempty(GuiInfo)
    disp(' ');
else
    GuiInfo.Value{end+1} = '';
end