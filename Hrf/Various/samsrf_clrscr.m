function samsrf_clrscr(str)
% Clears the screen - if SamSrfAnalysis was run this is in the GUI
% otherwise it is in the Matlab command window.

global GuiInfo

% Command window or GUI?
if isempty(GuiInfo)
    clc;
else
    GuiInfo.Value = {''};
end