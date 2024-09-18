function samsrf_disp(str)
% Displays the text in str - if SamSrfAnalysis was run this is in the GUI
% otherwise it is in the Matlab command window.

global GuiInfo

% Command window or GUI?
if isempty(GuiInfo)
    disp(str);
else
    if ischar(str)
        GuiInfo.Value{end+1} = str;
    elseif iscell(str)
        for i = 1:length(str)
            GuiInfo.Value{end+1} = str{i};
        end
    end
    if ~isstruct(GuiInfo)
        scroll(GuiInfo, 'bottom'); % Ensure latest output visible
    end
    pause(.1); % Otherwise updating can get stuck
end