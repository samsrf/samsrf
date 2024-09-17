function samsrf_progbar(P)
%
% samsrf_progbar(P)
%
% Displays a progress bar in the command window with the proportion P being complete. 
%
% 16/09/2024 - Added GUI support for SamSrf X (DSS)
%

global GuiInfo wb

% Empty progress bar
B = [' Progress: ' repmat('-',1,40)];

% Just started?
if P == 0
    if isempty(GuiInfo)
        fprintf(B);
    else
        % GuiInfo.Value{end+1} = B;
        wb = waitbar(0,'SamSrf analysis...');
    end
else
    % Draw progress bar
    x = 11 + round(P*40);
    if x > 11
        B(x) = '*';
    end
    if isempty(GuiInfo)
        fprintf([repmat('\b',1,51) B]);
    else
        % GuiInfo.Value{end} = B;
        waitbar(P,wb);
    end
end

% Is complete?
if P == 1
    if ~isempty(GuiInfo)
        close(wb);
    end
    samsrf_newline;
end

