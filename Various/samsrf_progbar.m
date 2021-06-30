function samsrf_progbar(P)
%
% samsrf_progbar(P)
%
% Displays a progress bar in the command window with the proportion P being complete. 
%
% 30/06/2021 - Written (DSS)

% Empty progress bar
B = [' Progress: ' repmat('-',1,40)];

% Just started?
if P == 0
    fprintf(B);
else
    % Draw progress bar
    x = 11 + round(P*40);
    if x > 11
        B(x) = '*';
    end
    fprintf([repmat('\b',1,51) B]);
end

% Is complete?
if P == 1
    new_line;
end