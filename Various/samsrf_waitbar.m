function h = samsrf_waitbar(ValueOrMessage, Handle)
%
% h = samsrf_waitbar(ValueOrMessage, Handle)
%
% Opens, updates, or closes the SamSrf waitbar & returns the handle to it.
%
%   ValueOrMessage: A string opens a new waitbar with this message
%                   A number between 0-1 it upates the waitbar (Handle must be defined then)
%                   '' closes the waitbar (Handle must be defined then)
%
%   Handle:         The figure handle to the waitbar. Must be defined when updating or closing
% 
% The function checks whether the waitbar is turned on in SamSrf_defaults and only displays it if this is the case.
%
% 27/05/2020 - Written (DSS)
%

%% Do we want a waitbar at all?
if nargout > 0
    % Opening new waitbar so check if desired
    if exist('SamSrf_defaults.mat', 'file')
        % Load default parameters
        load('SamSrf_defaults.mat');
        if ~exist('def_wb', 'var')
            % Wait bar undefined in defaults
            def_wb = true;
        end
        % Display waiting bars?
        wb = def_wb;
    else
        % If no default parameters exist
        wb = true;
    end
else
    % Updating waitbar so only do if desired
    if Handle == 0 
        wb = false;
    else
        wb = true;
    end
end

%% If waitbar is desired
if wb 
    if ischar(ValueOrMessage)
        % Message defined
        if nargin < 2
            % Message defined and no handle so start new waitbar
            warning off
            h = waitbar(0, {ValueOrMessage; pwd});
            warning on
        else
            if isempty(ValueOrMessage)
                % Empty message means close this waitbar
                if nargin < 2
                    error('Waitbar handle must be defined!');
                end
                close(Handle);
            else
                % Message defined but handle defined so update old waitbar
                warning off
                waitbar(0, Handle, {ValueOrMessage; pwd});
                warning on
            end
        end
    else
        % Value defined
        if nargin < 2
            error('Waitbar handle must be defined!');
        end
        waitbar(ValueOrMessage, Handle);
    end
else % No waitbar desired 
    % No handle
    h = 0;
end

