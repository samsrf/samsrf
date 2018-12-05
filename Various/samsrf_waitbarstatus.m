function wb = samsrf_waitbarstatus
%
% Returns the status of the SamSrf waitbar.

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
