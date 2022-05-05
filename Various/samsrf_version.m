function [v,d] = samsrf_version
%
% [v,d] = samsrf_version
%
% Returns the current SamSrf version (Only exists since Version 4.0).
% The second output d contains a string with the date of the version.
%
% This function used to include a 3rd output for the version of Octave support.
% Unsurprisingly, Sam is way too busy & life is too short for this & so the
% plan to improve Octave compatibility has been abandoned. I am not sad
% about this except that this means the excellent name SamOaSrf will never
% be used now. Fingers crossed that MapSrf will instead become a reality...
%

v = 8.1; % SamSrf version number (3rd digit is minor change)
d = '06-05-2022'; % Release date 
