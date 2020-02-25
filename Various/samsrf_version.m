function [v,d,o] = samsrf_version
%
% [v,d,o] = samsrf_version
%
% Returns the current SamSrf version (Only exists since Version 4.0).
% The second output d contains a string with the date of the version.
% The third output o is the SamOaSrf version for Octave support.
%

v = 6.231; % SamSrf version number
d = '25-02-2020'; % Release date 
o = 0.1; % SamOaSrf version number