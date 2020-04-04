function [v,d,o] = samsrf_version
%
% [v,d,o] = samsrf_version
%
% Returns the current SamSrf version (Only exists since Version 4.0).
% The second output d contains a string with the date of the version.
% The third output o is the SamOaSrf version for Octave support... 
%  but be aware that Octave support may never be fully implemented
%

v = 6.29993; % SamSrf version number
d = '04-04-2020'; % Release date 
o = 0.1; % SamOaSrf version number