function [v,d,o] = samsrf_version
%
% [v,d,o] = samsrf_version
%
% Returns the current SamSrf version (Only exists since Version 4.0).
% The second output d contains a string with the date of the version.
%
% The third output o is the SamOaSrf version for Octave support but be
%  aware that Octave support may never be fully implemented.
%  Life is too short for this. It's more likely MapSrf for Python will be
%  completed - and that pigs fly and hell freezes over - before I manage 
%  to get all of this to work in Octave...
%

v = 7.21; % SamSrf version number
d = '09-07-2021'; % Release date 
o = 0.1; % SamOaSrf version number
