function cmap = fireice(res)
%
% cmap = fscol([res=360])
%
% Returns a fire & ice colour map with the resolution res, which must be divisible by 4 (default = 360).
%
% 25/08/2021 - Written (DSS)
%

if nargin == 0
    res = 360;  
end
s = res / 4; % Steps from one colour to the next

% Combine ice with fire colour map
cmap = [linspace(.25,0,s)' linspace(1,0,s)' linspace(.75,1,s)'; ...
        zeros(s,2) linspace(1,0,s)'; ...
        fire(2*s)]; 
    