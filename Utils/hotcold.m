function cmap = hotcold(res)
%
% cmap = hotcold([res=256])
%
% Returns a hot-cold colour map with the resolution res, which must be a 
% multiple of 4. By default res is 256. 
%
% 20/04/2022 - SamSrf 8 version (DSS)
%

if nargin == 0
    res = 256;
end

steps = res/4;
cmap = [linspace(0,0,steps)' linspace(1,0,steps)' linspace(1,1,steps)'; ...
        linspace(0,0,steps)' linspace(0,0,steps)' linspace(1,0,steps)'; ...
        linspace(0,1,steps)' linspace(0,0,steps)' linspace(0,0,steps)'; ...
        linspace(1,1,steps)' linspace(0,1,steps)' linspace(0,0,steps)'];
