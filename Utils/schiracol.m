function cmap = schiracol(res)
%
% cmap = schiracol([res=360])
%
% Returns an approximation of the Schira polar colour map (from ACNS2018 poster) 
% with the resolution res. By default res is 360. This is a simple scheme which 
% is basically just HSV doubled.
%
% 20/04/2022 - SamSrf 8 version (DSS)
%

if nargin == 0
    res = 360;   % steps from one colour to the next
end
steps = round(res/3);

% Create the colour map
cmap = flipud([hsv(res/2); hsv(res/2)]);
cmap = [cmap(steps+1:end,:); cmap(1:steps,:)];