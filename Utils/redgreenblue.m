function cmap = redgreenblue(res)
%
% cmap = redgreenblue([res=360])
%
% Returns a colour map going from hot colours (yellow-orange-red) to green 
% and from there to cold colours (blue-cyan). The input res defines the 
% number of steps, which must be a multiple of 4. By default res is 256. 

if nargin == 0
    res = 256;   % steps from one colour to the next
end
steps = res / 4;

% Colour peaks
Y = [1 1 0];
R = [1 0 0];
G = [0 1 0];
B = [0 0 1];
C = [0 1 1];

% Create the colour map
cmap = [];
cmap = [cmap; linspace(Y(1), R(1), steps)', linspace(Y(2), R(2), steps)', linspace(Y(3), R(3), steps)'];
cmap = [cmap; linspace(R(1), G(1), steps)', linspace(R(2), G(2), steps)', linspace(R(3), G(3), steps)'];
cmap = [cmap; linspace(G(1), B(1), steps)', linspace(G(2), B(2), steps)', linspace(G(3), B(3), steps)'];
cmap = [cmap; linspace(B(1), C(1), steps)', linspace(B(2), C(2), steps)', linspace(B(3), C(3), steps)'];

% Upside down & too lazy to fix it
cmap = flipud(cmap); 