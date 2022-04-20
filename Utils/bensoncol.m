function cmap = bensoncol(res)
%
% cmap = bensoncol([res=360])
%
% Returns an approximation of the Benson polar colour map (from Bayesian retinotopy study) 
% with the resolution res, which should be a multiple of 8. By default res is 360. 
%
% It is only an approximation because their scheme is a weird rendition of
% HSV but with some extra colours in the ipsilateral hemifield.
%
% 20/04/2022 - SamSrf 8 version (DSS)
%

if nargin == 0
    res = 360;   % steps from one colour to the next
end
steps = round(res / 8);

% Colour peaks
R = [1 0 0];
Y = [1 1 0];
G = [0 1 0];
C = [0 1 1];
B = [0 0 1];
P = [1 0 1];

% Create the colour map
cmap = [];
cmap = [cmap; linspace(R(1), Y(1), steps)', linspace(R(2), Y(2), steps)', linspace(R(3), Y(3), steps)'];
cmap = [cmap; linspace(Y(1), G(1), steps)', linspace(Y(2), G(2), steps)', linspace(Y(3), G(3), steps)'];
cmap = [cmap; linspace(G(1), C(1), steps)', linspace(G(2), C(2), steps)', linspace(G(3), C(3), steps)'];
cmap = [cmap; linspace(C(1), B(1), steps)', linspace(C(2), B(2), steps)', linspace(C(3), B(3), steps)'];
cmap = [cmap; linspace(B(1), P(1), 2*steps)', linspace(B(2), P(2), 2*steps)', linspace(B(3), P(3), 2*steps)'];
cmap = [cmap; linspace(P(1), R(1), 2*steps)', linspace(P(2), R(2), 2*steps)', linspace(P(3), R(3), 2*steps)'];
% Reorganise
cmap = [cmap(res/2+1:end,:); cmap(1:res/2,:)];