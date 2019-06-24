function y = round_decs(x,d)
%y = round_decs(x,d)
%
% Rounds x to d decimals.
%

f = 10^d;
y = round(x*f)/f;