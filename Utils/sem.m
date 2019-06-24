function e = sem(x)
%Calculates standard error of the mean for x.
%This is across rows in x. For columns transpose.

e = sqrt(nanvar(x)./sum(~isnan(x)));