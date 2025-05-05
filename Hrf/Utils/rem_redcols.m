function y = rem_redcols(x)
%
% y = rem_redcols(x)
%
% Removes redundant columns from a matrix of time courses
%
% 20/04/2022 - SamSrf 8 version (DSS)
%

y = x(:,1);

for i = 2:size(x,2)
    d = sum(abs(y-repmat(x(:,i), 1, size(y,2))));   % Residuals of difference
    r = d==0;  % Find redundant rows
    if sum(r) == 0
        y = [y x(:,i)];
    end
end

    
        