function R = circcorr(A,B)
% R = circcorr(A,B)
%
% Calculates the circular correlation between angles (in degrees) in A and B.
% Observations in A and B should be in rows, variables in columns.
% If only A is defined, it computes the correlation of A with A.

% If only one input is defined
if nargin == 1
    B = A;
end

% Convert degrees to radians
A = A / 180 * pi;
B = B / 180 * pi;

% Sample size
n = size(A,1);

% Initialise correlation matrix
R = NaN(size(A,2), size(B,2));

% Loop thru columns of A
for ac = 1:size(A,2)
    % Current data vector
    a = A(:,ac);
    % Loop thru columns of B
    for bc = 1:size(B,2)
        % Current data vector
        b = B(:,bc);
        % Calculate correlation
        r = sum(sin(a-circmeanrad(a)) .* sin(b-circmeanrad(b))) ...
          / sqrt(sum(sin(a-circmeanrad(a)).^2) .* sum(sin(b-circmeanrad(b)).^2));
        R(ac,bc) = r;
    end
end

%% Circular mean (input in radians)
function Alpha = circmeanrad(X)

n = length(X); % Number of observations
Alpha = atan2(sum(sin(X)) / n, sum(cos(X)) / n); % Circular mean
