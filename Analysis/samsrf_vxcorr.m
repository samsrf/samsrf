function [p, r, df] = samsrf_vxcorr(A, B, tcA, tcB, IsCirc)
%
% [p r df] = samsrf_vxcorr(A, B, tcA, tcB, [IsCirc=false])
%
% Conducts a Pearson correlation between vector A and B. Both A and B are 
% row vectors with each column containing the values from a set of vertices. 
%
% tcA and tcB are matrices with time courses for the same set of vertices 
% as are used in A and B. Each column is a vertex and rows are volumes. 
% This is used to calculate the inter-correlation between vertices in each 
% data vector. Summing the correlation coefficients across the correlation 
% matrix (treating r<0 as 0) gives the degrees of freedom for the test.
%
% IMPORTANT: This concept assumes that the inter-correlation of time
% courses in tcA and tcB and vertex data in A and B, respectively, is
% comparable. So you cannot use raw time courses (Srf.Y) to determine the
% correlatedness of smoothed map data!
%
% The optional boolean IsCirc toggles whether a circular correlation is
% calculated (default is false).
%
% Returns the p-value, correlation coefficient, and degrees of freedom.
%
% 28/11/2018 - Created this function (DSS)
%

if nargin < 5
    IsCirc = false;
end

%% Correlation matrices
% Correlation matrix for tcB
rA = corr(tcA);
% Remove redundant cells
for i = 1:size(rA,1)
    rA(1:i,i) = NaN;
end
% Negative correlations count as independent
rA(rA < 0) = 0;

% Correlation matrix for tcB
rB = corr(tcB);
% Remove redundant cells
for i = 1:size(rB,1)
    rB(1:i,1) = NaN;
end
% Negative correlations count as independent
rB(rB < 0) = 0;

%% Adjusted sample sizes
wA = nanmean(1-rA); % Average correlation per vertex
nA = nansum(wA) + 1; % Weighted sum of vertices
if ~isempty(tcB)
    wB = nanmean(1-rB); % Average correlation per vertex
    nB = nansum(wB) + 1; % Weighted sum of vertices
end
n = (nA + nB)/2; % Average weighted number
df = n - 2; % Degrees of freedom 

%% Inferential statistics
% Degrees of freedom cannot be below 1
if df < 1
    df = 1;
end
% Correlation
if IsCirc
    r = circcorr(A',B');
else
    r = corr(A',B');
end
t = sign(r) .* (abs(r) .* sqrt(df ./ (1-r.^2)));

%% Calculate significance
p = 2 * (1 - tcdf(abs(t), df));
