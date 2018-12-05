function [p, t, df] = samsrf_vxttest(A, B, tcA, tcB)
%
% [p t df] = samsrf_vxttest(A, B, tcA, [tcB=[]])
%
% Conducts a t-test of whether the mean of vector A is significantly 
% different from the mean of vector B. Both A and B are row vectors 
% with each column containing the values from a set of vertices. 
%
% tcA and tcB are matrices with time courses for the same set of vertices 
% as are used in A and B. Each column is a vertex and rows are volumes. 
% This is used to calculate the inter-correlation between vertices in each 
% data vector. Summing the correlation coefficients across the correlation 
% matrix (treaing r<0 as 0) gives the degrees of freedom for the t-test.
%
% IMPORTANT: This concept assumes that the inter-correlation of time
% courses in tcA and tcB and vertex data in A and B, respectively, is
% comparable. So you cannot use raw time courses (Srf.Y) to determine the
% correlatedness of smoothed map data!
%
% There are three options:
%   1. One-sample t-test
%       B must be a scalar. tcB must be empty.
%   2. Paired t-test
%       A and B must be vectors. tcB must be empty.
%   3. Two-sample t-test
%       A and B must be vectors. tcA and tcB must be matrices.
%
% Returns the p-value, t-statistic, and degrees of freedom of the test.
%
% 31/10/2018 - Created this function (DSS)
%              Fixed issue with too low degrees of freedom (DSS)
% 01/11/2018 - Sample size is now plus 1 for mathematical reasons (DSS)
%

if nargin < 4
    tcB = [];
end

%% What test is this?
if length(B) == 1
    % One-sample t-test
    Type = 0; 
    if ~isempty(tcB)
        error('TcB must be empty for one-sample t-test!');
    end   
else 
    if length(A) ~= size(tcA,2)
        error('TcA must have the same number of vertices as A!');
    end
    if isempty(tcB)
        % Paired t-test
        Type = 1;
    else
        % Two-sample t-test
        Type = 2;
        if length(B) ~= size(tcB,2)
            error('TcB must have the same number of vertices as B!');
        end
    end
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
if ~isempty(tcB)
    rB = corr(tcB);
    % Remove redundant cells
    for i = 1:size(rB,1)
        rB(1:i,1) = NaN;
    end
    % Negative correlations count as independent
    rB(rB < 0) = 0;
end

%% Adjusted sample sizes
wA = nanmean(1-rA); % Average correlation per vertex
nA = nansum(wA) + 1; % Weighted sum of vertices
if ~isempty(tcB)
    wB = nanmean(1-rB); % Average correlation per vertex
    nB = nansum(wB) + 1; % Weighted sum of vertices
end

%% Inferential statistics
% Which test?
switch Type
    case 0
        % One-sample t-test
        disp('One-sample t-test');
        [~,~,~,S] = ttest(A,B); % Run t-test
        t = S.tstat; % t-statistic
        df = nA-1; % Degrees of freedom
        
    case 1
        % Paired t-test
        disp('Paired t-test');
        [~,~,~,S] = ttest(A,B); % Run t-test
        t = S.tstat; % t-statistic
        df = nA-1; % Degrees of freedom
        
    case 2
        % Two-sample t-test
        disp('Two-sample t-test');
        [~,~,~,S] = ttest2(A,B); % Run t-test
        t = S.tstat; % t-statistic
        df = nA+nB-2; % Degrees of freedom
        
    otherwise
        error('Something went wrong - not sure what test this is...');
end
% Degrees of freedom cannot be below 1
if df < 1
    df = 1;
end

%% Calculate significance
p = 2 * (1 - tcdf(abs(t), df));
