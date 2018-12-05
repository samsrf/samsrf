function Vs = samsrf_find_redundancy(Data, v)
%
% Vs = samsrf_find_redundancy(Data, v)
%
% Returns the vertex indices in Data which contain the same data 
% (e.g. time courses) as vertex v. This help reduce redundancy 
% in the analysis.
%

% Calculates the residuals between column v in Data and every column in the data set.
D = sum(abs(Data-repmat(Data(:,v),1,size(Data,2))));
% Finds columns where the difference is zero (i.e. redundant time series).
Vs = find(D==0);