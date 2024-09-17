function [Srf, vx] = samsrf_denoisemap(Srf, BetaThr, RemoveOrigins, SigmaThr)
%
% [Srf, vx] = samsrf_denoisemap(Srf, [BetaThr=[.01 3], RemoveOrigins=true, SigmaThr=0])
%
% Removes bad vertices in a pRF map in InSrf by setting their R^2 to zero.
% It removes any pRFs with Betas <= BetaThr(1) or greater than BetaThr(2). 
% By default it also removes pRFs whose x0 and y0 position is exactly zero 
% because these are typical artifacts (only works if both x0 and y0 exist).
% Also removes pRFs whose size Sigma <= SigmaThrsh.
%
%   BetaThr defaults to [.01 3] but you can tweak it.
%   RemoveOrigins is a boolean that toggles whether x0=y0=0 is removed. 
%   SigmaThr defaults to 0.
%
% It works by finding the relevant data rows for x0, y0, and Beta but it 
% asssumes that the first row in Srf.Data is R^2 because that's the convention.
%
% IMPORTANT: this function only works on Srf.Data. If Srf.Raw_Data exists it
% will not be affected. You would first need to assign that to Srf.Data.
%
% Also note that this function assumes single subject data. So if you have a
% multi-subject Srf, you need to do this separately for each subject.
%
% Returns the denoised Srf. The second output vx contains the removed vertex indeces.
%
% 20/04/2022 - SamSrf 8 version (DSS)
%

if nargin < 2
    BetaThr = [.01 3];
end
if nargin < 3
    RemoveOrigins = true;
end
if length(BetaThr) == 1
    BetaThr = [BetaThr 3];
end
if nargin < 4
    SigmaThr = 0;
end

% Relevant data rows
x0 = Srf.Data(strcmpi(Srf.Values, 'x0'),:);
y0 = Srf.Data(strcmpi(Srf.Values, 'y0'),:);
Betas = Srf.Data(strcmpi(Srf.Values, 'Beta'),:);
Sigma = Srf.Data(strcmpi(Srf.Values, 'Sigma'),:);
samsrf_disp(['Removing Betas <= ' num2str(BetaThr(1)) ' or > ' num2str(BetaThr(2))]);
samsrf_disp(['Removing Sigmas <= ' num2str(SigmaThr(1))]);

% If any data don't exist
if isempty(x0)
    samsrf_disp('No x0 data in this map');
    RemoveOrigins = false;
end
if isempty(y0)
    samsrf_disp('No y0 data in this map');
    RemoveOrigins = false;
end
if isempty(Betas)
    samsrf_error('No Beta data in this map!');
end

% Fine bad vertices
if RemoveOrigins
    % Remove pRFs perfectly at origin
    samsrf_disp('Removing pRFs located perfectly at origin');
    ovx = (x0 == 0 & y0 == 0);
else
    ovx = false(1,size(Srf.Data,2));
end
vx = Betas <= BetaThr(1) | Betas > BetaThr(2) | ovx | Sigma <= SigmaThr;
% Remove bad vertices
Srf.Data(1,vx) = 0;
% Denoised vertex indeces
vx = find(vx);
