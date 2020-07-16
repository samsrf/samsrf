function [Srf, vx] = samsrf_denoisemap(Srf)
%
% [Srf, vx] = samsrf_denoisemap(Srf)
%
% Removes bad vertices in a pRF map in InSrf by setting their R^2 to zero.
% It removes any pRFs with Betas <= 0.01 or greater than 3. It also removes pRFs
% whose x0 and y0 position is exactly zero because these are typical artifacts.
% It finds the relevant data rows for x0, y0, and Beta but it asssumes that the 
% first row in Srf.Data is R^2 because that's the convention.
%
% IMPORTANT: this function only works on Srf.Data. If Srf.Raw_Data exists it
% will not be affected. You would first need to assign that to Srf.Data.
%
% Also note that this function assumes single subject data. So if you have a
% multi-subject Srf, you need to do this separately for each subject.
%
% Returns the denoised Srf. The second output vx contains the removed vertex indeces.
%
% 16/07/2020 - SamSrf 7 version (DSS)
%

% Relevant data rows
x0 = Srf.Data(find(strcmpi(Srf.Values, 'x0')),:);
y0 = Srf.Data(find(strcmpi(Srf.Values, 'y0')),:);
Betas = Srf.Data(find(strcmpi(Srf.Values, 'Beta')),:);

% Fine bad vertices
vx = Betas <= 0.01 | Betas > 3 | (x0 == 0 & y0 == 0);
% Remove bad vertices
Srf.Data(1,vx) = 0;
% Denoised vertex indeces
vx = find(vx);
