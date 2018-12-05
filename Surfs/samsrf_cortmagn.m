function Srf = samsrf_cortmagn(Srf, Roi)
%
% Srf = samsrf_cortmagn(Srf, Roi])
%
% Uses the pRF parameters in Srf to calculate the local cortical magnification 
% factor for each vertex. This is done by dividing the average distance from 
% the vertex to its neighbours on the cortical surface by the equivalent distance 
% in visual space. To ensure linearity the square root of these measures is used.
% (see Harvey & Dumoulin 2011 J Neurosci for details).
%
% The input Srf must contain surface area measures. If it doesn't, the function
% simply ends without calculating anything.
%
% The optional Roi defines a ROI label to restrict the analysis. 
%
% Adds a 2*n matrix with the CMF and visual area values to the Srf.
% (Surface area is already stored in the anatomical data)
%
% 09/08/2018 - SamSrf 6 version (DSS)
%

%% Default parameters
if nargin < 2
    Roi = '';
end
% Check if waitbar is to be used
wb = samsrf_waitbarstatus;

% Expand Srf if necessary
Srf = samsrf_expand_srf(Srf);

% Surface area data 
Ctx = Srf.Area;

% Load ROI file
if ~isempty(Roi)
    mver = samsrf_loadlabel(Roi);
else
    mver = 1:size(Srf.Vertices,1);
end

% Determine visual area for each vertex
Vis = zeros(size(Srf.Data,2), 1);
Cmf = zeros(size(Srf.Data,2), 1);
if wb h = waitbar(0, 'Calculating local CMF...', 'Units', 'pixels', 'Position', [100 100 360 70]); end
for v = mver'
    if Srf.Data(1,v) >= 0.01  % No point doing this for crappy vertices
        % Calculate visual area (density)
        Vis(v) = samsrf_vertexarea(v, Srf.Data(2:3,:)', Srf.Faces);
        % Fix artifacts
        if Vis(v) == 0
            Vis(v) = 0.0001;
            Ctx(v) = 0;
        end
        % Calculate cortical magnification factors
        Cmf(v) = sqrt(Ctx(v)) / sqrt(Vis(v));
    end
    if wb waitbar(v/length(mver), h); end
end
if wb close(h); end

% Store in structure
Srf.Data = [Srf.Data; Cmf'; Vis'];
Srf.Values{end+1} = 'Cmf';
Srf.Values{end+1} = 'Visual area';
