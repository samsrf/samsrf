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
% 17/07/2020 - SamSrf 7 version (DSS)
% 23/10/2020 - Added support for parallel computing (DSS)
% 12/07/2021 - Added stand-by message since parallel progress reports are a pain (DSS)
%

%% Default parameters
if nargin < 2
    Roi = '';
end

% Expand Srf if necessary
Srf = samsrf_expand_srf(Srf);

% Surface area data 
Ctx = Srf.Area;

% Load ROI file
if ~isempty(Roi)
    mver = samsrf_loadlabel(Roi);
else
    mver = (1:size(Srf.Vertices,1))';
end

% Determine visual area for each vertex
Ctx = Ctx(mver); % Limit to ROI
Vis = zeros(1,length(mver));
Cmf = zeros(1,length(mver));
disp('Calculating cortical magnification factors...');
disp(' Please stand by...');
t0 = tic;
parfor v = 1:length(mver)
    if Srf.Data(1,mver(v)) >= 0.01  % No point doing this for crappy vertices
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
end
disp([' CMF computation completed in ' num2str(toc(t0)/60) ' minutes.']);

% Store in structure
Srf.Data = [Srf.Data; zeros(2,size(Srf.Data,2))];
Srf.Data(end-1,mver) = Cmf;
Srf.Data(end,mver) = Vis;
Srf.Values{end+1} = 'Cmf';
Srf.Values{end+1} = 'Visual Area';
