function [Visual_Space, Timecourses, Xysb] = samsrf_backproj_prf(Response, pRF_Data, PrfFcn, Eccentricity, Threshold, Stimulus, ColourMap, NormaliseByDensity)
%
% [Visual_Space, Timecourses, Xysb] = 
%    samsrf_backproj_prf(Response, pRF_Data, PrcFcn, Eccentricity, [Threshold=[0 1 0 Inf], Stimulus=[], ColourMap='berlin', NormaliseByDensity=true])
%
% Projects the activity values in Response back into visual space using the
% pRF parameters in pRF_Data. Effectively, this is a sum of all pRF profiles 
% weighted by their response. pRF_Data is taken from a Srf.Data field with 
% the vertices you want and using the first rows produced by samsrf_fit_prf. 
% You will need to include as many rows as parameters that your pRF model needs
% plus the R^2 field in the first row for filtering.
%
% PrfFcn defines the pRF profile function but assuming 100 px apertures.
% Through this you can choose to back-project any kind of pRF model. 
%
% Eccentricity defines the maximum eccentricity of the mapping stimulus.
% This is only used for defining the extent of the graph for example when 
% mapping the full square area of the screen. As discussed below, you can 
% further restrict the eccentricity of included vertices in Threshold(4).
%
% The optional input Threshold defines the minimal R^2 of the pRFs to be 
% projected and how the responses are clipped. The third and fourth value 
% defined the inner and outer eccentricity to be used. This defaults to the 
% most inclusive criteria, i.e. all R^2 and all eccentricities.
%
% The optional input argument Stimulus defines the stimulus you may want to 
% superimpose on this image. This must be a contrast image (with values 0-1) 
% and it will be added to the output images at a 20% contrast. It must be a 
% matrix of 100 x 100 x NumberOfVolumes.
%
% The optional input ColourMap defines the colour map to be used. This defaults 
% to 'berlin'. If the string is prefixed by '-' it inverts the scale.
%
% The optional NormaliseByDensity toggles density normalisation on or off: 
% In this normalisation, the final back-projection image is divided by a 
% back-projection of the pRF density (all pRFs summed with equal weights). 
%
% Returns a movie of the time series plotted back into visual space as a
% matrix of 100 x 100 x NumberOfVolumes. The second output contains the
% same movie but in intensity image format. The third output contains in 
% rows the X and Y coordinates, the Sigma and the Response for each vertex.
%
% 29/06/2020 - SamSrf 7 version (DSS)
% 20/09/2021 - Changed default colour map to berlin (DSS)
%              Fixed bug with default inputs (DSS)
%              Default clipping level is now 1 (DSS)
%

if nargin < 5
    Threshold = [];
end
if nargin < 6
    Stimulus = [];
end  
if nargin < 7
    ColourMap = 'berlin';
end
if nargin < 8
    NormaliseByDensity = true;
end

if isempty(Threshold)
    Threshold = [0 1 0 Inf];
elseif length(Threshold) == 1
    Threshold = [Threshold 1 0 Inf];
elseif length(Threshold) == 2
    Threshold = [Threshold 0 Inf];
elseif length(Threshold) == 3
    Threshold = [Threshold Inf];
end

% Open working figure
wf = figure;

% Whole time course
Timecourses = zeros(200, 200, size(Response,1)); % Backprojected responses
Ds = zeros(200, 200, size(Response,1)); % pRF densities

% Filter vertices
gof = pRF_Data(1,:); % Goodness of fit
ecc = sqrt(pRF_Data(2,:).^2 + pRF_Data(3,:).^2); % Eccentricity
Response = Response(:, gof > Threshold(1) & ecc > Threshold(3) & ecc < Threshold(4));
pRF_Data = pRF_Data(:, gof > Threshold(1) & ecc > Threshold(3) & ecc < Threshold(4));
Response(isnan(Response)) = 0;

% Parameter output
Xysb = [pRF_Data(2:4,:); Response];

% Normalise pRF position
pRF_Data(2,:) = pRF_Data(2,:) / Eccentricity;
pRF_Data(3,:) = pRF_Data(3,:) / Eccentricity;
pRF_Data(4,:) = pRF_Data(4,:) / Eccentricity;

% Threshold by response intensity
Response = Response / max(abs(Response(:)));

% Loop through volumes
for t = 1:size(Response,1)
% parfor t = 1:size(Response,1)
    % Loop through vertices
    Curr = zeros(200,200); % Current response
    Dcur = zeros(200,200);
    for v = 1:size(pRF_Data,2)
        Rfp = PrfFcn(pRF_Data(2:end,v)');
        Curr = Curr + Rfp * Response(t,v);
        Dcur = Dcur + Rfp;
    end
    % Store current frame
    Timecourses(:,:,t) = Curr;
    Ds(:,:,t) = Dcur;
end

% Normalize backprojected response by pRF density?
if NormaliseByDensity
    Timecourses = Timecourses ./ Ds;
end

% Clip images
Timecourses = Timecourses / max(Timecourses(:));
Timecourses(abs(Timecourses(:)) >= Threshold(2)) = sign(Timecourses(abs(Timecourses(:)) >= Threshold(2))) * Threshold(2);
Timecourses = Timecourses / max(abs(Timecourses(:)));

% Pseudocolour code
if ColourMap(1) == '-'
    ColourMap = ColourMap(2:end);
    InvCmap = true;
else
    InvCmap = false;
end
cm = colormap([ColourMap '(200)']);
if InvCmap
    cm = flipud(cm);
end
fTc = round(Timecourses / (max(abs(Timecourses(:)))*0.5) * 50) + 100;
fTc(fTc==0) = 1;
fTc(isnan(fTc)) = 100; 
Visual_Space = zeros(200, 200, 3, size(Response,1));
for t = 1:size(Response,1)
    Visual_Space(:,:,:,t) = reshape(cm(fTc(:,:,t),:), [200 200 3]);
end

% Superimpose stimulus?
if ~isempty(Stimulus)
    for t = 1:size(Visual_Space,4) 
        for c = 1:3 
            Visual_Space(50:149,50:149,c,t) = Visual_Space(50:149,50:149,c,t) + Stimulus(:,:,t) * 0.2; 
        end
    end
end

% Remove outside region
Visual_Space = Visual_Space(50:149,50:149,:,:);
Timecourses = Timecourses(50:149,50:149,:);

% Close working figure
close(wf); 
