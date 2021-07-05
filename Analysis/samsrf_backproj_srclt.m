function [Backprojections, X, Y, Weights, Numbers, Good_Vertices, Used_pRFs, SrclID, VtxInSrcl] = samsrf_backproj_srclt(Response, pRF_Data, Eccentricity, Threshold, Resolution, Mode, SkipVtxInSrcl)
%
% [Backprojections, X, Y, Weights, Numbers, Good_Vertices, Used_pRFs, SrclID, VtxInSrcl] = samsrf_backproj_srclt(Response, pRF_Data, Eccentricity, Threshold=[0, 0, Inf, 0, 360], Resolution=[0.1 1], Mode='Mean', SkipVtxInSrcl=true)
%
% Projects the activity values in Response back into visual space using the
%  pRF parameters in pRF_Data (a Srf.Data field with the vertices you want
%  and assuming the first four values produced by samsrf_fit_prf).
%  It then uses a searchlight procedure to calculate a summary statistic
%  within a searchlight (circular or square) at each visual field position.
%
% If the rows in Response arise from different experiments or GLMs, you can 
%  define different summary statistics for each (see below).
%
% Eccentricity defines the maximum eccentricity of the mapping stimulus.
%  This is used for defining the extent of the graph for example when mapping
%  the full square area of the screen and also for normalizing summary weights
%  in the nearest neighbour mode (see below). As discussed below, you can
%  further restrict the eccentricity range of included *data* in Threshold(3).
%
% The 1st input in Threshold defines the minimal R^2 of the pRFs to be projected.
%  The 2nd and 3rd input define the inner and outer eccentricity to be used.
%  The 4th and 5th input define the lower and upper polar angle to be used.
%  Defaults to the most inclusive criteria, i.e. all R^2, all eccentricities, 
%  and all polar angles.
%
% Resolution defines the granularity of the searchlight grid, i.e. the
%  inter-grid-point distance (1st value), & the size of the searchlight(2nd value).
%  If the 2nd value is
%     - positive, a circular searchlight is used and the 2nd value
%        represents the circle's radius.
%     - is negative, the searchlight instead determines the pRFs which are
%        the nearest N Euclidean neighbours to the searchlight centre. The
%        2nd value therefore represents N. (A separate function trimmed to only 
%        this functionality may be added in the future as this produces the 
%        nicest backprojections)
%     - is a complex number (i.e., 0.5i), a square searchlight is used and
%        the absolute value (i.e., the modulus) of the 2nd value (i.e., 0.5)
%        represents half of the square's side length.
%
%  Note suboptimal choices for the granularity of the searchlight grid and the
%   searchlight size can result in missing pRFs (see below for Used_pRFs).
%   Moreover, since the size of the searchlight grid is determined by the
%   maximum eccentrity of the mapping stimulus, the maximum eccentricity should
%   be divisible by the chosen granularity value. With suboptimal granularities,
%   the serachlight grid might not fully cover the eccentricity range.
%
% Mode defines the summary statistic for the searchlight:
%  Central tendency estimates
%   'Mean':         Arithmetic mean (default)
%   'Wmean':        Arithmetic mean weighted by distance (smoother)
%   'Median':       Median
%   'Mode':         Mode
%   'Maximum':      Maximum
%   'Minimum':      Minimum
%   'GeoMean':      Geometric mean
%   'Sum':          Sum
%   't-test':       t-test vs 0
%  Dispersion estimates
%   'StdDev'        Standard deviation
%   'MeanAbsDev'    Mean absolute deviation
%   'MedAbsDev'     Median absolute deviation
%   'CircMean'      Mean direction for circular data 
%      (Note that previous versions used an external toolbox for this. 
%       For consistency with other functions in SamSrf, the internal function 
%       circmean is used instead now (should give identical results). 
%       Moreover, we don't support the circular median any longer because 
%       calculating this is non-trivial & the algorithm for that should be 
%       chosen deliberately & declared explicitly in your Methods section)
%
%   You can apply a different mode for each row of Response, using a cell array 
%   as input, e.g. {'Mean' 'Median' 'Maximum'} 
%
% SkipVtxInSrcl toggles whether indidces flagging the vertices falling within
%  a given searchlight shall be produced. Defaults to false.
%
%
% *** OUTPUTS ***
%
% Returns a matrix in intensity image format where each level along the 3rd
%  dimension) is the backprojection for one row of Response. This intensity
%  image can be plotted as a contour graph or a 3D surface.
%
% X and Y contain the coordinate matrix (for contour or surf plot etc).
%
% Weights contains the sum of distances of pRFs in the searchlight to its centre.
%  The 1st level along the 3rd dimension contains a matrix filtered in a manner
%  identical to the chosen summary statistic (e.g., if the t-statistic is Inf
%  in a given searchlight, the corresponding weight will be set to 0). The
%  second level along the 3rd dimension contains an unfiltered matrix (i.e.,
%  if the t-statistic is Inf in a given searchlight, the corresponding
%  weight will be positive and non-zero).
%
% Numbers contains the number of pRFs in a given searchlight (could be used
%  as weights). As with Weights, the 1st level along the 3rd dimension contains
%  the filtered matrix and the 2nd level the unfiltered one. When using the
%  nearest neighbour mode, this has a 3rd level with a matrix containing the
%  maximum distance from the searchlight centre. This matrix is filtered in
%  the same fashion as the chosen summary statistic. An unfiltered matrix is
%  currently not available.
%
% Good_Vertices contains an index of good vertices allowing you to filter
%  out pRFs with a poor R^2 or outside the eccentricity range, and after
%  removing NaN values in pRF_Data and Response. This is helpful for subsequent
%  analyses and if you would like to apply the filtering here to another
%  experimental condition or when projecting data on the cortical surface.
%
% Used_pRFs counts for each pRF in pRF_Data (i.e., each column) how often it
%  was included in a searchlight. Similar to Wt and Nr, the first row contains
%  a filtered and the second row an unfiltered count.
% 
% SrclID contains an identifier for each searchlight.
% 
% VtxInSrcl contains indices flagging vertices falling in a specific
%  searchlight. The row number corresponds to the identifier in SrclID.
%
%
% 19/07/2020 - SamSrf 7 version (DSS & SuSt)
% 03/08/2020 - Removed dependency on external circular statistics toolbox (DSS)
% 05/07/2021 - Fixed bug when Response has only one row (DSS)
%

%% Default inputs
if nargin == 3
    Threshold = [];
end

if nargin <= 4
    Resolution = [0.1 1];
end

if nargin <= 5
    Mode = 'Mean';
end

if nargin <= 6
    SkipVtxInSrcl = true;
end

if ~iscell(Mode)
    Mode = repmat({Mode}, 1, size(Response,1));
end

InputThreshold = Threshold;
Threshold = [0 0 Inf 0 360];
Threshold(1:length(InputThreshold)) = InputThreshold;

% pRF mapping data
gof = pRF_Data(1,:); % Goodness of fit
ecc = sqrt(pRF_Data(2,:).^2 + pRF_Data(3,:).^2); % Eccentricity
sigma = pRF_Data(4,:); % pRF size

pol_theta = cart2pol(pRF_Data(2,:),pRF_Data(3,:));
pol_theta = rad2deg(pol_theta);
pol_theta(pol_theta < 0) = pol_theta(pol_theta < 0)+360;

%% Whole time course
% Grid points (i.e. searchlight midpoints)
[X,Y] = meshgrid(-Eccentricity:Resolution(1):Eccentricity, -Eccentricity:Resolution(1):Eccentricity);
Y = flipud(Y);

% If square searchlight selected
if ~isreal(Resolution(2))
    % Generate edges of searchlight square
    XEdgesLo = X-abs(Resolution(2));
    XEdgesUp = X+abs(Resolution(2));
    YEdgesLo = Y-abs(Resolution(2));
    YEdgesUp = Y+abs(Resolution(2));
end

Backprojections = zeros(size(X,1), size(X,2), size(Response,1));
Weights = zeros(size(X,1), size(X,2),2); % Weights based on number & distance to searchlight centre
if Resolution(2) < 0
    % When using nearest neighbour mode
    Nr_Dim = 3;
else
    Nr_Dim = 2;
end

Numbers = zeros(size(X,1), size(X,2), Nr_Dim); % Number of pRFs within the searchlight
SrclID = Numbers(:, :, 1);

%% Filter vertices
nanprf = sum(isnan(pRF_Data),1) > 0; % Shouldn't happen but just in case...
nanresp = sum(isnan(Response),1) > 0; % Determine NaNs in response
% Retrieve only good vertices
Good_Vertices = gof > Threshold(1) & ecc > Threshold(2) & ecc < Threshold(3) & sigma > 0 & nanresp == 0 & nanprf == 0 & ...
    (pol_theta >= Threshold(4) & pol_theta <= Threshold(5));
Cleaned_Response = Response(:, Good_Vertices);
Cleaned_pRF_Data = pRF_Data(:, Good_Vertices);
Used_pRFs = zeros(2,size(Cleaned_Response,2));

%% Vertices in searchlights
if SkipVtxInSrcl
    VtxInSrcl = [];
else
    % Matrix with searchlight-vertex index
    VtxInSrcl = nan(sum(sum(sqrt(X.^2 + Y.^2) <= Eccentricity + Resolution(1))), size(Cleaned_Response,2));
end

%% Feasibility check-ups
IdxGeomean = strcmpi(Mode, 'Geomean');
if sum(IdxGeomean) > 0 && sum(sum(Cleaned_Response(IdxGeomean,:) <= 0)) > 0
    error('Geometric mean won''t be defined for non-positive numbers.')
end

%% Backprojection with searchlight
disp('Backprojecting data into visual space...');
% Id counter
c_id = 1;

%% Loop through X-coordinates
for x = 1:size(X,2)
    %% Loop through Y-coordinates
    for y = 1:size(Y,1)
        % Only searchlight within eccentricity range (expansion by 1 inter-grid-point distance to omit missing vertices because of edge
        % effects (i.e. interaction of a square-shaped grid with a circular region of interest)
        if sqrt(X(y,x)^2 + Y(y,x)^2) <= Eccentricity + Resolution(1)
            % Euclidean distances from current searchlight centre
            ed = sqrt((Cleaned_pRF_Data(2,:)-X(y,x)).^2  + (Cleaned_pRF_Data(3,:)-Y(y,x)).^2);
            % If searchlight size is valid
            if abs(Resolution(2)) ~= 0
                % If circular searchlight selected
                if isreal(Resolution(2))
                    % If nearest neighbour mode selected
                    if Resolution(2) < 0
                        % Find nearest fixed number of pRFs from searchlight centre
                        [sed,sx] = sort(ed); % Sort by distance to centre
                        svx = find(sed) <= -Resolution(2); % Find nearest pRFs
                        
                        % Check whether enough pRFs have been identified
                        if sum(svx) < -Resolution(2)
                            error('Could not find enough pRFs. Adjust Resolution(2) or skip analysis.') % Really bad visual fied coverage
                        end
                        
                        vx = sx(svx); % Sequential indices for nearest neighbours
                        vx = ismember(ed, ed(vx)); % Logical indices for nearest neighbours
                        % Inverse distances from searchlight centre
                        dc = 1 - ed(vx)/(2*Eccentricity); % Normalised to 2*eccentricity (i.e. mapping area)
                        
                    elseif Resolution(2) > 0 % If standard mode selected
                        % Find all pRFs within searchlight circle
                        vx = ed < Resolution(2);
                        % Inverse distances from searchlight centre
                        dc = 1 - ed(vx)/Resolution(2); % Normalised to searchlight radius
                    end
                    
                else % If square searchlight selected
                    % Find all pRFs within searchlight square
                    vx = (Cleaned_pRF_Data(2,:) >= XEdgesLo(y,x) & Cleaned_pRF_Data(2,:) < XEdgesUp(y,x)) & ...
                        (Cleaned_pRF_Data(3,:) >= YEdgesLo(y,x) & Cleaned_pRF_Data(3,:) < YEdgesUp(y,x));
                    % Inverse distances from searchlight centre
                    dc = 1 - ed(vx)/sqrt(abs(Resolution(2))^2+abs(Resolution(2))^2); % Normalised to half of square diagonal (aka square radius)
                end
            else % If searchlight size invalid
                error('Resolution(2) must not be zero!');
            end
            
            %% Loop through volumes
            for v = 1:size(Cleaned_Response,1)
                
                CurrMode = Mode(1, v);
                
                %% Central tendency estimates
                if strcmpi(CurrMode, 'Mean')
                    % Arithmetic mean of values
                    curstat = mean(Cleaned_Response(v,vx));
                elseif strcmpi(CurrMode, 'Wmean')
                    % Distance weighted arithmetic mean of values
                    curstat = sum(Cleaned_Response(v,vx).*dc) / sum(dc);
                elseif strcmpi(CurrMode, 'Median')
                    % Median of values
                    curstat = median(Cleaned_Response(v,vx));
                elseif strcmpi(CurrMode, 'Mode')
                    % Median of values
                    curstat = mode(Cleaned_Response(v,vx));
                elseif strcmpi(CurrMode, 'Maximum')
                    % Maximum of values
                    curstat = max(Cleaned_Response(v,vx));
                elseif strcmpi(CurrMode, 'Minimum')
                    % Minimum of values
                    curstat = min(Cleaned_Response(v,vx));
                elseif strcmpi(CurrMode, 'GeoMean')
                    % Geometric mean of values
                    curstat = geomean(Cleaned_Response(v,vx));
                elseif strcmpi(CurrMode, 'Sum')
                    % Sum of values
                    curstat = sum(Cleaned_Response(v,vx));
                elseif strcmpi(CurrMode, 't-test')
                    % Calculate t-test vs zero
                    [~,~,~,TestStats] = ttest(Cleaned_Response(v,vx), 0);
                    curstat = TestStats.tstat;
                elseif strcmpi(CurrMode, 'CircMean')
                    % Circular mean
                    curstat = circmean(Cleaned_Response(v,vx));
                    
                 %% Dispersion estimates
                elseif strcmpi(CurrMode, 'StdDev')
                    % Standard deviation
                    curstat = std(Cleaned_Response(v,vx));
                elseif strcmpi(CurrMode, 'MeanAbsDev')
                    % Mean absolute deviation
                    curstat = mad(Cleaned_Response(v,vx), 0);
                elseif strcmpi(CurrMode, 'MedAbsDev')
                    % Median absolute deviation
                    curstat = mad(Cleaned_Response(v,vx), 1);
                
                 %% Unknown mode defined
                else
                    error(['Unknown summary statistic ' CurrMode ' specified.']);
                end
                
                %% Store current stat in pixel
                if ~isempty(curstat) && ~isinf(curstat) && ~isnan(curstat)
                    Backprojections(y,x,v) = curstat;
                end
                
                %% Store distances & numbers (only need to do once)
                if v == 1
                    % Calculated across all pRFs inside the searchlight
                    if ~isempty(curstat) && ~isinf(curstat) && ~isnan(curstat)
                        Weights(y,x,1) = nansum(dc); % Filtered weights based on inverse distances to searchlight centre
                        Numbers(y,x,1) = nansum(vx); % Filtered count of pRFs in searchlight
                        Used_pRFs(1,vx) = Used_pRFs(1,vx) + 1; % Filtered count of how often pRFs were used
                        
                        if SkipVtxInSrcl == false
                            VtxInSrcl(c_id,:) = vx; % Store searchlight-vertex index
                        end
                        
                        % When using nearest neighbour mode
                        if Resolution(2) < 0
                            % Filtered max distance of vertices from searchlight center
                            Numbers(y,x,3) = nanmax(ed(vx));
                        end
                    else
                        if ~SkipVtxInSrcl
                            VtxInSrcl(c_id,:) = zeros(1, size(VtxInSrcl,2)); % Zero out searchlight-vertex index
                        end
                    end
                    
                    SrclID(y,x)       = c_id; % Store srcl id
                    c_id              = c_id +1; % Update id counter
                    
                    Weights(y,x,2) = nansum(dc); % Unfiltered weights based on inverse distances to searchlight centre
                    Numbers(y,x,2) = nansum(vx); % Unfiltered count of pRFs in searchlight
                    Used_pRFs(2,vx) = Used_pRFs(2,vx) + 1; % Unfiltered count of how often pRFs were used
                    
                    if sum(isnan(dc)) > 0 || sum(isnan(vx)) > 0
                        warning('NaN detected. This should not happen. Please request bugfix.')
                    end
                end
            end
        end
    end
end
