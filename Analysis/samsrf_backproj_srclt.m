function [Bp, X, Y, Wt, Nr, Cleaned_pRF_Data, Cleaned_Response, Used_pRFs] = samsrf_backproj_srclt(Response, pRF_Data, Eccentricity, Threshold, Resolution, Mode)
%
% [Bp, X, Y, Wt, Nr, Cleaned_pRF_Data, Cleaned_Response, Used_pRFs] = samsrf_backproj_srclt(Response, pRF_Data, Eccentricity, [Threshold=[0, 0, Inf], Resolution=[0.1 1], Mode=['Mean'])
%
% Projects the activity values in Response back into visual space using the
%  pRF parameters in pRF_Data (a Srf.Data field with the vertices you want
%  and assuming the first four values produced by samsrf_fit_prf).
%  It then uses a searchlight procedure to calculate a summary statistic
%  within the search light radius at each visual field position.
%
% IMPORTANTLY, if the rows in Response arise from different experiments or
%  GLMs, these should be passed through the function iteratively. The reason
%  for this is that the function does not expect that the filter criteria
%  (for e.g., summary statistics) differ across rows.
%
% Eccentricity defines the maximum eccentricity of the mapping stimulus.
%  This is used for defining the extent of the graph for example when mapping
%  the full square area of the screen and also for normalizing summary weights 
%  in the nearest neighbour mode (see below). As discussed below, you can 
%  further restrict the eccentricity range of included *data* in Threshold(3).
%
% The 1st input in Threshold defines the minimal R^2 of the pRFs to be projected.
%  The 2nd and 3rd input define the inner and outer eccentricity to be used.
%  Defaults to the most inclusive criteria, i.e. all R^2 and all eccentricities.
%
% Resolution defines the granularity of the searchlight grid (1st value)
%  & the searchlight radius (2nd value). If the second value is negative
%  the searchlight instead determines the pRFs which are the nearest N
%  Euclidean neighbours to the searchlight centre (Note that this can
%  result in some pRFs not being used at all if the granularity of the
%  searchlight grid is not fine enough).
%
% Mode defines the summary statistic for the searchlight:
%  Central tendency estimates
%   'Mean':         Arithmetic mean (default)
%   'Median':       Median
%   'Mode':         Mode
%   'Maximum':      Maximum
%   'Minimum':      Minimum
%   'Geomean':      Geometric mean
%   'Sum':          Sum
%   't-test':       t-test vs 0
% Dispersion estimates
%   'StdDev'        Standard deviation
%   'MeanAbsDev'    Mean absolute deviation
%   'MedAbsDev'     Median absolutre deviation
%
% Returns a matrix in intensity image format where each level (that is,
%  along the 3rd dimension) is the backprojection for one row of Response.
%  This intensity image can be plotted as a contour graph or a 3D surface.
%
% X and Y contain the coordinate matrix (for contour or surf plot etc).
%
% Wt contains the sum of distances of pRFs in the searchlight to its centre.
%  The 1st level along the 3rd dimension contains a matrix filtered in a manner
%  identical to the chosen summary statistic (e.g., if the t-statistic is Inf
%  in a given searchlight, the corresponding weight will be set to 0). The
%  second level along the 3rd dimension contains an unfiltered matrix (i.e., 
%  if the t-statistic is Inf in a given searchlight, the corresponding 
%  weight will be positive and non-zero).
%
% Nr contains the number of pRFs in a given searchlight (could be used as
%  weights). As with Wt, the 1st level along the 3rd dimension contains the
%  filtered matrix and the 2nd level the unfiltered one. When using the
%  nearest neighbour mode, this has a 3rd level with a matrix containing the
%  maximum distance from the searchlight centre. This martix is filtered in
%  the same fashion as the chosen summary statistic. An unfilterd matrix is
%  currently not available.
%
% Cleaned_pRF_Data and Cleaned_Response contain the pRF_Data and Response inputs
%  after filtering pRFs with a poor R^2 or outside the eccentricity range, and
%  after removing rubbish NaN values.
%
% Used_pRFs count for each pRF in pRF_Data (i.e. each column) how often it
%  was included in a searchlight. This is more useful for nearest neighbour
%  mode (when Resolution(2)<0) than it is for the standard searchlight mode.
%  Similar to Wt and Nr, the first row contains a filtered and the
%  second row an unfiltered count. 

% 10/08/2018 - SamSrf 6 version (DSS & SuSt)
% 11/08/2018 - Fixed bug with method switch (DSS)
% 24/11/2018 - Added filtered weight matrices
%            - Removed NaN substitution outside loop as loop itself checks for this (SuSt)
% 26/11/2018 - Added warning message if NaN detected when summing up weights 
%              Simplified code for filtered weight matrices & storage of current summary statistic
%              Adapted calculation of used pRFs (filtered and unfiltered)
%              Maximum distance from searchlight centre now filtered like summary statistic
%              Fixed bug when calculating count weight matrix
%              Fixed bug when checking for standard mode (SuSt)
% 27/11/2018 - Nearest neighbour mode: distances are now normalized using 2*Eccentricity (SuSt)
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

if isempty(Threshold)
    Threshold = [0 0 Inf];
elseif length(Threshold) == 1
    Threshold = [Threshold 0 Inf];
elseif length(Threshold) == 2
    Threshold = [Threshold Inf];
end

% pRF map data
gof = pRF_Data(1,:); % Goodness of fit
ecc = sqrt(pRF_Data(2,:).^2 + pRF_Data(3,:).^2); % Eccentricity
sigma = pRF_Data(4,:); % pRF size

%% Whole time course
[X,Y] = meshgrid(-Eccentricity:Resolution(1):Eccentricity, -Eccentricity:Resolution(1):Eccentricity);
Y = flipud(Y);
Bp = zeros(size(X,1), size(X,2), size(Response,1));
Wt = zeros(size(X,1), size(X,2),2); % Weights based on number & distance to searchlight centre
if Resolution(2) < 0
    % When using nearest neighbour mode
    Nr_Dim = 3;
else
    Nr_Dim = 2;
end

Nr = zeros(size(X,1), size(X,2), Nr_Dim); % Number of pRFs within the searchlight

%% Filter vertices
nanprf = sum(isnan(pRF_Data)) > 0; % Shouldn't happen but just in case...
nanresp = sum(isnan(Response)) > 0; % Determine NaNs in response
% Retrieve only good vertices
gvx = gof > Threshold(1) & ecc > Threshold(2) & ecc < Threshold(3) & sigma > 0 & nanresp == 0 & nanprf == 0;
Cleaned_Response = Response(:, gvx);
Cleaned_pRF_Data = pRF_Data(:, gvx);
Used_pRFs = zeros(2,size(Cleaned_Response,2));

%% Feasibility check-ups
if strcmpi(Mode, 'Geomean') && sum(Cleaned_Response(:) <= 0) > 0
    error('Geometric mean won''t be defined for non-positive numbers.')
end

% Backprojection with searchlight
h = waitbar(0, 'Backprojecting volumes...');
% Loop through X-coordinates
for x = 1:size(X,2)
    % Loop through Y-coordinates
    for y = 1:size(Y,1)
        % Only searchlight within eccentricity range
        if sqrt(X(y,x)^2 + Y(y,x)^2) <= Eccentricity
            % Euclidean distances from current searchlight centre
            ed = sqrt((Cleaned_pRF_Data(2,:)-X(y,x)).^2  + (Cleaned_pRF_Data(3,:)-Y(y,x)).^2);
            % When using nearest neighbour mode
            if Resolution(2) < 0
                % Find nearest fixed number of pRFs from searchlight centre
                [sed,sx] = sort(ed); % Sort by distance to centre
                svx = find(sed) <= -Resolution(2); % Find nearest pRFs
                
                % Check whether enough pRFs have been identified
                if sum(svx) < -Resolution(2)
                    error('Could not find enough pRFs. Adjust Resolution(2) or skip analysis.') % Really bad visual fied coverage
                end
                
                vx = sx(svx); % Sequential indeces for nearest neighbours
                vx = ismember(ed, ed(vx)); % Logical indeces for nearest neighbours
                % Inverse distances from searchlight centre
                dc = 1 - ed(vx)/(2*Eccentricity); % Normalised to 2*eccentricity (i.e. mapping area)
            elseif Resolution(2) > 0
                % Find all pRFs within searchlight radius
                vx = ed < Resolution(2);
                % Inverse distances from searchlight centre
                dc = 1 - ed(vx)/Resolution(2); % Normalised to searchlight radius
            else
                error('Resolution(2) must not be zero!');
            end
            
            %% Loop through volumes
            for v = 1:size(Cleaned_Response,1)
                
                %% Central tendency estimates
                if strcmpi(Mode, 'Mean')
                    %% Arithmetic mean of values
                    curstat = mean(Cleaned_Response(v,vx));
                elseif strcmpi(Mode, 'Median')
                    %% Median of values
                    curstat = median(Cleaned_Response(v,vx));
                elseif strcmpi(Mode, 'Mode')
                    %% Median of values
                    curstat = mode(Cleaned_Response(v,vx));
                elseif strcmpi(Mode, 'Maximum')
                    %% Maximum of values
                    curstat = max(Cleaned_Response(v,vx));
                elseif strcmpi(Mode, 'Minimum')
                    %% Minimum of values
                    curstat = min(Cleaned_Response(v,vx));
                elseif strcmpi(Mode, 'Geomean')
                    %% Geometric mean of values
                    curstat = geomean(Cleaned_Response(v,vx));
                elseif strcmpi(Mode, 'Sum')
                    %% Geometric mean of values
                    curstat = sum(Cleaned_Response(v,vx));
                elseif strcmpi(Mode, 't-test')
                    %% Calculate t-test vs zero
                    [~,~,~,TestStats] = ttest(Cleaned_Response(v,vx), 0);
                    curstat = TestStats.tstat;
                    
                    %% Dispersion estimates
                elseif strcmpi(Mode, 'StdDev')
                    %% Standard deviation
                    curstat = std(Cleaned_Response(v,vx));
                elseif strcmpi(Mode, 'MeanAbsDev')
                    %% Mean absolute deviation
                    curstat = mad(Cleaned_Response(v,vx), 0);
                elseif strcmpi(Mode, 'MedAbsDev')
                    %% Median absolute deviation
                    curstat = mad(Cleaned_Response(v,vx), 1);
                else
                    %% Unknown mode defined
                    error(['Unknown summary statistic ' Mode ' specified.']);
                end
                %% Store current stat in pixel
                if ~isempty(curstat) && ~isinf(curstat) && ~isnan(curstat)
                    Bp(y,x,v) = curstat;
                end
                %% Store distances & numbers (only need to do once)
                if v == 1
                    % Calculated across all pRFs inside the searchlight
                    if ~isempty(curstat) && ~isinf(curstat) && ~isnan(curstat)
                        Wt(y,x,1) = nansum(dc); % Filtered weights based on inverse distances to searchlight centre
                        Nr(y,x,1) = nansum(vx); % Filtered count of pRFs in searchlight
                        Used_pRFs(1,vx) = Used_pRFs(1,vx) + 1; % Filtered count of how often pRFs were used
                        
                        % When using nearest neighbour mode
                        if Resolution(2) < 0
                            % Filtered max distance of vertices from searchlight center
                            Nr(y,x,3) = nanmax(ed(vx));
                        end
                        
                    end
                    
                    Wt(y,x,2) = nansum(dc); % Unfiltered weights based on inverse distances to searchlight centre
                    Nr(y,x,2) = nansum(vx); % Unfiltered count of pRFs in searchlight
                    Used_pRFs(2,vx) = Used_pRFs(2,vx) + 1; % Unfiltered count of how often pRFs were used
                    
                    if sum(isnan(dc)) > 0 || sum(isnan(vx)) > 0
                        warning('NaN detected. This should not happen. Please request bugfix.')
                    end
                end
            end
        end
    end
    waitbar(x/size(X,2),h);
end
close(h);
