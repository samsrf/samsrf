function Res = samsrf_wedgeplot(SrfDv, Value, SrfIv, Wedges, Rings,  Roi, Threshold, Mode, Cmap)
%
% Res = samsrf_wedgeplot(SrfDv, Value, SrfIv, Wedges, Rings, Roi, Threshold, Mode, Cmap)
%
% Plots the data defined by Value in SrfDv against the data in ValIv from SrfIv
%  (thus, SrfDv is the dependent variable, SrfIv the independent variable).
% It plots the summary statistic for each wedge segment of the visual field as defined.
%  The colour code plots the dependent variable.
%
% Note that this function tacitly assumes you are using pRF data for the 
%  independent variable but this is -not- mandatory. It will throw up a
%  warning if the first row of SrfIv.Values is not 'R^2'. In that case, 
%  you need to think carefully about how to filter your SrfIv.Data.
%
% The output Res contains the following columns:
%  Xs, Ys, Summary Statistics, Number of Vertices
%  
%
% SrfDv/SrfIv:  Srf structures. When plotting data within a map, use the -same- Srf.
%
% Value:        The value name as given by Srf.Values. Returns an error if 
%                more than one entry of that name exists.
%               This can also be 'Eccentricity' or 'Polar', in which case it 
%                may throw up an error if Srf.Data(2:3,:) are not X0 and Y0.
%               Prefixing Value by ':' uses unsmoothed data (if possible).
%                In some scenarios you might want to plot smoothed map
%                positions against unsmoothed dependent variables. In that case,
%                you will need to have a smoothed SrfIv and unsmoothed SrfDv!
%               Finally, this can also refer to anatomical statistics such as
%                'Area', 'Thickness', or 'Curvature'. (These will always be
%                anatomical so don't have any Srf.Values of the same names!)
%
% Wedges/Rings: Vectors with the polar and eccentricity boundaries of the segments.
%
% Roi:          ROI label (without the file extension).
%                Defaults to '' and data isn't restricted to any ROI.
%
% Threshold:    The R^2 threshold for the IV to be included (defaults to 0.01).
%                When plotting normalised R^2 this threshold is first applied 
%                to the noise ceiling and then again to the normalised R^2)
%               If this is a vector, the second & third entry define the
%                range of data values in SrfDv to be included (excluding
%                the values given). Defaults to [-Inf Inf] unless the data
%                you are plotting cannot be <=0 (e.g. Sigma) in which case
%                it defaults to [0 Inf].
%
% Mode:         How bins are statistically summarised:
%                'Mean':        Arithmetic mean (default)
%                'Median':      Median
%                'Sum':         Sum
%                'Geomean':     Geometric mean
%
% Cmap:         Colour map as string. Defaults to 'jet'.
%
% 28/10/2020 - Written (DSS)
% 26/11/2020 - Fixed bug with default inputs (DSS)
% 09/12/2020 - Changed output so it now contains NaN data for 0 vertex ROIs (DSS)
% 18/12/2020 - Now produces a proper wedge plot rather than scatter plot (EA)
%              Adjusted alpha scaling based on number of vertices per segment (DSS)
%

%% Expand Srfs if necessary
SrfDv = samsrf_expand_srf(SrfDv);    
SrfIv = samsrf_expand_srf(SrfIv);

%% Check compatibility
if size(SrfDv.Data,2) ~= size(SrfIv.Data,2)
    error('SrfDv & SrfIv are not the same mesh!');
end

%% Default inputs
if nargin == 5  Roi = ''; end % ROI defaults to all

if nargin <= 6  Threshold = [0.01 -Inf Inf]; end % Threshold defaults to R^2>0.01 and all values
% If Threshold not fully defined
if length(Threshold) == 1
    Threshold = [Threshold -Inf Inf];
elseif length(Threshold) == 2
    Threshold = [Threshold Inf];
end
% For measures where DV cannot be zero
if strcmpi(Value, 'Sigma') || strcmpi(Value, 'Centre') || strcmpi(Value, 'Surround') || strcmpi(Value, 'Visual area') ... 
    || strcmpi(Value, 'Fwhm') || strcmpi(Value, 'Spread') || strcmpi(Value, 'Sigma1') || strcmpi(Value, 'Sigma2') 
    if Threshold(2) < 0
        Threshold(2) = 0;
    end
end
 
if nargin <= 7  Mode = 'Mean'; end % Summary stat defaults to arithmetic mean
if nargin <= 8  Cmap = 'jet'; end % Colour scheme defaults to jet

%% Raw or smoothed data?
if Value(1) == ':'
    IsRaw = true;
    Value = Value(2:end);
else
    IsRaw = false;
end

%% Dependent variable
if IsRaw && isfield(SrfDv, 'Raw_Data')
    SrfDv.Data = SrfDv.Raw_Data;
end
% Which data to use?
if strcmpi(Value, 'Eccentricity')
    % Eccentricity
    Data = sqrt(SrfDv.Data(2,:).^2 + SrfDv.Data(3,:).^2);
elseif strcmpi(Value, 'Polar')
    % Polar angle
    Data = atan2(SrfDv.Data(3,:), SrfDv.Data(2,:)) / pi * 180;
elseif strcmpi(Value, 'Area')
    % Surface area
    Data = SrfDv.Area;
elseif strcmpi(Value, 'Thickness')
    % Cortical thickness
    Data = SrfDv.Thickness;
elseif strcmpi(Value, 'Curvature')
    % Curvature
    Data = SrfDv.Curvature;
elseif strcmpi(Value, 'nR^2')
    if sum(strcmpi(SrfDv.Values, 'Noise Ceiling'))
        % Normalised goodness-of-fit
        Nc = SrfDv.Data(strcmpi(SrfDv.Values, 'Noise Ceiling'),:); % Extract noise ceiling
        Data = SrfDv.Data(1,:) ./ Nc; % R^2 relative to noise ceiling
        Data(Nc <= Threshold(1)) = 0; % Threshold vertices under a noise ceiling level
    else
        warning('No Noise Ceiling in input Srf! Using R^2...');
        Data = SrfDv.Data(1,:);
    end
else
    %% Anything else
    loc = strcmpi(SrfDv.Values, Value);
    if sum(loc) > 1
        error([ValLab ' ' Value ' is ambiguous!']);
    elseif sum(loc) == 0
        error([ValLab ' ' Value ' does not exist!'])
    end
    Data = SrfDv.Data(loc,:);
end

%% Independent variable
if IsRaw && isfield(SrfIv, 'Raw_Data')
    SrfIv.Data = SrfIv.Raw_Data;
end
% Visual field coordinates 
X = SrfIv.Data(2,:);
Y = SrfIv.Data(3,:);
[Theta, Rho] = cart2pol(X,Y); % Polar angle & eccentricity
Theta = Theta / pi * 180; % Polar angles in degrees
Theta = mod(Theta, 360); % Angles between 0-360 

%% Filter data & restrict to ROI
GoF = true(1,size(SrfIv.Data,2)); % Goodness-of-fit placeholder
% Is 1st row R^2 (usually in pRF maps only)?
if strcmpi(SrfIv.Values{1}, 'R^2') || strcmpi(SrfIv.Values{1}, 'nR^2')
    GoF(SrfIv.Data(1,:) <= Threshold(1)) = false; % Unlabel bad fits
else
    % Warn user if 1st row isn't R^2
    warning('1st data row isn''t R^2 so using all good data... Do you really want this?');
end
% Limit data range & remove other rubbish
GoF(Data <= Threshold(2) | Data >= Threshold(3) | isnan(Data) | isnan(GoF) | isinf(Data)) = false;

% Load ROI label
if ~isempty(Roi)
    % ROI vertex indeces
    RoiVx = samsrf_loadlabel(Roi);
    if isnan(RoiVx)
        error(['Could not load ROI ' Roi '!']);
    end
    % Label ROI vertices
    RoiLab = false(1,size(SrfIv.Data,2));
    RoiLab(RoiVx) = true;
else
    % No ROI so label all
    RoiLab = true(1,size(SrfIv.Data,2));
end
% Filtered & restricted vectors
Data = Data(GoF & RoiLab); % Dependent variable
Theta = Theta(GoF & RoiLab); % Polar angle
Rho = Rho(GoF & RoiLab); % Eccentricity

%% Segment binning analysis
Res = [];
for R = 1:length(Rings)-1 % Loop thru eccentricity rings
    for t = 1:length(Wedges)-1 % Loop thru eccentricity rings
        % Select data
        CurDat = Data(Rho >= Rings(R) & Rho < Rings(R+1) & Theta >= Wedges(t) & Theta < Wedges(t+1));
        
        % Wedge location & number of vertices
        [CurX, CurY] = pol2cart(mean(Wedges(t:t+1))/180*pi, mean(Rings(R:R+1)));
        CurN = length(CurDat); 
        
        % Summary statistic 
        if strcmpi(Mode, 'Mean')
            CurDat = nanmean(CurDat); % Mean
        elseif strcmpi(Mode, 'Median')
            CurDat = nanmedian(CurDat); % Median
        elseif strcmpi(Mode, 'Sum')
            CurDat = nansum(CurDat); % Sum
        elseif strcmpi(Mode, 'Geomean')
            CurDat = exp(nanmean(log(CurDat))); % Mean
        else
            error('Invalid summary statistic specified!');
        end
        
        % Output results
        Res = [Res; CurX CurY CurDat CurN];
    end
end

%% Plot results
A = Res(:,4);
A(A > median(A)) = median(A);
A = log(A);
A(isinf(A)) = -0.001;
A = A - min(A);
A = A / max(A);
D = Res(:,3);
D(isnan(D)) = 0;
w = Wedges(2)-Wedges(1);
r = Rings(2)-Rings(1);

figure; axis square; hold on;
colors = colormap(Cmap);

for S = 1:length(D)
    [Th, R] = cart2pol(Res(S,1),Res(S,2));
    CurrD = D(S);
    cIdx = round(((CurrD - min(D))/ (max(D)-min(D)) * 255))+1;
    c = colors(cIdx, :);
    alpha = A(S);
    
    t = linspace(Th-deg2rad(w/2), Th+deg2rad(w/2));
    x1 = (R+r/2)*cos(t);
    x2 = fliplr((R-r/2)*cos(t));
    y1 = (R+r/2)*sin(t);
    y2 = fliplr((R-r/2)*sin(t));

    segm = fill([x1 x2], [y1 y2], c);
    set(segm, 'facealpha', alpha)
end

scatter(Res(:,1), Res(:,2), 30, D, 'filled');
colormap(Cmap);
colorbar
ax = gca;
text(ax.XLim(1), ax.YLim(1)+0.6, ['Alpha [0 1] = # of Vertices [' num2str(min(Res(:,4))) ' ' num2str(max(Res(:,4))) ']'])
xlabel('Horizontal position (deg)');
ylabel('Vertical position (deg)');
