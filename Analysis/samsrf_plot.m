function [Res, FigHdl] = samsrf_plot(SrfDv, ValDv, SrfIv, ValIv, Bins, Roi, Threshold, Mode, Colour, BootParams)
%
% [Res, FigHdl] = samsrf_plot(SrfDv, ValDv, SrfIv, ValIv, Bins, [Roi='', Threshold=[0.01 -Inf/0 Inf], Mode='Scatter', Colour='k', BootParams=[1000, 2.5, 97.5]])
%
% Plots the data defined by ValDv in SrfDv against the data in ValIv from SrfIv
%  (thus, SrfDv is the dependent variable, SrfIv the independent variable).
% A different set of plots can be made, including binned summary statistics,
%  scatter plots, or sliding window averages.
%
% Note that this function tacitly assumes you are using pRF data for the 
%  independent variable but this is -not- mandatory. It will throw up a
%  warning if the first row of SrfIv.Values is not 'R^2'. In that case, 
%  you need to think carefully about how to filter your SrfIv.Data.
%
% The output Res can be used for curve fitting (e.g. samsrf_indiecurve) 
%  or any other further analysis. It usually contains the following columns:
%  Bins, Summary Statistics, Confidence Intervals, Number of Vertices
%  (If generating a scatter plot, only the first two columns are created)
%
% The second output FigHdl is a figure handle to the main part of the plot
%  (i.e. the summary curve rather than the confidence intervals).
%
% IMPORTANT NOTE: Binning analyses like this can suffer from regression artifacts!
%   Results can be misleading & should be interpreted with caution.
%
%
% SrfDv/SrfIv:  Srf structures. When plotting data within a map, use the -same- Srf.
%
% ValDv/ValIv:  The value name as given by Srf.Values. Returns an error if 
%                more than one entry of that name exists.
%
%               This can also be 'Eccentricity' or 'Polar', in which case it 
%                may throw up an error if SrfIv.Data(2:3,:) are not X0 and Y0.
%
%               Prefixing ValDv/ValIv by ':' uses unsmoothed data (if possible).
%
%               It can also refer to anatomical statistics such as
%                'Area', 'Thickness', or 'Curvature'. (These will always be
%                anatomical so don't have any Srf.Values of the same names!)
%
%               Finally, if using M/EEG data you can also use 
%                'Sensor' or 'Time', equivalent to Srf.SensorNums or Srf.TimePts. 
%
% Bins:         Vector with the boundaries of the bins to be analysed.
%                0:8 would return bins centred on 0.5:7.5, respectively.
%                If the last component in this vector is NaN, the mean of the 
%                independent variable per bin is calculated. Otherwise, the 
%                conventional bin centre is plotted.
%
%               You can also decide to use a sliding window approach by
%                giving three entries: Starting-point, End-point, Window-width
%                In that case, Window-width must be negative so e.g.:
%                [0 8 -1] will produce sliding window with width 1 from 0 to 8.
%
%               When creating a scatter plot instead of binned summaries,
%                only the lowest and highest values are used to restrict the range.
%                If Bins is undefined this will then default to [-Inf Inf].
%
% Roi:          ROI label (without the file extension).
%                Defaults to '' and data isn't restricted to any ROI.
%               You can also define a vector of ROI vertex indeces directly.              
%
% Threshold:    The R^2 threshold for the IV to be included (defaults to 0.01).
%                When plotting normalised R^2 this threshold is first applied 
%                to the noise ceiling and then again to the normalised R^2)
%
%               If this is a vector, the second & third entry define the
%                range of data values in SrfDv to be included (excluding
%                the values given). Defaults to [-Inf Inf] unless the data
%                you are plotting cannot be <=0 (e.g. Sigma) in which case
%                it defaults to [0 Inf].
%
% Mode:         How bins are statistically summarised:
%                'Scatter':     Scatter plot without summary statistics (default)
%                'Mean':        Arithmetic mean 
%                'Median':      Median
%                'Sum':         Sum
%                'Geomean':     Geometric mean
%
% Colour:       Colour of the plot this produces (defaults to black = 'k').
%                Can also be a 1x3 vector with RGB values.
%
% BootParams:   1x3 vector with bootstrapping parameters for confidence intervals.
%                The first input defines the number of bootstrap iterations.
%                The second and third input define the percentiles for the CI.
%                If only two inputs are defined, the second input defines
%                the percentage of the interval (e.g. 95 for 95% CI).
%
% 20/09/2024 - Adjusted scaling of X-axis (DSS)
%

%% Expand Srfs if necessary
SrfDv = samsrf_expand_srf(SrfDv);    
SrfIv = samsrf_expand_srf(SrfIv);

%% Check compatibility
if size(SrfDv.Data,2) ~= size(SrfIv.Data,2)
    samsrf_error('SrfDv & SrfIv are not the same mesh!');
end

%% Is this M/EEG data?
if strcmpi(SrfDv.Hemisphere, 'eeg')
    SrfDv.Data(end+1:end+2,:) = [SrfDv.SensorNums; SrfDv.TimePts];
    SrfIv.Data(end+1:end+2,:) = [SrfIv.SensorNums; SrfIv.TimePts];
    SrfDv.Values{end+1} = 'Sensor';
    SrfIv.Values{end+1} = 'Sensor';
    SrfDv.Values{end+1} = 'Time';
    SrfIv.Values{end+1} = 'Time';
end

%% Plotting bin summary stat?
if isnan(Bins(end))
    IvSumStats = true;
    Bins = Bins(1:end-1);
else
    IvSumStats = false;
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
if strcmpi(ValDv, 'Sigma') || strcmpi(ValDv, 'Centre') || strcmpi(ValDv, 'Surround') || strcmpi(ValDv, 'Visual area') ... 
    || strcmpi(ValDv, 'Fwhm') || strcmpi(ValDv, 'Spread') || strcmpi(ValDv, 'Sigma1') || strcmpi(ValDv, 'Sigma2') 
    if Threshold(2) < 0
        Threshold(2) = 0;
    end
end
 
if nargin <= 7  Mode = 'Scatter'; end % Default mode is scatter plot
if nargin <= 8  Colour = 'k'; end % Colour defaults to black

if nargin <= 9  BootParams = [1000, 2.5, 97.5]; end % Bootstrapping defaults to 1000 iterations for 95%CI
% If BootParams not fully defined
if length(BootParams) == 1
    BootParams = [BootParams, 2.5, 97.5];
elseif length(BootParams) == 2
    % Only 2 parameters, so BootParams(2) defined CI%
    Pct = (100-BootParams(2)) / 2;
    BootParams = [BootParams(1) Pct 100-Pct];
end

%% Retrieve data
% Loop thru dependent & independent variable
for i_var = 1:2
    if i_var == 1
        Srf = SrfDv;
        Val = ValDv;
        ValLab = 'ValDv';
    elseif i_var == 2
        Srf = SrfIv;
        Val = ValIv;
        ValLab = 'ValIv';
    end       
    
    % Use unsmoothed data?
    RawLabel = '';
    if Val(1) == ':'
        Val = Val(2:end);
        if isfield(Srf, 'Raw_Data')
            Srf.Data = Srf.Raw_Data;
            RawLabel = ' (unsmoothed)';
        end
    end

    % Which data to use?
    if strcmpi(Val, 'Eccentricity')
        % Eccentricity
        Data = sqrt(Srf.Data(2,:).^2 + Srf.Data(3,:).^2);
    elseif strcmpi(Val, 'Polar')
        % Polar angle
        Data = atan2(Srf.Data(3,:), Srf.Data(2,:)) / pi * 180;
    elseif strcmpi(Val, 'Area')
        % Surface area
        Data = Srf.Area;
    elseif strcmpi(Val, 'Thickness')
        % Cortical thickness
        Data = Srf.Thickness;
    elseif strcmpi(Val, 'Curvature')
        % Curvature
        Data = Srf.Curvature;
    else
        % Anything else
        loc = strcmpi(Srf.Values, Val);
        if sum(loc) > 1
            samsrf_error([ValLab ' ' Val ' is ambiguous!']);
        elseif sum(loc) == 0
            samsrf_error([ValLab ' ' Val ' does not exist!'])
        end
        Data = Srf.Data(loc,:);
    end
        
    if i_var == 1
        DataDv = Data;
        ValDv = Val;
        RawLabelDv = RawLabel;
    elseif i_var == 2
        DataIv = Data;
        ValIv = Val;
        RawLabelIv = RawLabel;
    end
end

%% Filter data & restrict to ROI
GoF = true(1,size(SrfIv.Data,2)); % Goodness-of-fit placeholder
% Is 1st row R^2 (usually in pRF maps only)?
if strcmpi(SrfIv.Values{1}, 'R^2') || strcmpi(SrfIv.Values{1}, 'nR^2')
    GoF(SrfIv.Data(1,:) <= Threshold(1)) = false; % Unlabel bad fits
else
    % Warn user if 1st row isn't R^2
    samsrf_disp('WARNING: 1st data row isn''t R^2 so using all good data... Do you really want this?');
end
% Limit data range & remove other rubbish
GoF(DataDv <= Threshold(2) | DataDv >= Threshold(3) ...
    | isnan(DataDv) | isnan(DataIv) | isinf(DataDv) | isinf(DataIv)) = false;

% Load ROI label
if ~isempty(Roi)
    % ROI vertex indeces
    if ischar(Roi)
        RoiVx = samsrf_loadlabel(Roi); % Load from file
    else
        RoiVx = Roi; % Define indeces directly
        Roi = 'List of data points';
    end
    % Label ROI vertices
    RoiLab = false(1,size(SrfIv.Data,2));
    RoiLab(RoiVx) = true;
else
    % No ROI so label all
    RoiLab = true(1,size(SrfIv.Data,2));
end
% Filtered & restricted vectors
DataDv = DataDv(GoF & RoiLab); % Dependent variable
DataIv = DataIv(GoF & RoiLab); % Independent variable
Theta = mod(DataIv, 360); % Angles between 0-360 for circular plots

%% Produce outputs
Res = [];
if strcmpi(Mode, 'Scatter')
    if ~isempty(Bins)
        samsrf_disp('Scatter plot: lowest & highest bin are used as bounds');
    end
    % Scatter plot
    Res = [DataIv' DataDv'];
else % Binning analysis
    % Is sliding window desired?
    if length(Bins) == 3 && sign(Bins(3)) == -1
        SliWin = abs(Bins(3));
        Bins = Bins(1):(Bins(2)-Bins(1))/100:Bins(2);
    else
        SliWin = NaN;
    end
    % Binning analysis
    for b = 1:length(Bins)-1
        if isnan(SliWin)
            % Fixed binning analysis
            if strcmpi(ValIv,'Polar') || strcmpi(ValIv,'Phi') || strcmpi(ValIv,'Theta') || strcmpi(ValIv,'Phase')
                % Circular measures
                CurAng = mod(Bins(b),360); % Current bin start (0-360)
                CurRng = mod(Bins(b+1),360) - CurAng; % Range of angles in bin
                if CurRng < 0
                    CurRng = CurRng + 360; % Ensure range is positive
                end
                RelAng = Theta - CurAng; % Relative angle to bin start
                RelAng(RelAng < 0) = RelAng(RelAng < 0) + 360; % Anything negative must become positive
                CurDv = DataDv(RelAng <= CurRng); % Current DV data
                CurIv = DataIv(RelAng <= CurRng); % Current IV data
            else
                % Linear measures
                CurDv = DataDv(DataIv >= Bins(b) & DataIv < Bins(b+1)); % Current DV data
                CurIv = DataIv(DataIv >= Bins(b) & DataIv < Bins(b+1)); % Current IV data
            end
            if IvSumStats
                CurBin = nanmean(CurIv); % Mean 
            else
                CurBin = mean(Bins(b:b+1)); % Middle of bin
            end
        else
            % Sliding window analysis
            if strcmpi(ValIv,'Polar') || strcmpi(ValIv,'Phi') || strcmpi(ValIv,'Theta') || strcmpi(ValIv,'Phase')
                % Circular measures
                CurAng = mod(Bins(b)-SliWin/2,360); % Current bin start (0-360)
                CurRng = mod(Bins(b)+SliWin/2,360) - CurAng; % Range of angles in bin
                if CurRng < 0
                    CurRng = CurRng + 360; % Ensure range is positive
                end
                RelAng = Theta - CurAng; % Relative angle to bin start
                RelAng(RelAng < 0) = RelAng(RelAng < 0) + 360; % Anything negative must become positive
                CurDv = DataDv(RelAng <= CurRng); % Current DV data
                CurUv = DataIv(RelAng <= CurRng); % Current IV data
            else
                % Linear measures
                CurDv = DataDv(DataIv >= Bins(b)-SliWin/2 & DataIv < Bins(b)+SliWin/2); % Current DV data
                CurIv = DataIv(DataIv >= Bins(b)-SliWin/2 & DataIv < Bins(b)+SliWin/2); % Current IV data
            end
            if IvSumStats
                CurBin = nanmean(CurIv); % Mean 
            else
                CurBin = Bins(b); % Current bin
            end
        end
        CurN = length(CurDv); % Number of vertices
        
        %% Summary statistic & confidence interval
        if CurN > 1 
            if strcmpi(Mode, 'Mean')
                Bs = bootstrp(BootParams(1), @nanmean, CurDv); % Bootstrap distribution
                CurDv = nanmean(CurDv); % Mean
            elseif strcmpi(Mode, 'Median')
                Bs = bootstrp(BootParams(1), @nanmedian, CurDv); % Bootstrap distribution
                CurDv = nanmedian(CurDv); % Median
            elseif strcmpi(Mode, 'Sum')
                Bs = bootstrp(BootParams(1), @nansum, CurDv); % Bootstrap distribution
                CurDv = nansum(CurDv); % Sum
            elseif strcmpi(Mode, 'Geomean')
                Bs = exp(bootstrp(BootParams(1), @nanmean, log(CurDv))); % Bootstrap distribution
                CurDv = exp(nanmean(log(CurDv))); % Geometric mean
            else
                samsrf_error('Invalid summary statistic specified!');
            end
            CurCi = prctile(Bs, BootParams(2:3)); % 95% confidence interval
            % In case it is column vector
            if size(CurCi,1) > 1
                CurCi = CurCi'; % Was largely for Octave 4 compatibility & probably irrelevant
            end
            CurCi = CurCi - CurDv; % CIs relative to mean
        else
            CurCi = [NaN NaN]; %% No CI for single or zero entry
        end
        % If bin is empty
        if isempty(CurDv)
            CurDv = NaN;
        end        
        % Store bin results
        Res = [Res; CurBin CurDv CurCi CurN];        
    end
end

%% Plot data
if strcmpi(Mode, 'Scatter')
    % Scatter plot
    if nargin < 5
        Bins = [-Inf Inf];
    end
    Res(Res(:,1)<Bins(1) | Res(:,1)>Bins(end),:) = []; % Restrict range 
    hold on
    FigHdl = scatter(Res(:,1), Res(:,2), 'filled', 'markerfacecolor', Colour);
    alpha(FigHdl, 0.1);
    ylabel([ValDv RawLabelDv]);
    xlabel([ValIv RawLabelIv]);
else
    if strcmpi(ValIv,'Polar') || strcmpi(ValIv,'Phi') || strcmpi(ValIv,'Theta') || strcmpi(ValIv,'Phase')
        % Bin polar plot
        x = [1:size(Res,1) 1];
        FigHdl = polar(Res(x,1)/180*pi, Res(x,2), ['-o' Colour]); % Summary curve);
        hold on
        polar(Res(x,1)/180*pi, Res(x,2)+Res(x,3), [':' Colour]); % CI lower limit
        polar(Res(x,1)/180*pi, Res(x,2)+Res(x,4), [':' Colour]); % CI upper limit
        axis square
        s = nanmax(Res(x,2)+Res(x,3)*3);
        if ~isnan(s)
            axis([-s s -s s]);
        end
    else
        % Bin plot
        FigHdl = plot(Res(:,1), Res(:,2), 'linewidth', 2, 'color', Colour, 'Marker', 'o'); % Summary curve
        hold on
        plot(Res(:,1), Res(:,2)+Res(:,3), 'linestyle',':', 'color', Colour); % CI lower limit
        plot(Res(:,1), Res(:,2)+Res(:,4), 'linestyle',':', 'color', Colour); % CI upper limit
        ylabel([Mode ' ' ValDv RawLabelDv]);
        xlabel([ValIv RawLabelIv]);
        if length(Bins) == 3 && Bins(end) < 0
            % Sliding windows
            xlim(Bins(1:2));
        else
            % Fixed bins 
            xlim([Bins(1) Bins(end)]);
        end
    end
end

% Reformat title
Roi(strfind(Roi,'_')) = '-';
Roi([strfind(Roi,'/'),  strfind(Roi,'\')]) = ' ';
title(Roi);
