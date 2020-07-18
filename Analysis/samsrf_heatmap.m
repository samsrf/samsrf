function Img = samsrf_heatmap(X, Y, Data, Wts, Clipping, WtSatur, Cmap, Interpolate)
%
% Img = samsrf_heatmap(X, Y, Data, [Wts, Clipping=[0 Inf 0], WtSatur=50, Cmap='hotcold', Interpolate=true])
%
% Plots a heat-color map of the image in Data using the coordinate matrices
% X and Y (from samsrf_backproj_srclt). It further uses the weights in Wts 
% to determine the saturation of the colors. Wts defaults to all 1s. 
%
% Clipping is a vector defining above and below which value data are clipped.
%  The 3rd value toggles whether clipping should take sign into account.
%  Use this when you don't want the colour scale to be centred on 0. 
%  By default this is false, and in that case all values <= Clipping(1)
%  will not be displayed.
%
% WtSatur defines the percentile beyond which weights are fully saturated.
%  Use 0 for this when you turn weighting off but want a grey (unsaturated) background.
%  If you Wts is undefined or you provide an empty Wts matrix, the background 
%  is set to whatever colour 0 is (this is usually fine but may not always be desirable).
%
% Cmap, the colour map, defaults to 'hotcold'. Prefixing this with '-' inverts
%  the colour map. You can also provide a 256x3 RGB colour map matrix yourself.
%
% Interpolate is a boolean that toggles face interpolation (default = true).
%  Use this for smoother back-projection images.
%
% Optionally returns the image matrix Img, in which case plotting is turned off.
%
% 19/07/2020 - SamSrf 7 version (DSS)
%

%% Default inputs
if nargin == 3 || isempty(Wts)
    Wts = ones(size(Data));
end
if nargin <= 4
    Clipping = [0 Inf 0];
end
if nargin <= 5
    WtSatur = 50;
end
if nargin <= 6
    Cmap = 'hotcold';
end
if nargin <= 7 
    Interpolate = true;
end
if Interpolate
    Interpolate = 'interp';
else
    Interpolate = 'flat';
end

if length(Clipping) == 1
    Clipping = [Clipping Inf 0];
elseif length(Clipping) == 2
    Clipping = [Clipping 0];
end

%% Clipping data
if Clipping(3)
    % If retaining sign for clipping
    ClipData = Data;
    
    % Set out-of-range values to clipping thresholds
    Data(ClipData <= Clipping(1)) = Clipping(1);
    Data(ClipData >= Clipping(2)) = Clipping(2);
    % Minimum is 0
    Data = Data - Clipping(1); 
else    
    % If ignoring sign for clipping
    ClipData = abs(Data);
    Clipping(1:2) = abs(Clipping(1:2)); % Ensure no negative thresholds
    if Clipping(1) == Clipping(2)
        % Both clipping thresholds cannot be the same
        error('When ignoring sign, lower & upper clipping level cannot be identical!');
    elseif Clipping(1) > Clipping(2)
        % Lower clipping level must be below upper level
        error('When ignoring sign, lower clipping level must be below upper clipping level!');
    end
    
    % Set out-of-range values to clipping thresholds
    if mean(Wts(:)) ~=0 && var(Wts(:)) ~= 0 % If all weights are 1, weighting is off
        Wts(ClipData <= Clipping(1)) = 0; % Unweight values below clipping level
    end
    Data(ClipData <= Clipping(1)) = sign(Data(ClipData <= Clipping(1))) * Clipping(1); % Only appears if no weighting is used
    Data(ClipData >= Clipping(2)) = sign(Data(ClipData >= Clipping(2))) * Clipping(2);
    % Minimum is -1
    Data(Data > 0) = Data(Data > 0) - Clipping(1); 
    Data(Data < 0) = Data(Data < 0) + Clipping(1);
end
% Normalise by the clipping range
Data = Data / range(Clipping(1:2));

%% Weight saturation
Wts(isnan(Wts)) = 0; % Just in case
Wts(isinf(Wts)) = 0;
nzw = Wts ~= 0; % Non-zero weights
PctWt = prctile(Wts(nzw), WtSatur); % Percentile of weights
% If all weights empty
if ~isnan(PctWt)
    Wts = Wts / PctWt;
    Wts(Wts > 1) = 1;
end

%% Determine colours
if ischar(Cmap)
    % Create colour map
    if Cmap(1) == '-'
        InvertColours = true;
        Cmap = Cmap(2:end);
    else
        InvertColours = false;
    end
    Rgb = colormap(eval([Cmap '(256)']));
    if InvertColours
        Rgb = flipud(Rgb);
    end    
else
    % Colour map provided
    Rgb = Cmap;
end
colormap(Rgb);
% Grey image
Grey = ones(size(Data)) / 2;
% Convert data to color indeces
if Clipping(3)
    % If retaining sign for clipping
    Idx = ceil(Data * 255) + 1;
else
    % If ignoring sign for clipping
    Idx = ceil((Data+1) * 127) + 1;
    Idx(Idx > 256) = 256;
end

% Create color channels
Img = NaN(size(Data,1), size(Data,2), 3);
for c = 1:3
    Img(:,:,c) = (reshape(Rgb(Idx,c), size(Data,1), size(Data,2)) .* Wts) + (Grey .* (1-Wts));
end

%% Plot data 
if nargout == 0
    surf(X, Y, -Data + min(Data(:)), Img, 'EdgeColor', 'none', 'FaceColor', Interpolate);
    cb = colorbar;
    if Clipping(3)
        % If retaining sign for clipping
        set(cb, 'xtick', 0:.5:1, 'xticklabel', round([Clipping(1) (Clipping(1)+Clipping(2))/2 Clipping(2)],2));
    else
        % If ignoring sign for clipping
        if Clipping(1) == 0
            % Centre of colour scale is zero
            set(cb, 'xtick', 0:.5:1, 'xticklabel', round([-Clipping(2) 0 +Clipping(2)],2));
        else
            % Colour scale starts at lower clipping level
            set(cb, 'xtick', [0 .47 .53 1], 'xticklabel', round([-Clipping(2) -Clipping(1) +Clipping(1) +Clipping(2)],2));
        end
    end
    axis square
    xlabel('Horizontal position (deg)');
    ylabel('Vertical position (deg)');
    axis([-1 1 -1 1]*max(X(:)));
    hold on
    caxis([0 1]);
    view(0,90); % Ensure camera views from the top down
else
    close gcf
end