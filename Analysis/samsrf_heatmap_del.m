function samsrf_heatmap_del(Tri, X, Y, Data, Clipping, Cmap)
%
% samsrf_heatmap_del(Tri, X, Y, Data, [Clipping=[0 Inf 0], Cmap='hotcold'])
%
% Plots a heat-color map of Data using the Delaunay tesselation Tri, X and Y (from samsrf_backproj_del).  
%
% Clipping is a vector defining above and below which value data are clipped.
%  The 3rd value toggles whether clipping should take sign into account.
%  Use this when you don't want the colour scale to be centred on 0. 
%  By default this is false, and in that case all values <= Clipping(1)
%  will not be displayed.
%
% Cmap, the colour map, defaults to 'hotcold'. Prefixing this with '-' inverts
%  the colour map. You can also provide a 256x3 RGB colour map matrix yourself.
%
% 21/05/2019 - Created this function (DSS)
%

%% Default inputs
if nargin <= 4
    Clipping = [0 Inf 0];
end
if nargin <= 5
    Cmap = 'hotcold';
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
    
    Data(ClipData <= Clipping(1)) = sign(Data(ClipData <= Clipping(1))) * Clipping(1); % Only appears if no weighting is used
    Data(ClipData >= Clipping(2)) = sign(Data(ClipData >= Clipping(2))) * Clipping(2);
    % Minimum is -1
    Data(Data > 0) = Data(Data > 0) - Clipping(1); 
    Data(Data < 0) = Data(Data < 0) + Clipping(1);
end
% Normalise by the clipping range
Data = Data / range(Clipping(1:2));

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
Colours = NaN(size(Data,1), 3, size(Data,2));
for c = 1:size(Data,2)
    Colours(:,:,c) = Rgb(Idx(:,c),:);
end

%% Plot data 
trisurf(Tri, X, Y, zeros(size(X)), 'FaceVertexCData', Colours, 'EdgeColor', 'none', 'FaceColor', 'interp');
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
view([0 90]);
%     axis square
xlabel('Horizontal position (deg)');
ylabel('Vertical position (deg)');
axis([-1 1 -1 1]*max(X(:)));
hold on
caxis([0 1]);
