function samsrf_glm(SrfCell, X, Xnames, Roi, GlmFile)
%
% samsrf_glm(SrfCell, X, Xnames, [Roi='', GlmFile='glm'])
%
% Runs a GLM analysis on the surface data files in the cell array SrfCell. 
% The matrix X contains the design matrix. You may choose to convolve this with 
% a HRF (canonical or subject-specific). You may also add in covariates of 
% no interest. The function will add the constant term in the right hand 
% columns but everything else must be defined by you. Xnames is a cell 
% string defining the names of the regressors. Roi is the filename of a ROI
% label if you want to restrict to analysis to only that ROI.
%
% It should go without saying that all Srf's in SrfCell should contain
% equal numbers of vertices and they should have been expanded if
% necessary.
%
% The function saves a new surface data file called lh/rh_glm (or any name
% you define in GlmFile). This contains in each row the beta estimates of
% each regressor in the design matrix. The final row contains the residuals.
%
% 20/04/2022 - SamSrf 8 version (DSS)
%

if length(Xnames) ~= size(X,2)
    error('Number of names must match number of regressors in X!');
end

new_line; 
disp('Running general linear model analysis...');

%% Default parameters
if nargin < 4
    Roi = '';
    GlmFile = 'glm';
elseif nargin < 5
    GlmFile = 'glm';
end

%% If SrfCell is a string
if isa(SrfCell, 'char')
    SrfCell = {SrfCell};
end

%% Create data matrix
nRuns = length(SrfCell); % Number of runs
Y = []; % Data matrix
Ct = []; % Columns for constant term
% Add each run to data matrix (& add constant terms)
for r = 1:nRuns
    load(EnsurePath(SrfCell{r}));
    Srf = samsrf_expand_srf(Srf);
    Y = [Y; Srf.Data];
    % Add constant term
    ct = [];
    for t = 1:nRuns
        if t == r
            % Add ones for this run
            ct = [ct ones(size(Srf.Data,1),1)];
        else
            % Add zeros for other runs
            ct = [ct zeros(size(Srf.Data,1),1)];
        end
    end
    Ct = [Ct; ct];
end
X = [X Ct];

% Plot design matrix
figure('name', 'Design matrix');
hold on
cm = berlin(size(X,2)); % Colour scheme
% Plot each regressor
for i = 1:size(X,2)
    plot(X(:,i)/10 + i, 'linewidth', 2, 'color', cm(i,:));
end
% Regressor names
for r = 1:nRuns
    Xnames{end+1} = ['Constant #' num2str(r)];
end
set(gca, 'ytick', 1:size(X,2), 'yticklabel', Xnames);
xlabel('Volume (#)');
saveas(gcf, [GlmFile '.fig']);

%% Add version number
Srf.Version = samsrf_version;

%% Load ROI mask
if isempty(Roi)
    new_line; disp('Using all vertices...');
    mver = 1:size(Srf.Vertices,1);
else
    new_line; disp('Reading ROI mask...')
    mver = samsrf_loadlabel(Roi);
    disp([' Loading ' Roi ': ' num2str(size(mver,1)) ' vertices']);
end

%% Run linear regression 
disp('Running GLM on vertices...'); 
disp(' Please stand by...');
B = NaN(size(X,2)+1, length(mver));
parfor v = 1:length(mver)
    if ~isnan(sum(Y(:,mver(v))))
        [cb,bint,cr,rint,st] = regress(Y(:,mver(v)),X);
        B(:,v) = [cb; st(4)];  % Betas for current vertex & error variance
    end
end

%% Save as file
Srf.Functional = 'GLM coefficients';
Srf.Data = NaN(size(X,2)+1, size(Srf.Data,2));
Srf.Data(:,mver) = B; % Add betas into full surface map
% Ensure Xnames is a column
if size(Xnames,1) == 1
    Xnames = Xnames';
end
% Add regressor names to structure
Srf.Values = Xnames;
Srf.Values{end+1} = 'Error';
% Save design matrix & observed data
Srf.X_glm = X;
Srf.Y = Y;
% Compress Srf if necessary
Srf = samsrf_compress_srf(Srf, mver);
% Save Srf file
save([Srf.Hemisphere '_' GlmFile], 'Srf');
disp(['Saved ' Srf.Hemisphere '_' GlmFile '.mat.']);
new_line;


