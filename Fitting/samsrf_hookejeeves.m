function  [FinalParams, FinalErr] = samsrf_hookejeeves(ErrFcn, SeedParams, InitStep, OnlyPos, NumMaxIter, NumShrinks)
%
% [FinalParams, FinalErr] = samsrf_hookejeeves(ErrFcn, SeedParams, InitStep, [OnlyPos="All-False", NumMaxIter=100, NumShrinks=10])
%
% Optimises parameters for a function using the Hooke-Jeeves (pattern-search) algorithm.
% This method is in many situations faster than fminsearch or fminunc
% although this may come at the cost of losing some precision & those other
% algorithms probably outperform it when there are many free parameters.
% However, this function has some additional functionality for specifying
% search granularity & positive-only parameters useful for pRF fitting.
%
%   ErrFcn:         Error function (requires parameters; cf. fminsearch usage)
%   SeedParams:     SeedParameters (same as with fminsearch usage)
%   InitStep:       Initial step size for pattern search per parameter
%   OnlyPos:        Boolean per parameter to define if it must be > 0 
%                     (defaults to all false)
%   NumMaxIter:     Maximal number of iterations 
%                     (defaults to 100 but note that in pRF fitting this is only 15!) 
%   NumShrinks:     Maximal number of shrinking steps 
%                     (defaults to 10 but note that in pRF fitting this is only 3!)
%
% 13/04/2022 - Written (DSS)
%

%% Ensure sound inputs
NumParams = length(SeedParams); % Number of free parameters
if length(InitStep) ~= NumParams
    error('Mismatch between number of seed parameters & step sizes!');
end
if nargin < 4
    OnlyPos = false(1,NumParams);
end
if length(OnlyPos) ~= NumParams
    error('Mismatch between number of seed parameters & positive flags!');
end
OnlyPos = logical(OnlyPos); % In case doubles 
% Default number of maximum iterations
if nargin < 5
    NumMaxIter = 100;
end
% Default number of shrinking steps allowed
if nargin < 6
    NumShrinks = 10;
end

%% Initialise
CurParams = SeedParams; % Current parameters  

%% Offset matrix
Offsets = [eye(NumParams); -eye(NumParams)];
for p = 1:NumParams
    Offsets(:,p) = Offsets(:,p) * InitStep(p);
end
NumOfs = 2*NumParams; % Number of offsets

%% Loop thru iterations
CurErr = ErrFcn(CurParams); % Current error level
CountIter = 0; % Count iterations
CountNoBetterFit = 0; % Count times no better fit obtained 
while CountIter < NumMaxIter && CountNoBetterFit < NumShrinks
    NewParams = CurParams + Offsets; % Parameters with offsets
    NewErr = Inf(NumOfs,1); % Errors per offset
    % Loop thru offsets
    for i = 1:NumOfs
        if sum(NewParams(i,OnlyPos) <= 0) == 0
            NewErr(i,:) = ErrFcn(NewParams(i,:));
        end
    end
    % Check fit
    NewMinErr = min(NewErr);
    if NewMinErr < CurErr
        % New best error, so move current parameters
        CurParams = NewParams(find(NewErr == NewMinErr,1),:); % Update current parameters
        CurErr = NewMinErr; % Update current error
    else
        % Current error still best, so shrink offsets        
        Offsets = Offsets / 2; % Shrink offsets by half
        CountNoBetterFit = CountNoBetterFit + 1; % Count that no better fit found       
    end   
    CountIter = CountIter + 1; % Count iterations
end

disp(CountIter);

%% Finalise
FinalParams = CurParams;
FinalErr = CurErr;