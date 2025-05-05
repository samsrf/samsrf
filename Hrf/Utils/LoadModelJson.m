function Model = LoadModelJson(Filename)
%
% Model = LoadModelJson(Filename)
%
% Reads the model specification in the JSON file <Filename>.json into the struct Model.
% This is used by the SamSrf GUI but you could also use this to save your model
% without relying on binary Matlab format.
%
% Internally converts Matlab strings back into character arrays, and also
% converts the string containing the functions into function handles.
% But you need not worry about this.
%
% Note: This function only works from Matlab R2023b onwards! 
%       (Although you earlier versions allowed you to save XML files 
%        so you could changes it to that instead)
%
% 08/09/2024 - Written (DSS)
% 11/09/2024 - Added note about backwards compatibility (DSS)
%

% Read as a file
Model = readstruct([Filename '.json']);

% Convert function back from char array
if isfield(Model, 'Prf_Function') && ~isscalar(Model.Prf_Function)
    Model.Prf_Function = str2func(Model.Prf_Function);
end

% Convert seed-parameter function back from char array
if isfield(Model, 'SeedPar_Function') 
    Model.SeedPar_Function = str2func(Model.SeedPar_Function);
end

% Convert parameter names back from strings
if isfield(Model, 'Param_Names')
    Model.Param_Names = cellstr(Model.Param_Names);
end

% Convert any strings into character arrays
Fs = fieldnames(Model);
for i = 1:length(Fs)
    if isstring(Model.(Fs{i}))
        Model.(Fs{i}) = char(Model.(Fs{i}));
    end
end