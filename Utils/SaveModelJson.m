function SaveModelJson(Model, Filename)
%
% SaveModelJson(Model, Filename)
%
% Saves the model specification in struct Model as the JSON file <Filename>.json.
% This can be read back in using ReadModelJson.
%
% Internally converts function handles into character strings but you need
% not worry about this as it will be converted back when reading in.
%
% Note: This function only works from Matlab R2023b onwards! 
%       (Although you earlier versions allowed you to save XML files 
%        so you could changes it to that instead)
%
% 08/09/2024 - Written (DSS)
% 11/09/2024 - Added note about backwards compatibility (DSS)
%

% Convert function to char array
if isfield(Model, 'Prf_Function') && strcmpi(class(Model.Prf_Function),'function_handle')
    Model.Prf_Function = func2str(Model.Prf_Function);
end

% Convert seed-parameter function to char array
if isfield(Model, 'SeedPar_Function') 
    Model.SeedPar_Function = func2str(Model.SeedPar_Function);
end

% Now save as a file
writestruct(Model, [Filename '.json']);