function ApName = CreateApertures

% 20/04/2022 - SamSrf 8 version (DSS)

[f,p] = uigetfile('*.mat', 'Select results files', 'MultiSelect', 'on');
if isnumeric(f) && f == 0
    ApName = 0;
    return;
end
if ischar(f)
    f = {f};
end
ApFrm = [];

ApName = f{1};
usc = strfind(ApName, '_');
ApName = ApName(usc(1)+1:usc(2)-1);

nvols = 0;
for i = 1:length(f)
    Res = load([p f{i}]);
    if ~isfield(Res, 'ApFrm')
        samsrf_error('File contains no apertures!');
    end
    ApFrm(:,:,nvols+1:nvols+size(Res.ApFrm,3)) = Res.ApFrm;
    nvols = nvols + size(Res.ApFrm,3);
end

ApName = ['aps_' ApName];
save(ApName, 'ApFrm');