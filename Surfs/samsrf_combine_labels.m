function Status = samsrf_combine_labels(OldLabels, NewLabel)
%
% Status = samsrf_combine_labels(OldLabels, NewLabel)
%
% Combines the labels in the cell array OldLabels into a new label named NewLabel.
% Returns the Status of the operation: true = no error, false = error
%
% Note: This function assumes labels are from the same hemisphere! You can use
%       samsrf_bilat_label to create labels combining left & right hemispheres.
%
% 20/04/2022 - SamSrf 8 version (DSS)
%

Vs = [];
for i = 1:length(OldLabels)
    CurrLabel = Read_FreeSurfer([OldLabels{i} '.label']);
    if isnan(CurrLabel)
        Vs = NaN;
    else
        Vs = [Vs; CurrLabel];
    end
end
nver = size(Vs,1);

if sum(isnan(Vs(:,1))) == 0
    fid = fopen([NewLabel '.label'], 'w');
    fprintf(fid, '#! Converted from SamSurfer.\n');
    fprintf(fid, '%d\n', nver);
    for v = 1:nver
        fprintf(fid, '%d %5.3f %5.3f %5.3f %f\n', Vs(v,1), Vs(v,2), Vs(v,3), Vs(v,4), Vs(v,5));
    end
    fclose(fid);
    disp(['Saved ' NewLabel '.label']);
    Status = true;
else
    Status = false;
end