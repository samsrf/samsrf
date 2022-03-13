function AnonymiseLabels(Sure)
% 
% AnonymiseLabels([Sure=false])
%
% Loads all label files in a folder and anonymises them by redacting the
% anatomical coordinates stored in it (setting coordinates to zero). 
% They may not be compatible with FreeSurfer anymore (I haven't tried)...
%
% WARNING: Overwrites the label file so only do this if you have a copy &
%          are sure that you want to wipe the anatomical info!
%
% The optional input Sure toggles whether it asks if you're sure.
%
% 13/03/2022 - Written (DSS)

if nargin == 0
    Sure = false;
end

% List label files
Labels = dir('*.label');
Labels = {Labels.name}';
disp(Labels);

% Are we sure?
if ~Sure
    disp('This will wipe anatomical data from all these labels.');
    q = input('Are you sure you want to continue? (Y/N) ', 's');
    if upper(q) == 'Y'
        Sure = true;
    end
end

% We are sure
if Sure
    % Loop thru label files
    for i = 1:length(Labels)
        % Load old label file
        Data = Read_FreeSurfer(Labels{i});
        nver = size(Data,1); % Number of vertices

        % Save new label file
        fid = fopen(Labels{i}, 'w');
        fprintf(fid, '#! Converted from SamSurfer (Anonymised)\n');
        fprintf(fid, '%d\n', nver);
        for vi = 1:nver
            v = Data(vi,1);
            if size(Data,2) == 5
                fprintf(fid, '%d %5.3f %5.3f %5.3f %f\n', Data(vi,1), 0, 0, 0, Data(vi,5));
            elseif size(Data,2) == 6
                fprintf(fid, '%d %5.3f %5.3f %5.3f %f %f\n', Data(vi,1), 0, 0, 0, Data(vi,5), Data(vi,6));
            end
        end
        fclose(fid);
        disp(['Saved ' Labels{i}]);
    end
end