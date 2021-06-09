function [SrfL, SrfR] = samsrf_hemi_srfs(Srf)
%
% [SrfL, SrfR] = samsrf_hemi_srfs(Srf)
%
% Returns two Srf structures by separating the data from a bilateral Srf 
%  into its component hemispheres again. This is useful when you created
%  a bilateral Srf and deleted the hemisphere files to save space but then
%  realise you need them again for some reason...
%
%  Srf is automatically expanded. Returns SrfL & SrfR, respectively. 
%  These are not compressed or processed in any other way. Also, if you want 
%  to save them to a file you need to rename each as Srf before saving.
%
%
% Warning: May fail with very large data files due to lack of memory.
%
% 09/06/2021 - Written (DSS)
%

% Expand Srf
Srf = samsrf_expand_srf(Srf);

% Loop thru hemispheres
hem = {'lh' 'rh'};
for h = 1:2
    % Indeces for hemispheres
    if h == 1
        % Left hemisphere
        V = 1:Srf.Nvert_Lhem;
        F = sum(Srf.Faces <= Srf.Nvert_Lhem, 2) == 3;
    else
        % Right hemisphere
        V = Srf.Nvert_Lhem+1:size(Srf.Vertices,1);
        F = sum(Srf.Faces > Srf.Nvert_Lhem, 2) == 3;
    end    
    
    % Split hemispheres
    Cur = Srf;
    Cur.Hemisphere = hem{h};
    Cur = rmfield(Cur, 'Nvert_Lhem');
    % Subtract left hemisphere vertex indeces if right hemisphere
    Cur.Faces = Srf.Faces(F,:) - Srf.Nvert_Lhem*(h==2); 

    Cur.Vertices = Srf.Vertices(V,:);
    Cur.Normals = Srf.Normals(V,:);

    if isfield(Cur, 'Pial')
        Cur.Pial = Srf.Pial(V,:);
        Cur.Inflated = Srf.Inflated(V,:);
        Cur.Sphere = Srf.Sphere(V,:);
    end

    if isfield(Cur, 'Curvature')
        Cur.Curvature = Srf.Curvature(1,V); 
    else
        warning('Curvature data not in Srf...');
    end
    if isfield(Cur, 'Area')
        Cur.Area = Srf.Area(1,V); 
        Cur.Thickness = Srf.Thickness(1,V); 
    else
        warning('Cortical dimensions not in Srf...');
    end
    Cur.Data = Srf.Data(:,V); 

    if isfield(Cur, 'Raw_Data')
        Cur.Raw_Data = Srf.Raw_Data(:,V); 
    end
    if isfield(Cur, 'X')
        Cur.X = Srf.X(:,V); 
    end
    if isfield(Cur, 'Y')
        Cur.Y = Srf.Y(:,V); 
    end
    if isfield(Cur, 'Rmaps')
        Cur.Rmaps = Srf.Rmaps(:,V); 
    end
    if isfield(Cur, 'ConFlds')
        Cur.ConFlds = Srf.ConFlds(:,V); 
    end
    if isfield(Cur, 'Noise_Ceiling')
        Cur.Noise_Ceiling = Srf.Noise_Ceiling(:,V); 
    end

    % Remove mesh field if it exists
    if isfield(Cur, 'Meshes')
        Cur = rmfield(Cur, 'Meshes');
    end
    
    % Store in relevant structure
    if h == 1
        SrfL = Cur;
    else
        SrfR = Cur;
    end
end

% Shift hemispheres
if isfield(Cur, 'Sphere')
    SrfL.Pial(:,1) = SrfL.Pial(:,1)+3;
    SrfR.Pial(:,1) = SrfR.Pial(:,1)-3;
    SrfL.Inflated(:,1) = SrfL.Inflated(:,1)+50;
    SrfR.Inflated(:,1) = SrfR.Inflated(:,1)-50;
    SrfL.Sphere(:,1) = SrfL.Sphere(:,1)+105;
    SrfR.Sphere(:,1) = SrfR.Sphere(:,1)-105;
else
    warning('Anatomical surfaces not in Srf...');
end
