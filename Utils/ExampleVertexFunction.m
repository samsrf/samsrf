function ExampleVertexFunction
%
% Adapt this example script to change the default vertex selection function in samsrf_surf.
% This works on the current figure handle, so your surface map must be the active figure.
%
% Our suggestion is you make a copy of this script somewhere on your Matlab path and change 
% the name to something descriptive and adapt it according to your needs.
%
% 11/09/2021 - Written (DSS)
% 20/04/2022 - SamSrf 8 version (DSS)

rotate3d; % Ensure rotation mode since MatLab 2020 figure defaults to auto-selection
dcm_obj = datacursormode(gcf); % Data cursor mode object
set(dcm_obj, 'UpdateFcn', @gvf); % Change selection function

%% Vertex selection function 
% Update the various parts below to meet your needs & wishes
function txt = gvf(empt, event_obj)
    % Import vertex coordinates in mesh space 
    global Vertices 

    % Selected position in mesh space coordinates
    pos = get(event_obj, 'Position'); 
    
    % Selected vertex number
    v = find(Vertices(:,1) == pos(1) & Vertices(:,2) == pos(2) & Vertices(:,3) == pos(3), 1); 

    % Change this string to whatever you want the datatip to say
    txt = ''; 

    % Simple example of something you could do 
    disp(['Selected vertex: #' num2str(v)]);
end

end