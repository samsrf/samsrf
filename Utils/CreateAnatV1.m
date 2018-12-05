function CreateAnatV1(Srf, hemi, crit)
%
% CreateAnatV1(Srf, hemi, crit)
%
% Creates a V1 label from the anatomical prediction based on Hinds method.
% The first input Srf is the surface data needed to create the label.
% The second input hemi defines the hemisphere ('lh' or 'rh'). It can also
% define the direct path to the label folder in addition. The third input
% crit defines the minimum probability value to include.
%

% Split path & hemisphere
[p h] = fileparts(hemi);

% Expand Srf if needed
Srf = samsrf_expand_srf(Srf); 

% Read label data
Label = Read_FreeSurfer([hemi '.v1.prob.label']); 

% Create empty vertex data & fill with probabities
Srf.Data = zeros(1,size(Srf.Data,2)); 
Srf.Data(Label(:,1)+1) = Label(:,5); 

% Indeces of vertices inside ROI
Vx = find(Srf.Data >= crit); 

% Save label
samsrf_srf2label(Srf, [h '_v1_hinds' num2str(crit*100) '%'], 1, Vx); 
