function samsrf_error(str)
% Displays the error message in str - if SamSrfAnalysis was run this is 
% in the GUI otherwise it is in the Matlab command window.

samsrf_disp(['ERROR: 'srf]);
error('SamSrf Fatal Error!');