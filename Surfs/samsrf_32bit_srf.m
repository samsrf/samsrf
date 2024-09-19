function Srf = samsrf_32bit_srf(Srf)
%
% Srf = samsrf_32bit_srf(Srf)
%
% Converts all the 64 bit data (doubles) in Srf into 32 bit data (singles). 
% 32 bit data saves space & processing time during model fitting functions. 
% Note, however, that fitting results will not be identical to using 64 bit.
%
% As of SamSrf v9.6 using 32 bit is default but it can be turned back to 
% 64 bit by setting def_64bit in SamSrf_defaults to true. This function
% checks the default setting. If it is set to 64 bit, it doesn't convert 
% but gives a warning message about that.
%
% This function is called -automatically- by the functions creating Srf files
% (samsrf_gii2srf, etc.) as well as by samsrf_expand_srf. However, you can 
% force the latter to override this if you want to read in 64 bit data.
%
% 29/06/2023 - Written (DSS)
% 03/10/2023 - Removed overly verbose defaults messages (DSS)
%

SamSrfDefs = LoadSamSrfDefaults;
if exist('SamSrfDefs.defs_64bit', 'var') && SamSrfDefs.defs_64bit
    % Using 64 bit data (longer processing time but can avoid some errors!)
else
    % Using 32 bit data (default & saves disc space but could introduce errors!)
    fn = fieldnames(Srf);
    
    % Loop thru fields
    for i = 1:length(fn)
        if isa(Srf.(fn{i}), 'double')
            Srf.(fn{i}) = single(Srf.(fn{i}));
        end
    end
end
    