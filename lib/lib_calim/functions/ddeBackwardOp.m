function [DDE_FourierTIME]=ddeBackwardOp(DDE_FourierSPACE,H,param)
nTimeSlots      = param.nTimeSlots;
ddeSpatialDim1D = param.ddeSpatialDim^2;
ddeTemporalDim  = param.ddeTemporalDim;
 
%% ----------------------------------------------------
ProjectedDDE      = zeros(nTimeSlots,ddeSpatialDim1D);

for itimeSlot = 1:nTimeSlots
    if ~isempty(H{itimeSlot})
        ProjectedDDE(itimeSlot,:) =  H{itimeSlot}'*(DDE_FourierSPACE{itimeSlot});
    end
end

if param.flag_time_smoothness
    DDE_FourierTIME = zeros(ddeTemporalDim,ddeSpatialDim1D);
    for ddePIXEL = 1 : ddeSpatialDim1D
        DDE_FourierTIME(:,ddePIXEL) =  param.FWTime(ProjectedDDE(:,ddePIXEL));
    end
else
    DDE_FourierTIME =   ProjectedDDE;
end

end
