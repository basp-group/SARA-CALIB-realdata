function DDE_FourierTIME = ddeForwardBackwardOp(DDE_FourierTIME,H,param)
nTimeSlots      = param.nTimeSlots;
ddeSpatialDim1D    = param.ddeSpatialDim^2;
 
%% ----------------------------------------------------

if param.flag_time_smoothness
    BWTime           = param.BWTime;
    FWTime           = param.FWTime;
end
%% ----------------------------------------------------

DDE_FourierSPACE        = zeros(nTimeSlots,ddeSpatialDim1D);
if  param.flag_time_smoothness
    for ddePIXEL  = 1:ddeSpatialDim1D
        DDE_FourierSPACE(:,ddePIXEL) = BWTime(DDE_FourierTIME(:,ddePIXEL));
    end
else,  DDE_FourierSPACE = DDE_FourierTIME;
end
ProjectedDDE = zeros(nTimeSlots,ddeSpatialDim1D);
for itimeSlot = 1 : nTimeSlots
    if ~isempty(H{itimeSlot})
        ProjectedDDE(itimeSlot,:) = H{itimeSlot}'*( H{itimeSlot}*col(DDE_FourierSPACE(itimeSlot,:)));
    end
end
clear H;
if  param.flag_time_smoothness
    for ddePIXEL = 1 : ddeSpatialDim1D
        DDE_FourierTIME(:,ddePIXEL) = FWTime(ProjectedDDE(:,ddePIXEL));
    end
else, DDE_FourierTIME = ProjectedDDE;
end