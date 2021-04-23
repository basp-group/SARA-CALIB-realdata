function [DDE_SPACE]=ddeForwardOp(DDE_FourierTIME,H,param)
nTimeSlots      = param.nTimeSlots;
ddeSpatialDim1D    = param.ddeSpatialDim^2;

%% ------------------------------------------------
DDE_SPACE = cell(nTimeSlots,1);   
if param.flag_time_smoothness   
    BWTime= param.BWTime;
    DDE_ = zeros(nTimeSlots,ddeSpatialDim1D);
    for ddePIXEL  = 1:ddeSpatialDim1D
        DDE_(:,ddePIXEL)= BWTime(DDE_FourierTIME(:,ddePIXEL));
    end
else, DDE_ = DDE_FourierTIME;
end
for itimeSlot = 1 : nTimeSlots
    DDE_SPACE{itimeSlot}=[];
    if ~isempty(H{itimeSlot})
        DDE_SPACE{itimeSlot}= H{itimeSlot}*col(DDE_(itimeSlot,:));
    end
end
end
