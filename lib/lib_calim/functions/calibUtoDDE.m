function DDE = calibUtoDDE(U,param)
nTimeSlots         = param.nTimeSlots;
ddeSpatialDim1D     = param.ddeSpatialDim^2;


DDE   = zeros(nTimeSlots,ddeSpatialDim1D);
for colDDE  =1:ddeSpatialDim1D
    DDE(:,colDDE)= param.BWTime(U(:,colDDE));
end
end