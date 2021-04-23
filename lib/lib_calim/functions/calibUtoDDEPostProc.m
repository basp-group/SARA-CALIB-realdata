function [DDE1, DDE2] = calibUtoDDEPostProc(U1,U2,param)
nAntennas  = param.nAntennas;
nTimeSlots  = param.nTimeSlots;
ddeSpatialDim1D = param.ddeSpatialDim^2;
%%--------------------------------------------------
DDE1 = zeros(nAntennas,nTimeSlots,ddeSpatialDim1D);
DDE2 = zeros(nAntennas,nTimeSlots,ddeSpatialDim1D);
if param.flag_time_smoothness
    fprintf('\nApplying the temporal Fourier inverse of calibration kernels ..')
end
for AIdx = 1:nAntennas
    if param.flag_time_smoothness
        DDE1(AIdx,:,:) = fliplr(calibUtoDDE(U1{AIdx},param));
        DDE2(AIdx,:,:) = conj(calibUtoDDE(U2{AIdx},param));
    else
        DDE1(AIdx,:,:) = fliplr(U1{AIdx});
        DDE2(AIdx,:,:) = conj(U2{AIdx});
    end
end
end

