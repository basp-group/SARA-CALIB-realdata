function [GMAT,imFourierIndices] =calibComputeFourierIndicesSliceG(G0,SupD1,SupD2,Flag,param)
%% compute indices of FX and get coefficients of G involved in convolutions
Nxo          = param.Nxo;
Nyo          = param.Nyo;
nTimeSlots   = param.nTimeSlots;
nufftKernelDim  = param.nufftKernelDim;
ddeSpatialDim = param.ddeSpatialDim;
timeSlotStart    = param.timeSlotStart ;
%%
[Gsch_pos,Gsch_v1n,Gsch_v3n]=find((G0).'); clear G0;
Gsch_pos     = index1d(shift_ind(index2d(Gsch_pos,Nxo),Nxo,Nyo),Nxo);
snap_nszs    = 1;
indices1d    = cell(nTimeSlots,1);
Gmat_values  = cell(nTimeSlots,1);
Gmat_pos     = cell(nTimeSlots,1);
%% 
tstart1  = tic;
for itimeSlot = 1:nTimeSlots
    Gmat_values{itimeSlot} = [];
    Gmat_pos{itimeSlot}    = [];
    indices1d{itimeSlot}   = [];
    if nnz(Flag{itimeSlot})
        snap_nsze      = nufftKernelDim*(nnz(vertcat(Flag{1:itimeSlot})));
        if itimeSlot>1 
            snap_nszs = nufftKernelDim*nnz(vertcat(Flag{1:itimeSlot-1}))+1;
        end
        Gmat_pos{itimeSlot}    =  Gsch_pos(snap_nszs:snap_nsze);
        Gmat_values{itimeSlot} =  Gsch_v3n(snap_nszs:snap_nsze);
        indices1d{itimeSlot}   =  timeSlotStart(itimeSlot)-1 + ...
            vertcat(vertcat(SupD2{:,itimeSlot}),vertcat(SupD1{:,itimeSlot}));
    end
end
clear  Gsch_v3n;
[~,indices1d] = ismember(Gsch_v1n,vertcat(indices1d{:}));
indices1d     = unique(Gsch_pos((indices1d>0)));
indices2d     = index2d(indices1d,Nxo);
clear Gsch_v1n Gsch_pos indices1d;%save memory

imFourierIndices = util_findIdxConv(indices2d,Nxo,Nyo,ddeSpatialDim);

tend1 = floor(toc(tstart1));
fprintf('  %d s ',tend1);
%%
GMAT.Gmat_values = Gmat_values; clear Gmat_values 
GMAT.Gmat_pos    = Gmat_pos;    clear Gmat_pos 
end
