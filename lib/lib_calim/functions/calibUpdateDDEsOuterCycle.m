function [U1,U2,DiffDDEs] = calibUpdateDDEsOuterCycle(U1,U2,DATA,...
    ConvImFourier,MCv,MCvc,A1S20,A2S10, param)

ddeSpatialDim1D  = param.ddeSpatialDim^2;
nTimeSlots = param.nTimeSlots;
nAntennas  = param.nAntennas;
JUtot   = param.JUtot;

%% Init
ErrDATA2      = zeros(nAntennas,1);
DiffDDEAnt2   = zeros(nAntennas,1);
DiffDDEAnt1   = zeros(nAntennas,1);
DiffDDEs = inf ;
sizeDDE   = size(U1{1});
%% Calib Cycle
tStart = tic;
nbit_d        = 0;
while (DiffDDEs > param.tol_norm && nbit_d < 3 ) || ( nbit_d <= JUtot )
    nbit_d = nbit_d +1;
    %% U1  Cycle
    %compute DDE from U2
    SpatialDdeKernels = zeros(nAntennas,nTimeSlots,ddeSpatialDim1D);
    H = cell(nAntennas,1);
    for iAnt =1 :nAntennas
        if param.flag_time_smoothness
            SpatialDdeKernels(iAnt,:,:) = (calibUtoDDE(U2{iAnt},param));
        else, SpatialDdeKernels(iAnt,:,:) = U2{iAnt};
        end
        H{iAnt} = cell(nTimeSlots,1);
    end
    SpatialDdeKernels = permute(SpatialDdeKernels,[1 3 2]);
    % compute H operator
    for itimeSlot =1:nTimeSlots
        cDDESlice   = SpatialDdeKernels(:,:,itimeSlot);%reshape(CDDE(:,idSnap,:),[AntsNum,DSZ2]);
        cDDESliceFC = conj(fliplr(cDDESlice));
        convMatSlice = ConvImFourier{itimeSlot};
        for iAnt = 1:nAntennas
            Ant1Supp2 = A1S20{iAnt,itimeSlot}; A1S2L  = length(Ant1Supp2); HSlice1Conj = [];
            Ant2Supp1 = A2S10{iAnt,itimeSlot}; A2S1L  = length(Ant2Supp1); HSlice1 = [];
            if A1S2L, HSlice1Conj = conj(reshape(sum(cDDESlice(Ant1Supp2,:).* convMatSlice(MCv{iAnt,itimeSlot},:,:),2),A1S2L,ddeSpatialDim1D));  end
            if A2S1L, HSlice1 = fliplr(reshape(sum(cDDESliceFC(Ant2Supp1,:).* convMatSlice(MCvc{iAnt,itimeSlot},:,:),2),A2S1L,ddeSpatialDim1D));   end
            H{iAnt}{itimeSlot} = [HSlice1;HSlice1Conj];
        end
    end
    % FB update U1
    parfor iAnt =  1:nAntennas
        DDEold= U1{iAnt};
        param_fb   = param; %init
        DataSlice  = cell(1,nTimeSlots);
        for itimeSlot = 1:nTimeSlots
            DataSlice{itimeSlot} = (nonzeros(DATA(iAnt,:,itimeSlot)));
        end
        %-
        param_fb.Lips = op_norm_single_op(@(U) ddeForwardBackwardOp(U,H{iAnt},param), sizeDDE,1e-5,1000,0);
        BCst   = ddeBackwardOp(DataSlice,H{iAnt},param) + (U2{iAnt})*param_fb.nuo;
        GradOp = @(U)  ddeForwardBackwardOp(U,H{iAnt},param)+ param_fb.nuo*(U);
        %-
        U1{iAnt} = calibDDEsFwBwInnerCycle(U1{iAnt},BCst,GradOp,param_fb);
        %-
        DiffDDEAnt1(iAnt) = norm(col(abs(U1{iAnt} - DDEold)))./norm(col(abs(U1{iAnt})));
    end
    clear H ;
    %%
    %% U2  Cycle
    %% GET DATA, DDE & PARAMS% configure DATA for U2 update
    SpatialDdeKernels = zeros(nAntennas,nTimeSlots,ddeSpatialDim1D);
    %compute H operator
    H = cell(nAntennas,1);
    for iAnt  = 1 : nAntennas
        if param.flag_time_smoothness
            SpatialDdeKernels(iAnt,:,:) = (calibUtoDDE(U1{iAnt},param));
        else, SpatialDdeKernels(iAnt,:,:) = U1{iAnt};
        end
        H{iAnt} = cell(nTimeSlots,1);
    end
    SpatialDdeKernels = permute(SpatialDdeKernels,[1 3 2]);
    for itimeSlot = 1:nTimeSlots
        convMatSlice = ConvImFourier{itimeSlot};
        cDDESlice   = SpatialDdeKernels(:,:,itimeSlot);
        cDDESliceFC = conj(fliplr(cDDESlice));
        for iAnt = 1:nAntennas
            Ant1Supp2 = A1S20{iAnt,itimeSlot};  A1S2L = length(Ant1Supp2); HSlice2ConjM = [];
            Ant2Supp1 = A2S10{iAnt,itimeSlot};  A2S1L = length(Ant2Supp1); HSlice2M = [];
            if A1S2L, HSlice2ConjM =conj((reshape(sum(cDDESlice(Ant1Supp2,:).* convMatSlice(MCv{iAnt,itimeSlot} ,:,:),2),A1S2L,ddeSpatialDim1D))); end % meas. conj
            if A2S1L, HSlice2M = fliplr(reshape(sum(cDDESliceFC(Ant2Supp1,:).* convMatSlice(MCvc{iAnt,itimeSlot},:,:),2),A2S1L,ddeSpatialDim1D)); end%  meas.
            H{iAnt}{itimeSlot} =[HSlice2M ; HSlice2ConjM];
        end
    end
    % FB update U2
    parfor iAnt = 1:nAntennas
        DDEold = U2{iAnt};
        param_fb   = param;
        DataSlice  =  cell(1,nTimeSlots);
        for itimeSlot = 1:nTimeSlots
            DataSlice{itimeSlot} = (nonzeros(DATA(iAnt,:,itimeSlot)));
        end
        %-
        param_fb.Lips =  op_norm_single_op(@(U) ddeForwardBackwardOp(U,H{iAnt},param),sizeDDE,1e-5,1000,0);
        BCst   = ddeBackwardOp(DataSlice,H{iAnt},param) + (U1{iAnt})*param_fb.nuo;
        GradOp = @(U)  ddeForwardBackwardOp(U,H{iAnt},param)+ param_fb.nuo*(U);
        %-
        U2{iAnt}     = calibDDEsFwBwInnerCycle(U2{iAnt},BCst,GradOp,param_fb);
        %-
        DiffDDEAnt2(iAnt) = norm(col(abs(U2{iAnt} - DDEold)))./norm(col(abs(U2{iAnt})));
        DataSlice    = vertcat(DataSlice{:}); 
        DATASliceNl2 = norm(DataSlice);
        DATAModelU2  = ddeForwardOp(U2{iAnt},H{iAnt},param);
        ErrDATA2(iAnt) = sqrt(sum(abs(vertcat(DATAModelU2{:})-(DataSlice)).^2))./DATASliceNl2;
    end
    %
    DiffDDEs =  max(max(DiffDDEAnt2,DiffDDEAnt1));
    clear H ;
end
%%
clear   ConvImFourier DATA param  CDDE* ;%SAVE MEMEORY
%%
tempcal2 = toc(tStart);
fprintf('\nmse per antenna:%2.3f, mse: %2.3f ... %ds',max(ErrDATA2)*100,mean(ErrDATA2)*100,floor(tempcal2));
% stopping crit
end