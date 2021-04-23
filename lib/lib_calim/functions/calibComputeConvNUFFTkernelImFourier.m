function [ConvMat]= calibComputeConvNUFFTkernelImFourier(imFourierVal,imFourierIndices,GMAT,param)
%% Init
ddeSpatialDim1D = param.ddeSpatialDim^2;
nufftKernelDim  = param.nufftKernelDim;
nTimeSlots  = param.nTimeSlots;
Nxo  = param.Nxo;
Nyo  = param.Nyo;
%% G mat slices
GVals = GMAT.Gmat_values; 
GPos  = GMAT.Gmat_pos;
clear GMAT;
%% parpool constant
imFourierValPool = parallel.pool.Constant(imFourierVal); 
imFourierIndicesPool = parallel.pool.Constant(imFourierIndices); 
clear imFourierVal  imFourierIndices;%  save memory
%% computing matrices
tinit  = tic;
ConvMat = cell(nTimeSlots,1);
parfor itimeSlot =1:nTimeSlots
    LM      = index2d(1:ddeSpatialDim1D,sqrt(ddeSpatialDim1D));
    %% get Gridding kernels
    GValsT = GVals{itimeSlot};GVals{itimeSlot}=[];
    A12L   = floor(numel(GValsT)/nufftKernelDim);
    GPosT  = index2d(GPos{itimeSlot},Nxo);
    GPos{itimeSlot}  = [];%save memory
    GValsT = reshape(GValsT,[nufftKernelDim,A12L]);
    %% compute mats
    indicesB     = sub2ind([Nxo,Nyo],(col(GPosT(:,1)-LM(:,1).')+LM(:,1).'),(col(GPosT(:,2)-LM(:,2).')+LM(:,2).'));
    GPosT  = [];%save memory
    [~,indicesB] = ismember(indicesB,imFourierIndicesPool.Value);
    indicesB     = reshape(indicesB,nufftKernelDim,A12L,ddeSpatialDim1D,ddeSpatialDim1D) ;
    ConvMat{itimeSlot}= reshape((sum(GValsT.*imFourierValPool.Value(indicesB),1)),A12L,ddeSpatialDim1D,ddeSpatialDim1D);
    indicesB = []; GValsT = [];%save memory
end
clear  G*   imFourier*  ;% clear unnecessary vars to save memory
tEnd = floor(toc(tinit));
fprintf(" .. %ds  \n",tEnd);
end
