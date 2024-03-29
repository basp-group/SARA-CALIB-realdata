function RESULTS = solver_JointCalibImaging_BCFB(param_im, param_dde, param_algo, param_data)
%% ************************************************************************
% *************************************************************************
% Reconstruction algorithm for joint DDE calibration and imaging in RI:
% Algorithm based on a block-coordinate forward-backward algorithm
% Related paper: Dabbech et al. 2021
% *************************************************************************
%% Initialization 0
NumWorkers = feature('numcores');
% initial image model
MODEL = full(param_im.MODEL);
param_im  = rmfield(param_im,'MODEL');
imDims = size(MODEL);

%Number  of data sets
nDataSets = numel(param_data);

%Number of image facets
nFacets = param_algo.nFacets;

% flags
try flag_calib = param_algo.flag_calib;
catch,   flag_calib =1;
end
try flag_reweighting = param_algo.flag_reweighting;
catch,   flag_reweighting =1;
end
% Init, DDES, DATA, ANTENNA, FLAG variables
U1  = cell(nDataSets,1);
U2  = cell(nDataSets,1);
Y   = cell(nDataSets,1);
DATA  = cell(nDataSets,1);
Antennas = cell(nDataSets,1);
SupD1    = cell(nDataSets,1);
SupD2    = cell(nDataSets,1);
Flag     = cell(nDataSets,1);


%--------------------------------------------------------------------------
%% Init calib, imaging, param structures
fprintf('\nSetting up usefull params for calib and imaging, Init Vars and Funs');

scriptInSolverReadInputInitStruct;

%--------------------------------------------------------------------------
%% Compute Crit 0
RegDDE       = zeros(param_algo.max_it,1);
DataFidelity = zeros(param_algo.max_it,1);
RegImage     = zeros(param_algo.max_it,1);
CRIT_INNER   = zeros(param_algo.max_it,1);
% dde prior
for idataSet=1:nDataSets
for iAnt = 1:nAntennas(idataSet)
    RegDDE(1) = RegDDE(1) +...
        param_dde{idataSet}.nu*0.5*sum(sum(abs(U1{idataSet}{iAnt}-U2{idataSet}{iAnt}).^2));
end
end
% Data fidelity
deGriddingMat  = cell(nDataSets,1);
for idataSet =1 :nDataSets,  deGriddingMat{idataSet}=param_data{idataSet}.deGriddingMat;
    param_data{idataSet}.G0 =[];
end

% Measurment operator & its adjoint
FWOp_tmp = @(x) meas_op(deGriddingMat,x, param_im.tf); %forward op.
BWOp_tmp = @(y) adjoint_meas_op(deGriddingMat,y,param_im.tfadj); %backward op.
%Spectral norm of the measurement op.
MeasOpNorm = op_norm(FWOp_tmp,BWOp_tmp, imDims,1e-5,200,0);

% data fidelity term
ModelData =  FWOp_tmp(MODEL);
data_l2norm=0;
for idataSet =1 :nDataSets
    DataFidelity(1) =  DataFidelity(1) +sum(abs(DATA{idataSet}-ModelData{idataSet} ).^2);
    data_l2norm = data_l2norm + sum(abs(DATA{idataSet}).^2);
end

%  fidelity of the initial model image, in both image & data spaces
ResidualIm =BWOp_tmp(DATA)- BWOp_tmp(ModelData);
RESULTS.snr0Im  = SNR((BWOp_tmp(DATA)), (BWOp_tmp(ModelData)));
RESULTS.snr0Vis = 10*log10(data_l2norm./DataFidelity(1) );
clear ModelData  ;

% image regularisation
param_algo.eta  = param_algo.lambda*MeasOpNorm ;
RegImage(1)  =    param_algo.eta  * sum(param_algo.Psitw(MODEL)) ;

% objective function
CRIT_INNER(1) = DataFidelity(1) +  RegDDE(1) + RegImage(1) ;
% outer & inner iterations
nitr_inner = 0;
nitr_outer = 0;
critDiff_outer =1;
CRIT_OUTER = inf;
doReweight = 0;
% imaging params
param_fb_im.MeasOpNorm  = MeasOpNorm;
param_fb_im.eta = param_algo.eta;

%% splitting degridding matrix, computing Fourier indices ..
%--------------------------------------------------------------------------
fprintf('\nSlice Gridding matrix  & compute Fourier indices ')
tStart = tic;
GmatTimeSlices = cell(nDataSets,1);
imFourierIndices = cell(nDataSets,1);
calibMCvS2 = cell(nDataSets,1);
calibMCvS1 = cell(nDataSets,1);
calibA1S2  = cell(nDataSets,1);
calibA2S1  = cell(nDataSets,1);
parpool(NumWorkers);
for idataSet = 1:nDataSets
    fprintf('\n!! Data set %d !! ',idataSet);
    % Computing FX indices-------------------------------------------------
    fprintf('\nComput. indices of the image Fourier transform..')
    fprintf('\nSlice the degridding matrix in time slots ..')
    [GmatTimeSlices{idataSet},imFourierIndices{idataSet}]= ...
        calibComputeFourierIndicesSliceG(deGriddingMat{idataSet},...
        SupD1{idataSet}, SupD2{idataSet},Flag{idataSet},gen_param(idataSet));
    
    deGriddingMat{idataSet} = [];
    % flag re-organise-----------------------------------------------------
    fprintf('\nRe-organising flags');
    for itimeSlot = 1:nTimeSlots(idataSet)
        dummy = Flag{idataSet}{itimeSlot};
        Flag{idataSet}{itimeSlot}(1) = dummy(1);
        for i = 2:length(dummy)
            Flag{idataSet}{itimeSlot}(i) = dummy(i)*sum(dummy(1:i));
        end
    end
    % Computing Conv. indices----------------------------------------------
    if flag_calib
    fprintf('\nGather antenna indices for DDE calibration ');
    end
    [calibMCvS2{idataSet},calibMCvS1{idataSet},calibA1S2{idataSet},calibA2S1{idataSet}] = ...
        calibComputeAntennaIndices(SupD2{idataSet},SupD1{idataSet},...
        Antennas{idataSet},Flag{idataSet}, gen_param(idataSet));
end
clear SupD2 SupD1 *dummy  param_ComputeIndFXGDDEsU;%save memory

tEnd = floor( toc(tStart) );
fprintf('\nDone! in %d seconds\n',tEnd)

%% load results from previous run
startItr =1;
if param_algo.load
    fprintf('\nLoading previous results ..')
    
    if ~isempty( param_algo.results_file  )
        load(param_algo.results_file,'RESULTS')
    else, load(ResultsFile,'RESULTS');
    end
    MODEL = full(RESULTS.MODEL);
    U1      = RESULTS.U1;
    U2      = RESULTS.U2;
    imWeights = RESULTS.imWeights;
    
    try     critDiff_inner  = RESULTS.critDiff_inner;
    catch,  critDiff_inner = 1;
    end
    % stopping crit & iteration indices
    CRIT_INNER  = RESULTS.CRIT_INNER;
    nitr_inner = RESULTS.nitr_inner;
    nitr_outer = RESULTS.nitr_outer;
    CRIT_OUTER   = RESULTS.CRIT_OUTER;
    try imReweightItr =RESULTS.imReweightItr;end
    try imReweightItrs= RESULTS.imReweightItr;end
    try doReweight =RESULTS.doReweight;end
    startItr = RESULTS.currItr +1;
end

%-----------------------------------------------------
%% ALGORITHM
fprintf('\n**************************************\n')
fprintf('********* STARTING ALGORITHM *********')
fprintf('\n**************************************\n')
runTime = 0;


for nit_tot = startItr:param_algo.max_it
    tItrI =tic;
    %% ************
    %  Calibration
    % *************
    if nit_tot > 1 && flag_calib
        tCalib =tic;
        %%************************************
        fprintf(  '\n\n\n* Calibration   *\n ')
        %*************************************
        calibCrit = ones(nDataSets,1);
        imFourier = param_im.tf_shift(MODEL); %%fftshifted
        delete(gcp('nocreate'));
        parpool(NumWorkers);
        for idataSet = 1:nDataSets
            % Compute comvolution matrices between the NUFFT kernels & the
            % Fourier of the image
            fprintf('\n!! Data set %d !!',idataSet);
            fprintf('\nCompute convolutions between NUFFT kernels and the Fourier of image ..\n ')
            ConvMat = calibComputeConvNUFFTkernelImFourier(imFourier(imFourierIndices{idataSet}),imFourierIndices{idataSet},...
                GmatTimeSlices{idataSet},gen_param(idataSet));%
            % Update DDEs: outer cycle
            fprintf('\nUpdate DDEs ..');
            [U1{idataSet},U2{idataSet},calibCrit(idataSet)] = calibUpdateDDEsOuterCycle(...
                U1{idataSet},U2{idataSet},Y{idataSet},ConvMat,calibMCvS2{idataSet},calibMCvS1{idataSet},...
                calibA1S2{idataSet},calibA2S1{idataSet},param_calib(idataSet));%...
            clear ConvMat; %save memory
        end
        tCalibE = floor(toc(tCalib));
        fprintf(" ->Calib done in %d",tCalibE)
    end
    clear imFourier;
    
    %% ***************************
    %  Update measurement operator
    % ****************************
    tStart = tic;
    if  flag_calib && nit_tot>1
        fprintf('\n\n\n* Update measurement operator * \n');
        deGriddingMat =  cell(nDataSets,1);
        % update de-gridding matrices
        for idataSet  = 1:nDataSets
            fprintf('\n!! Data set %d !! ',idataSet);
            [RESULTS.DDE1{idataSet},RESULTS.DDE2{idataSet}] =...
                calibUtoDDEPostProc(U1{idataSet},U2{idataSet},param_calib(idataSet));
            fprintf('\nUpdating de-gridding matrix  for imaging .. ')
            deGriddingMat{idataSet} = calibUpdateDeGriddingMatrix(RESULTS.DDE1{idataSet},...
                RESULTS.DDE2{idataSet},Antennas{idataSet},Flag{idataSet},...
                GmatTimeSlices{idataSet},gen_param(idataSet));
            
        end
        %  Update Operators
        FWOp_tmp = @(x) meas_op(deGriddingMat,x, param_im.tf);
        BWOp_tmp = @(y) adjoint_meas_op(deGriddingMat,y,param_im.tfadj);
        fprintf('\nComputing operator''s spectral norm .. ')
        param_fb_im.MeasOpNorm = op_norm(FWOp_tmp,BWOp_tmp,imDims,1e-6,1000,0);
        
        
        tEnd = floor(toc(tStart))+1;
        fprintf('\nDone! in %d s\n',tEnd)
    end
    %% **************
    %  Update Image
    % ***************
    if nit_tot >1
        fprintf('\n*      IMAGING      * \n')
    end
    tStartIm = tic;
    delete(gcp('nocreate'))
    parpool(nFacets);
    
    scriptInSolverImaging;
    delete(gcp('nocreate'))

    tEndIm = floor(toc(tStartIm))+1;
    fprintf('--> Done! in %ds\n',tEndIm);
    tItrE   = toc(tItrI);
    runTime = runTime+tItrE;
    %% **********************
    % Gather current resutls
    % ***********************
    if nit_tot >1
        fprintf('\n\n*   Gather results of itr. %d   *\n\n',nit_tot)
    end
    scriptInSolverGetResults;
    
    %display numbers
    fprintf("\n_________________________________________________________________\n")
    fprintf('Itr GLOBAL = %d, CRIT = %1.4f, TIME:%1.1f s, RW:%d\n',...
        nit_tot,critDiff_inner,tItrE,imReweightItr-1);
    if nit_tot>1
        fprintf('Relative variation: \nData_fid: %f \nReg_DDE: %f \nReg_Im: %f',...
            dataFid_rel_var,regDDE_rel_var,regIm_rel_var);
    end
    fprintf('\nSNR IM EST: %2.3f \nSNR VIS: %2.3f (Init %2.3f)',...
        RESULTS.snrEIm(end),RESULTS.snrEVis(end),RESULTS.snr0Vis(1))
    fprintf("\n________________________________________________________________\n")
    
    
    % ****************************************************************
    % Global stopping crit.
    if nit_tot >1
        fprintf(  '* Checking stopping criteria * ')
    end
    if (nitr_inner >= param_algo.min_it_crit_inner && critDiff_inner <= param_algo.tol_crit_inner)...
            || (nitr_inner >= param_algo.max_it_crit_inner)
        nitr_inner = 1;
        if flag_reweighting
            fprintf('\n!! Reweighting in the next iteration !!\n');
            doReweight = 1;
        end
        nitr_outer = nitr_outer +1;
        CRIT_OUTER(nitr_outer) =CRIT_OUTER_TMP;
        if nitr_outer>1
            critDiff_outer= abs(((CRIT_OUTER(nitr_outer) - CRIT_OUTER(nitr_outer-1))/CRIT_OUTER(nitr_outer-1) ));
            fprintf('\nCRIT global rel. var. : %f \n', critDiff_outer);
        end
    else
        doReweight = 0;
        nitr_inner = nitr_inner+1;
    end
    if ((nit_tot>2) && abs(critDiff_outer)< param_algo.tol_crit_outer) ||...
            ( nitr_outer> param_algo.max_it_crit_outer)
        fprintf("\nGlobal Stop Crit. reached, it glob %d, critDiff %f\n",nit_tot,critDiff_inner)
        break;
    end
    fprintf("\n________________________________________________________________\n")
    % save current results
    RESULTS.critDiff_inner = critDiff_inner;
    RESULTS.CRIT_INNER    = CRIT_INNER;
    RESULTS.nitr_inner = nitr_inner;
    RESULTS.nitr_outer = nitr_outer;
    RESULTS.CRIT_OUTER   = CRIT_OUTER;
    RESULTS.doReweight   = doReweight;
    if  param_algo.save && nit_tot>1
        fprintf('\n* Saving intermittent results *\n');
        save(ResultsFile,'RESULTS')
    end
    fprintf("\n________________________________________________________________\n")
    % clear memory
    RESULTS.imWeights =[];
    RESULTS.DDE1=[];
    RESULTS.DDE2=[];
    RESULTS.U1=[];
    RESULTS.U2=[];
    RESULTS.MODEL=[];
end
fprintf("\n\nJoint Imaging and Calibration finished in %d\n\n", floor(runTime));
fprintf('\n**************************************\n')
fprintf('********* END OF ALGORITHM *********')
fprintf('\n**************************************\n')
%% Final variables
RESULTS.MODEL_IMAGE = MODEL; %reconstructed image
RESULTS.ESTIMATED_U1 = U1; %ddes
RESULTS.ESTIMATED_U2 = U2; %ddes
end
