calib_param_dde =cell(nDataSets,1);
calib_param_data=cell(nDataSets,1);
calib_param_im = [];
calib_param_algo = [];
%% Compute noise map
fprintf('\nBuilding noise map to determine floor level ..')
FW_OP = @(x) meas_op(PreCalib.Gw,x,param_nufft.A);
BW_OP = @(y) adjoint_meas_op(PreCalib.Gw,y,param_nufft.At);
dirac = zeros(ImDims);
dirac(ImDims(1)./2 +1,ImDims(2)./2 +1) =1;
PSF = BW_OP(FW_OP(dirac));
PSFPeak = max(max(PSF));
clear dirac;
noise = cell(nDataSets,1);
for idataSet =1 :nDataSets
    noise{idataSet} = (randn(data{idataSet}.DataDims,1)+1i*(randn(data{idataSet}.DataDims,1)))./sqrt(2);
end

imNoiseMap = BW_OP(noise)./PSFPeak;

clear noise PSF BW_OP FW_OP ;
%% temporal settings for calib
fprintf("\nSetting up temporal params. for calibration ...")
temporalNUFFTKernelDim    = 8;
temporalNUFFTOversampling = 2;
ddesTimeDim = zeros(nDataSets,1);
ddesTimeFourierDim  = zeros(nDataSets,1);
ntimeSlots = zeros(nDataSets,1);
timeFlagVect = cell(nDataSets,1);
flag_time_smoothness = ones(nDataSets,1);

for idataSet=1:nDataSets
    ctimeObsTimes = unique(cell2mat(PreCalib.timeslot{idataSet}.time));  % time slots
    ntimeSlots(idataSet) = numel(ctimeObsTimes);                         % Nber of time slots
    ctimeDiff    = ctimeObsTimes(2:end) -  ctimeObsTimes(1:end-1);
    
    % get integration time: time step:
    ctimeStep = floor(median(nonzeros(ctimeDiff)));
    if median(nonzeros(ctimeDiff))>ctimeStep+0.5
        ctimeStep =ctimeStep+1;
    end
    
    % compute DDE kernel temporal dimension and set NUFFT operators
    ctimeUniformSamples = floor(ctimeObsTimes/ctimeStep);
    ctimeDim = 2*floor((ctimeUniformSamples(end) - ctimeUniformSamples(1))/2);
    ddesTimeFourierDim(idataSet) =  floor(ctimeDim/ddeTemporalDimDataSets(idataSet))+...
        (~mod(floor(ctimeDim/ddeTemporalDimDataSets(idataSet)),2)); %odd size
    
    if ddesTimeFourierDim(idataSet)>= ntimeSlots(idataSet)
        flag_time_smoothness(idataSet) = 0;
        ddesTimeDim(idataSet) = ntimeSlots(idataSet);
        if flag_calibration
            fprintf('\nINFO: no temporal smoothness is adopted for data set %d due to sparse temporal sampling!',idataSet)
        end
        
    else
        ddesTimeDim(idataSet)= ddesTimeFourierDim(idataSet) ;
        ctimeContSamples = ctimeObsTimes/ctimeStep;
        ctimeContSamplesShifted = (ctimeContSamples   - ctimeContSamples(floor(ntimeSlots(idataSet)*0.5)+1))./ctimeStep;
        ctimeContSmaplesShiftedRadians = ctimeContSamplesShifted*pi./(max(abs(ctimeContSamplesShifted)));
        [calib_param_dde{idataSet}.ATime , calib_param_dde{idataSet}.AtTime, calib_param_dde{idataSet}.GwTime,~] = op_p_nufft_time_calib(ctimeContSmaplesShiftedRadians(:), ...
            ddesTimeDim(idataSet),temporalNUFFTKernelDim,temporalNUFFTOversampling*ddesTimeDim(idataSet),(ddesTimeDim(idataSet)-1)*0.5);
        if flag_calibration
        fprintf('\nINFO: Temporal smoothness is adopted for data set %d, DDE temporal dim %d.',idataSet,ddesTimeDim(idataSet))
        end
    end
    clear ctime* ;
    % store flagged time slots
    timeFlagVect{idataSet}  = find(~cell2mat(PreCalib.timeslot{idataSet}.time));
    PreCalib.timeslot{idataSet}.time=[];
end

%% Resturcture data for calib & other vars
nActiveAntennas = zeros(nDataSets,1);
uniqueAntennas = cell(nDataSets,1);
fprintf("\nRe-structuring the data for calibration ...\n")

for idataSet =1:nDataSets
    % get number of antennas, and indices of active antennas
    uniqueAntennas{idataSet} = unique(([unique(vertcat(PreCalib.timeslot{idataSet}.ant2{:})); ...
        unique(vertcat(PreCalib.timeslot{idataSet}.ant1{:}))]));
    nActiveAntennas(idataSet) = length(uniqueAntennas{idataSet});
    % init data and antenna support structures
    Y_Curr  =  zeros(nActiveAntennas(idataSet),nActiveAntennas(idataSet),ntimeSlots(idataSet));
    SupportAntAlpha1Curr =  cell(nActiveAntennas(idataSet),ntimeSlots(idataSet));
    SupportAntAlpha2Curr =  cell(nActiveAntennas(idataSet),ntimeSlots(idataSet));
    %
    for itimeSlot = 1:ntimeSlots(idataSet)
        if ~ismember(itimeSlot,timeFlagVect{idataSet})
            timeSlotDataCurr = zeros(nActiveAntennas(idataSet));
            na2Curr = PreCalib.timeslot{idataSet}.ant2{itimeSlot};
            na1Curr = PreCalib.timeslot{idataSet}.ant1{itimeSlot};
            for iAnt=1:nActiveAntennas(idataSet)
                iAntPos  = 0;
                try  iAntPos = PreCalib.timeslot{idataSet}.involved_ants_idx{itimeSlot}(uniqueAntennas{idataSet}(iAnt));
                end
                if iAntPos
                    alpha_id = uniqueAntennas{idataSet}(iAnt);
                    index1id   = find(na1Curr==alpha_id);
                    index2id   = find(na2Curr==alpha_id);
                    
                    SupportAntAlpha1Curr{iAnt,itimeSlot}    = flipud(index1id);
                    SupportAntAlpha2Curr{iAnt,itimeSlot}    = flipud(index2id);
                    
                    if ~isempty(index1id)
                        [~,ants2Curr]   = ismember(na2Curr(index1id),uniqueAntennas{idataSet});
                        timeSlotDataCurr(iAnt,ants2Curr) = data{idataSet}.y{itimeSlot}(index1id);
                    end
                    clear  index2id index1id alpha_id;
                end
            end
            if PreCalib.timeslot{idataSet}.sizes(itimeSlot)~= nnz(timeSlotDataCurr)
                msg='FATAL ERROR! Check the data structure';
                error(msg)
            end
            timeSlotDataCurr = timeSlotDataCurr+conj(transpose(timeSlotDataCurr));
            if nnz(timeSlotDataCurr)
                for iAnt = 1:nActiveAntennas(idataSet)
                    iAntPos =  PreCalib.timeslot{idataSet}.involved_ants_idx{itimeSlot}(uniqueAntennas{idataSet}(iAnt));
                    if iAntPos
                        Y_Curr(iAnt,:,itimeSlot) = (transpose(timeSlotDataCurr(iAnt, :))) ;
                    end
                end
            end
            clear timeSlotDataCurr;
        end
        
        % data: antennas for each time slot
        calib_param_data{idataSet}.antennas.ant1{itimeSlot} =  int32(PreCalib.timeslot{idataSet}.ant1{itimeSlot}); PreCalib.timeslot{idataSet}.ant1{itimeSlot}=[];
        calib_param_data{idataSet}.antennas.ant2{itimeSlot} =  int32(PreCalib.timeslot{idataSet}.ant2{itimeSlot});PreCalib.timeslot{idataSet}.ant2{itimeSlot}=[];
        
    end
    
    %calib: antenna supports 
    calib_param_dde{idataSet}.Support_Dalpha1 = SupportAntAlpha1Curr; % Antenna 1 support
    calib_param_dde{idataSet}.Support_Dalpha2 = SupportAntAlpha2Curr; % Antenna 2 support
    
    calib_param_data{idataSet}.YY =Y_Curr; %data: structured for calib
    Y_Curr =[];
     
    calib_param_data{idataSet}.deGriddingMat  = PreCalib.Gw{idataSet}; %data :  degridding matrix
    PreCalib.Gw{idataSet}=[];
    
    calib_param_data{idataSet}.y   = cell2mat(data{idataSet}.y);    %data :structured for imaging 
    data{idataSet}.y =[];
    
    calib_param_data{idataSet}.flag = (data{idataSet}.flag);  % flags
    data{idataSet}.flag =[];
    
    clear *Curr;
end


%% Init params
fprintf('\nSetting up the different structures of imaging and calibration parameters ..')
%Data & specs
for idataSet = 1:nDataSets    
    %% DDEs - calib - params 
    calib_param_dde{idataSet}.ddeSpatialDim    = ddeSpacialDimDataSets(idataSet); % ddes: kernel spatial size along 1 dim.
    calib_param_dde{idataSet}.ddeTemporalDim    = ddesTimeDim(idataSet);%ddesTimeFourierDim(idataSet); % ddes: temporal size (odd)
    calib_param_dde{idataSet}.dde_theta_max     = ddeAmplitudeBounds; % ddes: l_inf radius on the dde components 
    calib_param_dde{idataSet}.die_theta_max     = ddeAmplitudeBounds; % ddes: l_inf radius on the dde components 
    calib_param_dde{idataSet}.flag_time_smoothness = flag_time_smoothness(idataSet); % smoothness in time flag
    calib_param_dde{idataSet}.JUtot             = 5;%FW-BW ddes updates: max number iter in outer cycle
    calib_param_dde{idataSet}.JU1               = 5;%FW-BW ddes updates:  max number iter in inner cycle (U1)
    calib_param_dde{idataSet}.JU2               = 5;%FW-BW ddes updates:  max number iter in inner cycle (U2)
    calib_param_dde{idataSet}.tol_norm          = 1e-4;%FW-BW ddes updates: relative variation lower bound
    calib_param_dde{idataSet}.nu   = data{idataSet}.MeasOpNorm/nActiveAntennas(idataSet); %FW-BW ddes updates: regularisation param
    calib_param_dde{idataSet}.init_random = 1;
    calib_param_dde{idataSet}.init = 0;
    
    %% data - measurement operator
    calib_param_data{idataSet}.nufftKernelDim = param_nufft.Ky*param_nufft.Kx; % NUFFT: kernel support
    calib_param_data{idataSet}.nMeas = numel( calib_param_data{idataSet}.y );  % Nber of meas.
    calib_param_data{idataSet}.G_row_size = (2*ddeSpacialDimDataSets(idataSet)+param_nufft.Ky-2)^2; % deGridding matrix : row support size
    calib_param_data{idataSet}.nTimeSlots      = ntimeSlots(idataSet);  % Nber of time slots
    calib_param_data{idataSet}.nAntennas      = nActiveAntennas(idataSet); %Nber of active antennas
    calib_param_data{idataSet}.timeSlotStart    = PreCalib.timeslot{idataSet}.start; %Index of first data point
    calib_param_data{idataSet}.timeSlotEnd      = PreCalib.timeslot{idataSet}.end; % Index of last data point 
    calib_param_data{idataSet}.unique_ants      = uniqueAntennas{idataSet}; %Antennas ids
    PreCalib.timeslot{idataSet}.end=[];    PreCalib.timeslot{idataSet}.start =[];   
end

%% algorithmic parameters
param_algo.results_file =[];
param_algo.max_it   = 200;
param_algo.imJmax    = 150; % imaging cycle: max nber of itr 
param_algo.imJmin    = 10;  % imaging cycle: max nber of itr 
param_algo.tol_im    = 5e-5;% imaging cycle: lower bound on the rel var of the image

if flag_reweighting
    param_algo.tol_crit_inner = 5e-3;  % inner algo: lower bound on the rel var of the obj. funct
    param_algo.max_it_crit_inner = 10; % inner algo: max nber of iter
    param_algo.min_it_crit_inner = 3;  % inner algo: min nber of iter
    param_algo.tol_crit_outer = 1.5e-3;  % outer algo: lower bound on the rel var of the obj. funct.
    param_algo.max_it_crit_outer = 40; % outer algo: max nber of weighted pbs to be solved.
else %Solving the l1 problem
    param_algo.tol_crit_inner = 3e-4;  % inner algo: lower bound on the rel var of the obj. funct
    param_algo.max_it_crit_inner = param_algo.max_it ; % inner algo: max nber of iter
    param_algo.min_it_crit_inner = 3;  % inner algo: min nber of iter
    param_algo.tol_crit_outer = 1.5e-3;  % outer algo: lower bound on the rel var of the obj. funct.
    param_algo.max_it_crit_outer = 0;  %outer algo: max nber of weighted pbs to be solved.
end
param_algo.lambda   = imLambda;
param_algo.load     = flag_load_previous_results;
param_algo.flag_calib  = flag_calibration;
param_algo.flag_reweighting = flag_reweighting;
if flag_load_previous_results
    fprintf('\nINFO: loading the results of a previous run');
    param_algo.results_file = path_previous_results;
end

param_algo.save = flag_save_tmp_results;
if flag_save_tmp_results
    fprintf('\nINFO: will be saving results of the current itr of BCFB ..\n');
end 

%% imaging parameters

param_imaging.NumWorkers = NumWorkers; %matlab available workers
param_imaging.thres     = imInitThres; % hard threshold level to be applied on  init. image
param_imaging.Ni        = [imDimy imDimx]; %image dim.
param_imaging.Nt        = [param_nufft.oy param_nufft.ox].* [imDimy imDimx];
param_imaging.tf        = param_nufft.A;
param_imaging.tf_shift  = @(X)so_fft2_shift(X,param_imaging.Nt,param_nufft.Scale);
param_imaging.tfadj     = param_nufft.At;
param_imaging.MODEL     = (imInitModel); clear imInitModel;
param_imaging.xA        = sparse(imDimy, imDimx);


param_imaging.nFacetsPerDim  = param_dict.nFacetsPerDim; % number of facets along each dimension
% facet params
delete(gcp('nocreate'))
[param_algo.Psiw, param_algo.Psitw] = op_sp_wlt_basis(param_dict.basis, param_dict.nlevel, imDimy, imDimx);
param_algo.nFacets = prod(param_imaging.nFacetsPerDim); %number of facets
parpool(param_algo.nFacets);
% Build faceted wavelet operators
[paramPsit,paramPsi,param_algo.imFacets]=facetInitParams(...
    param_imaging.Ni, param_dict.nlevel,param_imaging.nFacetsPerDim,param_dict.basis);
param_algo.facetPsit = @(x) facet_sdwt2(x,paramPsit);
param_algo.facetPsi = @(x)  facet_isdwt2(x,paramPsi);

% Get floor level in wavelet domain from noise map
imFacets = param_algo.imFacets;
Psit  = param_algo.facetPsit;
spmd(param_algo.nFacets)
    if labindex==1
        for iFacet = 2 : param_algo.nFacets
            labSend(imNoiseMap(imFacets(iFacet).XI(1):imFacets(iFacet).XI(2),...
                imFacets(iFacet).XJ(1):imFacets(iFacet).XJ(2)),iFacet,iFacet);
        end
        FacetsCmpst = imNoiseMap(imFacets(1).XI(1):imFacets(1).XI(2),...
            imFacets(1).XJ(1):imFacets(1).XJ(2));
    elseif labindex<= param_algo.nFacets
        FacetsCmpst =  labReceive(1,labindex);
    end
end
PsitE = Psit(FacetsCmpst);
spmd(param_algo.nFacets)
    imReweightBounds= std(PsitE(:));
end
for iFacet =1:param_algo.nFacets
    param_imaging.imReweightBounds{iFacet}  = imReweightBounds{iFacet};
end
clear  FacetsCmpst  PsitE imNoiseMap imFacets imReweightBounds Psit;

%% clear unnecessary vars
clear PreCalib data  temporal* tmp uv* time* nFacets  InAntAll  dir* ;
clear param_dict  uniqueAntennas flag*  ms_specific_param  *OP;

%% path files
if param_algo.flag_calib
    param_imaging.path=[DirResults,'/Calim_S',num2str(ddeSpacialDim),'_ID',num2str(RunID),...
    '_temporalRatio',num2str(ddeTemporalRatio),...
    '_regL1',num2str(imLambda),...
    '_thres',num2str(param_imaging.thres(1)),'/'];
else
    param_imaging.path=[DirResults,'/Imaging_ID',num2str(RunID),...
    '_regL1',num2str(imLambda),...
    '_thres',num2str(param_imaging.thres(1)),'/'];
    
end

mkdir(param_imaging.path)



%%
try
    delete(gcp('nocreate'))
end
