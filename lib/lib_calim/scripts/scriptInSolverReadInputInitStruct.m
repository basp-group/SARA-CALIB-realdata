%% Joint Calib and Imaging: reading input vars
Nxo = param_im.Nt(1); % dim in Fourier
Nyo = param_im.Nt(2); % dim in Fourier
imDimFourier = prod(param_im.Nt);
%
%
nTimeSlots = zeros(nDataSets,1);
nAntennas  = zeros(nDataSets,1);
ddeSpatialDim     = zeros(nDataSets,1);
ddeTemporalDim    = zeros(nDataSets,1);
dieSpatialFreq    = zeros(nDataSets,1);
dieTempralFreq    = zeros(nDataSets,1);
nufftKernelDim    = zeros(nDataSets,1);
GmatRowSZ  = zeros(nDataSets,1);
timeSlotStartRow  = cell(nDataSets,1);
timeSlotEndRow    = cell(nDataSets,1);
for idataSet = 1:nDataSets
    
    nAntennas(idataSet)    = param_data{idataSet}.nAntennas;% nber of antennas
    nTimeSlots(idataSet)   = param_data{idataSet}.nTimeSlots;%nber of time slots
    ddeSpatialDim(idataSet)   = param_dde{idataSet}.ddeSpatialDim; %dde spatial dim
    ddeTemporalDim(idataSet)  = param_dde{idataSet}.ddeTemporalDim;% dde temporal dim
    
    nufftKernelDim(idataSet) = param_data{idataSet}.nufftKernelDim;%nufft kernel dim
    GmatRowSZ(idataSet)      = param_data{idataSet}.G_row_size; % dim of the de-gridding matrix rows
    dieSpatialFreq(idataSet) = (ddeSpatialDim(idataSet)*ddeSpatialDim(idataSet)+1)*0.5; % dde kernel: zero spatial freq
    dieTempralFreq(idataSet) = (ddeTemporalDim(idataSet)+1)/2;%dde kernel: zero temporal freq
    %
    DATA{idataSet}         = param_data{idataSet}.y;
    param_data{idataSet}   = rmfield(param_data{idataSet},'y');
    %
    Antennas{idataSet}     = param_data{idataSet}.antennas;
    param_data{idataSet}   = rmfield(param_data{idataSet},'antennas');
    %
    SupD2{idataSet}        = param_dde{idataSet}.Support_Dalpha2;
    param_dde{idataSet}    = rmfield(param_dde{idataSet},'Support_Dalpha2');
    %
    SupD1{idataSet}        = param_dde{idataSet}.Support_Dalpha1;
    param_dde{idataSet}    = rmfield(param_dde{idataSet},'Support_Dalpha1');
    %
    AntennasIndex{idataSet} = param_data{idataSet}.unique_ants;
    param_data{idataSet}    = rmfield(param_data{idataSet},'unique_ants');
    %
    timeSlotStartRow{idataSet} = (param_data{idataSet}.timeSlotStart);
    timeSlotEndRow{idataSet}   = (param_data{idataSet}.timeSlotEnd);
    param_data{idataSet}    = rmfield(param_data{idataSet},'timeSlotStart');
    param_data{idataSet}    = rmfield(param_data{idataSet},'timeSlotEnd');
    %
    Y{idataSet}             = param_data{idataSet}.YY;
    param_data{idataSet}    = rmfield(param_data{idataSet},'YY');
    %
    Flag{idataSet} = param_data{idataSet}.flag;
    param_data{idataSet} = rmfield(param_data{idataSet},'flag');
end

%% -----------------------------------------------------------------------
%% Joint Calibration and Imaging: Structures, params involved in calib and imaging
%% structs calib
for idataSet=1:nDataSets
    InactiveTimeSlots  = zeros(1,nTimeSlots(idataSet));
    for itimeSlot = 1:nTimeSlots(idataSet)
        InactiveTimeSlots(itimeSlot)= itimeSlot;
        if nnz(Flag{idataSet}{itimeSlot})
            InactiveTimeSlots(itimeSlot)= 0;
        end
    end
    InactiveTimeSlots = nonzeros(InactiveTimeSlots);
    %
    gen_param(idataSet).Nxo = Nxo; % Dim 1 of the oversampled Fourier domain
    gen_param(idataSet).Nyo = Nyo; % Dim 2 of the oversampled Fourier domain
    gen_param(idataSet).nMeas =  param_data{idataSet}.nMeas; % Nber of measurements
    gen_param(idataSet).ddeSpatialDim   = ddeSpatialDim(idataSet) ; % DDE: dim 1 in the spatial Fourier domain
    gen_param(idataSet).nufftKernelDim  = nufftKernelDim(idataSet);% support size of the NUFFT kernels
    gen_param(idataSet).GmatRowSZ  = GmatRowSZ(idataSet);% support size of each row of the de-griddinng matrix
    gen_param(idataSet).nTimeSlots = nTimeSlots(idataSet); % Nber of time slots
    gen_param(idataSet).nAntennas  = nAntennas(idataSet); %Nber of active antennas
    gen_param(idataSet).AntennasIndex = AntennasIndex{idataSet}; %Antennasa ID
    gen_param(idataSet).InactiveTimeSlot = InactiveTimeSlots;% flagged time slots
    gen_param(idataSet).nMeasPerTimeSlot = timeSlotEndRow{idataSet} - timeSlotStartRow{idataSet} +1; % Nber of mesurements per time slot
    gen_param(idataSet).timeSlotStart   = timeSlotStartRow{idataSet}; % indices of the rows marking the time slots
    %
    param_calib(idataSet).nAntennas      = nAntennas(idataSet);
    param_calib(idataSet).nuo            = param_dde{idataSet}.nu;
    param_calib(idataSet).JUo            = param_dde{idataSet}.JU1;
    param_calib(idataSet).JUtot          = param_dde{idataSet}.JUtot;
    param_calib(idataSet).tol_norm       = param_dde{idataSet}.tol_norm;
    param_calib(idataSet).dde_theta_maxo = param_dde{idataSet}.dde_theta_max;
    param_calib(idataSet).die_theta_maxo = param_dde{idataSet}.die_theta_max;
    param_calib(idataSet).dieSpatialFreq  = dieSpatialFreq(idataSet);
    param_calib(idataSet).dieTempralFreq  = dieTempralFreq(idataSet);
    param_calib(idataSet).ddeSpatialDim   = ddeSpatialDim(idataSet);
    param_calib(idataSet).ddeTemporalDim  = ddeTemporalDim(idataSet);
    param_calib(idataSet).nTimeSlots      = nTimeSlots(idataSet);
    param_calib(idataSet).InactiveTimeSlot = InactiveTimeSlots;
    param_calib(idataSet).flag_time_smoothness = param_dde{idataSet}.flag_time_smoothness;
    % DDEs: Temporal Fourier operators
    if param_dde{idataSet}.flag_time_smoothness
        param_calib(idataSet).BWTime =  @(U) (param_dde{idataSet}.GwTime*param_dde{idataSet}.ATime((flipud(U)))); %backward operator: inverse temporal Fourier transform
        param_calib(idataSet).FWTime =  @(D) flipud((param_dde{idataSet}.AtTime(param_dde{idataSet}.GwTime'*(D))));% forward operator: temporal Fourier transform
    end
end

% clear vars
clear timeSlotEndRow  timeSlotStartRow  InactiveTimeSlot;

%% init U1 & U2
if flag_calib
fprintf('\nInitialisation of the DDEs ..')
end
for idataSet=1:nDataSets
    try  eval(param_dde{idataSet}.init_random);
        maxDDEU = 1e-6;
    catch,  param_dde{idataSet}.init_random=0;
    end
    U1{idataSet}=cell(nAntennas(idataSet),1); %#ok<SAGROW>
    UDim = [ddeTemporalDim(idataSet),ddeSpatialDim(idataSet)*ddeSpatialDim(idataSet)];

    if param_dde{idataSet}.flag_time_smoothness
        % temporal smoothness
        if  param_dde{idataSet}.init_random
            for iAnt=1:nAntennas(idataSet)
                U1{idataSet}{iAnt}= maxDDEU*(randn(UDim) + 1i*randn(UDim));
                U1{idataSet}{iAnt}(dieTempralFreq(idataSet),dieSpatialFreq(idataSet)) = 1 + U1{idataSet}{iAnt}(dieTempralFreq,dieSpatialFreq(idataSet));
            end
        else
            for iAnt=1:nAntennas(idataSet)
                U1{idataSet}{iAnt}=zeros(UDim);
                U1{idataSet}{iAnt}(dieTempralFreq(idataSet),dieSpatialFreq(idataSet)) =1.0 ;
            end
        end
    else
        % No temporal smoothness
        if  param_dde{idataSet}.init_random
            for iAnt=1:nAntennas(idataSet)
                U1{idataSet}{iAnt}= maxDDEU*(randn(UDim) + 1i*randn(UDim));
                U1{idataSet}{iAnt}(:,dieSpatialFreq(idataSet))=1.0+U1{idataSet}{iAnt}(:,dieSpatialFreq(idataSet));
            end
        else
            for iAnt=1:nAntennas(idataSet)
                U1{idataSet}{iAnt}=zeros(UDim);
                U1{idataSet}{iAnt}(:,dieSpatialFreq(idataSet)) =1.0 ;
            end
        end
    end
    if flag_calib
        fprintf('\nINFO: Data set %d: DDEs spatial dimension: %d x %d, DDEs temporal dimension %d  ',...
        idataSet,ddeSpatialDim(idataSet),ddeSpatialDim(idataSet),UDim(1))
    end
end
U2  = U1;

%% imaging params
imFacets   = param_algo.imFacets;  %facets positions
imReweightFloorLevel = param_im.imReweightBounds;    % floor level in the weighting scheme

% Proximal step params
param_l1.LowerBound =0; %positivity
param_l1.verbose    = 1 ;
param_l1.max_itr    = 200;  % max number of iterations
param_l1.xA         = param_im.xA;
param_l1.nFacets    = nFacets;              % number of facets
param_l1.facetPsi   = param_algo.facetPsi;  % Psit
param_l1.facetPsit  = param_algo.facetPsit; % Psi
param_l1.imFacets   = param_algo.imFacets;  % facets positions
param_l1.floorLevel = imReweightFloorLevel; % floor level in reweighting scheme

% image cycle params
param_fb_im.imJmin    = param_algo.imJmin;
param_fb_im.imJmax    = param_algo.imJmax; % maximum number of itr
param_fb_im.tol_im   = param_algo.tol_im; % lower bound on the relative variation
param_fb_im.lambda  = param_algo.lambda; %regularisation parameter
param_fb_im.nFacets = nFacets;    % number of facets
param_fb_im.imFacets  =imFacets;  % facets positions
param_fb_im.param_l1= param_l1;   % parameters of the proximal step

% usefull function
SNR = @(x, xtrue) 20 * log10( sqrt( sum( abs(xtrue(:)).^2 ) / sum( abs(x(:)-xtrue(:)).^2 ) ) );
%%  files
if param_algo.save
    FileNameSlice= ['S',num2str(max(ddeSpatialDim)),'x',num2str(max(ddeSpatialDim)),'_ImReg',num2str(param_algo.lambda)];
    ResFileInit    = [param_im.path,'InitResidual0_',FileNameSlice,'.fits'];
    ImFileInit      = [param_im.path,'InitSol_',FileNameSlice,'.fits'];
    ImFile       = [param_im.path,'tmpSol_',FileNameSlice,'.fits'];
    ResFile      = [param_im.path,'tmpResidual_',FileNameSlice,'.fits'];
    ResultsFile  = [param_im.path,'tmpResults_',FileNameSlice,'.mat'];
end
