%% ************************************************************************
% *************************************************************************
% Reconstruction algorithm for joint DDE calibration and imaging in RI:
% Algorithm based on a block-coordinate forward-backward algorithm
% Related paper: Dabbech et al. 2021 arXiv:2102.00065
% authors: Arwa Dabbech, Audrey Repetti,Rick Perley, Oleg Smirnow & Yves
% Wiaux
% *************************************************************************
clear; clc; close all;
%% General paths
fprintf('\nSetting project path .. ')
DirProject = pwd;
DirProject = [DirProject,filesep];% main DIR

%% General params
RunID      = 10;   % Job id
srcName    = '3c391'; % source name

% Imaging
imPixelSize= 2.5; % cellsize in asec
imDimx     = 512; % image dim 1
imDimy     = 512; % image dim 2
imLambda     = 2e-6 ; % image regularisation parameter
imInitThres  = 5e-2 ; % hard threshold to be applied to the initialised image in order to keep true signal only

% Calibration
ddeSpacialDim   = 3  ; % square root of the DDEs spatial dimension: odd integer is required !!  imaging only: 0|  die: 1 | dde: 3,5,7 ..
ddeTemporalDim  = 16  ; % DDEs temporal raduction ratio
ddeAmplitudeBounds   = 0.01 ; % radius of the l_inf ball on the die and dde components.

% ms specific params & flags
param.ms_dataSetsNameTag={'ch1'}; %  name tags of the multiple  data sets
%
param.DirProject = DirProject; % optional: default: current folder
%
param.nFacetsPerDim = [2 2]; % optional: num of facets in the image along each dim: 
%default: will consider facets of minimum dim 256x256 and using available num of cores
%total number of facets should <= number of matlab workers
%
param.ms_scrFieldId = 2; %optional : field id of the source: default 0
param.ms_freqId =1;      %optional: observation freq ID (if wideband data):default 1
%
param.flag_calibration =1; %optional: default 1
param.flag_reweighting =1;  %optional: default 1
param.flag_load_previous_results = 0;  %optional: resume from a previous run: default 0
param.path_previous_results=[];  %optional: default empty
param.flag_save_tmp_results =1 ; %optional: save the latest itr results of BCFB :default 1
%% info
fprintf("\n________________________________________________________________\n")
fprintf('\nJoint calibration and imaging of %s ',srcName)
fprintf('\nData sets name tags: ')
disp(param.ms_dataSetsNameTag)
fprintf('\nINFO: pixel size %f asec',imPixelSize)
fprintf('\nINFO: image dimensions %d x %d',imDimy,imDimx)
fprintf("\n________________________________________________________________\n")

%% main function
RI_JointCalibImaging(srcName,imPixelSize,imDimx,imDimy, imLambda,imInitThres,...
    ddeSpacialDim,ddeTemporalDim,ddeAmplitudeBounds,param,RunID);
