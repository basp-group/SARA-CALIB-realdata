function [MODEL,output] = imUpdateImage(DATA, MODEL,FW,BW,weightsCmpst,param)
%% Init / prox params
Gamma     = 1.98/ (param.MeasOpNorm); % factor 2 because of the criterion ||y-Phix || ^2 step verified paper Chouzenoux et al. 2016
Lambda    = Gamma * param.eta;

imJmax  = param.imJmax;
imJmin   = param.imJmin;
tol_im  = param.tol_im;

iter_x = 1;
Dist   = 1;
param_prox = param.param_l1;
param_prox.imInit  = sparse(MODEL);

%%  imaging Cycle :Forward-backward steps
DirtyIm = BW(DATA);
while  iter_x <=imJmax
    tst      = tic;
    prevEstimate    = MODEL;
    %Gradient step
    PProx = MODEL - Gamma*(BW(FW(MODEL))-DirtyIm);
   
    %Proximal step
    [MODEL,imL1]  = solver_prox_L1_faceted_dict(PProx,Lambda,param_prox,weightsCmpst,param.imFacets) ;
    iter_x   = iter_x +1;

    Dist = sqrt(sum(sum((MODEL-prevEstimate ).^2)))/sqrt(sum(sum(MODEL.^2))) ;
    tend = toc(tst); 
    if ~mod(iter_x,10)
        fprintf("\nIter %d, dist %f ",iter_x,Dist)
    end
    % ---------------------------------------------------------------------
    % Stopping criterion
    % ---------------------------------------------------------------------
    if ( iter_x > imJmin && ( Dist < tol_im ) )
        fprintf('\nCrit: %f Num Itr: %d \n',Dist,iter_x)
        break;
    end
    % ---------------------------------------------------------------------
end
clear PProx  param_prox weightsCmpst;
%% imaging output
ModelVis = FW(MODEL);
nL2DirtyIm     = sqrt(sum(sum(DirtyIm.^2))); %  l2 norm of the dirty image
nL2ResidualIm  =  sqrt(sum(sum(abs(DirtyIm - BW(ModelVis)).^2)));%  l2 norm of the residual image
clear  DirtyIm;

if iscell(DATA)
    visL2Res = 0; 
    visL2Data = 0;
    for i =1:numel(DATA)
        visL2Res = visL2Res +sum(sum(abs(DATA{i} - ModelVis{i}).^2));
        visL2Data = visL2Data + sum(sum(abs(DATA{i} ).^2));
    end
    clear  ModelVis;
else 
    visL2Res = norm(DATA-FW(MODEL)).^2;
    visL2Data = norm(abs(DATA)).^2;
end
clear DATA;

output.visL2Res  = visL2Res;
output.visSNR    = 10 * log10(visL2Data /output.visL2Res );
output.imSNR     = 20 * log10(nL2DirtyIm/nL2ResidualIm);
output.imL1Prior       = imL1.L1;
output.imLogPrior      = imL1.Log;

