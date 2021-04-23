function [MODEL,nrmL1] = solver_prox_L1_faceted_dict(Y , lambda, param,weightsCmpst,imFacets)
% PROJ_L1 - Proximal operator with L1 norm
% sol = solver_prox_L1_full_image(x, lambda, param) solves:
%   min_{z} 0.5*||x - z||_2^2 + lambda * ||Psit (xA + z)||_1
% References:
% [1] M.J. Fadili and J-L. Starck, "Monotone operator splitting for
% optimization problems in sparse recovery" , IEEE ICIP, Cairo,
% Egypt, 2009.
% [2] Amir Beck and Marc Teboulle, "A Fast Iterative Shrinkage-Thresholding
% Algorithm for Linear Inverse Problems",  SIAM Journal on Imaging Sciences
% 2 (2009), no. 1, 183--202.
%% Optional input arguments -- >  to be updated
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'nu'), param.nu = 1; end
if ~isfield(param, 'rel_obj'), param.rel_obj = 1e-4; end
if ~isfield(param, 'max_iter'), param.max_iter = 100; end
if ~isfield(param, 'weights'), param.weights = 1; end
if ~isfield(param, 'xA'), param.xA = zeros(param.Nyo,param.Nxo); end
if ~isfield(param, 'nFacets'), param.nFacets = 1; end
%%  Initializations
% algo params
nu   = param.nu;
if nu~=1, error("check value of nu");end
seuil = lambda ;
relative_obj = param.rel_obj;
max_iter     = param.max_iter ;
% sparsity ops
Psit  = param.facetPsit;
Psi   = param.facetPsi;
% num of facets
nFacets     = param.nFacets;
% positivity
try  LowerBound   = param.LowerBound;
catch, LowerBound = 0;
end
% image initialisation
try imInit  = full(param.imInit);
catch,    imInit = full(param.xA);
end
% log prior param (nominator of the weights)
try     rho_ = max(1e-8,param.floorLevel);
catch , rho_ = ones(nFacets,1)*1e-8;
end
% clear  memory
clear param
%
itrL1 = 1;obj=-1;
%% Useful functions
soft_star  = @(z, T) z - (sign(z).*max(abs(z)-T, 0));
%% init vars
tau          = 2/3*0.9;
%%
spmd(nFacets)
    if labindex==nFacets
        imCmpst = imInit;
        
        for iFacet = 1 : nFacets -1
            sliceD = imCmpst(imFacets(iFacet).XI(1):imFacets(iFacet).XI(2),imFacets(iFacet).XJ(1):imFacets(iFacet).XJ(2));
            labSend(sliceD,iFacet,iFacet);
        end
        imFacetedCmpst = imCmpst(imFacets(nFacets).XI(1):imFacets(nFacets).XI(2),imFacets(nFacets).XJ(1):imFacets(nFacets).XJ(2));
        
        YCmpst = Y;
        
    else
        imFacetedCmpst =  labReceive(nFacets,labindex);
    end
end

dualL1Cmpst  = Psit(imFacetedCmpst);
PsitImFacetedCmpst  = dualL1Cmpst;
clear  imFacetCmpst imInit Y sliceD;
%%
while 1
    
    PsiDualCmpst = Psi(dualL1Cmpst);
    
    %%    update dual
    spmd(nFacets)
  
        if labindex==nFacets
            % send old facet slices
            %update primal
            imCmpst  = max((1-tau)*imCmpst + tau* (YCmpst - PsiDualCmpst) , LowerBound); % FW
            %  send new facet slices
            for iFacet = 1 : nFacets -1
                labSend(imCmpst(imFacets(iFacet).XI(1):imFacets(iFacet).XI(2),imFacets(iFacet).XJ(1):imFacets(iFacet).XJ(2)) , iFacet, iFacet);
            end
            
            % facet slices
            imFacetedCmpst     = imCmpst(imFacets(nFacets).XI(1):imFacets(nFacets).XI(2),imFacets(nFacets).XJ(1):imFacets(nFacets).XJ(2));            
            dualL1Cmpst = dualL1Cmpst - PsitImFacetedCmpst;
            
        else
            % pre-update dual
            dualL1Cmpst = dualL1Cmpst - PsitImFacetedCmpst;
            % facet slices
            imFacetedCmpst =  labReceive(nFacets,labindex);
        end
        
    end
    
    PsitImFacetedCmpst = Psit(imFacetedCmpst);
    
    clear  XcCmpst;
    
    spmd(nFacets)
        
        if labindex ==nFacets
            nrmL2 =  0.5 *(sum(sum((imCmpst -YCmpst).^2)));
        end
        dualL1Cmpst = soft_star(dualL1Cmpst + 2 * PsitImFacetedCmpst ,  seuil.*weightsCmpst) ;
        dummyL1 = sum(abs(weightsCmpst.*PsitImFacetedCmpst));
        dummyL1 = gplus(dummyL1,nFacets);
      
    end
    
    %% crit
    
    nrmL1.L1   =  dummyL1{nFacets};
    prev_obj =  obj(itrL1);
    itrL1    =  itrL1 + 1;
    obj(itrL1) = nrmL2{nFacets} + lambda * nrmL1.L1 ;
    rel_obj = abs(obj(itrL1)-prev_obj)/obj(itrL1);
    
    % Stopping criterion
    if (rel_obj < relative_obj) || itrL1 >= max_iter
        % fprintf('\nProx Iter %i, prox_fval = %e, rel_fval = %e\n',iter_L1-1, obj(end), rel_obj);
        break;
    end
    
    
end

spmd(nFacets)
    dummyLog      =  rho_(labindex)*sum(log(rho_(labindex) + abs(PsitImFacetedCmpst)));
    dummyLog      = gplus(dummyLog,nFacets);
end
nrmL1.Log = dummyLog{nFacets};

MODEL = imCmpst{nFacets};

clear   im*  *Cmpst Xprev d*;
end


