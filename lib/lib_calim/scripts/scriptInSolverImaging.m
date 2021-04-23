%% Update Image - script
if nit_tot ==1
    %% weights initialisation
    imWeights = cell(nFacets,1); 
    imReweightItr  = 1 ;
    for iFacet =1:nFacets, imWeights{iFacet} = 1;
    end    
    imReweightItrs =0;
else
    %% Reweighinting 
    imWeightsCmpst =Composite();
    if doReweight && flag_reweighting
        fprintf('\n \t!!!! Itr:%d  Reweighting %d !!!!! \n',nit_tot,imReweightItr)
        
        imReweightItrs(imReweightItr+1)=  nit_tot;
          imReweightItr   =  1+ imReweightItr;
        Psit  = param_algo.facetPsit;
        spmd(nFacets)
            if labindex==1
                for iFacet = 2 : nFacets
                    labSend(MODEL(imFacets(iFacet).XI(1):imFacets(iFacet).XI(2),imFacets(iFacet).XJ(1):imFacets(iFacet).XJ(2)),iFacet,iFacet);
                end
                MODEL_DIST = MODEL(imFacets(1).XI(1):imFacets(1).XI(2),imFacets(1).XJ(1):imFacets(1).XJ(2));
            elseif labindex<=nFacets, MODEL_DIST =  labReceive(1,labindex);
            end
        end
        
        PsitModel = Psit(MODEL_DIST);
        spmd(nFacets)
            imWeightsCmpst = (imReweightFloorLevel{labindex}./(imReweightFloorLevel{labindex} + abs(PsitModel)));
        end
        clear PsitModel MODEL_DIST;
        
        imWeights =cell(1,nFacets);
        for iFacet =1:nFacets
            imWeights{iFacet} = imWeightsCmpst{iFacet};
        end
        RESULTS.imWeights = imWeights;
    else
        spmd(nFacets)
            imWeightsCmpst = imWeights{labindex};
        end
    end    
  
    %% Imaging
    
    [MODEL,ImagingOutput]   = imUpdateImage(DATA,MODEL,FWOp_tmp,BWOp_tmp,imWeightsCmpst,param_fb_im);
    
    %% RESULTS
    RESULTS.imReweightItr  = imReweightItr-1;
    RESULTS.imWeights = imWeights;
    RESULTS.imReweightItrs = imReweightItrs; 
    RESULTS.logPrior(nit_tot)=ImagingOutput.imLogPrior ;
    
    
end


