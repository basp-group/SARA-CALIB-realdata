%--------------------------------------------------------------------------

fprintf('\nCompute residual''s normalisation factor  ..')
dirac = zeros(imDims);
dirac(imDims(1)/2+1,imDims(2)/2+1) = 1;
psf = BWOp_tmp(FWOp_tmp(dirac));
clear dirac;
RESULTS.maxpsf(nit_tot) = max(max(psf));
clear  psf;

%--------------------------------------------------------------------------
fprintf('\nCompute objective function & stats of the residuals in image and data spaces  ..')


if nit_tot ==1
    if param_algo.save
        fitswrite(single(MODEL),ImFileInit);
        fitswrite(single(ResidualIm)./RESULTS.maxpsf(nit_tot) ,ResFileInit);
    end
    clear ResidualIm  ;
    critDiff_inner = inf;
    RESULTS.snrEIm(nit_tot)  = RESULTS.snr0Im; %image fidelity in image space
    RESULTS.snrEVis(nit_tot) = RESULTS.snr0Vis;
    
else
    
    RESULTS.currItr = nit_tot;
    %residuals in both spaces
    ResidualVis = cell(nDataSets,1);
    ModelVis = FWOp_tmp(MODEL);
    for idataSet =1 :nDataSets
        ResidualVis{idataSet} =  DATA{idataSet} - ModelVis{idataSet}  ;
    end
    RESULTS.ResidualVis = cell2mat(ResidualVis); % residual vis.
    RESULTS.ResidualIm  = BWOp_tmp(ResidualVis)./RESULTS.maxpsf(nit_tot); % residual image
    clear ResidualVis;
    
    %  dde prior val
    RegDDE(nit_tot)  = 0;
    for idataSet=1:nDataSets
        for iAnt = 1:nAntennas(idataSet)
            RegDDE(nit_tot) = RegDDE(nit_tot) + ...
                param_dde{idataSet}.nu*0.5*sum(sum(abs(U1{idataSet}{iAnt}-U2{idataSet}{iAnt}).^2));
        end
    end 
    % fidelity to data term
    DataFidelity(nit_tot)  = ImagingOutput.visL2Res;
  
    %  image prior val
    RegImage(nit_tot)  = ImagingOutput.imL1Prior * param_algo.eta;
    
    % Inner objective function
    CRIT_INNER(nit_tot) = DataFidelity(nit_tot) + RegDDE(nit_tot) + RegImage(nit_tot) ;
    %Outer objective function
    CRIT_OUTER_TMP = DataFidelity(nit_tot) + RegDDE(nit_tot) + param_algo.eta*RESULTS.logPrior(nit_tot);
    
    % relative variation on the different terms

    critDiff_inner = (((CRIT_INNER(nit_tot-1) - CRIT_INNER(nit_tot))/CRIT_INNER(nit_tot-1) ));
    dataFid_rel_var = (DataFidelity(nit_tot-1)-DataFidelity(nit_tot))/(DataFidelity(nit_tot-1));
    regDDE_rel_var = (RegDDE(nit_tot-1)-RegDDE(nit_tot))/(RegDDE(nit_tot-1));
    regIm_rel_var = (RegImage(nit_tot-1)-RegImage(nit_tot))/(RegImage(nit_tot-1));
    
    % add results and stats to structure
    RESULTS.snrEIm(nit_tot)  = ImagingOutput.imSNR; %image fidelity in image space
    RESULTS.snrEVis(nit_tot) = ImagingOutput.visSNR;%image fidelity in vis. space
    clear ImagingOutput;
    
    % stats of residuals (STD, SKEWNESS, KURTOSIS)
    RESULTS.ResImSTD(nit_tot)  = std((RESULTS.ResidualIm(:)));
    RESULTS.ResImskewness(nit_tot)  = skewness((RESULTS.ResidualIm(:)));
    RESULTS.ResImkurtosis(nit_tot)  = kurtosis((RESULTS.ResidualIm(:)));
    RESULTS.ResVisSTD(nit_tot)  = std([imag(RESULTS.ResidualVis ) ;real(RESULTS.ResidualVis )]);
    RESULTS.ResVisKurtosis(nit_tot) = kurtosis([imag(RESULTS.ResidualVis ) ;real(RESULTS.ResidualVis )]);
    RESULTS.ResVisSkewness(nit_tot) = skewness([imag(RESULTS.ResidualVis ) ;real(RESULTS.ResidualVis )]);
    
    % save fits files
    if param_algo.save
        fitswrite(single(MODEL),ImFile);
        fitswrite(single(RESULTS.ResidualIm),ResFile);
    end     
end



RESULTS.U1      = U1;
RESULTS.U2      = U2;
RESULTS.MODEL = MODEL;
RESULTS.MeasOpNorm = MeasOpNorm;
RESULTS.runTime = runTime;


