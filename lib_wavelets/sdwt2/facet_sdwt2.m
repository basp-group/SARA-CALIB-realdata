function [SPsitLx]=facet_sdwt2(ImX,param)
% sdwt2
Q =(length(param));
spmd 
    if labindex<=Q
    zerosNum = param(labindex).offsetL + ...
               param(labindex).offsetR + ...
               param(labindex).dims_overlap_ref;    
    x_overlap= zeros([zerosNum(:)',1]);
    x_overlap(param(labindex).offsetL(1)+1:end-param(labindex).offsetR(1),...
              param(labindex).offsetL(2)+1:end-param(labindex).offsetR(2))...
        = ImX;    
    % forward operator [put the following instructions into a parfeval for parallelisation]
    SPsitLx = sdwt21_sara(x_overlap,...
        param(labindex).I,...
        param(labindex).dims,...
        param(labindex).offset, ...
        param(labindex).status,...
        param(labindex).nlevel,...
        param(labindex).wavelet);
    end
end

end

