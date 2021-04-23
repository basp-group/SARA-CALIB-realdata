function [param_facetIn, param_facetOut]=facet_init_params(N,nlevel,QDim,wavelet)
Qy =  QDim(1);% 4;    % number of facets along the dimension y
Qx =  QDim(2);%8;    % number of facets along the dimension y
Q = Qx*Qy; % total number of facets
rg_y = domain_decomposition(Qy, N(1));
rg_x = domain_decomposition(Qx, N(2));
segDims = zeros(Q, 4);
for qx = 1:Qx
    for qy = 1:Qy
        q = (qx-1)*Qy+qy;
        segDims(q, :) = [rg_y(qy, 1)-1, rg_x(qx, 1)-1, rg_y(qy,2)-rg_y(qy,1)+1, rg_x(qx,2)-rg_x(qx,1)+1];
    end
end
param_facet.I = segDims(:, 1:2);    % starting index of each facet
param_facet.dims = segDims(:, 3:4); % facet dimensions
param_facet.Q=Q;
param_facet.wavelet =wavelet;
param_facet.QDim=QDim;
param_facet.nlevel = nlevel;
param_facet.N=N;

% Wavelet parameters
n = 1:(length(wavelet)-1); % index of daubechies wavelet
% dwtmode('zpd','nodisp');
param_facet.L = [2*n,0].'; % filter length

[param_facet.I_overlap_ref_nc, param_facet.dims_overlap_ref_nc,...
    param_facet.I_overlap_ref,param_facet.dims_overlap_ref, ...
    param_facet.I_overlap, param_facet.dims_overlap, ...
    param_facet.I_overlap_nc, param_facet.dims_overlap_nc, ...
    param_facet.status, param_facet.offset, ~, ~, ...
    param_facet.offsetL, param_facet.offsetR] = generate_segdwt_indices(N,...
    param_facet.I, param_facet.dims, nlevel, wavelet, param_facet.L);

xsol=zeros(N);
SPsitLx =cell(Q,1);
for q = 1:Q
    offsetLbis = param_facet.offsetL(q,:);
    offsetRbis = param_facet.offsetR(q,:);
    dimsOverlapRefbis = param_facet.dims_overlap_ref(q,:) ;
    IOverlapRefbis = param_facet.I_overlap_ref(q,:);
    
    zerosNum = dimsOverlapRefbis + offsetLbis + offsetRbis;
    x_overlap = zeros([zerosNum(:)',1]);
    
    x_overlap(offsetLbis(1)+1:end-offsetRbis(1),...
        offsetLbis(2)+1:end-offsetRbis(2))...
        = xsol(IOverlapRefbis(1)+1:IOverlapRefbis(1)+dimsOverlapRefbis(1), ...
        IOverlapRefbis(2)+1:IOverlapRefbis(2)+dimsOverlapRefbis(2));
    % forward operator [put the following instructions into a parfeval for parallelisation]
    [SPsitLx{q}, ~, param_facet.dims_PsitLx{q}, param_facet.Ncoefs{q}] = sdwt2_sara_init(x_overlap, param_facet.I(q, :),...
        param_facet.dims(q, :), param_facet.offset, param_facet.status(q, :), nlevel, wavelet);
    
    
    [~, param_facet.I_overlap{q}, param_facet.dims_overlap{q}] = isdwt2_sara_init(SPsitLx{q}, param_facet.I(q, :),...
        param_facet.dims(q, :), param_facet.I_overlap_nc{q}, param_facet.dims_overlap_nc{q}, param_facet.Ncoefs{q}, N, nlevel, wavelet);
    
end

param_facetIn.I = param_facet.I;
param_facetIn.dims =  param_facet.dims;
param_facetIn.offset=param_facet.offset;
param_facetIn.status = param_facet.status;
param_facetIn.Q=Q;
param_facetIn.nlevel =nlevel;
param_facetIn.wavelet =wavelet;
param_facetIn.offsetL= param_facet.offsetL;
param_facetIn.offsetR = param_facet.offsetR;
param_facetIn.dims_overlap_ref = param_facet.dims_overlap_ref ;
param_facetIn.I_overlap_ref = param_facet.I_overlap_ref;
param_facetIn.N =N;
param_facetIn.L= param_facet.L;

param_facetOut.N=N;
param_facetOut.Q=Q;
param_facetOut.nlevel =nlevel;
param_facetOut.wavelet =wavelet;
param_facetOut.I = param_facet.I;
param_facetOut.dims   =    param_facet.dims;
param_facetOut.I_overlap_nc = param_facet.I_overlap_nc;
param_facetOut.dims_overlap_nc = param_facet.dims_overlap_nc;
param_facetOut.Ncoefs = param_facet.Ncoefs;
param_facetOut.I_overlap =param_facet.I_overlap;
param_facetOut.dims_overlap =param_facet.dims_overlap;
param_facetOut.L = param_facet.L;
end
