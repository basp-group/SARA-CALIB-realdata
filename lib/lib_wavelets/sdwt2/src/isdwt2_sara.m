function PsiSty  = isdwt2_sara(SPsitLx, I, dims, I_overlap_nc, ...
    dims_overlap_nc, Ncoefs, N, J, wavelet)
% Inverse operator to compute the contribution of the SARA prior resulting
% from a single facet.
%-------------------------------------------------------------------------%
%%
% Input:
%
% > SPsitLx          wavelet coefficients obtained from the facet of
%                    interest
% > I                starting index of the facet (without overlap) [1, 2]
% > dims             facet dimension (without overlap) [1, 2]
% > I_overlap_nc     [M, 2]
% > dims_overlap_nc  [M, 2]
% > Ncoefs           [M, 2]
% > N                image dimension [1, 2]
% > J                number of scales considered
% > wavelet          name of the wavelets considered {1, M}
%
% Output:
%
% < PsiSty  inverse transform
%
%-------------------------------------------------------------------------%
%%
% Debug:
% [23/10/18] ok.
% [19/11/18] code acceleration.
%-------------------------------------------------------------------------%
%%
M = numel(wavelet);
PsiSty = []; %cell(M, 1);
start = 1;
start_coefs = 1;
for w = 1:M
    % inverse transform
    if ~strcmp(wavelet{w}, 'self')
        [lo_r, hi_r] = wfilters(wavelet{w}, 'r'); % reconstruction filters
        Ncoefs_m = Ncoefs(start_coefs:start_coefs+(J+1)-1,:);
        s = 3*sum(prod(Ncoefs_m(1:end-1, :), 2)) + prod(Ncoefs_m(end,:));
        SPsitLx_m = SPsitLx(start:start+s-1);
        [PsiSty_m,~,~] = isdwt2(SPsitLx_m, I, I_overlap_nc(w, :), ...
            dims, dims_overlap_nc(w, :), Ncoefs_m, N,...
            lo_r, hi_r, J);
        start = start + s;
        start_coefs = start_coefs + (J+1);
    else
        s = prod(dims);
        PsiSty_m = SPsitLx(start:start+s-1);
        start = start + s;
        start_coefs = start_coefs + 1;
    end
    PsiSty = [PsiSty; PsiSty_m(:)];
end
PsiSty =PsiSty./sqrt(M);
end
