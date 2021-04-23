function SPsitLx = sdwt21_sara(Slice, I, dims, offset, status, J, wavelet)
% Forward operator to compute the contribution of the SARA prior relative
% to a single facet.
%%
%-------------------------------------------------------------------------%
% Input:
%
% > x_overlap  ...
% > I          starting index of the facet (without overlap) [1, 2]
% > dims       facet dimension (without overlap) [1, 2]
% > offset     offset to be considered for each dictionary w.r.t. the largest overlapping
%              facet x_overlap [M, 1]
% > status     indicates whether the facet is the first (or the last)
%              along each dimension {2, 1}
% > J          number of scales considered
% > wavelet    name of the wavelets considered {1, M}
%
% Output:
%
% < SPsitLx
% < Ij
% < dims_PsitLx
% < Ncoefs
%
%-------------------------------------------------------------------------%
%%
% Debug:
% [23/10/18] ok.
%modified by Arwa [12/03/19]
%-------------------------------------------------------------------------%
%%
M = numel(wavelet);
Slice = Slice./sqrt(M);
szS   = size(Slice);

LId   = zeros(2, M);
for i = 1:2
if I(i) > 0
    LId(i,:) = offset;% no offset if I(q, i) == 0
    LId(i,M) = LId(i,M) + mod(I(i), 2^J);
end
end
LId = LId +1;

SPsitLx = [];
for w = 1:M-1       % forward transform
    dummy = sdwt21(Slice(LId(1,w):szS(1),LId(2,w):szS(2)),I,dims,status,wavelet{w},J);
    SPsitLx = [SPsitLx ; dummy];
end

SPsitLx = [SPsitLx ;col(Slice(LId(1,M):end, LId(2,M):end))];
end

