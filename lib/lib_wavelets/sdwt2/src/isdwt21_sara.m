function Facet  = isdwt21_sara(SPsitLx, I, dims, IOverlapNc, ...
    dimsOverlapNc, Ncoefs, N, J, wavelet,dims_o,Io)
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
% modified by Arwa [12/03/19]
%-------------------------------------------------------------------------%
%% dimensions
M       = numel(wavelet); % number of tranforms
% szFacet = max(dims_o) ;% Facet size
idS     = min(Io);% Facet position in the full image
idE     = N - max(Io + dims_o);% Facet position in the full image
Io      = Io - idS  +1;  % Slices position in the Facet
Oo      = Io + dims_o -1 ; % Slices position in the Facet

corner = IOverlapNc;
dimensions =dimsOverlapNc;
temLIdxs =  zeros(size(IOverlapNc));
temRIdxs =  zeros(size(IOverlapNc));
for j =1:size(IOverlapNc,1)
    for i=1:2
        if corner(j,i)<0
            temLIdxs(j,i)= -corner(j,i);
            dimensions(j,i) = dimensions(j,i) - temLIdxs(j,i);
            corner(j,i)=0;
        end
        if corner(j,i)+dimensions(j,i)>=N(i)
            temRIdxs(j,i) = corner(j,i)+dimensions(j,i) - N(i);
            dimensions(j,i) = N(i)-corner(j,i);
        end
    end
end
LIs     = 1 + temLIdxs;
RIs     = temRIdxs;
% folorwing calculation is wrong for Nx~=Ny
%LIs     = 1 +max(0 ,-IOverlapNc);% Facet  cropping positions
%RIs     = max(0,(max(0 , IOverlapNc)+max(dimsOverlapNc,dimsOverlapNc - LIs) -N));


Ncoefs   = reshape(Ncoefs,J+1,M,2);
NcoefsRE = reshape((prod(Ncoefs,3)),[J+1,M]);  
sE = col(3*sum(NcoefsRE(1:J,:),1)) + col(NcoefsRE(J+1,:)) ;
  
%% inverse transform
Facet = reshape(SPsitLx(end -prod(dims_o(M,:))+1:end),dims_o(M,:));% last coeffs --> dirac basis
Facet = padarray(padarray(Facet,Io(M,:)-1,'pre'),max(dims_o) - Oo(M,:),'post'); %dirac basis
sS  = 1;
for w =1: M-1
% dummy  =  Facet(Io(w,1):Oo(w,1),Io(w,2):Oo(w,2));
Facet(Io(w,1):Oo(w,1),Io(w,2):Oo(w,2))= Facet(Io(w,1):Oo(w,1),Io(w,2):Oo(w,2))+...
    isdwt21(SPsitLx(sS:sS+sE(w)),I,dims,Ncoefs(:,w,:),wavelet{w},J,LIs(w,:),RIs(w,:));
sS = sS+sE(w);
end

%% full image size
Facet = padarray(padarray(Facet./sqrt(M),idS,'pre'),idE,'post');

end



%%
% Facet = [Facet; SPsitLx(end -prod(dims_o(M,:))+1:end)]./sqrt(M);
% FacetO.data = Facet./sqrt(M);%
% FacetO.idS = idS;
% FacetO.idE = idE;
% FacetO.N = N;
% tend(M)=toc(ts);
% disp([sum(tend)])
