function [G] = calibUpdateDeGriddingMatrix(DDE1,DDE2,Antennas,flag,GMAT,param)
Nxo = param.Nxo;
Nyo = param.Nyo;
N   = Nxo*Nyo;
nMeas  = param.nMeas; 
ddeSpatialDim = param.ddeSpatialDim;
nTimeSlots    = param.nTimeSlots;
rowSupportSize = param.GmatRowSZ;
nufftKernelDim = param.nufftKernelDim;
InactiveTimeSlot = param.InactiveTimeSlot;
AntennasIndex = param.AntennasIndex;
nMeasPerTimeSlot = param.nMeasPerTimeSlot;

cAnt1 = Antennas.ant1;
cAnt2 = Antennas.ant2;
clear Antennas;

GmatRows = GMAT.Gmat_pos;
GmatVals = GMAT.Gmat_values;
clear param GMAT;

%%
deGriddingMatRows = cell(nTimeSlots,1);
deGriddingMatVals = cell(nTimeSlots,1);

parfor itimeSlot = 1:nTimeSlots
    if ~ismember(itimeSlot,InactiveTimeSlot)
        nRowsPerTimeSlot  = nMeasPerTimeSlot(itimeSlot);
        flagTimeSlot      = flag{itimeSlot};
        [~,ants1itimeSlot]  = ismember(nonzeros(cAnt1{itimeSlot}),AntennasIndex);
        [~,ants2itimeSlot]  = ismember(nonzeros(cAnt2{itimeSlot}),AntennasIndex);  
  
        DDE1itimeSlot  = squeeze(DDE1(:,itimeSlot,:));  
        DDE1itimeSlot  = (DDE1itimeSlot(ants1itimeSlot,:)); %d1o already flipped

        DDE2itimeSlot  = squeeze(DDE2(:,itimeSlot,:));   
        DDE2itimeSlot  = DDE2itimeSlot(ants2itimeSlot,:);  % d2o is already conj.
      
        grows   = GmatRows{itimeSlot};  
        GmatRows{itimeSlot}=[];
        gvals   = GmatVals{itimeSlot};
        GmatVals{itimeSlot}=[];
       
        GmatUpdatedRows   = zeros(rowSupportSize,nnz(flagTimeSlot));
        GmatUpdatedValues = zeros(rowSupportSize,nnz(flagTimeSlot));
      
        for iRow = 1 : nRowsPerTimeSlot
            nnziRow = flagTimeSlot(iRow);
            if nnziRow
                indices = col((nnziRow-1)*nufftKernelDim +1 : nnziRow*nufftKernelDim);
                DDEsConvKernel = conv2(reshape(DDE2itimeSlot(nnziRow,:),[ddeSpatialDim,ddeSpatialDim]),reshape(DDE1itimeSlot(nnziRow,:),[ddeSpatialDim,ddeSpatialDim]),'full');
                [GmatUpdatedRows(:,nnziRow),GmatUpdatedValues(:,nnziRow)]=  sconv2_mod(grows(indices),gvals(indices),DDEsConvKernel,Nxo,Nyo);
            end
        end
        grows =[]; gvals =[];
        deGriddingMatVals{itimeSlot} = GmatUpdatedValues(:);
        GmatUpdatedValues=[];
        deGriddingMatRows{itimeSlot} = index1d(shift_ind(index2d(GmatUpdatedRows(:),Nxo),Nxo,Nyo),Nxo); 
        GmatUpdatedRows =[]; 
    end
end
clear Gmat*  DDE*; %save memory
rows = find(vertcat(flag{:}));
clear flag;
G    = sparse(col(repelem(rows,rowSupportSize)),vertcat(deGriddingMatRows{:}),vertcat(deGriddingMatVals{:}),nMeas,N);
clear deGriddingMat* rows* ;
end

%%
function [pos,val] = sconv2_mod(jo,a, B,m,n)
% C = sconv21d(A, B)
% Like conv2 but suitable for convolution of sparse matrices
% Author: Bruno Luong <brunoluong@yahoo.com>
% Date creation: 15/April/2013
% See also: sconv2
idx= index2d(jo,m);
i=idx(:,1);
j=idx(:,2);
[p, q] = size(B);
[k, l, b] = find(B);
[I, K] = ndgrid(i, k);
[J, L] = ndgrid(j, l);
C = a(:)*b(:).';

i = I(:)+K(:)-ceil((p+1)/2);
j = J(:)+L(:)-ceil((q+1)/2);
b = i > 0 & i <= m & ...
    j > 0 & j <= n;
C = sparse(i(b),j(b), C(b), m, n);
[i,j,val]= find(C);
pos= m.*(j-1)+i;
end
