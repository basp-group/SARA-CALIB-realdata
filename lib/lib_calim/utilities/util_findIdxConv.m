function [idxOut] = util_findIdxConv(GIdx,Nxo,Nyo,S)

spmd
    gId = codistributed(GIdx,codistributor('1d',1));
    LId = getLocalPart(gId);
end
clear GIdx gId;

funG = @(M) getInd(M,Nxo,Nyo,S);
spmd
    idC = funG(LId);
end
idxOut = unique(vertcat(idC{:}));
end


%%
function  idU = getInd(Mat,Nxo,Nyo,S)
lm = index2d(1:S^2,S);
rowN = size(Mat,1);
blkSz = 1000;
blk  = floor(rowN/blkSz); 
if blk ==0; blk = 1;end
idxS = 1;
idU  = [];
for blkIdx = 1:blk
    if blkIdx ==blk
        M = Mat(idxS:end,:);
        idxE = size(Mat,1);
    else
        idxE = min((blkIdx*blkSz),rowN);
        M = Mat(idxS:idxE,:);
    end
    v1 = col(col(M(:,1) - lm(:,1).') + lm(:,1).');
    v2 = col(col(M(:,2) - lm(:,2).') + lm(:,2).');
    v1 = unique(sub2ind([Nxo,Nyo],v1,v2));v2=[];
    idU = unique([idU(:); v1(:)]);v1=[];
    idxS = idxE +1 ;
end
end
