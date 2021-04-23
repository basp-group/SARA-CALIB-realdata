function [PsiSty] = isdwt21(SPsitLx,I,SzEnd,Ncoefs,wavelet,J,temLIdxs,temRIdxs)
%% inverse wavelet transform
[lo, hi]= wfilters(wavelet, 'r');
SzdJ = Ncoefs(1:J,:)+ [(2.^(J-(1:J-1).')-1)*(length(lo)-2)+floor(mod(I,2^J)./(2.^(1:J-1).'));zeros(1,2)];
dS3  = 3*prod(Ncoefs(1:J,:),2);
% case J+1
Pos = sum(dS3) ;
PsiSty = reshape(SPsitLx(1+Pos : prod(Ncoefs(J+1,:))+Pos), Ncoefs(J+1,:));
%case j<J+1
for j = J:-1:1
    dSPsitLx = reshape(SPsitLx(Pos - dS3(j)+1 : Pos),[Ncoefs(j,:),3]);
    Pos  = Pos - dS3(j);
    if j < J, dSPsitLx = padarray(dSPsitLx,SzdJ(j,:)-Ncoefs(j,:),'pre');
    end
    if j == 1,  idx = SzEnd + (2^J-1)*(length(lo)-2) + mod(I,2^J);
    else, idx = SzdJ(j-1,:);
    end
    % upsampling, convolution and cropping    
    PsiSty = upConv2r(upConv2c(PsiSty,lo,idx(1)) + upConv2c(dSPsitLx(:,:,1),hi,idx(1))   ,lo, idx(2)) + ...
        upConv2r(upConv2c(dSPsitLx(:,:,2),lo,idx(1))+upConv2c(dSPsitLx(:,:,3),hi,idx(1)) ,hi, idx(2));   
      
end
%% Crop

PsiSty = PsiSty(temLIdxs:size(PsiSty,1)-temRIdxs(1),temLIdxs(2):size(PsiSty,2)-temRIdxs(2));
end

%% internal function
function y = upConv2c(x, w, rows)
z = zeros(2*size(x, 1), size(x,2));
z(1:2:end, :) = x;
y = conv2(z, w.');
y = y(1:rows, :);
end
function y = upConv2r(x, w, cols)
z = zeros(size(x,1), 2*size(x, 2));
z(:, 1:2:end) = x;
y = conv2(z, w);
y = y(:, 1:cols);
end



% %     PsiSty = upConv2D(PsiSty,lo,lo,idr,idc)+ ...
% %         upConv2D(dSPsitLx(:,:,1),hi,lo,idr,idc)+ ...
% %         upConv2D(dSPsitLx(:,:,2),lo,hi,idr,idc)+ ...
% %         upConv2D(dSPsitLx(:,:,3),hi,hi,idr,idc);
% % function x = upConv2D(x, wc, wr, rows, cols)
% % z = zeros(2*size(x,1),size(x,2));%dyadup(x,0,'r');
% % z(1:2:end, :) = x;
% % y = conv2(z, wc.');
% % z  = zeros(rows,2*size(y,2));
% % z(:,1:2:end) =  y(1:rows, :);
% % y = conv2(z, wr);
% % x = y(:,1:cols);
% % end
