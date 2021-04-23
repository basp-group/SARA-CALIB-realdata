function idx=index2d(k,Nx)

j=floor(k/Nx);
i=mod(k,Nx);

j(i>0)=j(i>0)+1;
i(i==0)=Nx;
    
idx=[i(:) j(:)];


end
