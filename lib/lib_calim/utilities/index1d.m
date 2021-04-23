function idx=index1d(k,Nx)

i=k(:,1);
j=k(:,2);
idx= Nx.*(j-1)+i;

end
