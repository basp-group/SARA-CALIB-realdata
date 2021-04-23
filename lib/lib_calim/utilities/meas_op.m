function y = meas_op(G,x,A)

y = cell(size(G));
x = A(x);
clear param_nufft.A;
for i = 1:size(G,1)
    for j=1:size(G,2)
        y{i,j}= G{i,j}*x;
        G{i,j}=[];
    end
end
end