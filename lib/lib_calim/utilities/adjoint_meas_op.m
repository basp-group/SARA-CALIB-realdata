function x = adjoint_meas_op(G,y,At)
for i = 1:size(G,1)
    for j=1:size(G,2)
        if i==1 && j==1
            x = (G{i,j}' * y{i,j});
        else
            x = x+(G{i,j}' * y{i,j});
        end
        G{i,j} =[]; y{i,j} =[];
    end
end
x = real(At(full(x)));


end
