function val = op_norm(A, At, im_size, tol, max_iter, verbose)
% computes the maximum eigen value of the compund 
% operator AtA

x = randn(im_size);
x = x./(sqrt(sum(sum(x.^2))));
init_val = 1;
for k = 1:(max_iter)
    x = At(A(x));
    val = sqrt(sum(sum(x.^2)));
    rel_var = abs(val-init_val)/init_val;
    if (verbose > 1),   fprintf('Iter = %i, norm = %e \n',k,val);
    end
    if (rel_var < max(2e-6,tol)),  break;
    end
    init_val = val;
    x = x./val;
    
end

if (verbose > 0),fprintf('Norm = %e \n\n', val);
end

end

