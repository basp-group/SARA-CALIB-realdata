function val = op_norm_single_op(A, im_size, tol, max_iter, verbose)
% computes the maximum eigen value of the compund 
% operator AtA
x = randn(im_size);
x = x./sqrt(sum(sum(abs(x.^2))));
init_val = 1;
for k = 1:max_iter
    x = A(x);
    val = sqrt(sum(sum(abs(x.^2))));
    rel_var = abs(val-init_val)/init_val;
    if (verbose > 1)
        fprintf('Iter = %i, norm = %e \n',k,val);
    end
    if (rel_var < tol)
       break;
    end
    init_val = val;
    x = x./val;    
end

if (verbose > 0)
    fprintf('Itr %d Norm = %e \n\n',k, val);
end

end

