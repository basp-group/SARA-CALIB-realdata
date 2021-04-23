function val = op_norm_mod(A, im_size, tol, max_iter, verbose)
% computes the maximum eigen value of the compund 
% from Alex's onose lib
% operator A -->AtA
x = randn(im_size)+1i*randn(im_size);
x = x./norm(x(:));
init_val = 1;

for k = 1:max_iter
    x = (A(x));
    val = norm(x(:));
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
    fprintf('Norm = %e \n\n', val);
end

end

