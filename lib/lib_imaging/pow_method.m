function val = pow_method(A, At, im_size)
%Computes the maximum eigen value of the compund 
%operator AtA
%   


x=randn(im_size);
x=x/norm(x(:));


p = 1 + 10^(-5) ;
pnew = 1 ;

n = 1 ;

epsilon = 10^(-6) ;

nmax = 200;

cond = abs( pnew-p ) / pnew ;


while ( cond >= epsilon && n < nmax)
    xnew=At(A(x));
    if iscell(xnew)        
        xnew = cell2mat(xnew);
        xnew =reshape(xnew,im_size);
    end
    p=pnew;
    
    pnew = norm(xnew(:)) / norm(x(:));
    cond = abs(  pnew-p ) / pnew ;
    
    x = xnew;
    n = n+1 ;
    
end
val = p ;


end

