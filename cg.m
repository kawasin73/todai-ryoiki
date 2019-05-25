#
# Conjugate Gradient Method
#

# implementation of Conjugate Gradient Method
function [x, iter, errors] = cg(A, b, x0, nmax, tol)
    # setup
    iter = 0;
    x = x0;
    r = b - (A * x);
    p = r;
    rho2 = 0;
    rho = r' * r;

    # preallocate errors array
    errors = NaN(1, nmax);
    errors(1) = rho;

    while iter < nmax-1 && rho > tol
        iter++;
        if iter != 1
            beta = rho / rho2;
            p = r + beta * p;
        endif

        alpha = rho / (p' * A * p);
        x = x + alpha * p;
        rho2 = rho;
        r = r - alpha * (A * p);
        rho = r' * r;
        errors(iter + 1) = rho;
        x;
    endwhile

    return
endfunction
