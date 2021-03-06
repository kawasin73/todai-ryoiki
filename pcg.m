%
% Preconditioned Conjugate Gradient Method
%

% implementation of Preconditioned Conjugate Gradient Method
function [x, iter, errors] = pcg(A, b, x0, nmax, tol, apply_P)
    % setup
    iter = 0;
    x = x0;
    r = b - (A * x);
    z = apply_P(r);
    p = z;

    rho2 = 0;
    rho = r' * z;

    % preallocate errors array
    errors = NaN(1, nmax);
    errors(1) = r' * r;

    while iter < nmax-1 && errors(iter+1) > tol
        iter = iter + 1;
        if iter ~= 1
            beta = rho / rho2;
            p = z + beta * p;
        end

        alpha = rho / (p' * A * p);
        x = x + alpha * p;
        rho2 = rho;
        r = r - alpha * (A * p);
        z = apply_P(r);
        rho = r' * z;
        errors(iter + 1) = r' * r;
    end

    return
end
