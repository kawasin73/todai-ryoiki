#
# Setup Functions
#

# build Kg matrix size (nx * ny) from Ke (matrix of 1 element entry 8 * 8)
function result = build_global_matrix(ematrix, nx, ny, is_s)
    # validate
    if nx < 2 || ny < 2
        error("nx and ny must greater than or equals to 2")
    elseif size(ematrix, 1) != 8 || size(ematrix, 2) != 8
        error("matrix of element must be 8 * 8 matrix")
    endif

    # build base matrix
    points = nx * ny * 2;

    if is_s
        result = sparse(points, points);
    else
        result = zeros(points, points);
    endif

    # add element entry to each point
    for j = 1:(ny - 1)
        for i = 1:(nx - 1)
            base = (i - 1) * 2;
            base_lower = 2 * nx * (j - 1);
            base_upper = base_lower + 2 * nx;
            con = [base+1, base+2, base+3, base+4, base+3, base+4, base+1, base+2];
            con(1:4) += base_lower;
            con(5:8) += base_upper;
            # add each element entry
            result(con, con) += ematrix;
        endfor
    endfor

    # return result matrix
    return;
endfunction

# build free point entry index array
function idx = build_free_idx(nx, ny)
    idx = [];
    for j = 1:ny
        startidx = 3 + (j-1)*2*nx;
        endidx = 2*nx + (j-1)*2*nx;
        idx = horzcat(idx, startidx:endidx);
    endfor

    # return free index array
    return
endfunction

#
# Element Matrix
#
e = [ 494505.49450549431 178571.42857142852 -302197.80219780206 -13736.263736263727 -247252.74725274718 -178571.42857142849 54945.054945054952 13736.263736263734
 178571.42857142852 494505.49450549431 13736.263736263731 54945.054945054952 -178571.42857142852 -247252.74725274718 -13736.263736263727 -302197.80219780206
 -302197.80219780206 13736.263736263727 494505.49450549425 -178571.42857142852 54945.054945054944 -13736.263736263732 -247252.74725274718 178571.42857142849
 -13736.263736263736 54945.054945054952 -178571.42857142852 494505.49450549431 13736.263736263732 -302197.802197802 178571.42857142849 -247252.74725274718
 -247252.74725274718 -178571.42857142852 54945.054945054944 13736.263736263732 494505.49450549425 178571.42857142852 -302197.802197802 -13736.263736263725
 -178571.42857142852 -247252.74725274718 -13736.263736263729 -302197.802197802 178571.42857142852 494505.49450549425 13736.263736263732 54945.054945054944
 54945.054945054952 -13736.263736263736 -247252.74725274718 178571.42857142849 -302197.802197802 13736.263736263725 494505.49450549431 -178571.42857142852
 13736.263736263729 -302197.80219780206 178571.42857142852 -247252.74725274718 -13736.263736263732 54945.054945054944 -178571.42857142852 494505.49450549431];
 size(e)

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

    while iter < nmax && rho > tol
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

#
# Preconditioned Conjugate Gradient Method
#

# implementation of Preconditioned Conjugate Gradient Method
function [x, iter, errors] = pcg(A, b, x0, nmax, tol, P)
    # setup
    iter = 0;
    x = x0;
    r = b - (A * x);
    z = P \ r;
    p = z;

    rho2 = 0;
    rho = r' * z;

    # preallocate errors array
    errors = NaN(1, nmax);
    errors(1) = r' * r;

    while iter < nmax && errors(iter+1) > tol
        iter++;
        if iter != 1
            beta = rho / rho2;
            p = z + beta * p;
        endif

        alpha = rho / (p' * A * p);
        x = x + alpha * p;
        rho2 = rho;
        r = r - alpha * (A * p);
        z = P \ r;
        rho = r' * z;
        errors(iter + 1) = r' * r;
        x;
    endwhile

    return
endfunction

#
# THIS IS MAIN TEST
#

# parameters
cases = [2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300];
nmax = 20000;
tol = 10 ^ -6;

# result containers
cerrors = NaN(size(cases, 2), nmax);
citers = NaN(1, size(cases, 2));
ctimes = NaN(1, size(cases, 2));
perrors = NaN(size(cases, 2), nmax);
piters = NaN(1, size(cases, 2));
ptimes = NaN(1, size(cases, 2));

# Execute test
k = 0;
for n = cases
    k++
    # setup to create matrix and force
    K = build_global_matrix(e, n, n, true);
    free_idx = build_free_idx(n, n);
    f = zeros(1, 2 * n * n);
    f(2 * n) = -1;

    # Ku = f
    Ku = K(free_idx, free_idx);
    fu = f(free_idx);
    x0 = zeros(1, size(fu, 2));

    # execute CG method
    disp(sprintf("start CG method n == %d", n));
    [_, start, _] = cputime();
    [_, iter, errors] = cg(Ku, fu', x0', nmax, tol);
    [_, finish, _] = cputime();
    total = finish - start;
    disp(sprintf("finish CG method %d iteration, %d seconds", iter, total));
    # save CG method result
    cerrors(k, :) = errors;
    citers(k) = iter;
    ctimes(k) = total;

    P = diag(diag(Ku));

    # execute PCG method
    disp(sprintf("start PCG method n == %d", n));
    [_, start, _] = cputime();
    [_, iter, errors] = pcg(Ku, fu', x0', nmax, tol, P);
    [_, finish, _] = cputime();
    total = finish - start;
    disp(sprintf("finish PCG method %d iteration, %d seconds", iter, total));

    # save PCG method result
    perrors(k, :) = errors;
    piters(k) = iter;
    ptimes(k) = total;
endfor

#
# Print and Save results
#

# save cputime result
loglog(cases, ctimes', ';CG method;', cases, ptimes', ';PCG method;');
title("calculation time");
ylabel("calculation time");
xlabel("size of object");
print("cputime.png", "-dpng")
csvwrite("cputime.csv", [cases', ctimes', ptimes'])

# save iterations result
loglog(cases, citers', ';CG method;', cases, piters', ';PCG method;');
title("count of iteration");
ylabel("iterations");
xlabel("size of object");
print("iterations.png", "-dpng")
csvwrite("iterations.csv", [cases', citers', piters'])

# save errors result
maxn = max(citers(end), piters(end))
semilogy(1:maxn, cerrors(end, 1:maxn), ';CG method;', 1:maxn, perrors(end, 1:maxn), ';PCG method;');
title("size of error");
ylabel("error");
xlabel("iteration");
print("errors.png", "-dpng")
csvwrite("errors.csv", [(1:maxn)', cerrors(end, 1:maxn)', perrors(end, 1:maxn)'])
