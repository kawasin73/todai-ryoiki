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

function A = get_element_matrix(i, j, nx, ny, As)
    eidx = 0;
    switch (j)
    case 1
        eidx = 0;
    case ny
        eidx = 6;
    otherwise
        eidx = 3;
    end
    switch (i)
    case 1
        eidx += 1;
    case nx
        eidx += 3;
    otherwise
        eidx += 2;
    end
    A = As(:,:,eidx);
    return
endfunction

function idx = get_new_idx(nx, i, j)
    idx = build_index_for_element(nx, i, j);
    idx(1:4) -= 2 * j;
    idx(5:8) -= 2 * (j+1);
    if i == 1
        idx = idx(3:6);
    endif
    return
endfunction

# build preconditioning matrix by EBE (Element By Element Method)
function result = apply_ebe_prematrix(r, Adis, Ls, Ds, Us)
    # r = Ads \ r
    r = Adis * r;
    r = Ls \ r;
    r = Ds \ r;
    r = Us \ r;
    r = Adis * r;

    result = r;
    return;
endfunction

function apply_P = build_ebe_prematrix(Au, ematrix, nx, ny)
    # validate
    if nx < 2 || ny < 2
        error("nx and ny must greater than or equals to 2")
    elseif size(ematrix, 1) != 8 || size(ematrix, 2) != 8
        error("matrix of element must be 8 * 8 matrix")
    endif
    
    # diagonal matrix of A
    Ad = diag(diag(Au));
    # Ad ^ 1/2
    Adis = inv(sqrt(Ad));
    
    A9 = build_global_matrix(ematrix, 3, 3, false);
    Ad9 = inv(sqrt(diag(diag(A9))));
    
    L9s = zeros(8, 8, 9);
    D9s = zeros(8, 8, 9);
    U9s = zeros(8, 8, 9);
    
    eidx = 1;
    # caluculate 9 element matrix
    for j = 1:3
        for i = 1:3
            idx = build_index_for_element(3, i, j);
            Ae = eye(8,8) + Ad9(idx, idx) * (ematrix - diag(diag(ematrix))) * Ad9(idx, idx);
            
            if i == 1
                tmp = zeros(8, 8);
                tmp(1:4, 1:4) = Ae(3:6, 3:6);
                Ae = tmp;
            endif
        
            # https://octave.sourceforge.io/octave/function/lu.html
            # When called with two or three output arguments and a sparse input matrix, lu does not attempt to perform sparsity preserving column permutations
            [L, U] = lu(Ae);
            D = diag(diag(U));
            
            # mutiply L, D, L'
            L9s(:, :, eidx) = L;
            D9s(:, :, eidx) = D;
            U9s(:, :, eidx) = L';
            eidx += 1;
        endfor
    endfor
    
    points = nx * (ny+1) * 2;
    
    Us = speye(points, points);
    Ds = speye(points, points);
    Ls = speye(points, points);

    for j = 1:ny
        for i = 1:nx
            idx = get_new_idx(nx, i, j);
            
            # calc Us
            U = get_element_matrix(i, j, nx, ny, U9s);
            if i == 1
                U = U(1:4, 1:4);
            endif
            Us(idx, idx) = Us(idx, idx) * U;

            # calc Ds
            D = get_element_matrix(i, j, nx, ny, D9s);
            if i == 1
                D = D(1:4, 1:4);
            endif
            Ds(idx, idx) = D * Ds(idx, idx);
            
            # calc Ls
            L = get_element_matrix(i, j, nx, ny, L9s);
            if i == 1
                L = L(1:4, 1:4);
            endif
            Ls(idx, idx) = L * Ls(idx, idx);
        endfor
    endfor
    
    apply_P = @(r) apply_ebe_prematrix(r, Adis, Ls, Ds, Us);
    return;
endfunction


#
# THIS IS MAIN TEST
#

# suppress warning
warning('off', 'all');

# parameters
cases = [2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]; # , 200, 300
nmax = 20000;
tol = 10 ^ -6;

# result containers
cerrors = NaN(size(cases, 2), nmax);
citers = NaN(1, size(cases, 2));
ctimes = NaN(1, size(cases, 2));
perrors = NaN(size(cases, 2), nmax);
piters = NaN(1, size(cases, 2));
ptimes = NaN(1, size(cases, 2));
ebeerrors = NaN(size(cases, 2), nmax);
ebeiters = NaN(1, size(cases, 2));
ebetimes = NaN(1, size(cases, 2));
ebecgtimes = NaN(1, size(cases, 2));


# Execute test
k = 0;
for n = cases
    k++
    # setup to create matrix and force
    K = build_global_matrix(e, n, n, true);
    free_idx = build_free_idx(n, n);
    f = zeros(1, 2 * (n+1) * (n+1));
    f(2 * (n+1)) = -1;

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
    [_, iter, errors] = pcg(Ku, fu', x0', nmax, tol, @(r) P \ r);
    [_, finish, _] = cputime();
    total = finish - start;
    disp(sprintf("finish PCG method %d iteration, %d seconds", iter, total));
    # save PCG method result
    perrors(k, :) = errors;
    piters(k) = iter;
    ptimes(k) = total;

    # execute PCG method with EBE
    disp(sprintf("start PCG method with EBE n == %d", n));
    [_, startpre, _] = cputime();
    apply_Pu = build_ebe_prematrix(Ku, e, n, n, free_idx);
    [_, start, _] = cputime();
    [_, iter, errors] = pcg(Ku, fu', x0', nmax, tol, apply_Pu);
    [_, finish, _] = cputime();
    total = finish - startpre;
    disp(sprintf("finish PCG method %d iteration, %d seconds", iter, total));
    # save PCG method result
    ebeerrors(k, :) = errors;
    ebeiters(k) = iter;
    ebetimes(k) = total;
    ebecgtimes(k) = finish - start;
endfor



#
# Print and Save results
#

# save cputime result
loglog(cases, ctimes', ';CG method;', cases, ptimes', ';PCG method;', cases, ebetimes', ';PCG with EBE;', cases, ebecgtimes', ';EBE without P;', cases, ebetimes' - ebecgtimes', ';EBE;');
title("calculation time");
ylabel("calculation time");
xlabel("size of object");
print("cputime.png", "-dpng")
csvwrite("cputime.csv", [cases', ctimes', ptimes', ebetimes', ebecgtimes'])

# save iterations result
loglog(cases, citers', ';CG method;', cases, piters', ';PCG method;', cases, ebeiters', ';PCG with EBE;');
title("count of iteration");
ylabel("iterations");
xlabel("size of object");
print("iterations.png", "-dpng")
csvwrite("iterations.csv", [cases', citers', piters', ebeiters'])

# save errors result
maxn = max(citers(end), piters(end))
semilogy(1:maxn, cerrors(end, 1:maxn), ';CG method;', 1:maxn, perrors(end, 1:maxn), ';PCG method;', 1:maxn, ebeerrors(end, 1:maxn), ';PCG with EBE;');
title("size of error");
ylabel("error");
xlabel("iteration");
print("errors.png", "-dpng")
csvwrite("errors.csv", [(1:maxn)', cerrors(end, 1:maxn)', perrors(end, 1:maxn)', ebeerrors(end, 1:maxn)'])
