%
% Element Matrix
%
e = [ 494505.49450549431 178571.42857142852 -302197.80219780206 -13736.263736263727 -247252.74725274718 -178571.42857142849 54945.054945054952 13736.263736263734
 178571.42857142852 494505.49450549431 13736.263736263731 54945.054945054952 -178571.42857142852 -247252.74725274718 -13736.263736263727 -302197.80219780206
 -302197.80219780206 13736.263736263727 494505.49450549425 -178571.42857142852 54945.054945054944 -13736.263736263732 -247252.74725274718 178571.42857142849
 -13736.263736263736 54945.054945054952 -178571.42857142852 494505.49450549431 13736.263736263732 -302197.802197802 178571.42857142849 -247252.74725274718
 -247252.74725274718 -178571.42857142852 54945.054945054944 13736.263736263732 494505.49450549425 178571.42857142852 -302197.802197802 -13736.263736263725
 -178571.42857142852 -247252.74725274718 -13736.263736263729 -302197.802197802 178571.42857142852 494505.49450549425 13736.263736263732 54945.054945054944
 54945.054945054952 -13736.263736263736 -247252.74725274718 178571.42857142849 -302197.802197802 13736.263736263725 494505.49450549431 -178571.42857142852
 13736.263736263729 -302197.80219780206 178571.42857142852 -247252.74725274718 -13736.263736263732 54945.054945054944 -178571.42857142852 494505.49450549431];

%
% THIS IS MAIN TEST
%

% suppress warning
warning('off', 'all');

% parameters
cases = [2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60];
nmax = 20000;
tol = 10 ^ -10;

% dof array
dofs = NaN(1, size(cases, 2));
for i = 1:size(cases,2)
    dofs(i) = dof(cases(i) * 10, cases(i));
end

% result containers
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

% Execute test
k = 0;
for n = cases
    k = k + 1;
    x = n * 10;
    y = n;
    % setup to create matrix and force
    K = build_global_matrix(e, x, y, true);
    free_idx = build_free_idx(x, y);
    f = zeros(1, 2 * (x+1) * (y+1));
    f(2 * (x+1)) = -1;

    % Ku = f
    Ku = K(free_idx, free_idx);
    fu = f(free_idx);
    x0 = zeros(1, size(fu, 2));

    % execute CG method
    disp(sprintf("start CG method n == %d", n));
    tic();
    [~, iter, errors] = cg(Ku, fu', x0', nmax, tol);
    total = toc();
    disp(sprintf("finish CG method %d iteration, %d seconds", iter, total));
    % save CG method result
    cerrors(k, :) = errors;
    citers(k) = iter;
    ctimes(k) = total;

    P = diag(diag(Ku));

    % execute PCG method
    disp(sprintf("start PCG method n == %d", n));
    tic();
    [~, iter, errors] = pcg(Ku, fu', x0', nmax, tol, @(r) P \ r);
    total = toc();
    disp(sprintf("finish PCG method %d iteration, %d seconds", iter, total));
    % save PCG method result
    perrors(k, :) = errors;
    piters(k) = iter;
    ptimes(k) = total;

    % execute PCG method with EBE
    disp(sprintf("start PCG method with EBE n == %d", n));
    tid = tic();
    
    % ~ 100 elements
    % apply_Pu = build_ebe_from_big_element(K, e, n, n, free_idx);

    % non-parallel mode
    % ~ 50 elements
    % apply_Pu = build_ebe_from_each_element(Ku, e, x, y, false);

    % parallel mode
    apply_Pu = build_ebe_from_each_element(Ku, e, x, y, true);
    
    % ~ 100 elements
    % apply_Pu = build_ebe_from_each_element_compact(Ku, e, n, n);
    
    
    setuptime = toc(tid);
    [~, iter, errors] = pcg(Ku, fu', x0', nmax, tol, apply_Pu);
    total = toc(tid);
    disp(sprintf("finish PCG method %d iteration, %d seconds", iter, total));
    % save PCG method result
    ebeerrors(k, :) = errors;
    ebeiters(k) = iter;
    ebetimes(k) = total;
    ebecgtimes(k) = total - setuptime;
end



%
% Print and Save results
%

% save cputime result
loglog(dofs, ctimes', dofs, ptimes', dofs, ebetimes', dofs, ebecgtimes', dofs, ebetimes' - ebecgtimes');
legend('CG method', 'PCG method', 'PCG with EBE', 'EBE without P', 'EBE');
title("calculation time");
ylabel("calculation time");
xlabel("degrees of freedom");
legend("location", "northwest");
print("cputime-a2.png", "-dpng")
csvwrite("cputime-a2.csv", [dofs', ctimes', ptimes', ebetimes', ebecgtimes'])

% save iterations result
plot(dofs, citers', dofs, piters', dofs, ebeiters');
legend('CG method', 'PCG method', 'PCG with EBE');
title("count of iteration");
ylabel("iterations");
xlabel("degrees of freedom");
legend("location", "northwest");
print("iterations-a2.png", "-dpng")
csvwrite("iterations-a2.csv", [dofs', citers', piters', ebeiters'])

% save errors result
maxn = max(citers(end), piters(end))
semilogy(1:maxn, cerrors(end, 1:maxn), 1:maxn, perrors(end, 1:maxn), 1:maxn, ebeerrors(end, 1:maxn));
legend('CG method', 'PCG method', 'PCG with EBE', 'EBE without P', 'EBE');
title("size of error");
ylabel("error");
xlabel("iteration");
print("errors-a2.png", "-dpng")
csvwrite("errors-a2.csv", [(1:maxn)', cerrors(end, 1:maxn)', perrors(end, 1:maxn)', ebeerrors(end, 1:maxn)'])
