function apply_P = build_ebe_from_each_element(Au, ematrix, nx, ny)
    % validate
    if nx < 2 || ny < 2
        error("nx and ny must greater than or equals to 2")
    elseif size(ematrix, 1) != 8 || size(ematrix, 2) != 8
        error("matrix of element must be 8 * 8 matrix")
    endif
    
    % diagonal matrix of A
    Ad = diag(diag(Au));
    % Ad ^ 1/2
    Adis = inv(sqrt(Ad));
    
    A9 = build_global_matrix(ematrix, 3, 3, false);
    Ad9 = inv(sqrt(diag(diag(A9))));
    
    Ls = zeros(8, 8, 9);
    Ds = zeros(8, 8, 9);
    Us = zeros(8, 8, 9);
    
    eidx = 1;
    % caluculate each element matrix
    for j = 1:3
        for i = 1:3
            idx = build_index_for_element(3, i, j);
            Ae = eye(8,8) + Ad9(idx, idx) * (ematrix - diag(diag(ematrix))) * Ad9(idx, idx);
            
            if i == 1
                tmp = zeros(8, 8);
                tmp(1:4, 1:4) = Ae(3:6, 3:6);
                Ae = tmp;
            endif
        
            % https://octave.sourceforge.io/octave/function/lu.html
            % When called with two or three output arguments and a sparse input matrix, lu does not attempt to perform sparsity preserving column permutations
            [L, U] = lu(Ae);
            D = diag(diag(U));
            
            % mutiply L, D, L'
            Ls(:, :, eidx) = L;
            Ds(:, :, eidx) = D;
            Us(:, :, eidx) = L';
            eidx += 1;
        endfor
    endfor
    
    apply_P = @(r) apply_ebe_from_each_element(r, Adis, Ls, Ds, Us, nx, ny);
    return;
endfunction

function result = apply_ebe_from_each_element(r, Adis, Ls, Ds, Us, nx, ny)
    % r = Ads \ r
    r = Adis * r;
    
    % r = Ls \ r
    for j = 1:ny
        for i = 1:nx
            L = get_element_matrix(i, j, nx, ny, Ls);
            if i == 1
                L = L(1:4, 1:4);
            endif
            idx = build_new_index_for_element(nx, i, j);
            
            r(idx) = L \ r(idx);
        endfor
    endfor

    % r = Ds \ r
    for j = 1:ny
        for i = 1:nx
            D = get_element_matrix(i, j, nx, ny, Ds);
            if i == 1
                D = D(1:4, 1:4);
            endif
            idx = build_new_index_for_element(nx, i, j);
            
            r(idx) = D \ r(idx);
        endfor
    endfor

    % r = Us \ r
    for j = ny:-1:1
        for i = nx:-1:1
            U = get_element_matrix(i, j, nx, ny, Us);
            if i == 1
                U = U(1:4, 1:4);
            endif
            idx = build_new_index_for_element(nx, i, j);
            
            r(idx) = U \ r(idx);
        endfor
    endfor

    % r = Ads \ r
    r = Adis * r;

    result = r;
    return;
endfunction
