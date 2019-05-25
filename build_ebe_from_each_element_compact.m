function apply_P = build_ebe_from_each_element_compact(Au, ematrix, nx, ny)
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
            idx = build_new_index_for_element(nx, i, j);
            
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
    
    apply_P = @(r) apply_ebe_from_each_element_compact(r, Adis, Ls, Ds, Us);
    return;
endfunction

function result = apply_ebe_from_each_element_compact(r, Adis, Ls, Ds, Us)
    r = Adis * r;
    r = Ls \ r;
    r = Ds \ r;
    r = Us \ r;
    r = Adis * r;

    result = r;
    return;
endfunction