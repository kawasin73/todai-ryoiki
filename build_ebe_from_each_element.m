function apply_P = build_ebe_from_each_element(Au, ematrix, nx, ny, is_parallel)
    % validate
    if nx < 2 || ny < 2
        error("nx and ny must greater than or equals to 2")
    elseif size(ematrix, 1) ~= 8 || size(ematrix, 2) ~= 8
        error("matrix of element must be 8 * 8 matrix")
    end
    
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
            end
        
            % https://octave.sourceforge.io/octave/function/lu.html
            % When called with two or three output arguments and a sparse input matrix, lu does not attempt to perform sparsity preserving column permutations
            [L, U] = lu(Ae);
            D = diag(diag(U));
            
            % mutiply L, D, L'
            Ls(:, :, eidx) = L;
            Ds(:, :, eidx) = D;
            Us(:, :, eidx) = L';
            eidx = eidx + 1;
        end
    end
    
    if is_parallel
        apply_P = @(r) apply_ebe_from_each_element_parallel(r, Adis, Ls, Ds, Us, nx, ny);
    else
        apply_P = @(r) apply_ebe_from_each_element(r, Adis, Ls, Ds, Us, nx, ny);
    end
    return;
end

function r = apply_ebe_parallel(r, As, nx, ny, is_reverse)
    if is_reverse
        rows = 1:-1:0;
        cols = 2:-1:1;
    else
        rows = 0:1;
        cols = 1:2;
    end
    for row = rows
        if mod(ny, 2) == 0
            loop = ny / 2;
        else
            loop = (ny + 1 - 2 * row) / 2;
        end
        % prefix shifts segment for row.
        % if row is 0 then each segment is row 1, 3, 5, ...
        % if row is 1 then each segment is row 2, 4, 6, ...
        prefix = row * 2 * nx;
        % r is vector size of (2 * nx * (ny + 1), 1).
        % it is not allowed to access (and update) segment of r by indexing
        % like `r( 2*li+1:2*li ) = [x, x]` in `parfor`, except indexing simply
        % using loop variable like `r(li) = x`.
        % in order to update segment of r, reshape r to `rmat` which size is
        % (segment_size:loop_size) and access `rmat` by `rmat(:,li) = [...]`.
        %
        % segment_size is nx * 4 which is entries for 1 row of elements
        rmat = reshape(r(prefix+1:prefix+4*nx*loop), 4*nx, loop);
        parfor li = 1:loop
            rseg = rmat(:, li);
            j = 2 * li - 1 + row;
            % update entries in 2 colors of 1 row of elements sequentially.
            for s = cols
                % update entries in 1 color in this loop.
                for i = s:2:nx
                    A = get_element_matrix(i, j, nx, ny, As);
                    base_upper = 2 * (i-2) + 2*nx;
                    base_lower = 2 * (i-2);
                    idx = [base_lower+1, base_lower+2, base_lower+3, base_lower+4 base_upper+3, base_upper+4, base_upper+1, base_upper+2];

                    if i == 1
                        A = A(1:4, 1:4);
                        idx = idx(3:6);
                    end
                    rseg(idx) = A \ rseg(idx);
                end
            end
            rmat(:, li) = rseg;
        end
        r(prefix+1:prefix+4*nx*loop) = reshape(rmat, 4*nx*loop, 1);
    end
    return
end

function result = apply_ebe_from_each_element_parallel(r, Adis, Ls, Ds, Us, nx, ny)
    % r = Ads \ r
    r = Adis * r;

    % r = Ls \ r
    r = apply_ebe_parallel(r, Ls, nx, ny, false);

    % r = Ds \ r
    r = apply_ebe_parallel(r, Ds, nx, ny, false);

    % r = Us \ r
    r = apply_ebe_parallel(r, Us, nx, ny, true);

    % r = Ads \ r
    r = Adis * r;

    result = r;
    return;
end

function result = apply_ebe_from_each_element(r, Adis, Ls, Ds, Us, nx, ny)
    % r = Ads \ r
    r = Adis * r;
    
    % r = Ls \ r
    for j = 1:ny
        for i = 1:nx
            L = get_element_matrix(i, j, nx, ny, Ls);
            if i == 1
                L = L(1:4, 1:4);
            end
            idx = build_new_index_for_element(nx, i, j);
            
            r(idx) = L \ r(idx);
        end
    end

    % r = Ds \ r
    for j = 1:ny
        for i = 1:nx
            D = get_element_matrix(i, j, nx, ny, Ds);
            if i == 1
                D = D(1:4, 1:4);
            end
            idx = build_new_index_for_element(nx, i, j);
            
            r(idx) = D \ r(idx);
        end
    end

    % r = Us \ r
    for j = ny:-1:1
        for i = nx:-1:1
            U = get_element_matrix(i, j, nx, ny, Us);
            if i == 1
                U = U(1:4, 1:4);
            end
            idx = build_new_index_for_element(nx, i, j);
            
            r(idx) = U \ r(idx);
        end
    end

    % r = Ads \ r
    r = Adis * r;

    result = r;
    return;
end
