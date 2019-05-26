% build preconditioning matrix by EBE (Element By Element Method)
function apply_P = build_ebe_from_big_element(A, ematrix, nx, ny, free_idx)
    points = (nx+1) * (ny+1) * 2;

    % validate
    if nx < 1 || ny < 1
        error("nx and ny must greater than or equals to 2")
    else size(ematrix, 1) != 8 || size(ematrix, 2) != 8
        error("matrix of element must be 8 * 8 matrix")
    else size(A, 1) != points || size(A, 2) != points
        error("matrix A is not match to nx, ny")
    end

    % diagonal matrix of A
    Ad = diag(diag(A));

    Eye = speye(points, points);
    Ls = speye(points, points);
    Ds = speye(points, points);
    Ad = sparse(Ad);

    % Ad ^ 1/2
    Ads = sqrt(Ad);
    % Ad ^ -1/2
    Adsi = inv(Ads);

    % caluculate each element matrix
    for j = 1:(ny)
        for i = 1:(nx)
            % create element matrix
            idx = build_index_for_element(nx, i, j);
            Ae = sparse(points, points);
            Ae(idx, idx) += ematrix;

            % scaled, regularzied element array
            rAe = Eye + Adsi * (Ae - diag(diag(Ae))) * Adsi;

            % https://octave.sourceforge.io/octave/function/lu.html
            % When called with two or three output arguments and a sparse input matrix, lu does not attempt to perform sparsity preserving column permutations
            [L, U] = lu(rAe);
            D = diag(diag(U));

            % mutiply L, D, L'
            Ls = Ls * L;
            Ds = Ds * D;
        end
    end

    apply_P = @(r) apply_ebe_from_big_element(r, Adsi(free_idx, free_idx), Ls(free_idx, free_idx), Ds(free_idx, free_idx), Ls'(free_idx, free_idx));
    return;
end

function result = apply_ebe_from_big_element(r, Adsi, Ls, Ds, Us)
    r = Adsi * r;
    r = Ls \ r;
    r = Ds \ r;
    r = Us \ r;
    r = Adsi * r;
    result = r;
    return;
end
