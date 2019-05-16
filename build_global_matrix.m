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
            idx = build_index_for_element(nx, i, j);
            # add each element entry
            result(idx, idx) += ematrix;
        endfor
    endfor

    # return result matrix
    return;
endfunction