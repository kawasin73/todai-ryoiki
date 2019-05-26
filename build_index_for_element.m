% build index in matrix for element
% returns 8 size vector
function idx = build_index_for_element(nx, i, j)
    base = (i - 1) * 2;
    base_lower = 2 * (nx + 1) * (j - 1);
    base_upper = base_lower + 2 * (nx + 1);
    idx = [base+1, base+2, base+3, base+4, base+3, base+4, base+1, base+2];
    idx(1:4) += base_lower;
    idx(5:8) += base_upper;

    return;
endfunction
