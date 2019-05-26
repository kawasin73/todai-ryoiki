% build index in unknown matrix for element
% returns 8 size vector or 4 size vector (4 size when left fixed element)
function idx = build_new_index_for_element(nx, i, j)
    idx = build_index_for_element(nx, i, j);
    idx(1:4) -= 2 * j;
    idx(5:8) -= 2 * (j+1);
    if i == 1
        idx = idx(3:6);
    end
    return
end
