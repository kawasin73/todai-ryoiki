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
