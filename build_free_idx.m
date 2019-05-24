# build free point entry index array
function idx = build_free_idx(nx, ny)
    idx = [];
    for j = 1:(ny+1)
        startidx = 3 + (j-1)*2*(nx+1);
        endidx = 2*(nx+1) + (j-1)*2*(nx+1);
        idx = horzcat(idx, startidx:endidx);
    endfor

    # return free index array
    return
endfunction
