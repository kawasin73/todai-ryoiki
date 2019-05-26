% get matrix from As which match type 1 ~ 9.
function A = get_element_matrix(i, j, nx, ny, As)
    eidx = 0;
    switch (j)
    case 1
        eidx = 0;
    case ny
        eidx = 6;
    otherwise
        eidx = 3;
    end
    switch (i)
    case 1
        eidx += 1;
    case nx
        eidx += 3;
    otherwise
        eidx += 2;
    end
    A = As(:,:,eidx);
    return
end
