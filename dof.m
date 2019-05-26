% dof returns degree of freedom
function n = dof(nx, ny)
    n = (nx) * (ny + 1) * 2 - 3;
    return
end
