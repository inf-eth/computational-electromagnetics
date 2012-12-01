% wpmsquared calculation for cloaking problem.
function return_val = wpmsquared(i, j, w)

[pI pJ pPMLw pra prb pdt pdelta px0 py0] = Parameters;

return_val = w^2*(uinf(i,j)-uphi(i,j));