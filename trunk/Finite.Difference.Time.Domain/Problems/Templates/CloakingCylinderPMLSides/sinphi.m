% sinphi calculation for 2D cloaking problem.
function return_val = sinphi(i, j)

[pI pJ pPMLw pra prb pdt pdelta px0 py0] = Parameters;

x = (i*delta) - px0;
y = (j*delta) - py0;
r = sqrt(x^2 + y^2);

return_val = y/r;