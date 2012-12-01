% sinphi calculation for 2D cloaking problem.
function return_val = sinphi ( i, j )

[pI pJ pPMLw pra prb pdt pdelta px0 py0] = Parameters;

x = (i*pdelta) - px0;
y = (j*pdelta) - py0;
r = sqrt(x^2 + y^2);

if r == 0
    return_val = 0;
else
    return_val = y/r;
end