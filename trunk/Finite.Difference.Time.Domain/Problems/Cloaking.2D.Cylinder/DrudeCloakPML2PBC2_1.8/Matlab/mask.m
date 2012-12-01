% E-field mask calculation for cloaking problem.
function return_val = mask(i, j)

[pI pJ pPMLw pra prb pdt pdelta px0 py0] = Parameters;

a = pra;
x = (i*pdelta) - px0;
y = (j*pdelta) - py0;
r = sqrt(x^2 + y^2);

if r^2 <= a^2
    return_val = 0;
else
    return_val = 1;
end
