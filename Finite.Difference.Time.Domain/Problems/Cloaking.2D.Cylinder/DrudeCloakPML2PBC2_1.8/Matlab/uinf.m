% uinf calculation for cloaking problem.
function return_val = uinf(i, j)

[pI pJ pPMLw pra prb pdt pdelta px0 py0] = Parameters;

a = pra;
b = prb;
x = (i*pdelta) - px0;
y = (j*pdelta) - py0;
r = sqrt(x^2 + y^2);

uinf = 1;

if r^2 < b^2
    if r^2 > a^2
        return_val = uinf;
    else
        return_val = 1;
    end
else
    return_val = 1;
end