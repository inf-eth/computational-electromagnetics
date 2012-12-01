% ge calculation for cloaking problem.
function return_val = gammae(i, j, w)

[pI pJ pPMLw pra prb pdt pdelta px0 py0] = Parameters;

a = pra;
b = prb;
x = (i*pdelta) - px0;
y = (j*pdelta) - py0;
r = sqrt(x^2 + y^2);

g = 0;%w/16;

if r^2 < b^2
    if r^2 > a^2
        return_val = g;
    else
        return_val = 0;
    end
else
    return_val = 0;
end
