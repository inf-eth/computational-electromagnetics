% urr calculation for cloaking problem.
function return_val = urr(i, j)

[pI pJ pPMLw pra prb pdt pdelta px0 py0] = Parameters;

a = pra;
b = prb;
x = (i*pdelta) - px0;
y = (j*pdelta) - py0;
r = sqrt(x^2 + y^2);

urr = (r-a)/r;      % ideal
%urr = ((b/(b-1))*((r-a)/r))^2;  % Zhao reduced.
% urr = (r/(r-a))^2;  % Pendry reduced.
%urr = 1;

if r^2 < b^2
    if r^2 > a^2
        return_val = urr;
    else
        return_val = 1;
    end
else
    return_val = 1;
end
