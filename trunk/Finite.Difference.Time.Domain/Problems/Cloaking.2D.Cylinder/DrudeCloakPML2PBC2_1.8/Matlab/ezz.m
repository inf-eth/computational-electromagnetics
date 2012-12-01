% ezz calculation for cloaking problem.
function return_val = ezz (i, j)

[pI pJ pPMLw pra prb pdt pdelta px0 py0] = Parameters;

a = pra;
b = prb;
x = (i*pdelta) - px0;
y = (j*pdelta) - py0;
r = sqrt(x^2 + y^2);

ezz = ((b/(b-a))^2) * ((r-a)/r);    % Ideal.
%ezz = (b/(b-a))^2;  % Pendry reduced
% ezz = b/(b-a); % Zhao reduced.
%ezz = 1; % Free space?

if r^2 < b^2
    if r^2 > a^2
        return_val = ezz;
    else
        return_val = 1;
    end
else
    return_val = 1;
end
