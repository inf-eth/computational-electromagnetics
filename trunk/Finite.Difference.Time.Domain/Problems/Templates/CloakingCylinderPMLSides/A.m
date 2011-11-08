% A calculation for cloaking problem.
function return_val = A ( i, j )

[ISize JSize XCenter YCenter delta ra rb DT PMLw] = Parameters;

b = rb;
a = ra;
x = (i-XCenter) * delta;
y = (j-YCenter) * delta;
r = sqrt ( x^2 + y^2 );
offset = delta * 0;

A = b/(b-a); % Zhao reduced.
% A = 1;  % Free space simulation?

if (i-XCenter)^2+(j-YCenter)^2 < ((rb/delta))^2 
    if (i-XCenter)^2+(j-YCenter)^2 > ((ra/delta)+offset)^2
        return_val = A;
    else
        return_val = 1;
    end
else
    return_val = 1;
end