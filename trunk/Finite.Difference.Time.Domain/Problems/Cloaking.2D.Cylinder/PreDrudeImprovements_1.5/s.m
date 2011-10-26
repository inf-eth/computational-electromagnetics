% conductance calculation for Yee scattering problem.
function return_val = s ( i, j )

[Size XCenter YCenter delta ra rb] = Parameters;

b = rb;
a = ra;
x = (i-XCenter) * delta;
y = (j-YCenter) * delta;
r = sqrt ( x^2 + y^2 );

if (i-XCenter)^2+(j-YCenter)^2 <= (ra/delta)^2
    return_val = 0;
else
    return_val = 1;
end