% cosphi calculation for 2D cloaking problem.
function return_val = cosphi ( i, j )

[ISize JSize XCenter YCenter delta ra rb DT PMLw] = Parameters;

b = rb;
a = ra;
x = (i-XCenter) * delta;
y = (j-YCenter) * delta;
r = sqrt ( x^2 + y^2 );

return_val = x/r;