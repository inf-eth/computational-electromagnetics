% sinphi calculation for 2D cloaking problem.
function return_val = sinphi ( i, j )

[Size XCenter YCenter delta ra rb] = Parameters;

b = rb;
a = ra;
x = (i-XCenter) * delta;
y = (j-YCenter) * delta;
r = sqrt ( x^2 + y^2 );

return_val = y/r;