% wpmsquareduxx calculation for cloaking problem.
function return_val = wpmsquareduxx( i, j, w )

[Size XCenter YCenter delta ra rb DT] = Parameters;

sinwDT = sin ( w * DT / 2 );
coswDT = cos ( w * DT / 2 );

b = rb;
a = ra;
x = (i-XCenter) * delta;
y = (j-YCenter) * delta;
r = sqrt ( x^2 + y^2 );
if (i-XCenter)^2+(j-YCenter)^2 < (rb/delta)^2
    
    if  (i-XCenter)^2+(j-YCenter)^2 > (ra/delta)^2
       u = ((r-a)/r)^2;
       return_val = ( 2*sinwDT * ( -2*(u-1)*sinwDT - urrdp(i, j)*gamma(i, j)*DT*coswDT ) ) / ( (DT^2)*(coswDT^2) );
    else
        return_val = 1;
    end
else
    return_val = 1;
end