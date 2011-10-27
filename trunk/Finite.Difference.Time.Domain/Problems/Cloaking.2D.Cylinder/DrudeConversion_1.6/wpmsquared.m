% wpmsquared calculation for cloaking problem.
function return_val = wpmsquared( i, j, w )

[Size XCenter YCenter delta ra rb] = Parameters;

DT = delta / ( sqrt(2) * 3e8 );
sinwDT = sin ( w * DT / 2 );
coswDT = cos ( w * DT / 2 );

return_val = ( 2*sinwDT * ( -2*(urr(i, j)-1)*sinwDT - urrdp(i, j)*gamma(i, j)*DT*coswDT ) ) / ( (DT^2)*(coswDT^2) );
