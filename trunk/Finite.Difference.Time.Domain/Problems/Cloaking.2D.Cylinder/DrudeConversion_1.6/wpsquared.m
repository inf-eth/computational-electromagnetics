% wpsquared calculation for cloaking problem.
function return_val = wpsquared( i, j, w )

[Size XCenter YCenter delta ra rb] = Parameters;

A = rb/(rb-ra);

return_val = w^2*(1-(ezz(i, j)/A));
return_val = 0;