% wpsquared calculation for cloaking problem.
function return_val = wpsquared( i, j, w )

[Size XCenter YCenter delta ra rb] = Parameters;

return_val = w^2*(1-(ezz(i, j)/A(i, j)));