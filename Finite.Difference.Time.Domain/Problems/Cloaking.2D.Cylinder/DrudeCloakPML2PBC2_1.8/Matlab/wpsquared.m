% wpsquared calculation for cloaking problem.
function return_val = wpsquared(i, j, w)

[ISize JSize XCenter YCenter delta ra rb DT PMLw] = Parameters;

return_val = w^2*(einf(i,j)-ezz(i,j)/A(i,j));
