% er calculation for Yee scattering problem.
function return_val = er ( i, j )
delta = 3e-3;
b = 0.06;
a = 0.04;
x = (i-25.5) * delta;
y = (j-24.5) * delta;
r = sqrt ( x^2 + y^2 );

if (i-25.5)^2+(j-24.5)^2 < (b/delta)^2 & (i-25.5)^2+(j-24.5)^2 > (a/delta)^2
    
    exx = (r/(r-a)) + ( (a^2-2*a*r)/((r-a)*r^3) )*x^2;
    exy = ( (a^2-2*a*r)/((r-a)*r^3) )*x*y;
    eyx = exy;
    eyy = (r/(r-a)) + ( (a^2-2*a*r)/((r-a)*r^3) )*y^2;
    ezz = (b/(b-a))^2 * ((r-a)/r);
    
    return_val = ezz;
    %return_val = 4;
else
    return_val = 1;
end