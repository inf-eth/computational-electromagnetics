% ur calculation for Yee scattering problem.
function return_val = ur ( i, j )
delta = 3e-3;
b = 0.06;
a = 0.02;
x = (i-25.5) * delta;
y = (j-24.5) * delta;
r = sqrt ( x^2 + y^2 );

if (i-25.5)^2+(j-24.5)^2 < (b/delta)^2 & (i-25.5)^2+(j-24.5)^2 > (a/delta)^2
    
    uxx = (r/(r-a)) + ( (a^2-2*a*r)/((r-a)*r^3) )*x^2;
    uxy = ( (a^2-2*a*r)/((r-a)*r^3) )*x*y;
    uyx = uxy;
    uyy = (r/(r-a)) + ( (a^2-2*a*r)/((r-a)*r^3) )*y^2;
    uzz = (b/(b-a))^2 * ((r-a)/r);
    
    return_val = [uxx uxy; uyx uyy];
    %return_val = [1 0; 0 1];
else
    return_val = [1 0; 0 1];
end