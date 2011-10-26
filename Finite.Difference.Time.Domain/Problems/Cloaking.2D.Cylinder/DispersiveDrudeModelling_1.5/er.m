% er calculation for Yee scattering problem.
function return_val = er ( i, j )

[Size XCenter YCenter delta ra rb] = Parameters;

b = rb;
a = ra;
x = (i-XCenter) * delta;
y = (j-YCenter) * delta;
r = sqrt ( x^2 + y^2 );
offset = delta * 3;

if (i-XCenter)^2+(j-YCenter)^2 < ((rb/delta))^2 
    
%     %exx = (r/(r-a)) + ( (a^2-2*a*r)/((r-a)*r^3) )*x^2;
%     %exy = ( (a^2-2*a*r)/((r-a)*r^3) )*x*y;
%     %eyx = exy;
%     %eyy = (r/(r-a)) + ( (a^2-2*a*r)/((r-a)*r^3) )*y^2;
    if (i-XCenter)^2+(j-YCenter)^2 > ((ra/delta)+offset)^2
        ezz = (b/(b-a))^2;
        return_val = ezz;
%         return_val = 1;
    else
        return_val = 1;%e80;
    end
%     return_val = 1;
else
    return_val = 1;
end