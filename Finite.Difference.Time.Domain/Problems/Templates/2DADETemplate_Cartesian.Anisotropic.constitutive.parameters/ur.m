% ur tensor calculations for 2D ADE FDTD template.
function return_val = ur ( i, j )

[Size XCenter YCenter delta ra rb] = Parameters;

b = rb;
a = ra;
x = (i-XCenter) * delta;
y = (j-YCenter) * delta;
r = sqrt ( x^2 + y^2 );

sinphi = y/r;
cosphi = x/r;
ur = (r-a)/r;
uphi = r/(r-a);

if (i-XCenter)^2+(j-YCenter)^2 < (rb/delta)^2
    
    if  (i-XCenter)^2+(j-YCenter)^2 > (ra/delta)^2
    
%     %uxx = (r/(r-a)) + ( (a^2-2*a*r)/((r-a)*r^3) )*x^2;
%     %uxy = ( (a^2-2*a*r)/((r-a)*r^3) )*x*y;
%     %uyx = uxy;
%     %uyy = (r/(r-a)) + ( (a^2-2*a*r)/((r-a)*r^3) )*y^2;
%         uzz = (b/(b-a))^2 * ((r-a)/r);
%        u = ((r-a)/r)^2;
%         uxx = ur * (cosphi)^2 + uphi * (sinphi)^2;
%         uxy = (ur - uphi) * sinphi * cosphi;
%         uyx = uxy;
%         uyy = ur * (sinphi)^2 + uphi * (cosphi)^2;
%         return_val = [uxx uxy; uyx uyy];
%     if ( u < 0.1 )
%         u = 0.1;
%     end
        u = 1;
        return_val = [u 0; 0 u];
    else
        return_val = [1 0; 0 1];
    end
else
    return_val = [1 0; 0 1];
end