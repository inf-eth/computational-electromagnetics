% cylinder points calculation for cloaking problem.
function return_val = iscylinder ( i, j )

[Size XCenter YCenter delta ra rb] = Parameters;

b = rb;
a = ra;
x = (i-XCenter) * delta;
y = (j-YCenter) * delta;
r = sqrt ( x^2 + y^2 );

if (i-XCenter)^2+(j-YCenter)^2 < ((rb/delta))^2 
    
    if (i-XCenter)^2+(j-YCenter)^2 > ((ra/delta))^2
         return_val = 1;
    else
        return_val = 0;
    end
else
    return_val = 0;
end