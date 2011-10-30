% er calculation for 2D ADE FDTD template.
function return_val = er ( i, j )

[Size XCenter YCenter delta ra rb] = Parameters;

b = rb;
a = ra;
x = (i-XCenter) * delta;
y = (j-YCenter) * delta;
r = sqrt ( x^2 + y^2 );
offset = delta * 0;

if (i-XCenter)^2+(j-YCenter)^2 < ((rb/delta))^2 
    if (i-XCenter)^2+(j-YCenter)^2 > ((ra/delta)+offset)^2
%         ezz = (b/(b-a))^2;
        ezz = 1;
        return_val = ezz;
    else
        return_val = 1;%e80;    % A high value of e can be used to simulate a PEC object.
    end
else
    return_val = 1;
end