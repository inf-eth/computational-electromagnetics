% urr calculation for cloaking problem.
function return_val = urr ( i, j )

[Size XCenter YCenter delta ra rb] = Parameters;

offset = delta * 0;
b = rb;
a = ra + offset;
x = (i-XCenter) * delta;
y = (j-YCenter) * delta;
r = sqrt ( x^2 + y^2 );

if (i-XCenter)^2+(j-YCenter)^2 < ((b/delta))^2 
    
    if (i-XCenter)^2+(j-YCenter)^2 > ((a/delta)+offset)^2
%         urr = (r-a)/r;      % ideal
        urr = ((b/(b-1))*((r-a)/r))^2;  % Zhao reduced.
%         urr = (r/(r-a))^2;  % Pendry reduced.
        return_val = urr;
        return_val = 1;
    else
        return_val = 1;%e80;
    end
else
    return_val = 1;
end