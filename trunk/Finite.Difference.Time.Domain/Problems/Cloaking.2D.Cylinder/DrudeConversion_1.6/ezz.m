% ezz calculation for cloaking problem.
function return_val = ezz ( i, j )

[Size XCenter YCenter delta ra rb] = Parameters;

b = rb;
a = ra;
x = (i-XCenter) * delta;
y = (j-YCenter) * delta;
r = sqrt ( x^2 + y^2 );
offset = delta * 0;

A = b/(b-a);
ezz = ((b/(b-a))^2) * ((r-a)/r);    % Ideal.
% ezz = (b/(b-a))^2;  % Pendry reduced
% ezz = b/(b-a); % Zhao reduced.

   
if (i-XCenter)^2+(j-YCenter)^2 < ((rb/delta))^2 
    
    if (i-XCenter)^2+(j-YCenter)^2 > ((ra/delta)+offset)^2
        return_val = ezz;        
    else
        return_val = A;%e80;
    end
else
    return_val = A;
end