% A calculation for cloaking problem.
function return_val = A ( i, j )

[Size XCenter YCenter delta ra rb] = Parameters;

b = rb;
a = ra;
x = (i-XCenter) * delta;
y = (j-YCenter) * delta;
r = sqrt ( x^2 + y^2 );
offset = delta * 0;

if (i-XCenter)^2+(j-YCenter)^2 < ((rb/delta))^2 
    
    if (i-XCenter)^2+(j-YCenter)^2 > ((ra/delta)+offset)^2
        A = b/(b-a); % Zhao reduced.
        return_val = A;
        return_val = 1;
    else
        return_val = 1;%e80;
    end
else
    return_val = 1;
end