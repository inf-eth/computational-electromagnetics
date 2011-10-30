% uphi calculation for cloaking problem.
function return_val = uphi ( i, j )

[Size XCenter YCenter delta ra rb] = Parameters;

offset = delta * 0;
b = rb;
a = ra + offset;
x = (i-XCenter) * delta;
y = (j-YCenter) * delta;
r = sqrt ( x^2 + y^2 );

uphi = r/(r-a); % ideal
% uphi = (r/(r-a))^2; % Zhao reduced
% uphi = 1; % Pendry reduced
% uphi = 1; % Free space?

if (i-XCenter)^2+(j-YCenter)^2 < ((b/delta)-offset)^2 
    
    if (i-XCenter)^2+(j-YCenter)^2 > ((a/delta)+offset)^2
        return_val = uphi;
    else
        return_val = 1;
    end
else
    return_val = 1;
end