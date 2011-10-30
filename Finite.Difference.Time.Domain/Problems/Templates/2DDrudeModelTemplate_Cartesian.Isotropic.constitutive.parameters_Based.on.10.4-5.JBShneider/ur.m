% ur calculation 2D FDTD Drude model template from 10.4-5 JB Shneider.
function return_val = ur ( i, j )

[Size XCenter YCenter delta ra rb] = Parameters;

if (i-XCenter)^2+(j-YCenter)^2 < (rb/delta)^2
    if  (i-XCenter)^2+(j-YCenter)^2 > (ra/delta)^2
        ur = 1;
        return_val = ur;
    else
        return_val = 1;
    end
else
    return_val = 1;
end