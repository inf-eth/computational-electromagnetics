% er calculation 2D FDTD Drude model template from 10.4-5 JB Shneider.
function return_val = er ( i, j )

[Size XCenter YCenter delta ISlab JSlab] = Parameters;

if i > (XCenter - ISlab/(2*delta)) && i < (XCenter + ISlab/(2*delta)) && j > (YCenter - JSlab/(2*delta)) && j < (YCenter + JSlab/(2*delta))
    return_val = -1;
else
    return_val = 1;
end