% conductance calculation for Yee scattering problem.
function return_val = s ( i, j )

[ISize JSize XCenter YCenter delta ra rb DT PMLw] = Parameters;

if (i-XCenter)^2+(j-YCenter)^2 <= (ra/delta)^2
    return_val = 0;
else
    return_val = 1;
end