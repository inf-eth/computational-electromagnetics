% conductance calculation for Yee scattering problem.
function return_val = s ( i, j )
delta = 3e-3;
rb = 0.06;
ra = 0.02;
if (i-25.5)^2+(j-24.5)^2 < (ra/delta)^2
    return_val = 1;
else
    return_val = 1;
end