% wpsquared calculation for 2D FDTD Drude model template from 10.4-5 JB Shneider.
function return_val = wpsquared( i, j, w )
return_val = w^2*(1-er(i, j));