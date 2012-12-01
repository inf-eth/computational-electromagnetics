% wpmsquared corrected calculation for cloaking problem.
function return_val = wpmsquaredc(i, j, w)

[pI pJ pPMLw pra prb pdt pdelta px0 py0] = Parameters;

sinwdt = sin(w*pdt/2);
coswdt = cos(w*pdt/2);

%return_val = (2*sinwdt*(-2*(urr(i, j)-1)*sinwdt))/((dt^2)*(coswdt^2));
return_val = 0;

