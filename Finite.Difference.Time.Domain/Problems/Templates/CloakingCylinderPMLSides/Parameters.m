% Simulation related parameters.

function [ISize JSize XCenter YCenter delta ra rb DT PMLw dtscalar] = Parameters

ISize = 200;
JSize = 200;
XCenter = (ISize+1)/2;
YCenter = (JSize-1)/2;
delta = 2.5e-3;

ra = 0.1;
rb = 0.2;
Cl = 299792458;
dtscalar = 2;
DT = delta / (sqrt(2) * Cl) / dtscalar;
PMLw = 50;   % PML width. It is automatically assumed that outermost PML layer is PEC. Effectively PML is (PMLw - 1) cells wide.