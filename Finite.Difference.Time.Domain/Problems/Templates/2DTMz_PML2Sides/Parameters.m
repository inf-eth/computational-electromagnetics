% Simulation related parameters.

function [ISize JSize XCenter YCenter delta ra rb DT PMLw] = Parameters

ISize = 100;
JSize = 100;
XCenter = (ISize+1)/2;
YCenter = (JSize-1)/2;
delta = 5.0e-3;

ra = 0.1;
rb = 0.2;
Cl = 299792458;
DT = delta / (sqrt(2) * Cl);
PMLw = 7;   % PML width. It is automatically assumed that outermost PML layer is PEC. Effectively PML is (PMLw - 1) cells wide.