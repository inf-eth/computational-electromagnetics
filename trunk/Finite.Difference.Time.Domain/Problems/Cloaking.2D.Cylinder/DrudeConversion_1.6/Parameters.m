function [Size XCenter YCenter delta ra rb DT PMLw] = Parameters
% Simulation related parameters.
Size = 250;
XCenter = (Size+1)/2;
YCenter = (Size-1)/2;
delta = 2.0e-3;

ra = 0.1;
rb = 0.2;
Cl = 299792458;
DT = delta / (sqrt(2) * Cl);

PMLw = 5;   % Width or thickness of PML layer
