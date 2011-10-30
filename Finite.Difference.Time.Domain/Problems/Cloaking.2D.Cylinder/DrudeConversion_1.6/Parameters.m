function [Size XCenter YCenter delta ra rb DT] = Parameters
% Simulation related parameters.
Size = 300;
XCenter = (Size+1)/2;
YCenter = (Size-1)/2;
delta = 2.5e-3;

ra = 0.1;
rb = 0.2;
Cl = 299792458;
DT = delta / (sqrt(2) * Cl);
