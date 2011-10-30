function [Size XCenter YCenter delta ISlab JSlab DT] = Parameters
% Simulation related parameters.
Size = 200;
XCenter = (Size+1)/2;
YCenter = (Size-1)/2;
delta = 2.5e-3;

ISlab = 0.6;
JSlab = 0.1;
Cl = 3e8;
DT = delta / (sqrt(2) * Cl);