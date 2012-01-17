% Simulation related parameters.

function [ISize JSize XCenter YCenter delta ra rb DT PMLw dtscalar Skin TissueIW TissueIIW PulseWidth TumourX TumourY TumourRadius] = Parameters

ISize = 200;
JSize = 300;
XCenter = (ISize+1)/2;
YCenter = (JSize-1)/2;
delta = 1.0e-3;

% Skin layers.
% |Skin |TissueI|TissueII|
% | 5mm | 40mm  | 60mm   |

Skin = 5e-3;
TissueIW = 40e-3;
TissueIIW = 60e-3;
TumourX = 0.0e-3;
TumourY = 50e-3;
TumourRadius = 4e-3;

ra = 0.1;
rb = 0.2;
Cl = 299792458;
dtscalar = 4;
DT = delta / (sqrt(2) * Cl);
PulseWidth = 15;
PMLw = 50;   % PML width. It is automatically assumed that outermost PML layer is PEC. Effectively PML is (PMLw - 1) cells wide.