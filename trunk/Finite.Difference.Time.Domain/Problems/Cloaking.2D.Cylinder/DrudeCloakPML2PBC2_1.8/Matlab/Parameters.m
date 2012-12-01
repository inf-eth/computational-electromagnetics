% Simulation related parameters.

function [I J PMLw ra rb dt delta x0 y0] = Parameters

I = 300;
J = 300;
PMLw = 50;

ra = 0.1;
rb = 0.2;

dt = 0.25e-11;
delta = 3e-3;% PML width.

x0 = delta*I/2; % x0.
y0 = delta*(PMLw+J/2); % y0.
