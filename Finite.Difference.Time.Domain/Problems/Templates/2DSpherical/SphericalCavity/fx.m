function [ output_args ] = fx( w, s0 )
%FX Summary of this function goes here
%   Detailed explanation goes here
Ra = 150e-6;

e0 = (1e-9)/(36*pi);
u0 = (1e-7)*4*pi;

ed = e0;
ud = u0;
nd = (ud/ed)^0.5;

erinf = 1;
sc = s0;
ec = e0*(erinf-1j*(sc/(w*e0)));
uc = u0;
nc = (uc/ec)^0.5;

Bd = w*(ud*ed)^0.5;
Bc = w*(uc*ec)^0.5;
output_args = nd*(jnu(0,Bd*Ra)/jnu(1,Bd*Ra)-1/Bd*Ra)-nc*(h2nu(0,Bc*Ra)/h2nu(1,Bc*Ra)-1/Bc*Ra);
end

