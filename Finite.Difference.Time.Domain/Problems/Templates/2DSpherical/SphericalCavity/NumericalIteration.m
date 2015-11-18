clc
clear all

Ra = 150e-6;
s0 = 1e2;

e0 = (1e-9)/(36*pi);
u0 = (1e-7)*4*pi;

ed = e0;
ud = u0;
nd = (ud/ed)^0.5;
wI = 2.74370/(Ra*(ud*ed)^0.5);

erinf = 1;
sc = s0;
w0 = wI
w = w0 % ?
ec = e0*(erinf-1j*(sc/(w*e0)))
uc = u0
nc = (uc/ec)^0.5;

Bd = w0*(ud*ed)^0.5
Bc = w0*(uc*ec)^0.5

%x = jnu(1, 2.74370)
%y = jnup(1, 2.74370)
x0 = 2.74370;

% Newton-Raphston method: x1 = x0-f(x0)/f'(x0);
n = 100;
wn = w0;
delta = 0.025;
tolerance = 1e-15;
t0 = 0;
tf = 0;
for i=1:n
    temp = wn;
    t0 = tf;
    wn = wn-fx(wn,s0)/fxp(wn,s0,delta);
    tf = abs(wn-temp)/abs(wn);
    if tf < tolerance
        break;
    end
    if tf > t0 && i ~= 1
        wn = temp;
        break;
    end
end
i
wn
abs(wn)
residual = fx(wn,s0)

% Plot data
ec = e0*(erinf-1j*(sc/(wn*e0)));
nc = (uc/ec)^0.5;

Bd = wn*(ud*ed)^0.5;
Bc = wn*(uc*ec)^0.5;

% ZTM
r1 = 0:1e-6:150e-6;
r2 = 150e-6:1e-6:450e-6;
ZTM1 = 1j*nd*(jnu(0,Bd*r1)./jnu(1,Bd*r1)-1./Bd*r1);
ZTM2 = 1j*nc*(h2nu(0,Bc*r2)./h2nu(1,Bc*r2)-1./Bc*r2);

figure(1)
subplot(121)
plot(r1/Ra,real(ZTM1))
hold on
plot(r2/Ra,real(ZTM2))
axis([0 inf 1 10000])
set(gca,'yscale','log');
subplot(122)
plot(r1/Ra,imag(ZTM1))
hold on
plot(r2/Ra,imag(ZTM2))
axis([0 inf 0.1 10000])
set(gca,'yscale','log');