clc
clear all

% Simulation parameters.
SIZE = 2024; % No. of spatial steps
SlabLeft = round(SIZE/3); % Location of left end of Slab.
SlabRight = round(2*SIZE/3); % Location of right end of Slab
MaxTime = 2*SIZE; % No. of time steps
PulseWidth = round(SIZE/8); % Controls width of Gaussian Pulse
td = PulseWidth; % Temporal delay in pulse.
source = 10; % Location of source

% Constants.
c = 3e8;
pi = 3.141592654;
e0 = (1e-9)/(36*pi);
u0 = (1e-7)*4*pi;

dt = 0.5e-11;
dz = 3e-3;
Sc = c * dt/dz

l = PulseWidth*dz;
f = c/l
w = 2*pi*f;
k0 = w/c; % Free space wave number.

% Choice of source.
% 1. Gaussian 2. Sine wave.
SourceChoice = 1;

% Initialization.
Ex = zeros(SIZE, 3); % x-component of E-field
Dx = zeros(SIZE, 3); % x-component of D
Hy = zeros(SIZE, 3); % y-component of H-field
By = zeros(SIZE, 3); % y-component of B

% Incident and Transmitted Fields.
Exi = zeros(MaxTime);
Ext = zeros(MaxTime);
x1 = SlabLeft+1; % Position of observation.

% Refractive Index calculations.
Z1 = SlabLeft + 50;
z1 = Z1*dz;
Z2 = SlabLeft + 60;
z2 = Z2*dz;
Exz1 = zeros(MaxTime);
Exz2 = zeros(MaxTime);

einf = ones(SIZE,1);
einf(SlabLeft:SlabRight) = 1; % einf(Drude) or er in slab.
uinf = ones(SIZE,1);
uinf(SlabLeft:SlabRight) = 1; % uinf(Drude) or ur in slab.
wpsq = zeros(SIZE,1);
wpsq(SlabLeft:SlabRight) = 2*w^2; % DNG(Drude) value of wp squared.
wpmsq = zeros(SIZE,1);
wpmsq(SlabLeft:SlabRight) = 2*w^2; % DNG(Drude) value of wpm squared.

a = 4./(e0*(4*einf+wpsq*dt^2));
b = -1*e0*a.*einf;
c = -1*(wpsq*dt^2)./(4*einf+(wpsq*dt^2));
am = 4./(u0*(4*uinf+wpmsq*dt^2));
bm = -1*u0*am.*uinf;
cm = -1*(wpmsq*dt^2)./(4*uinf+(wpmsq*dt^2));

ExSnapshots = zeros(SIZE, MaxTime); % Data for plotting.
frame = 1;

n1 = 1;
n2 = 2;
% Outer loop for time-stepping.
tic
% Test loop for incident field in free space.
for q = 0:MaxTime
   
    % Calculation of Hy using update difference equation for Hy. This is time step q.
    Hy(1:SIZE-1,n2) = Hy(1:SIZE-1,n1) + ( ( Ex(1:SIZE-1,n1) - Ex(2:SIZE,n1) ) * dt/(u0*dz) );
        
    % ABC for H at SIZE.
    Hy(SIZE,n2) = Hy(SIZE-1,n1) + (Sc-1)/(Sc+1)*(Hy(SIZE-1,n2) - Hy(SIZE,n1) );
 
    % Calculation of Ex using updated difference equation for Ex. This is time step q+1/2.
    Ex(2:SIZE,n2) = Ex(2:SIZE, n1) + ( dt/(e0*dz)*(Hy(1:SIZE-1, n2) - Hy(2:SIZE, n2)) );
    
    % ABC for E at 1.
    Ex(1,n2) = Ex(2,n1) + (Sc-1)/(Sc+1)*(Ex(2,n2) - Ex(2,n1));
    
    % Source.
    if SourceChoice == 1
    Ex(source,n2) = Ex(source,n2) + exp( -1*((q-td)/(PulseWidth/4))^2 ) * Sc;
    elseif SourceChoice == 2
    Ex(source,n2) = Ex(source,n2) + sin(2*pi*f*(q)*dt) * Sc;
    end
    
    Exi(q+1) = Ex(x1,n2); % Incident field is left of slab.
    
    temp = n1;
    n1 = n2;
    n2 = temp;
end
% Reinitialization of fields for actual simulation.
Ex = zeros(SIZE, 3); % x-component of E-field
Hy = zeros(SIZE, 3); % y-component of H-field
% Actual simulation with scatterer.
for q = 0:MaxTime
    
    % Storing past fields.
    Ex(:,3) = Ex(:,n2);
    Dx(:,3) = Dx(:,n2);
    Hy(:,3) = Hy(:,n2);
    By(:,3) = By(:,n2);
    
    % Calculation of Hy using update difference equation for Hy. This is time step q.
    By(1:SIZE-1,n2) = By(1:SIZE-1,n1) + ( ( Ex(1:SIZE-1,n1) - Ex(2:SIZE,n1) ) * dt/(dz) );
    Hy(:,n2) = am.*(By(:,n2)-2*By(:,n1)+By(:,3))+bm.*(-2*Hy(:,n1)+Hy(:,3))+cm.*(2*Hy(:,n1)+Hy(:,3));
    
    % ABC for H at SIZE.
    Hy(SIZE,n2) = Hy(SIZE-1,n1) + (Sc-1)/(Sc+1)*(Hy(SIZE-1,n2) - Hy(SIZE,n1) );
    By(SIZE,n2) = u0*Hy(SIZE,n2);

    % Calculation of Ex using updated difference equation for Ex. This is time step q+1/2.
    Dx(2:SIZE,n2) = Dx(2:SIZE, n1) + ( dt/(dz)*(Hy(1:SIZE-1, n2) - Hy(2:SIZE, n2)) );
    Ex(:,n2) = a.*(Dx(:,n2)-2*Dx(:,n1)+Dx(:,3))+b.*(-2*Ex(:,n1)+Ex(:,3))+c.*(2*Ex(:,n1)+Ex(:,3));
    
    % ABC for E at 1.
    Ex(1,n2) = Ex(2,n1) + (Sc-1)/(Sc+1)*(Ex(2,n2) - Ex(2,n1));
    Dx(1,n2) = e0*Ex(1,n2);
    
    % Source.
    if SourceChoice == 1
    Ex(source,n2) = Ex(source,n2) + exp( -1*((q-td)/(PulseWidth/4))^2 ) * Sc;
    elseif SourceChoice == 2
    Ex(source,n2) = Ex(source,n2) + sin(2*pi*f*(q)*dt) * Sc;
    end
    Dx(source,n2) = e0*Ex(source,n2);

    ExSnapshots(:,frame) = Ex(:,n2);
    frame=frame+1;
    
    Ext(q+1) = Ex(x1,n2);
    
    % Fields for calculation of refractive index.
    Exz1(q+1) = Ex(Z1, n2);
    Exz2(q+1) = Ex(Z2, n2);
    
    temp = n1;
    n1 = n2;
    n2 = temp;
end
toc
% Postprocessing.
Fs = 1/dt;                    % Sampling frequency
T = dt;                       % Sample time
L = length(Exi);              % Length of signal
t = (0:L-1)*T;                % Time vector

figure(1)
subplot(211)
plot(Fs*t,Exi)
title('Incident Field')
xlabel('time')
figure(2)
subplot(211)
plot(Fs*t,Ext)
title('Transmitted Field')
xlabel('time')

NFFT = 2^nextpow2(L); % Next power of 2 from length of Exi
EXI = fft(Exi,NFFT)/L;
EXT = fft(Ext,NFFT)/L;
EXZ1 = fft(Exz1,NFFT)/L;
EXZ2 = fft(Exz2,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);

% Plot single-sided amplitude spectrum.
figure(1)
subplot(212)
EXIp = 2*abs(EXI(1:NFFT/2+1));
plot(f(1:50),EXIp(1:50)) 
title('Single-Sided Amplitude Spectrum of Exi(t)')
xlabel('Frequency (Hz)')
ylabel('|EXI(f)|')
figure(2)
subplot(212)
EXTp = 2*abs(EXT(1:NFFT/2+1));
plot(f(1:50),EXTp(1:50)) 
title('Single-Sided Amplitude Spectrum of Ext(t)')
xlabel('Frequency (Hz)')
ylabel('|EXT(f)|')

% Transmission Coefficient.
figure(3)
subplot(211)
TAU = abs(EXT(1:NFFT/2+1)./EXI(1:NFFT/2+1));
plot(f(1:50),TAU(1:50))
title('Transmission Coefficient')
xlabel('Frequency (Hz)')
ylabel('|EXT(f)/EXI(f)|')
axis([-1 1 -2 2])
axis 'auto x'
subplot(212)
plot(f(1:50),1-TAU(1:50))
title('Reflection Coefficient')
xlabel('Frequency (Hz)')
ylabel('1-|EXT(f)/EXI(f)|')
axis([-1 1 -2 2])
axis 'auto x'

% Refractive Index calculations.
nFDTD = (1/(1i*k0*(z1-z2))).*log(EXZ2(1:NFFT/2+1)./EXZ1(1:NFFT/2+1));
figure(4)
subplot(211)
plot(f(1:50),real(nFDTD(1:50)));
title('Refractive index re(n)')
xlabel('Frequency (Hz)')
ylabel('re(n)')
subplot(212)
plot(f(1:50),imag(nFDTD(1:50)));
title('Refractive index im(n)')
xlabel('Frequency (Hz)')
ylabel('im(n)')

% Simulation animation.
for i=1:frame-1
    figure (5)
    plot ( ExSnapshots(:,i) )
    axis([0 SIZE -0.6 0.6])
    xlabel('Spatial step (k)')
    ylabel('Electric field (Ex)')
end