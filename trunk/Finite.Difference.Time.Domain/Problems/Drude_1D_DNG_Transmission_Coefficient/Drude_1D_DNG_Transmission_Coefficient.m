clc
clear all
close all

% Simulation parameters.
SIZE = 4*1024; % No. of spatial steps
SlabLeft = round(SIZE/3); % Location of left end of Slab.
SlabRight = round(2*SIZE/3); % Location of right end of Slab
MaxTime = 4*SIZE; % No. of time steps
PulseWidth = round(SIZE/8); % Controls width of Gaussian Pulse
td = PulseWidth; % Temporal delay in pulse.
source = 10; % Location of source
SnapshotInterval = 16; % Amount of time delay between snaps.

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
fmax = 1/(2*dt)
w = 2*pi*f;
k0 = w/c; % Free space wave number.

% Choice of source.
% 1. Gaussian 2. Sine wave 3. Ricker wavelet
SourceChoice = 1;
if SourceChoice == 3
    fp = f; % Peak frequency
    dr = PulseWidth*dt*2; % Delay
end

% Initialization.
Ex = zeros(SIZE, 3); % x-component of E-field
Dx = zeros(SIZE, 3); % x-component of D
Hy = zeros(SIZE, 3); % y-component of H-field
By = zeros(SIZE, 3); % y-component of B

% Incident and Transmitted Fields.
Exi = zeros(MaxTime, 1);
Ext = zeros(MaxTime, 1);
Extt = zeros(MaxTime, 1);
x1 = SlabLeft+1; % Position of observation.

% Refractive Index calculations.
Z1 = SlabLeft + 50;
z1 = Z1*dz;
Z2 = SlabLeft + 60;
z2 = Z2*dz;
Exz1 = zeros(MaxTime, 1);
Exz2 = zeros(MaxTime, 1);

% Power calculations.
Exip = zeros(MaxTime, 1);
Hyip = zeros(MaxTime, 1);
Extp = zeros(MaxTime, 1);
Hytp = zeros(MaxTime, 1);

einf = ones(SIZE,1);
einf(SlabLeft:SlabRight) = 1; % einf(Drude) or er in slab.
uinf = ones(SIZE,1);
uinf(SlabLeft:SlabRight) = 1; % uinf(Drude) or ur in slab.
wpesq = zeros(SIZE,1);
wpesq(SlabLeft:SlabRight) = 2*w^2; % DNG(Drude) value of wpe squared in slab.
wpmsq = zeros(SIZE,1);
wpmsq(SlabLeft:SlabRight) = 2*w^2; % DNG(Drude) value of wpm squared in slab.
ge = zeros(SIZE,1);
ge(SlabLeft:SlabRight) = w/16; % Electric collision frequency in slab.
gm = zeros(SIZE,1);
gm(SlabLeft:SlabRight) = w/16; % Magnetic collision frequency in slab.

a0 = (4*dt^2)./(e0*(4*einf+dt^2*wpesq+2*dt*einf.*ge));
a = (1/dt^2)*a0;
b = (1/(2*dt))*ge.*a0;
c = (e0/dt^2)*einf.*a0;
d = (-1*e0/4).*wpesq.*a0;
e = (1/(2*dt))*e0*einf.*ge.*a0;
am0 = (4*dt^2)./(u0*(4*uinf+dt^2*wpmsq+2*dt*uinf.*gm));
am = (1/dt^2)*am0;
bm = (1/(2*dt))*gm.*am0;
cm = (u0/dt^2)*uinf.*am0;
dm = (-1*u0/4).*wpmsq.*am0;
em = (1/(2*dt))*u0*uinf.*gm.*am0;

ExSnapshots = zeros(SIZE, MaxTime/SnapshotInterval); % Data for plotting.
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
    elseif SourceChoice == 3
    Ex(source,n2) = Ex(source,n2) + (1-2*(pi*fp*(q*dt-dr))^2)*exp(-1*(pi*fp*(q*dt-dr))^2) * Sc;
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
    Hy(:,n2) = am.*(By(:,n2)-2*By(:,n1)+By(:,3))+bm.*(By(:,n2)-By(:,3))+cm.*(2*Hy(:,n1)-Hy(:,3))+dm.*(2*Hy(:,n1)+Hy(:,3))+em.*(Hy(:,3));
    
    % ABC for H at SIZE.
    Hy(SIZE,n2) = Hy(SIZE-1,n1) + (Sc-1)/(Sc+1)*(Hy(SIZE-1,n2) - Hy(SIZE,n1) );
    By(SIZE,n2) = u0*Hy(SIZE,n2);

    % Calculation of Ex using updated difference equation for Ex. This is time step q+1/2.
    Dx(2:SIZE,n2) = Dx(2:SIZE, n1) + ( dt/(dz)*(Hy(1:SIZE-1, n2) - Hy(2:SIZE, n2)) );
    Ex(:,n2) = a.*(Dx(:,n2)-2*Dx(:,n1)+Dx(:,3))+b.*(Dx(:,n2)-Dx(:,3))+c.*(2*Ex(:,n1)-Ex(:,3))+d.*(2*Ex(:,n1)+Ex(:,3))+e.*(Ex(:,3));
    
    % ABC for E at 1.
    Ex(1,n2) = Ex(2,n1) + (Sc-1)/(Sc+1)*(Ex(2,n2) - Ex(2,n1));
    Dx(1,n2) = e0*Ex(1,n2);
    
    % Source.
    if SourceChoice == 1
    Ex(source,n2) = Ex(source,n2) + exp( -1*((q-td)/(PulseWidth/4))^2 ) * Sc;
    elseif SourceChoice == 2
    Ex(source,n2) = Ex(source,n2) + sin(2*pi*f*(q)*dt) * Sc;
    elseif SourceChoice == 3
    Ex(source,n2) = Ex(source,n2) + (1-2*(pi*fp*(q*dt-dr))^2)*exp(-1*(pi*fp*(q*dt-dr))^2) * Sc;
    end
    Dx(source,n2) = e0*Ex(source,n2);

    if mod(q,SnapshotInterval) == 0
        ExSnapshots(:,frame) = Ex(:,n2);
        frame=frame+1;
    end
    
    Ext(q+1) = Ex(x1,n2);
    Extt(q+1) = Ex(SlabRight+10,n2);
    
    % Fields for calculation of refractive index.
    Exz1(q+1) = Ex(Z1, n2);
    Exz2(q+1) = Ex(Z2, n2);
    % Fields for power calculations.
    Exip(q+1) = Ex(SlabLeft-1, n2);
    Hyip(q+1) = Hy(SlabLeft-1, n2);
    Extp(q+1) = Ex(SlabLeft+1, n2);
    Hytp(q+1) = Hy(SlabLeft+1, n2);
    
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
fspan = 100;                  % Points to plot in frequency domain

figure(1)
subplot(211)
plot(Fs*t, Exi, 'LineWidth', 2.0, 'Color', 'b')
set(gca, 'FontSize', 10, 'FontWeight', 'b')
title('Incident Field', 'FontSize', 12, 'FontWeight', 'b')
xlabel('time', 'FontSize', 11, 'FontWeight', 'b')
grid on
figure(2)
subplot(211)
plot(Fs*t, Ext, 'LineWidth', 2.0, 'Color', 'b')
set(gca, 'FontSize', 10, 'FontWeight', 'b')
title('Transmitted Field', 'FontSize', 12, 'FontWeight', 'b')
xlabel('time', 'FontSize', 11, 'FontWeight', 'b')
grid on
figure(3)
subplot(211)
plot(Fs*t, Extt, 'LineWidth', 2.0, 'Color', 'b')
set(gca, 'FontSize', 10, 'FontWeight', 'b')
title('Transmitted Field Beyond Slab', 'FontSize', 12, 'FontWeight', 'b')
xlabel('time', 'FontSize', 11, 'FontWeight', 'b')
grid on

NFFT = 2^nextpow2(L); % Next power of 2 from length of Exi
% Incident and Transmitted fields.
EXI = fft(Exi,NFFT)/L;
EXT = fft(Ext,NFFT)/L;
EXTT = fft(Extt,NFFT)/L;
% Refractive index calculations.
EXZ1 = fft(Exz1,NFFT)/L;
EXZ2 = fft(Exz2,NFFT)/L;
% Fields for power calculations.
EXIP = fft(Exip,NFFT)/L;
HYIP = fft(Hyip,NFFT)/L;
EXTP = fft(Extp,NFFT)/L;
HYTP = fft(Hytp,NFFT)/L;
SI = EXIP.*conj(HYIP);
ST = EXTP.*conj(HYTP);
SIave = 0.5*real(SI);
STave = 0.5*real(ST);
f = Fs/2*linspace(0,1,NFFT/2+1);

% Plot single-sided amplitude spectrum.
figure(1)
subplot(212)
EXIp = 2*abs(EXI(1:NFFT/2+1));
plot(f(1:fspan), EXIp(1:fspan), 'LineWidth', 2.0, 'Color', 'r')
set(gca, 'FontSize', 10, 'FontWeight', 'b')
title('Single-Sided Amplitude Spectrum of Exi(t)', 'FontSize', 12, 'FontWeight', 'b')
xlabel('Frequency (Hz)', 'FontSize', 11, 'FontWeight', 'b')
ylabel('|EXI(f)|', 'FontSize', 11, 'FontWeight', 'b')
grid on
figure(2)
subplot(212)
EXTp = 2*abs(EXT(1:NFFT/2+1));
plot(f(1:fspan), EXTp(1:fspan), 'LineWidth', 2.0, 'Color', 'r')
set(gca, 'FontSize', 10, 'FontWeight', 'b')
title('Single-Sided Amplitude Spectrum of Ext(t)', 'FontSize', 12, 'FontWeight', 'b')
xlabel('Frequency (Hz)', 'FontSize', 11, 'FontWeight', 'b')
ylabel('|EXT(f)|', 'FontSize', 11, 'FontWeight', 'b')
grid on
figure(3)
subplot(212)
EXTTp = 2*abs(EXTT(1:NFFT/2+1));
plot(f(1:fspan), EXTTp(1:fspan), 'LineWidth', 2.0, 'Color', 'r')
set(gca, 'FontSize', 10, 'FontWeight', 'b')
title('Single-Sided Amplitude Spectrum of Extt(t)', 'FontSize', 12, 'FontWeight', 'b')
xlabel('Frequency (Hz)', 'FontSize', 11, 'FontWeight', 'b')
ylabel('|EXT(f)|', 'FontSize', 11, 'FontWeight', 'b')
grid on

% Transmission Coefficient.
figure(4)
subplot(211)
TAU = abs(EXT(1:NFFT/2+1)./EXI(1:NFFT/2+1));
plot(f(1:fspan), TAU(1:fspan), 'LineWidth', 2.0, 'Color', 'b')
set(gca, 'FontSize', 10, 'FontWeight', 'b')
title('Transmission Coefficient', 'FontSize', 12, 'FontWeight', 'b')
xlabel('Frequency (Hz)', 'FontSize', 11, 'FontWeight', 'b')
ylabel('|EXT(f)/EXI(f)|', 'FontSize', 11, 'FontWeight', 'b')
axis([-1 1 -2 2])
axis 'auto x'
grid on
subplot(212)
plot(f(1:fspan), 1-TAU(1:fspan), 'LineWidth', 2.0, 'Color', 'b')
set(gca, 'FontSize', 10, 'FontWeight', 'b')
title('Reflection Coefficient', 'FontSize', 12, 'FontWeight', 'b')
xlabel('Frequency (Hz)', 'FontSize', 11, 'FontWeight', 'b')
ylabel('1-|EXT(f)/EXI(f)|', 'FontSize', 11, 'FontWeight', 'b')
axis([-1 1 -2 2])
axis 'auto x'
grid on

% Refractive Index calculations.
nFDTD = (1/(1i*k0*(z1-z2))).*log(EXZ2(1:NFFT/2+1)./EXZ1(1:NFFT/2+1));
figure(5)
subplot(211)
plot(f(1:fspan), real(nFDTD(1:fspan)), 'LineWidth', 2.0, 'Color', 'b');
set(gca, 'FontSize', 10, 'FontWeight', 'b')
title('Refractive index re(n)', 'FontSize', 12, 'FontWeight', 'b')
xlabel('Frequency (Hz)', 'FontSize', 11)
ylabel('re(n)', 'FontSize', 11)
grid on
subplot(212)
plot(f(1:fspan), imag(nFDTD(1:fspan)), 'LineWidth', 2.0, 'Color', 'r');
set(gca, 'FontSize', 10, 'FontWeight', 'b')
title('Refractive index im(n)', 'FontSize', 12, 'FontWeight', 'b')
xlabel('Frequency (Hz)', 'FontSize', 11, 'FontWeight', 'b')
ylabel('im(n)', 'FontSize', 11, 'FontWeight', 'b')
grid on

% Power calculations
figure(6)
subplot(211)
plot(f(1:fspan), SIave(1:fspan), 'LineWidth', 2.0, 'Color', 'b');
set(gca, 'FontSize', 10, 'FontWeight', 'b')
title('Incident average power', 'FontSize', 12, 'FontWeight', 'b')
xlabel('Frequency (Hz)', 'FontSize', 11)
ylabel('Power', 'FontSize', 11)
grid on
subplot(212)
plot(f(1:fspan), STave(1:fspan), 'LineWidth', 2.0, 'Color', 'r');
set(gca, 'FontSize', 10, 'FontWeight', 'b')
title('Transmitted average power', 'FontSize', 12, 'FontWeight', 'b')
xlabel('Frequency (Hz)', 'FontSize', 11, 'FontWeight', 'b')
ylabel('Power', 'FontSize', 11, 'FontWeight', 'b')
grid on

% Simulation animation.
for i=1:frame-1
    figure (7)
    % Scatterer boundaries.
    hold off
    plot([SlabLeft SlabLeft], [-1 1], 'Color', 'r');
    hold on
    plot([SlabRight SlabRight], [-1 1], 'Color', 'r');
    plot(ExSnapshots(:,i), 'LineWidth', 2.0, 'Color', 'b');
    set(gca, 'FontSize', 10, 'FontWeight', 'b')
    axis([0 SIZE -1 1])
    title('Time Domain Simulation', 'FontSize', 12, 'FontWeight', 'b')
    xlabel('Spatial step (k)', 'FontSize', 11, 'FontWeight', 'b')
    ylabel('Electric field (Ex)', 'FontSize', 11, 'FontWeight', 'b')
    grid on
end
