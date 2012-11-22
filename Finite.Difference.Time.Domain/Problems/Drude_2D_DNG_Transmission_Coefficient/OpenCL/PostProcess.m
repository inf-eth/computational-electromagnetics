% Plots field snapshots as a movie. Taken from appendix of
% Understanding FDTD, J. B Schneider.
% Input files are assumed to be in binary format.
clc
clear all
Precision = 'double';

fidp = fopen ('./FieldData/Parameters.smp', 'r', 'l');
if fidp == -1
    return;
end
dt = fread (fidp, 1, Precision);
k0 = fread (fidp, 1, Precision);
y1 = fread (fidp, 1, Precision);
y2 = fread (fidp, 1, Precision);

datap = fread (fidp, 9, 'uint');
fclose (fidp);
I = datap(1)
J = datap(2)
PMLw = datap(3)
MaxTime = datap(4)
SnapshotResolution = datap(5)
SnapshotInterval = datap(6)
SlabLeft = datap(7)
SlabRight = datap(8)
simTime = datap(9)
Size = [I J+2*PMLw];    % Spatial size or width w.

fid = fopen ('./FieldData/Ezi.fdt', 'r', 'l');
if fid == -1
    return;
end
Ezi = fread(fid, MaxTime, Precision);
fclose (fid);

fid = fopen ('./FieldData/Ezt.fdt', 'r', 'l');
if fid == -1
    return;
end
Ezt = fread(fid, MaxTime, Precision);
fclose (fid);

fid = fopen ('./FieldData/Eztt.fdt', 'r', 'l');
if fid == -1
    return;
end
Eztt = fread(fid, MaxTime, Precision);
fclose (fid);

fid = fopen ('./FieldData/Ezy1.fdt', 'r', 'l');
if fid == -1
    return;
end
Ezy1 = fread(fid, MaxTime, Precision);
fclose (fid);

fid = fopen ('./FieldData/Ezy2.fdt', 'r', 'l');
if fid == -1
    return;
end
Ezy2 = fread(fid, MaxTime, Precision);
fclose (fid);

% Postprocessing.
Fs = 1/dt;                    % Sampling frenuency
T = dt;                       % Sample time
L = length(Ezi);              % Length of signal
t = (0:L-1)*T;                % Time vector
fspan = 100;                  % Points to plot in frenuency domain

figure(1)
subplot(211)
plot(Fs*t, Ezi, 'LineWidth', 2.0, 'Color', 'b')
set(gca, 'FontSize', 10, 'FontWeight', 'b')
title('Incident Field', 'FontSize', 12, 'FontWeight', 'b')
xlabel('time step (n)', 'FontSize', 11, 'FontWeight', 'b')
grid on
figure(2)
subplot(211)
plot(Fs*t, Ezt, 'LineWidth', 2.0, 'Color', 'b')
set(gca, 'FontSize', 10, 'FontWeight', 'b')
title('Transmitted Field', 'FontSize', 12, 'FontWeight', 'b')
xlabel('time step (n)', 'FontSize', 11, 'FontWeight', 'b')
grid on
figure(3)
subplot(211)
plot(Fs*t, Eztt, 'LineWidth', 2.0, 'Color', 'b')
set(gca, 'FontSize', 10, 'FontWeight', 'b')
title('Transmitted Field Beyond Slab', 'FontSize', 12, 'FontWeight', 'b')
xlabel('time step (n)', 'FontSize', 11, 'FontWeight', 'b')
grid on

NFFT = 2^nextpow2(L); % Next power of 2 from length of Exi
% Incident and Transmitted fields.
EZI = fft(Ezi,NFFT)/L;
EZT = fft(Ezt,NFFT)/L;
EZTT = fft(Eztt,NFFT)/L;
% Refractive index calculations.
EZY1 = fft(Ezy1,NFFT)/L;
EZY2 = fft(Ezy2,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);

% Plot single-sided amplitude spectrum.
figure(1)
subplot(212)
EZIp = 2*abs(EZI(1:NFFT/2+1));
plot(f(1:fspan), EZIp(1:fspan), 'LineWidth', 2.0, 'Color', 'r')
set(gca, 'FontSize', 10, 'FontWeight', 'b')
title('Single-Sided Amplitude Spectrum of Ezi(t)', 'FontSize', 12, 'FontWeight', 'b')
xlabel('Frenuency (Hz)', 'FontSize', 11, 'FontWeight', 'b')
ylabel('|EZI(f)|', 'FontSize', 11, 'FontWeight', 'b')
grid on
figure(2)
subplot(212)
EZTp = 2*abs(EZT(1:NFFT/2+1));
plot(f(1:fspan), EZTp(1:fspan), 'LineWidth', 2.0, 'Color', 'r')
set(gca, 'FontSize', 10, 'FontWeight', 'b')
title('Single-Sided Amplitude Spectrum of Ezt(t)', 'FontSize', 12, 'FontWeight', 'b')
xlabel('Frenuency (Hz)', 'FontSize', 11, 'FontWeight', 'b')
ylabel('|EZT(f)|', 'FontSize', 11, 'FontWeight', 'b')
grid on
figure(3)
subplot(212)
EZTTp = 2*abs(EZTT(1:NFFT/2+1));
plot(f(1:fspan), EZTTp(1:fspan), 'LineWidth', 2.0, 'Color', 'r')
set(gca, 'FontSize', 10, 'FontWeight', 'b')
title('Single-Sided Amplitude Spectrum of Eztt(t)', 'FontSize', 12, 'FontWeight', 'b')
xlabel('Frenuency (Hz)', 'FontSize', 11, 'FontWeight', 'b')
ylabel('|EZT(f)|', 'FontSize', 11, 'FontWeight', 'b')
grid on

% Transmission Coefficient.
figure(4)
subplot(211)
TAU = abs(EZT(1:NFFT/2+1)./EZI(1:NFFT/2+1));
plot(f(1:fspan), TAU(1:fspan), 'LineWidth', 2.0, 'Color', 'b')
set(gca, 'FontSize', 10, 'FontWeight', 'b')
title('Transmission Coefficient', 'FontSize', 12, 'FontWeight', 'b')
xlabel('Frenuency (Hz)', 'FontSize', 11, 'FontWeight', 'b')
ylabel('|EZT(f)/EZI(f)|', 'FontSize', 11, 'FontWeight', 'b')
axis([-1 1 -2 2])
axis 'auto i'
grid on
subplot(212)
plot(f(1:fspan), 1-TAU(1:fspan), 'LineWidth', 2.0, 'Color', 'b')
set(gca, 'FontSize', 10, 'FontWeight', 'b')
title('Reflection Coefficient', 'FontSize', 12, 'FontWeight', 'b')
xlabel('Frenuency (Hz)', 'FontSize', 11, 'FontWeight', 'b')
ylabel('1-|EZT(f)/EZI(f)|', 'FontSize', 11, 'FontWeight', 'b')
axis([-1 1 -2 2])
axis 'auto i'
grid on

% Refractive Index calculations.
nFDTD = (1/(1i*k0*(y1-y2))).*log(EZY2(1:NFFT/2+1)./EZY1(1:NFFT/2+1));
figure(5)
subplot(211)
plot(f(1:fspan), real(nFDTD(1:fspan)), 'LineWidth', 2.0, 'Color', 'b');
set(gca, 'FontSize', 10, 'FontWeight', 'b')
title('Refractive index re(n)', 'FontSize', 12, 'FontWeight', 'b')
xlabel('Frenuency (Hz)', 'FontSize', 11)
ylabel('re(n)', 'FontSize', 11)
grid on
subplot(212)
plot(f(1:fspan), imag(nFDTD(1:fspan)), 'LineWidth', 2.0, 'Color', 'r');
set(gca, 'FontSize', 10, 'FontWeight', 'b')
title('Refractive index im(n)', 'FontSize', 12, 'FontWeight', 'b')
xlabel('Frenuency (Hz)', 'FontSize', 11, 'FontWeight', 'b')
ylabel('im(n)', 'FontSize', 11, 'FontWeight', 'b')
grid on

basename = './FieldData/Ez';
for i=1:simTime
    filename = sprintf ('%s%d.fdt', basename, i);
    fid = fopen (filename, 'r', 'l');
    if fid == -1
        return;
    end
    
    EzSnapshot = fread (fid, Size, Precision);
    SizeS=size(EzSnapshot);
    figure (6)
    mesh ( EzSnapshot (1:SnapshotResolution:SizeS(1), 1:SnapshotResolution:SizeS(2)) );
    view (4, 4)
    xlim([0 SizeS(2)])
    ylim([0 SizeS(1)])
    zlim([-1 1])
    %caxis([-0.1 0.6])
    xlabel ('j-axis')
    ylabel ('i-axis')
    %colorbar
    
    figure (7)
    mesh ( EzSnapshot (1:SnapshotResolution:SizeS(1), 1:SnapshotResolution:SizeS(2)) );
    view (0, 90)
    xlim([0 SizeS(2)])
    ylim([0 SizeS(1)])
    zlim([-10 10])
    %caxis([-0.1 0.6])
    xlabel ('j-axis')
    ylabel ('i-axis')
    %colorbar
    
    fclose (fid);
end
