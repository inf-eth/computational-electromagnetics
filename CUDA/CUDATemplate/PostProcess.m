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
z1 = fread (fidp, 1, Precision);
z2 = fread (fidp, 1, Precision);
datap = fread (fidp, 6, 'uint');
fclose (fidp);
Size = datap(1)
MaxTime = datap(2)
SnapshotInterval = datap(3)
SlabLeft = datap(4)
SlabRight = datap(5)
simTime = datap(6)
size = [Size 1];    % Spatial size or width w.

fid = fopen ('./FieldData/Exi.fdt', 'r', 'l');
if fid == -1
    return;
end
Exi = fread(fid, MaxTime, Precision);
fclose (fid);

fid = fopen ('./FieldData/Ext.fdt', 'r', 'l');
if fid == -1
    return;
end
Ext = fread(fid, MaxTime, Precision);
fclose (fid);

fid = fopen ('./FieldData/Extt.fdt', 'r', 'l');
if fid == -1
    return;
end
Extt = fread(fid, MaxTime, Precision);
fclose (fid);

fid = fopen ('./FieldData/Exz1.fdt', 'r', 'l');
if fid == -1
    return;
end
Exz1 = fread(fid, MaxTime, Precision);
fclose (fid);

fid = fopen ('./FieldData/Exz2.fdt', 'r', 'l');
if fid == -1
    return;
end
Exz2 = fread(fid, MaxTime, Precision);
fclose (fid);

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

basename = './FieldData/Ex';
for i=1:simTime
    filename = sprintf ('%s%d.fdt', basename, i);
    fid = fopen (filename, 'r', 'l');
    if fid == -1
        return;
    end

    ExSnapshot = fread (fid, size, Precision);
    figure (6)
    % Scatterer boundaries.
    hold off
    plot([SlabLeft SlabLeft], [-1 1], 'Color', 'r');
    hold on
    plot([SlabRight SlabRight], [-1 1], 'Color', 'r');
    plot(ExSnapshot, 'LineWidth', 2.0, 'Color', 'b');
    set(gca, 'FontSize', 10, 'FontWeight', 'b')
    axis([0 Size -1 1])
    title('Time Domain Simulation', 'FontSize', 12, 'FontWeight', 'b')
    xlabel('Spatial step (k)', 'FontSize', 11, 'FontWeight', 'b')
    ylabel('Electric field (Ex)', 'FontSize', 11, 'FontWeight', 'b')
    grid on
    fclose (fid);
end
