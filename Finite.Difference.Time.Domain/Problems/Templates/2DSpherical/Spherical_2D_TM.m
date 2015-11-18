clear all
clc
% Simulation parameters.
Scalar = 1;
PMLw = 30;
SIZEr = 256*Scalar+PMLw; % No. of spatial steps
SIZEphi = SIZEr-2*PMLw; % No. of spatial steps
MaxTime = SIZEr*4*Scalar; % No. of time steps
PulseWidth = round(SIZEr/8); % Controls width of Gaussian Pulse
td = PulseWidth; % Temporal delay in pulse.
imp0 = 377.0; % Impedence of free space
sourceX = 2; % Location of sourcer
sourceY = 1:SIZEphi; % Location of sourcephi
Ra = 150e-6;

% Choice of source.
% 1. Gaussian 2. Sine wave 3. Ricker wavelet
SourceChoice = 2;

SaveFields = 0; % 0. No, 1. Yes.
SnapshotResolution = 1; % Snapshot resolution. 1 is best.
SnapshotInterval = 8*Scalar; % Amount of time delay between snaps.

% Constants.
pi = 3.141592654;
e0 = (1e-9)/(36*pi);
u0 = (1e-7)*4*pi;
c = 1/sqrt(e0*u0);

dt = 0.3e-14;
dr = 1e-6;
dphi = 2*pi/SIZEphi;
Sc = c*dt/dr

l = PulseWidth*dr;
f = 703.6e9%0.8727e12%c/(1*l)
fmax = 1/(2*dt)
w = 2*pi*f;
k0 = w/c; % Free space wave number.
% Ricker wavelet parameters.
if SourceChoice == 3
    fp = f; % Peak frequency
    dr = PulseWidth*dt*2; % Delay
end

% Material
Partition = ceil(Ra/dr)%SIZE/4;
er1 = 1;
er2 = 1;
ur1 = 1;
ur2 = 1;
sig1 = 0;
sig2 = 65;
sigm1 = 0;
sigm2 = 0;

mu = u0*ones(SIZEr, SIZEphi);
epsEr = e0*ones(SIZEr, SIZEphi);
epsEphi = e0*ones(SIZEr, SIZEphi);
sigr = sig1+zeros(SIZEr, SIZEphi);
sigphi = sig1+zeros(SIZEr, SIZEphi);
sigmr = sigm1+zeros(SIZEr, SIZEphi);
sigmphi = sigm1+zeros(SIZEr, SIZEphi);

for i=1:SIZEr
    for j=1:SIZEphi
        if (i-1)*dr<=Ra
            epsEphi(i,j) = epsEphi(i,j)*er1;
            sigr(i,j) = sig1;
        else
            epsEphi(i,j) = epsEphi(i,j)*er2;
            sigr(i,j) = sig2;
        end

        if (i-0.5)*dr<=Ra
            mu(i,j) = mu(i,j)*ur1;
            epsEr(i,j) = epsEr(i,j)*er1;
            sigphi(i,j) = sig1;
            sigmr(i,j) = sigm1;
            sigmphi(i,j) = sigm1;
        else
            mu(i,j) = mu(i,j)*ur2;
            epsEr(i,j) = epsEr(i,j)*er2;
            sigphi(i,j) = sig2;
            sigmr(i,j) = sigm2;
            sigmphi(i,j) = sigm2;
        end
    end
end

% Order of polynomial
order=6;
% Required reflection co-efficient
gammap=1e-9;
% Polynomial model for sigma
sigmamax=(-log10(gammap)*(order+1)*e0*c)/(2*PMLw*dr);
Bound=((epsEphi(SIZEr-PMLw,1)/e0)*sigmamax)/((PMLw^order)*(order+1));
x=0:1:PMLw;
 for i=1:1:SIZEphi
%     sigr(PMLw+1:-1:1,i)=Lower*((x+0.5*ones(1,PMLw+1)).^(order+1)-(x-0.5*[0 ones(1,PMLw)]).^(order+1));
     sigr(SIZEr-PMLw:1:SIZEr,i)=Bound*((x+0.5*ones(1,PMLw+1)).^(order+1)-(x-0.5*[0 ones(1,PMLw)]).^(order+1));
 end

% Magnetic conductivity
sigmr=(sigr.*mu)./epsEphi;
sigmphi=(sigphi.*mu)./epsEr;
% ==================================================================

% Initialization.
Hthetar = zeros(SIZEr, SIZEphi); % x-component of E-field
Hthetaphi = zeros(SIZEr, SIZEphi); % y-component of E-field
Htheta = zeros(SIZEr, SIZEphi); % Total E-field
Er = zeros(SIZEr, SIZEphi); % x-component of H-field
Ephi = zeros(SIZEr, SIZEphi); % y-component of H-field

ri = zeros(SIZEr, SIZEphi);
rialt = zeros(SIZEr, SIZEphi);

for i=1:1:SIZEr
    ri(i,:) = i-1;
    rialt(i,:) = i-0.5;
end

% Amplitude and phase calculations.
t1 = floor(MaxTime/Scalar)
t2 = floor(t1+1/f/4/dt)
HthetaAbs = zeros(SIZEr, SIZEphi);
HthetaPhase = zeros(SIZEr, SIZEphi);
ErAbs = zeros(SIZEr, SIZEphi);
ErPhase = zeros(SIZEr, SIZEphi);
EphiAbs = zeros(SIZEr, SIZEphi);
EphiPhase = zeros(SIZEr, SIZEphi);
 
HthetaPhasor = zeros(SIZEr, SIZEphi);
ErPhasor = zeros(SIZEr, SIZEphi);
EphiPhasor = zeros(SIZEr, SIZEphi);

%HthetaSaved = zeros(SIZEr, MaxTime);
EphiSaved = zeros(SIZEr, MaxTime);

if SaveFields == 1
    EphiSnapshots = zeros(ceil(SIZEr/SnapshotResolution), ceil(SIZEphi/SnapshotResolution), ceil(MaxTime/SnapshotInterval)); % Data for plotting.
    frame = 1;
end

% Outer loop for time-stepping.
tic
for q = 2:MaxTime
    
    i=2:SIZEr;
    j=1:SIZEphi;
    
    % Calculation of Ephi using update difference equation for Ephi. This is time step q.
    Ephi(i,j) = (1-sigr(i,j)*dt./(2*epsEphi(i,j)))./(1+sigr(i,j)*dt./(2*epsEphi(i,j))).*Ephi(i,j) + ((Htheta(i,j).*rialt(i,j) - Htheta(i-1,j).*rialt(i-1,j)) .* ((dt./(dr*epsEphi(i,j).*ri(i,j)))./(1+sigr(i,j)*dt./(2*epsEphi(i,j)))));
    
    i=1:SIZEr;
    j=2:SIZEphi;
    % Calculation of Er using update difference equation for Er. This is time step q.
    Er(i,j) = (1-sigphi(i,j)*dt./(2*epsEr(i,j)))./(1+sigphi(i,j)*dt./(2*epsEr(i,j))).*Er(i,j) + ((Htheta(i,j-1).*rialt(i,j-1) - Htheta(i,j).*rialt(i,j)) .* ((dt./(dphi*epsEr(i,j).*rialt(i,j).^2))./(1+sigphi(i,j)*dt./(2*epsEr(i,j)))));
    % PBC for Er
    j=1;
    Er(i,j) = (1-sigphi(i,j)*dt./(2*epsEr(i,j)))./(1+sigphi(i,j)*dt./(2*epsEr(i,j))).*Er(i,j) + ((Htheta(i,SIZEphi).*rialt(i,SIZEphi) - Htheta(i,j).*rialt(i,j)) .* ((dt./(dphi*epsEr(i,j).*rialt(i,j).^2))./(1+sigphi(i,j)*dt./(2*epsEr(i,j)))));
        
    i=1:SIZEr-1;
    j=1:SIZEphi;
    % Calculation of Hthetar using updated difference equation. This is time step q+1/2.
    Hthetar(i,j) = (1-sigmr(i,j)*dt./(2*mu(i,j)))./(1+sigmr(i,j)*dt./(2*mu(i,j))).*Hthetar(i,j) + ((Ephi(i+1,j).*ri(i+1,j) - Ephi(i,j).*ri(i,j)) .* ((dt./(dr*mu(i,j).*rialt(i,j)))./(1+sigmr(i,j)*dt./(2*mu(i,j)))));
    
    i=1:SIZEr;
    j=1:SIZEphi-1;
    % Calculation of Hthetaphi using updated difference equation. This is time step q+1/2.
    Hthetaphi(i,j) = (1-sigmphi(i,j)*dt./(2*mu(i,j)))./(1+sigmphi(i,j)*dt./(2*mu(i,j))).*Hthetaphi(i,j) + ((Er(i,j+1) - Er(i,j)) .* ((dt./(dphi*mu(i,j).*rialt(i,j)))./(1+sigmphi(i,j)*dt./(2*mu(i,j)))));
    % PBC for Hthetaphi
    j=SIZEphi;
    Hthetaphi(i,j) = (1-sigmphi(i,j)*dt./(2*mu(i,j)))./(1+sigmphi(i,j)*dt./(2*mu(i,j))).*Hthetaphi(i,j) + ((Er(i,1) - Er(i,j)) .* ((dt./(dphi*mu(i,j).*rialt(i,j)))./(1+sigmphi(i,j)*dt./(2*mu(i,j)))));
    
    
    Htheta = Hthetar + Hthetaphi;
    EphiSaved(:,q) = Ephi(:,1);
    
    sourceX = 2;%Partition;
    % Recording absolute value
    if q > floor(MaxTime/Scalar/2)
    HthetaAbs = max(HthetaAbs, abs(Htheta(:,:)));
    ErAbs = max(ErAbs, abs(Er(:,:)));
    EphiAbs = max(EphiAbs, abs(Ephi(:,:)));
    %HthetaAbs(sourceX,:) = (HthetaAbs(sourceX+1,:) + HthetaAbs(sourceX-1,:))./2;
    EphiAbs(sourceX,:) = (EphiAbs(sourceX+1,:) + EphiAbs(sourceX-1,:))./2;
    %HrAbs(sourceX,:) = (HrAbs(sourceX+1,:) + HrAbs(sourceX-1,:))./2;
    %HphiAbs(sourceX,:) = (HphiAbs(sourceX+1,:) + HphiAbs(sourceX-1,:))./2;
    end
    
    % Source.
    if q < floor(MaxTime)
        if SourceChoice == 1
        Ephi(sourceX,sourceY) = Ephi(sourceX,sourceY) + exp( -1*((q-td)/(PulseWidth/4))^2 ) * Sc;
        elseif SourceChoice == 2
        Ephi(sourceX,sourceY) = Ephi(sourceX,sourceY) + 400*sin(2*pi*f*(q)*dt) * Sc;
        elseif SourceChoice == 3
        Ephi(sourceX,sourceY) = Ephi(sourceX,sourceY) + (1-2*(pi*fp*(q*dt-dr))^2)*exp(-1*(pi*fp*(q*dt-dr))^2) * Sc;
        end
    end

    if (SaveFields == 1 && mod(q, SnapshotInterval) == 0)
        EphiSnapshots(:,:,q/SnapshotInterval+1) = Ephi(1+(0:SnapshotResolution:(SIZEr-1)), 1+(0:SnapshotResolution:(SIZEphi-1)));
    end
    
end
toc
% Simulation animation.
if SaveFields == 1
    % Electric field snapshots.
    sizeS=size(EphiSnapshots);
    for i=1:(MaxTime/SnapshotInterval)-1
        PlotData = EphiSnapshots (:, :, i)/max(max(EphiSnapshots(:,:,i)));
        %PlotData = PlotData/max(max(PlotData));
        figure (6)
        mesh ( EphiSnapshots (:, :, i) );
        view (4, 4)
        xlim([0 sizeS(2)])
        ylim([0 sizeS(1)])
        zlim([-1.1 1.1])
        %caxis([-0.1 0.6])
        xlabel ('j-axis')
        ylabel ('i-axis')
        %colorbar

        figure (7)
        mesh ( EphiSnapshots (:, :, i) );
        view (0, 90)
        xlim([0 sizeS(2)])
        ylim([0 sizeS(1)])
        zlim([-1.1 1.1])
        %caxis([-0.1 0.6])
        xlabel ('j-axis')
        ylabel ('i-axis')
        %colorbar

    end
end
Partition = round(Ra/dr);
for i=1:Scalar:MaxTime
    figure (2)
    subplot(211)
    hold off
    plot([(Partition)*dr/Ra (Partition)*dr/Ra], [-1.1 1.1], 'Color', 'r');
    hold on
    %plot([(sourceX-Partition)*delta/Ra (sourceX-Partition)*delta/Ra], [-1.1 1.1], 'Color', 'r');
    plot((0:SIZEr-1)*dr/Ra,EphiSaved(:,i))
    xlim([0 (SIZEr-1)*dr/Ra])
    ylim([-1.1 1.1])
    xlabel('r/Ra')
    ylabel('Electric field (Ephi) at sourceY')
    
%     subplot(212)
%     hold off
%     plot([(Partition)*delta/Ra (Partition)*delta/Ra], [-1.1 1.1], 'Color', 'r');
%     hold on
%     plot((0:SIZE-1)*delta/Ra,Hy(:,i))
%     xlim([0 (SIZE-1)*dz/Ra])
%     ylim([-1.1/imp0 1.1/imp0])
%     xlabel('r/Ra')
%     ylabel('Magnetic field (Hy)')
end
% 
% figure(3)
% subplot(211)
% hold off
% plot([(Partition)*dz/Ra (Partition)*dz/Ra], [-1.1 1.1], 'Color', 'r');
% hold on
% plot((0:SIZE-1)*dz/Ra,ExAbs)
% xlim([0 (SIZE-1)*dz/Ra])
% ylim([-1.1 1.1])
% xlabel('r/Ra')
% ylabel('Electric field (Ex)')
% 
% subplot(212)
% hold off
% plot([(Partition)*dz/Ra (Partition)*dz/Ra], [-1.1 1.1], 'Color', 'r');
% hold on
% plot((0:SIZE-1)*dz/Ra,HyAbs)
% xlim([0 (SIZE-1)*dz/Ra])
% ylim([-1.1/imp0 1.1/imp0])
% xlabel('r/Ra')
% ylabel('Magnetic field (Hy)')
% 
figure(4)
hold off
plot([(Partition)*dr/Ra (Partition)*dr/Ra], [-1.1 1.1], 'Color', 'r');
hold on
plot((0:SIZEr-1)*dr/Ra,EphiAbs./HthetaAbs)
xlim([0 (SIZEr-1)*dr/Ra])
%ylim([-1.1/imp0 1.1/imp0])
xlabel('r/Ra')
ylabel('Magnitude of wave impedance')
% 
% % ================== Postprocessing =====================
% % Reference: Chapter 5 section 5.6 from Understanding FDTD.
% % for i=1:SIZE
% %     if abs(Ex(i,t1)) < 1e-12 && Ex(i,t2) > 0
% %         ExPhase(i) = -pi/2;
% %     elseif abs(Ex(i,t1)) < 1e-12 && Ex(i,t2) > 0
% %         ExPhase(i) = -pi/2;
% %     else
% %         ExPhase(i) = atan((cos(2*pi*f*t2*dt)-Ex(i,t2)/Ex(i,t1))/(sin(2*pi*f*t2*dt)));
% %     end
% %     if abs(Hy(i,t1)) < 1e-12 && Hy(i,t2) > 0
% %         HyPhase(i) = -pi/2;
% %     elseif abs(Hy(i,t1)) < 1e-12 && Hy(i,t2) > 0
% %         HyPhase(i) = -pi/2;
% %     else
% %         HyPhase(i) = atan((cos(2*pi*f*t2*dt)-Hy(i,t2)/Hy(i,t1))/(sin(2*pi*f*t2*dt)));
% %     end
% % end
% % 
% % for i=1:SIZE
% %     if abs(Ex(i,t1)) >= abs(Ex(i,t2))
% %         ExAbs(i) = Ex(i,t1)/cos(ExPhase(i));
% %     else
% %         ExAbs(i) = Ex(i,t2)/(cos(2*pi*f*t2*dt)*cos(ExPhase(i))-sin(2*pi*f*t2*dt)*sin(ExPhase(i)));
% %     end
% %     if abs(Hy(i,t1)) >= abs(Hy(i,t2))
% %         HyAbs(i) = Hy(i,t1)/cos(HyPhase(i));
% %     else
% %         HyAbs(i) = Hy(i,t2)/(cos(2*pi*f*t2*dt)*cos(HyPhase(i))-sin(2*pi*f*t2*dt)*sin(HyPhase(i)));
% %     end
% % end
% % 
% % for i=1:SIZE
% %     if ExAbs(i) < 0 && ExPhase(i) >= 0
% %         ExAbs(i) = -1*ExAbs(i);
% %         ExPhase(i) = ExPhase(i)-pi;
% %     elseif ExAbs(i) < 0 && ExPhase(i) < 0
% %         ExAbs(i) = -1*ExAbs(i);
% %         ExPhase(i) = ExPhase(i)+pi;
% %     end
% %     if HyAbs(i) < 0 && HyPhase(i) >= 0
% %         HyAbs(i) = -1*HyAbs(i);
% %         HyPhase(i) = HyPhase(i)-pi;
% %     elseif HyAbs(i) < 0 && HyPhase(i) < 0
% %         HyAbs(i) = -1*HyAbs(i);
% %         HyPhase(i) = HyPhase(i)+pi;
% %     end
% % end
% 
HthetaPhase = acos(Htheta./HthetaAbs);
EphiPhase = acos(Ephi./EphiAbs);
HthetaPhasor = HthetaAbs.*exp(1j.*HthetaPhase);
EphiPhasor = EphiAbs.*exp(1j.*EphiPhase);
% 
figure(5)
subplot(211)
hold off
plot([(Partition)*dr/Ra (Partition)*dr/Ra], [-1.1 1.1], 'Color', 'r');
hold on
plot((0:SIZEr-1)*dr/Ra,EphiAbs)
xlim([0 (SIZEr-1)*dr/Ra])
ylim([-1.1 1.1])
xlabel('r/Ra')
ylabel('Magnitude of Electric field (Ex)')

subplot(212)
hold off
plot([(Partition)*dr/Ra (Partition)*dr/Ra], [-1.1 1.1], 'Color', 'r');
hold on
plot((0:SIZEr-1)*dr/Ra,EphiPhase)
xlim([0 (SIZEr-1)*dr/Ra])
ylim([-pi pi])
xlabel('r/Ra')
ylabel('Phase')
% 
figure(6)
subplot(211)
hold off
plot([(Partition)*dr/Ra (Partition)*dr/Ra], [-1.1 1.1], 'Color', 'r');
hold on
plot((0:SIZEr-1)*dr/Ra,HthetaAbs)
xlim([0 (SIZEr-1)*dr/Ra])
ylim([-1.1/imp0 1.1/imp0])
xlabel('r/Ra')
ylabel('Magnitude of Magnetic field (Hy)')

subplot(212)
hold off
plot([(Partition)*dr/Ra (Partition)*dr/Ra], [-1.1 1.1], 'Color', 'r');
hold on
plot((0:SIZEr-1)*dr/Ra,HthetaPhase)
xlim([0 (SIZEr-1)*dr/Ra])
ylim([-pi pi])
xlabel('r/Ra')
ylabel('Phase')
% 
figure(7)
subplot(211)
hold off
plot([(Partition)*dr/Ra (Partition)*dr/Ra], [-1.1 1.1], 'Color', 'r');
hold on
plot((0:SIZEr-1)*dr/Ra,real(EphiPhasor./HthetaPhasor))
xlim([0 (SIZEr-1)*dr/Ra])
%ylim([-1.1/imp0 1.1/imp0])
xlabel('r/Ra')
ylabel('Real part of wave impedance')

subplot(212)
hold off
plot([(Partition)*dr/Ra (Partition)*dr/Ra], [-1.1 1.1], 'Color', 'r');
hold on
plot((0:SIZEr-1)*dr/Ra,imag(EphiPhasor./HthetaPhasor))
xlim([0 (SIZEr-1)*dr/Ra])
%ylim([-1.1/imp0 1.1/imp0])
xlabel('r/Ra')
ylabel('Imag part of wave impedance')