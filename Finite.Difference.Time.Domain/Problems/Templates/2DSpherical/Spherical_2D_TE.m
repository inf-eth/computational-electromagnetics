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
sig2 = 1e2;
sigm1 = 0;
sigm2 = 0;

eps = e0*ones(SIZEr, SIZEphi);
muHr = u0*ones(SIZEr, SIZEphi);
muHphi = u0*ones(SIZEr, SIZEphi);
sigr = sig1+zeros(SIZEr, SIZEphi);
sigphi = sig1+zeros(SIZEr, SIZEphi);
sigmr = sigm1+zeros(SIZEr, SIZEphi);
sigmphi = sigm1+zeros(SIZEr, SIZEphi);

for i=1:SIZEr
    for j=1:SIZEphi
        if (i-1)*dr<=Ra
            eps(i,j) = eps(i,j)*er1;
            sigr(i,j) = sig1;
            sigphi(i,j) = sig1;
            muHr(i,j) = muHr(i,j)*ur1;
            sigmphi(i,j) = sigm1;
        else
            eps(i,j) = eps(i,j)*er2;
            sigr(i,j) = sig2;
            sigphi(i,j) = sig2;
            muHr(i,j) = muHr(i,j)*ur2;
            sigmphi(i,j) = sigm2;
        end

        if (i-0.5)*dr<=Ra
            muHphi(i,j) = muHphi(i,j)*ur1;
            sigmr(i,j) = sigm1;
        else
            muHphi(i,j) = muHphi(i,j)*ur2;
            sigmr(i,j) = sigm2;
        end
    end
end

% Order of polynomial
order=6;
% Required reflection co-efficient
gammap=1e-9;
% Polynomial model for sigma
sigmamax=(-log10(gammap)*(order+1)*e0*c)/(2*PMLw*dr);
% Lower=((eps(sourceX,PMLw)/e0)*sigmamax)/((PMLw^order)*(order+1));
% Upper=((eps(sourceX,SIZE-PMLw)/e0)*sigmamax)/((PMLw^order)*(order+1));
% Left=((eps(PMLw,sourceY)/e0)*sigmamax)/((PMLw^order)*(order+1));
Bound=((eps(SIZEr-PMLw,1)/e0)*sigmamax)/((PMLw^order)*(order+1));
x=0:1:PMLw;
 for i=1:1:SIZEphi
%     sigr(PMLw+1:-1:1,i)=Lower*((x+0.5*ones(1,PMLw+1)).^(order+1)-(x-0.5*[0 ones(1,PMLw)]).^(order+1));
     sigr(SIZEr-PMLw:1:SIZEr,i)=Bound*((x+0.5*ones(1,PMLw+1)).^(order+1)-(x-0.5*[0 ones(1,PMLw)]).^(order+1));
 end
% for i=1:1:SIZEr
%     sigphi(i,PMLw+1:-1:1)=Left*((x+0.5*ones(1,PMLw+1)).^(order+1)-(x-0.5*[0 ones(1,PMLw)]).^(order+1))';
%    sigphi(i,SIZEr-PMLw:1:SIZEr)=Bound*((x+0.5*ones(1,PMLw+1)).^(order+1)-(x-0.5*[0 ones(1,PMLw)]).^(order+1))';
% end

% Magnetic conductivity
sigmr=(sigr.*muHphi)./eps;
sigmphi=(sigphi.*muHr)./eps;
% ==================================================================

% Initialization.
Ethetar = zeros(SIZEr, SIZEphi); % x-component of E-field
Ethetaphi = zeros(SIZEr, SIZEphi); % y-component of E-field
Etheta = zeros(SIZEr, SIZEphi); % Total E-field
Hr = zeros(SIZEr, SIZEphi); % x-component of H-field
Hphi = zeros(SIZEr, SIZEphi); % y-component of H-field

ri = zeros(SIZEr, SIZEphi);
rialt = zeros(SIZEr, SIZEphi);

for i=1:1:SIZEr
    ri(i,:) = i-1;
    rialt(i,:) = i-0.5;
end

% Amplitude and phase calculations.
t1 = floor(MaxTime/Scalar)
t2 = floor(t1+1/f/4/dt)
EthetaAbs = zeros(SIZEr, SIZEphi);
EthetaPhase = zeros(SIZEr, SIZEphi);
HrAbs = zeros(SIZEr, SIZEphi);
HrPhase = zeros(SIZEr, SIZEphi);
HphiAbs = zeros(SIZEr, SIZEphi);
HphiPhase = zeros(SIZEr, SIZEphi);
 
EthetaPhasor = zeros(SIZEr, SIZEphi);
HrPhasor = zeros(SIZEr, SIZEphi);
HphiPhasor = zeros(SIZEr, SIZEphi);

EthetaSaved = zeros(SIZEr, MaxTime);

if SaveFields == 1
    EthetaSnapshots = zeros(ceil(SIZEr/SnapshotResolution), ceil(SIZEphi/SnapshotResolution), ceil(MaxTime/SnapshotInterval)); % Data for plotting.
    frame = 1;
end

% Outer loop for time-stepping.
tic
for q = 2:MaxTime
    
    i=1:SIZEr-1;
    j=1:SIZEphi;
    
    % Calculation of Hphi using update difference equation for Hphi. This is time step q.
    Hphi(i,j) = (1-sigmr(i,j)*dt./(2*muHphi(i,j)))./(1+sigmr(i,j)*dt./(2*muHphi(i,j))).*Hphi(i,j) + ((Etheta(2:SIZEr,j).*ri(2:SIZEr,j) - Etheta(1:SIZEr-1,j).*ri(1:SIZEr-1,j)) .* ((dt./(dr*muHphi(i,j).*rialt(i,j)))./(1+sigmr(i,j)*dt./(2*muHphi(i,j)))));
    
    i=2:SIZEr;
    j=1:SIZEphi-1;
    % Calculation of Hr using update difference equation for Hr. This is time step q.
    Hr(i,j) = (1-sigmphi(i,j)*dt./(2*muHr(i,j)))./(1+sigmphi(i,j)*dt./(2*muHr(i,j))).*Hr(i,j) + ((Etheta(i,1:SIZEphi-1).*ri(i,1:SIZEphi-1) - Etheta(i,2:SIZEphi).*ri(i,2:SIZEphi)) .* ((dt./(dphi*muHr(i,j).*rialt(i,j).^2))./(1+sigmphi(i,j)*dt./(2*muHr(i,j)))));
    % PBC for Hr
    j=SIZEphi;
    Hr(i,j) = (1-sigmphi(i,j)*dt./(2*muHr(i,j)))./(1+sigmphi(i,j)*dt./(2*muHr(i,j))).*Hr(i,j) + ((Etheta(i,SIZEphi-1).*ri(i,SIZEphi-1) - Etheta(i,1).*ri(i,1)) .* ((dt./(dphi*muHr(i,j).*rialt(i,j).^2))./(1+sigmphi(i,j)*dt./(2*muHr(i,j)))));
        
    i=2:SIZEr;
    j=1:SIZEphi;
    % Calculation of Ethetar using updated difference equation for Ez. This is time step q+1/2.
    Ethetar(i,j) = (1-sigr(i,j)*dt./(2*eps(i,j)))./(1+sigr(i,j)*dt./(2*eps(i,j))).*Ethetar(i,j) + ((Hphi(2:SIZEr,j).*rialt(2:SIZEr,j) - Hphi(1:SIZEr-1,j).*rialt(1:SIZEr-1,j)) .* ((dt./(dr*eps(i,j).*ri(i,j)))./(1+sigr(i,j)*dt./(2*eps(i,j)))));
    
    i=2:SIZEr;
    j=2:SIZEphi;
    % Calculation of Ezy using updated difference equation for Ez. This is time step q+1/2.
    Ethetaphi(i,j) = (1-sigphi(i,j)*dt./(2*eps(i,j)))./(1+sigphi(i,j)*dt./(2*eps(i,j))).*Ethetaphi(i,j) + ((Hr(i,1:SIZEphi-1) - Hr(i,2:SIZEphi)) .* ((dt./(dphi*eps(i,j).*ri(i,j)))./(1+sigphi(i,j)*dt./(2*eps(i,j)))));
    % PBC for Ethetaphi
    j=1;
    Ethetaphi(i,j) = (1-sigphi(i,j)*dt./(2*eps(i,j)))./(1+sigphi(i,j)*dt./(2*eps(i,j))).*Ethetaphi(i,j) + ((Hr(i,SIZEphi) - Hr(i,1)) .* ((dt./(eps(i,j)*dr))./(1+sigphi(i,j)*dt./(2*eps(i,j)))));
    
    
    Etheta = Ethetar + Ethetaphi;
    EthetaSaved(:,q) = Etheta(:,1);
    
    sourceX = Partition-1;
    % Recording absolute value
    if q > floor(MaxTime/Scalar/2)
    EthetaAbs = max(EthetaAbs, abs(Etheta(:,:)));
    HrAbs = max(HrAbs, abs(Hr(:,:)));
    HphiAbs = max(HphiAbs, abs(Hphi(:,:)));
    EthetaAbs(sourceX,:) = (EthetaAbs(sourceX+1,:) + EthetaAbs(sourceX-1,:))./2;
    %HrAbs(sourceX,:) = (HrAbs(sourceX+1,:) + HrAbs(sourceX-1,:))./2;
    %HphiAbs(sourceX,:) = (HphiAbs(sourceX+1,:) + HphiAbs(sourceX-1,:))./2;
    end
    
    % Source.
    if q < floor(MaxTime)
        if SourceChoice == 1
        Etheta(sourceX,sourceY) = Etheta(sourceX,sourceY) + exp( -1*((q-td)/(PulseWidth/4))^2 ) * Sc;
        elseif SourceChoice == 2
        Etheta(sourceX,sourceY) = Etheta(sourceX,sourceY) + 1*sin(2*pi*f*(q)*dt) * Sc;
        elseif SourceChoice == 3
        Etheta(sourceX,sourceY) = Etheta(sourceX,sourceY) + (1-2*(pi*fp*(q*dt-dr))^2)*exp(-1*(pi*fp*(q*dt-dr))^2) * Sc;
        end
    end

    if (SaveFields == 1 && mod(q, SnapshotInterval) == 0)
        EthetaSnapshots(:,:,q/SnapshotInterval+1) = Etheta(1+(0:SnapshotResolution:(SIZEr-1)), 1+(0:SnapshotResolution:(SIZEphi-1)));
    end
    
end
toc
% Simulation animation.
if SaveFields == 1
    % Electric field snapshots.
    sizeS=size(EthetaSnapshots);
    for i=1:(MaxTime/SnapshotInterval)-1
        PlotData = EthetaSnapshots (:, :, i)/max(max(EthetaSnapshots(:,:,i)));
        %PlotData = PlotData/max(max(PlotData));
        figure (6)
        mesh ( EthetaSnapshots (:, :, i) );
        view (4, 4)
        xlim([0 sizeS(2)])
        ylim([0 sizeS(1)])
        zlim([-1.1 1.1])
        %caxis([-0.1 0.6])
        xlabel ('j-axis')
        ylabel ('i-axis')
        %colorbar

        figure (7)
        mesh ( EthetaSnapshots (:, :, i) );
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
    plot((0:SIZEr-1)*dr/Ra,EthetaSaved(:,i))
    xlim([0 (SIZEr-1)*dr/Ra])
    ylim([-1.1 1.1])
    xlabel('r/Ra')
    ylabel('Electric field (Etheta) at sourceY')
    
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
plot((0:SIZEr-1)*dr/Ra,EthetaAbs./HphiAbs)
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
EthetaPhase = acos(Etheta./EthetaAbs);
HphiPhase = acos(Hphi./HphiAbs);
EthetaPhasor = EthetaAbs.*exp(1j.*EthetaPhase);
HphiPhasor = HphiAbs.*exp(1j.*HphiPhase);
% 
% figure(5)
% subplot(211)
% hold off
% plot([(Partition)*dz/Ra (Partition)*dz/Ra], [-1.1 1.1], 'Color', 'r');
% hold on
% plot((0:SIZE-1)*dz/Ra,ExAbs)
% xlim([0 (SIZE-1)*dz/Ra])
% ylim([-1.1 1.1])
% xlabel('r/Ra')
% ylabel('Magnitude of Electric field (Ex)')
% 
% subplot(212)
% hold off
% plot([(Partition)*dz/Ra (Partition)*dz/Ra], [-1.1 1.1], 'Color', 'r');
% hold on
% plot((0:SIZE-1)*dz/Ra,ExPhase)
% xlim([0 (SIZE-1)*dz/Ra])
% ylim([-pi pi])
% xlabel('r/Ra')
% ylabel('Phase')
% 
% figure(6)
% subplot(211)
% hold off
% plot([(Partition)*dz/Ra (Partition)*dz/Ra], [-1.1 1.1], 'Color', 'r');
% hold on
% plot((0:SIZE-1)*dz/Ra,HyAbs)
% xlim([0 (SIZE-1)*dz/Ra])
% ylim([-1.1/imp0 1.1/imp0])
% xlabel('r/Ra')
% ylabel('Magnitude of Magnetic field (Hy)')
% 
% subplot(212)
% hold off
% plot([(Partition)*dz/Ra (Partition)*dz/Ra], [-1.1 1.1], 'Color', 'r');
% hold on
% plot((0:SIZE-1)*dz/Ra,HyPhase)
% xlim([0 (SIZE-1)*dz/Ra])
% ylim([-pi pi])
% xlabel('r/Ra')
% ylabel('Phase')
% 
% figure(7)
% subplot(211)
% hold off
% plot([(Partition)*dr/Ra (Partition)*dr/Ra], [-1.1 1.1], 'Color', 'r');
% hold on
% plot((0:SIZEr-1)*dr/Ra,real(EthetaPhasor./HphiPhasor))
% xlim([0 (SIZEr-1)*dr/Ra])
% %ylim([-1.1/imp0 1.1/imp0])
% xlabel('r/Ra')
% ylabel('Real part of wave impedance')
% 
% subplot(212)
% hold off
% plot([(Partition)*dr/Ra (Partition)*dr/Ra], [-1.1 1.1], 'Color', 'r');
% hold on
% plot((0:SIZEr-1)*dr/Ra,imag(EthetaPhasor./HphiPhasor))
% xlim([0 (SIZEr-1)*dr/Ra])
% %ylim([-1.1/imp0 1.1/imp0])
% xlabel('r/Ra')
% ylabel('Imag part of wave impedance')