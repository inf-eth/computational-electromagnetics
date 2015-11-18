clear all
clc
% Simulation parameters.
Scalar = 1;
PMLw = 10;
SIZE = 256*Scalar+2*PMLw; % No. of spatial steps
MaxTime = SIZE*6*Scalar; % No. of time steps
PulseWidth = round(SIZE/8); % Controls width of Gaussian Pulse
td = PulseWidth; % Temporal delay in pulse.
imp0 = 377.0; % Impedence of free space
sourceX = ceil((SIZE+1)/2); % Location of sourceX
sourceY = ceil((SIZE+1)/2); % Location of sourceX
Ra = 150e-6;

% Choice of source.
% 1. Gaussian 2. Sine wave 3. Ricker wavelet
SourceChoice = 2;

SaveFields = 1; % 0. No, 1. Yes.
SnapshotResolution = 1; % Snapshot resolution. 1 is best.
SnapshotInterval = 6*Scalar; % Amount of time delay between snaps.

% Constants.
pi = 3.141592654;
e0 = (1e-9)/(36*pi);
u0 = (1e-7)*4*pi;
c = 1/sqrt(e0*u0);

dt = 0.15e-14;
delta = 0.95e-6;
Sc = c*dt/delta

l = PulseWidth*delta;
f = 3.3e12%0.8727e13%c/(1*l)
fmax = 1/(2*dt)
w = 2*pi*f;
k0 = w/c; % Free space wave number.
% Ricker wavelet parameters.
if SourceChoice == 3
    fp = f; % Peak frequency
    dr = PulseWidth*dt*2; % Delay
end

% Material
Partition = round(Ra/delta)%SIZE/4;
er1 = 1;
er2 = 1;
ur1 = 1;
ur2 = 1;
sig1 = 0;
sig2 = 1e2;
sigm1 = 0;
sigm2 = 0;

c1 = 1/sqrt(er1*e0*ur1*u0);
c2 = 1/sqrt(er2*e0*ur2*u0);
Sc1 = c1*dt/delta
Sc2 = c2*dt/delta

eps = e0*ones(SIZE, SIZE);
muHx = u0*ones(SIZE, SIZE+1);
muHy = u0*ones(SIZE+1, SIZE);
sigx = sig1+zeros(SIZE, SIZE);
sigy = sig1+zeros(SIZE, SIZE);
sigmx = sigm1+zeros(SIZE+1, SIZE);
sigmy = sigm1+zeros(SIZE, SIZE+1);

for i=1:SIZE+1
    for j=1:SIZE+1
        x = (i-sourceX)*delta;
        y = (j-sourceY)*delta;
        xc = (i-sourceX)*delta;
        yc = (j-sourceY)*delta;
        r = sqrt(x^2 + y^2);
        rc1 = sqrt(xc^2 + y^2);
        rc2 = sqrt(x^2 + yc^2);
        if i<SIZE+1 && j<SIZE+1 
            if r<=Ra
                eps(i,j) = eps(i,j)*er1;
                sigx(i,j) = sig1;
                sigy(i,j) = sig1;
            else
                eps(i,j) = eps(i,j)*er2;
                sigx(i,j) = sig2;
                sigy(i,j) = sig2;
            end
        end
        if i<SIZE+1
            if rc2<=Ra
                muHx(i,j) = muHx(i,j)*ur1;
                sigmy(i,j) = sigm1;
            else
                muHx(i,j) = muHx(i,j)*ur2;
                sigmy(i,j) = sigm2;
            end
        end
        if j<SIZE+1
            if rc1<=Ra
                muHy(i,j) = muHy(i,j)*ur1;
                sigmx(i,j) = sigm1;
            else
                muHy(i,j) = muHy(i,j)*ur2;
                sigmx(i,j) = sigm2;
            end
        end
    end
end

% Order of polynomial
order=6;
% Required reflection co-efficient
gammap=1e-6;
% Polynomial model for sigma
sigmamax=(-log10(gammap)*(order+1)*e0*c)/(2*PMLw*delta);
Lower=((eps(sourceX,PMLw)/e0)*sigmamax)/((PMLw^order)*(order+1));
Upper=((eps(sourceX,SIZE-PMLw)/e0)*sigmamax)/((PMLw^order)*(order+1));
Left=((eps(PMLw,sourceY)/e0)*sigmamax)/((PMLw^order)*(order+1));
Right=((eps(SIZE-PMLw,sourceY)/e0)*sigmamax)/((PMLw^order)*(order+1));
x=0:1:PMLw;
for i=1:1:SIZE
    sigx(PMLw+1:-1:1,i)=Lower*((x+0.5*ones(1,PMLw+1)).^(order+1)-(x-0.5*[0 ones(1,PMLw)]).^(order+1));
    sigx(SIZE-PMLw:1:SIZE,i)=Upper*((x+0.5*ones(1,PMLw+1)).^(order+1)-(x-0.5*[0 ones(1,PMLw)]).^(order+1));
end
for i=1:1:SIZE
    sigy(i,PMLw+1:-1:1)=Left*((x+0.5*ones(1,PMLw+1)).^(order+1)-(x-0.5*[0 ones(1,PMLw)]).^(order+1))';
    sigy(i,SIZE-PMLw:1:SIZE)=Right*((x+0.5*ones(1,PMLw+1)).^(order+1)-(x-0.5*[0 ones(1,PMLw)]).^(order+1))';
end

% Magnetic conductivity
sigmx(1:SIZE,:)=(sigx.*muHy(1:SIZE,:))./eps;
sigmy(:,1:SIZE)=(sigy.*muHx(:,1:SIZE))./eps;
% ==================================================================

% Initialization.
Ezx = zeros(SIZE, SIZE); % x-component of E-field
Ezy = zeros(SIZE, SIZE); % y-component of E-field
Ez = zeros(SIZE, SIZE); % Total E-field
Hx = zeros(SIZE, SIZE+1); % x-component of H-field
Hy = zeros(SIZE+1, SIZE); % y-component of H-field

% Amplitude and phase calculations.
t1 = floor(MaxTime/Scalar)
t2 = floor(t1+1/f/4/dt)
EzAbs = zeros(SIZE, SIZE);
EzPhase = zeros(SIZE, SIZE);
HxAbs = zeros(SIZE, SIZE+1);
HxPhase = zeros(SIZE, SIZE+1);
HyAbs = zeros(SIZE+1, SIZE);
HyPhase = zeros(SIZE+1, SIZE);

EzPhasor = zeros(SIZE, SIZE);
HxPhasor = zeros(SIZE, SIZE+1);
HyPhasor = zeros(SIZE+1, SIZE);

EzSaved = zeros(SIZE, MaxTime);

if SaveFields == 1
    EzSnapshots = zeros(ceil(SIZE/SnapshotResolution), ceil(SIZE/SnapshotResolution), ceil(MaxTime/SnapshotInterval)); % Data for plotting.
    frame = 1;
end

% Outer loop for time-stepping.
tic
for q = 2:MaxTime
    
    i=2:SIZE;
    j=i;
    
    % Calculation of Hx using update difference equation for Hx. This is time step q.
    Hx(:,j) = (1-sigmy(:,j)*dt./(2*muHx(:,j)))./(1+sigmy(:,j)*dt./(2*muHx(:,j))).*Hx(:,j) + ((Ez(:,1:SIZE-1) - Ez(:,2:SIZE)) .* ((dt./(muHx(:,j)*delta))./(1+sigmy(:,j)*dt./(2*muHx(:,j)))));
    
    % Calculation of Hy using update difference equation for Hy. This is time step q.
    Hy(i,:) = (1-sigmx(i,:)*dt./(2*muHy(i,:)))./(1+sigmx(i,:)*dt./(2*muHy(i,:))).*Hy(i,:) + ((Ez(2:SIZE,:) - Ez(1:SIZE-1,:)) .* ((dt./(muHy(i,:)*delta))./(1+sigmx(i,:)*dt./(2*muHy(i,:)))));
    
    % Calculation of Ezx using updated difference equation for Ez. This is time step q+1/2.
    Ezx(:,:) = (1-sigx*dt./(2*eps))./(1+sigx*dt./(2*eps)).*Ezx(:,:) + ((Hy(2:SIZE+1,:) - Hy(1:SIZE,:)) .* ((dt./(eps*delta))./(1+sigx*dt./(2*eps))));
    % Calculation of Ezy using updated difference equation for Ez. This is time step q+1/2.
    Ezy(:,:) = (1-sigy*dt./(2*eps))./(1+sigy*dt./(2*eps)).*Ezy(:,:) + ((Hx(:,1:SIZE) - Hx(:,2:SIZE+1)) .* ((dt./(eps*delta))./(1+sigy*dt./(2*eps))));
    
    Ez = Ezx + Ezy;
    EzSaved(:,q) = Ez(:,sourceY);
    % Recording absolute value
    if q > floor(MaxTime/Scalar)
    EzAbs = max(EzAbs, abs(Ez(:,:)));
    HxAbs = max(HxAbs, abs(Hx(:,:)));
    HyAbs = max(HyAbs, abs(Hy(:,:)));
    end
    
    % Source.
    if q < floor(MaxTime)
        if SourceChoice == 1
        Ezy(sourceX,sourceY) = Ezy(sourceX,sourceY) + exp( -1*((q-td)/(PulseWidth/4))^2 ) * Sc;
        elseif SourceChoice == 2
        Ezy(sourceX,sourceY) = Ezy(sourceX,sourceY) + sin(2*pi*f*(q)*dt) * Sc;
        elseif SourceChoice == 3
        Ezy(sourceX,sourceY) = Ezy(sourceX,sourceY) + (1-2*(pi*fp*(q*dt-dr))^2)*exp(-1*(pi*fp*(q*dt-dr))^2) * Sc;
        end
    end

    if (SaveFields == 1 && mod(q, SnapshotInterval) == 0)
        EzSnapshots(:,:,q/SnapshotInterval+1) = Ez(1+(0:SnapshotResolution:(SIZE-1)), 1+(0:SnapshotResolution:(SIZE-1)));
    end
    
end
toc
% Simulation animation.
if SaveFields == 1
    % Electric field snapshots.
    sizeS=size(EzSnapshots);
    for i=1:(MaxTime/SnapshotInterval)-1
        PlotData = EzSnapshots (:, :, i)/max(max(EzSnapshots(:,:,i)));
        %PlotData = PlotData/max(max(PlotData));
        figure (6)
        mesh ( PlotData );
        view (4, 4)
        xlim([0 sizeS(2)])
        ylim([0 sizeS(1)])
        zlim([-1.1 1.1])
        %caxis([-0.1 0.6])
        xlabel ('j-axis')
        ylabel ('i-axis')
        %colorbar

        figure (7)
        mesh ( PlotData );
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
Partition = round(Ra/delta);
for i=1:Scalar:MaxTime
    figure (2)
    subplot(211)
    hold off
    plot([(sourceX+Partition)*delta/Ra (sourceX+Partition)*delta/Ra], [-1.1 1.1], 'Color', 'r');
    hold on
    plot([(sourceX-Partition)*delta/Ra (sourceX-Partition)*delta/Ra], [-1.1 1.1], 'Color', 'r');
    plot((0:SIZE-1)*delta/Ra,EzSaved(:,i))
    xlim([0 (SIZE-1)*delta/Ra])
    ylim([-0.1 0.1])
    xlabel('r/Ra')
    ylabel('Electric field (Ez) at sourceY')
    
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
% figure(4)
% hold off
% plot([(Partition)*dz/Ra (Partition)*dz/Ra], [-1.1 1.1], 'Color', 'r');
% hold on
% plot((0:SIZE-1)*dz/Ra,ExAbs./HyAbs)
% xlim([0 (SIZE-1)*dz/Ra])
% %ylim([-1.1/imp0 1.1/imp0])
% xlabel('r/Ra')
% ylabel('Magnitude of wave impedance')
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
% ExPhase = acos(Ex(:,MaxTime)./ExAbs);
% HyPhase = acos(Hy(:,MaxTime)./HyAbs);
% ExPhasor = ExAbs.*exp(1j.*ExPhase);
% HyPhasor = HyAbs.*exp(1j.*HyPhase);
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
% plot([(Partition)*dz/Ra (Partition)*dz/Ra], [-1.1 1.1], 'Color', 'r');
% hold on
% plot((0:SIZE-1)*dz/Ra,real(ExPhasor./HyPhasor))
% xlim([0 (SIZE-1)*dz/Ra])
% %ylim([-1.1/imp0 1.1/imp0])
% xlabel('r/Ra')
% ylabel('Real part of wave impedance')
% 
% subplot(212)
% hold off
% plot([(Partition)*dz/Ra (Partition)*dz/Ra], [-1.1 1.1], 'Color', 'r');
% hold on
% plot((0:SIZE-1)*dz/Ra,imag(ExPhasor./HyPhasor))
% xlim([0 (SIZE-1)*dz/Ra])
% %ylim([-1.1/imp0 1.1/imp0])
% xlabel('r/Ra')
% ylabel('Imag part of wave impedance')