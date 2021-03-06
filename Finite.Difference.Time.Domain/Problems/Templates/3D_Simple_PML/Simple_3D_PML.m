clear all
clc
% Simulation parameters.
Scalar = 1;
PMLw = 10;
SIZE = 70*Scalar+2*PMLw; % No. of spatial steps
MaxTime = SIZE*12*Scalar; % No. of time steps
PulseWidth = round(SIZE*8); % Controls width of Gaussian Pulse
td = PulseWidth/8; % Temporal delay in pulse.
imp0 = 377.0; % Impedence of free space
sourceX = ceil((SIZE+1)/2); % Location of sourceX
sourceY = ceil((SIZE+1)/2); % Location of sourceY
sourceZ = ceil((SIZE+1)/2); % Location of sourceZ
Ra = 150e-6;

% Choice of source.
% 1. Gaussian 2. Sine wave 3. Ricker wavelet
SourceChoice = 2;

SaveFields = 1; % 0. No, 1. Yes.
SnapshotResolution = 1; % Snapshot resolution. 1 is best.
SnapshotInterval = Scalar; % Amount of time delay between snaps.

% Constants.
pi = 3.141592654;
e0 = (1e-9)/(36*pi);
u0 = (1e-7)*4*pi;
c = 1/sqrt(e0*u0);

dt = 0.25e-14;
delta = 5.0e-6;
CourantNumber = c*dt/delta
Sc = c*dt/delta/(3^0.5)
1/sqrt(3)
l = PulseWidth*delta;
f = 5.3839e12%703.6e9%1.0e13%0.8727e13%c/(1*l)
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
sig2 = 100;
sigm1 = 0;
sigm2 = 0;

c1 = 1/sqrt(er1*e0*ur1*u0);
c2 = 1/sqrt(er2*e0*ur2*u0);
Sc1 = c1*dt/delta
Sc2 = c2*dt/delta

epsEx = e0*ones(SIZE, SIZE, SIZE);
epsEy = e0*ones(SIZE, SIZE, SIZE);
epsEz = e0*ones(SIZE, SIZE, SIZE);
muHx = u0*ones(SIZE, SIZE, SIZE);
muHy = u0*ones(SIZE, SIZE, SIZE);
muHz = u0*ones(SIZE, SIZE, SIZE);

sigx = sig1*ones(SIZE, SIZE, SIZE);
sigy = sig1*ones(SIZE, SIZE, SIZE);
sigz = sig1*ones(SIZE, SIZE, SIZE);
sigmx = sigm1*ones(SIZE, SIZE, SIZE);
sigmy = sigm1*ones(SIZE, SIZE, SIZE);
sigmz = sigm1*ones(SIZE, SIZE, SIZE);

for i=1:SIZE
    for j=1:SIZE
        for k=1:SIZE
            x = (i-sourceX)*delta;
            y = (j-sourceY)*delta;
            z = (k-sourceZ)*delta;
            xc = (i-sourceX)*delta;
            yc = (j-sourceY)*delta;
            zc = (k-sourceZ)*delta;
            r = sqrt(x^2 + y^2 + z^2);
            rc1 = sqrt(xc^2 + y^2 + z^2);
            rc2 = sqrt(x^2 + yc^2 + z^2);
            if r<=Ra
                epsEx(i,j,k) = epsEx(i,j,k)*er1;
                epsEy(i,j,k) = epsEy(i,j,k)*er1;
                epsEz(i,j,k) = epsEz(i,j,k)*er1;
                sigx(i,j,k) = sig1;
                sigy(i,j,k) = sig1;
                sigz(i,j,k) = sig1;
                
                muHx(i,j,k) = muHx(i,j,k)*ur1;
                muHy(i,j,k) = muHy(i,j,k)*ur1;
                muHz(i,j,k) = muHz(i,j,k)*ur1;
                sigmx(i,j,k) = sigm1;
                sigmy(i,j,k) = sigm1;
                sigmz(i,j,k) = sigm1;
            else
                epsEx(i,j,k) = epsEx(i,j,k)*er2;
                epsEy(i,j,k) = epsEy(i,j,k)*er2;
                epsEz(i,j,k) = epsEz(i,j,k)*er2;
                sigx(i,j,k) = sig2;
                sigy(i,j,k) = sig2;
                sigz(i,j,k) = sig2;
                
                muHx(i,j,k) = muHx(i,j,k)*ur2;
                muHy(i,j,k) = muHy(i,j,k)*ur2;
                muHz(i,j,k) = muHz(i,j,k)*ur2;
                sigmx(i,j,k) = sigm2;
                sigmy(i,j,k) = sigm2;
                sigmz(i,j,k) = sigm2;
            end
        end
    end
end

% =====================================================================
% Reference: http://www.mathworks.com/matlabcentral/fileexchange/35578-2d-fdtd-of-a-region-with-perfectly-matched-layer-boundary
% Conductivity continuation correction.
sigcorr=sig2;
% Permittivity at PML.
epscorr=e0;
% Order of polynomial
order=6;
% Required reflection co-efficient
gammap=1e-6;
% Polynomial model for sigma
sigmamax=(-log10(gammap)*(order+1)*e0*c)/(2*PMLw*delta);
Bound=((epscorr/e0)*sigmamax)/((PMLw^order)*(order+1));
x=0:1:PMLw;

for m=1:SIZE
    for i=1:1:SIZE
        sigx(PMLw+1:-1:1,i,m)=sigcorr+Bound*((x+0.5*ones(1,PMLw+1)).^(order+1)-(x-0.5*[0 ones(1,PMLw)]).^(order+1));
        sigx(SIZE-PMLw:1:SIZE,i,m)=sigcorr+Bound*((x+0.5*ones(1,PMLw+1)).^(order+1)-(x-0.5*[0 ones(1,PMLw)]).^(order+1));
        sigy(i,PMLw+1:-1:1,m)=sigcorr+Bound*((x+0.5*ones(1,PMLw+1)).^(order+1)-(x-0.5*[0 ones(1,PMLw)]).^(order+1))';
        sigy(i,SIZE-PMLw:1:SIZE,m)=sigcorr+Bound*((x+0.5*ones(1,PMLw+1)).^(order+1)-(x-0.5*[0 ones(1,PMLw)]).^(order+1))';
        sigz(m,i,PMLw+1:-1:1)=sigcorr+Bound*((x+0.5*ones(1,PMLw+1)).^(order+1)-(x-0.5*[0 ones(1,PMLw)]).^(order+1))';
        sigz(m,i,SIZE-PMLw:1:SIZE)=sigcorr+Bound*((x+0.5*ones(1,PMLw+1)).^(order+1)-(x-0.5*[0 ones(1,PMLw)]).^(order+1))';
    end
end
% for i=2:SIZE
%     sigx(:,:,i) = sigx(:,:,1);
%     sigy(:,:,i) = sigy(:,:,1);
%     sigz(i,:,:) = sigz(1,:,:);
% end

% Magnetic conductivity
sigmx=(sigx.*muHx)./epsEx;
sigmy=(sigy.*muHy)./epsEy;
sigmz=(sigz.*muHz)./epsEz;
% ==================================================================

Ex = zeros(SIZE-1, SIZE, SIZE);
Ey = zeros(SIZE, SIZE-1, SIZE);
Ez = zeros(SIZE, SIZE, SIZE-1);

Hx = zeros(SIZE, SIZE-1, SIZE-1);
Hy = zeros(SIZE-1, SIZE, SIZE-1);
Hz = zeros(SIZE-1, SIZE-1, SIZE);

Exy = zeros(SIZE-1, SIZE, SIZE);
Exz = zeros(SIZE-1, SIZE, SIZE);
Eyx = zeros(SIZE, SIZE-1, SIZE);
Eyz = zeros(SIZE, SIZE-1, SIZE);
Ezx = zeros(SIZE, SIZE, SIZE-1);
Ezy = zeros(SIZE, SIZE, SIZE-1);

Hxy = zeros(SIZE, SIZE-1, SIZE-1);
Hxz = zeros(SIZE, SIZE-1, SIZE-1);
Hyx = zeros(SIZE-1, SIZE, SIZE-1);
Hyz = zeros(SIZE-1, SIZE, SIZE-1);
Hzx = zeros(SIZE-1, SIZE-1, SIZE);
Hzy = zeros(SIZE-1, SIZE-1, SIZE);

CExe = (1-sigx*dt./(2*epsEx))./(1+sigx*dt./(2*epsEx)) .* ones(SIZE, SIZE, SIZE);
CEye = (1-sigy*dt./(2*epsEy))./(1+sigy*dt./(2*epsEy)) .* ones(SIZE, SIZE, SIZE);
CEze = (1-sigz*dt./(2*epsEz))./(1+sigz*dt./(2*epsEz)) .* ones(SIZE, SIZE, SIZE);

CHxm = (1-sigmx*dt./(2*muHx))./(1+sigmx*dt./(2*muHx)) .* ones(SIZE, SIZE, SIZE);
CHym = (1-sigmy*dt./(2*muHy))./(1+sigmy*dt./(2*muHy)) .* ones(SIZE, SIZE, SIZE);
CHzm = (1-sigmz*dt./(2*muHz))./(1+sigmz*dt./(2*muHz)) .* ones(SIZE, SIZE, SIZE);

CHxe = ((dt./(muHx*delta))./(1+sigmx*dt./(2*muHx))) .* ones(SIZE, SIZE, SIZE);
CHye = ((dt./(muHy*delta))./(1+sigmy*dt./(2*muHy))) .* ones(SIZE, SIZE, SIZE);
CHze = ((dt./(muHz*delta))./(1+sigmz*dt./(2*muHz))) .* ones(SIZE, SIZE, SIZE);

CExm = ((dt./(epsEx*delta))./(1+sigx*dt./(2*epsEx))) .* ones(SIZE, SIZE, SIZE);
CEym = ((dt./(epsEy*delta))./(1+sigy*dt./(2*epsEy))) .* ones(SIZE, SIZE, SIZE);
CEzm = ((dt./(epsEz*delta))./(1+sigz*dt./(2*epsEz))) .* ones(SIZE, SIZE, SIZE);

% Amplitude and phase calculations.
% t1 = floor(MaxTime/Scalar)
% t2 = floor(t1+1/f/4/dt)
% EzAbs = zeros(SIZE, SIZE);
% EzPhase = zeros(SIZE, SIZE);
% HxAbs = zeros(SIZE, SIZE+1);
% HxPhase = zeros(SIZE, SIZE+1);
% HyAbs = zeros(SIZE+1, SIZE);
% HyPhase = zeros(SIZE+1, SIZE);
% 
% EzPhasor = zeros(SIZE, SIZE);
% HxPhasor = zeros(SIZE, SIZE+1);
% HyPhasor = zeros(SIZE+1, SIZE);
% 
% EzSaved = zeros(SIZE, MaxTime);
% 
if SaveFields == 1
    EzSnapshots = zeros(ceil(SIZE/SnapshotResolution), ceil(SIZE/SnapshotResolution), ceil(MaxTime/SnapshotInterval)); % Data for plotting.
    frame = 1;
    Slice = floor(sourceZ);
end

EAbs = zeros(SIZE-1,SIZE-1);
HAbs = zeros(SIZE-1,SIZE-1);

ExAbs = zeros(SIZE-1,SIZE-1);
EyAbs = zeros(SIZE-1,SIZE-1);
EzAbs = zeros(SIZE-1,SIZE-1);

HxAbs = zeros(SIZE-1,SIZE-1);
HyAbs = zeros(SIZE-1,SIZE-1);
HzAbs = zeros(SIZE-1,SIZE-1);


MaxTime
% Outer loop for time-stepping.
tic
for q = 2:MaxTime
    
%     Ex = zeros(SIZE-1, SIZE, SIZE);
% Ey = zeros(SIZE, SIZE-1, SIZE);
% Ez = zeros(SIZE, SIZE, SIZE-1);
% 
% Hx = zeros(SIZE, SIZE-1, SIZE-1);
% Hy = zeros(SIZE-1, SIZE, SIZE-1);
% Hz = zeros(SIZE-1, SIZE-1, SIZE);

%     % Magnetic field components.
%     i = 1:SIZE; j = 1:SIZE-1; k = 1:SIZE-1;
%     Hx(i,j,k) = CHxm(i,j,k) .* Hx(i,j,k) + CHxe(i,j,k) .* ((Ey(i,j,k+1) - Ey(i,j,k)) - (Ez(i,j+1,k) - Ez(i,j,k)));
%     
%     i = 1:SIZE-1; j = 1:SIZE; k = 1:SIZE-1;
%     Hy(i,j,k) = CHym(i,j,k) .* Hy(i,j,k) + CHye(i,j,k) .* ((Ez(i+1,j,k) - Ez(i,j,k)) - (Ex(i,j,k+1) - Ex(i,j,k)));
%     
%     i = 1:SIZE-1; j = 1:SIZE-1; k = 1:SIZE;
%     Hz(i,j,k) = CHzm(i,j,k) .* Hz(i,j,k) + CHze(i,j,k) .* ((Ex(i,j+1,k) - Ex(i,j,k)) - (Ey(i+1,j,k) - Ey(i,j,k)));
%     
%     % Electric field components.
%     i = 1:SIZE-1; j = 2:SIZE-1; k = 2:SIZE-1;
%     Ex(i,j,k) = CExe(i,j,k) .* Ex(i,j,k) + CExm(i,j,k) .* ((Hz(i,j,k) - Hz(i,j-1,k)) - (Hy(i,j,k) - Hy(i,j,k-1)));
%     
%     i = 2:SIZE-1; j = 1:SIZE-1; k = 2:SIZE-1;
%     Ey(i,j,k) = CEye(i,j,k) .* Ey(i,j,k) + CEym(i,j,k) .* ((Hx(i,j,k) - Hx(i,j,k-1)) - (Hz(i,j,k) - Hz(i-1,j,k)));
%     
%     i = 2:SIZE-1; j = 2:SIZE-1; k = 1:SIZE-1;
%     Ez(i,j,k) = CEze(i,j,k) .* Ez(i,j,k) + CEzm(i,j,k) .* ((Hy(i,j,k) - Hy(i-1,j,k)) - (Hx(i,j,k) - Hx(i,j-1,k)));

    % Magnetic field components.
    i = 1:SIZE; j = 1:SIZE-1; k = 1:SIZE-1;
    Hxy(i,j,k) = CHym(i,j,k) .* Hxy(i,j,k) + CHye(i,j,k) .* (Ey(i,j,k+1) - Ey(i,j,k));
    Hxz(i,j,k) = CHzm(i,j,k) .* Hxz(i,j,k) + CHze(i,j,k) .* (-Ez(i,j+1,k) + Ez(i,j,k));
    
    i = 1:SIZE-1; j = 1:SIZE; k = 1:SIZE-1;
    Hyx(i,j,k) = CHxm(i,j,k) .* Hyx(i,j,k) + CHxe(i,j,k) .* (-Ex(i,j,k+1) + Ex(i,j,k));
    Hyz(i,j,k) = CHzm(i,j,k) .* Hyz(i,j,k) + CHze(i,j,k) .* (Ez(i+1,j,k) - Ez(i,j,k));
    
    i = 1:SIZE-1; j = 1:SIZE-1; k = 1:SIZE;
    Hzx(i,j,k) = CHxm(i,j,k) .* Hzx(i,j,k) + CHxe(i,j,k) .* (Ex(i,j+1,k) - Ex(i,j,k));
    Hzy(i,j,k) = CHym(i,j,k) .* Hzy(i,j,k) + CHye(i,j,k) .* (-Ey(i+1,j,k) + Ey(i,j,k));
    
    Hx = Hxy+Hxz;
    Hy = Hyx+Hyz;
    Hz = Hzx+Hzy;
    
    % Electric field components.
    i = 1:SIZE-1; j = 2:SIZE-1; k = 2:SIZE-1;
    Exy(i,j,k) = CEye(i,j,k) .* Exy(i,j,k) + CEym(i,j,k) .* (-Hy(i,j,k) + Hy(i,j,k-1));
    Exz(i,j,k) = CEze(i,j,k) .* Exz(i,j,k) + CEzm(i,j,k) .* (Hz(i,j,k) - Hz(i,j-1,k));
    
    i = 2:SIZE-1; j = 1:SIZE-1; k = 2:SIZE-1;
    Eyx(i,j,k) = CExe(i,j,k) .* Eyx(i,j,k) + CExm(i,j,k) .* (Hx(i,j,k) - Hx(i,j,k-1));
    Eyz(i,j,k) = CEze(i,j,k) .* Eyz(i,j,k) + CEzm(i,j,k) .* (-Hz(i,j,k) + Hz(i-1,j,k));
    
    i = 2:SIZE-1; j = 2:SIZE-1; k = 1:SIZE-1;
    Ezx(i,j,k) = CExe(i,j,k) .* Ezx(i,j,k) + CExm(i,j,k) .* (-Hx(i,j,k) + Hx(i,j-1,k));
    Ezy(i,j,k) = CEye(i,j,k) .* Ezy(i,j,k) + CEym(i,j,k) .* (Hy(i,j,k) - Hy(i-1,j,k));

    Ex = Exy+Exz;
    Ey = Eyx+Eyz;
    Ez = Ezx+Ezy;
    
    % Recording absolute value
    if q > floor(MaxTime/Scalar/2)
        EAbs = max(EAbs, sqrt(Ex(1:SIZE-1,1:SIZE-1,sourceZ).^2+Ey(1:SIZE-1,1:SIZE-1,sourceZ).^2+Ez(1:SIZE-1,1:SIZE-1,sourceZ).^2));
        HAbs = max(HAbs, sqrt(Hx(1:SIZE-1,1:SIZE-1,sourceZ).^2+Hy(1:SIZE-1,1:SIZE-1,sourceZ).^2+Hz(1:SIZE-1,1:SIZE-1,sourceZ).^2));
        
        ExAbs = max(ExAbs, abs(Ex(1:SIZE-1,1:SIZE-1,sourceZ)));
        EyAbs = max(EyAbs, abs(Ey(1:SIZE-1,1:SIZE-1,sourceZ)));
        EzAbs = max(EzAbs, abs(Ez(1:SIZE-1,1:SIZE-1,sourceZ)));
        
        HxAbs = max(HxAbs, abs(Hx(1:SIZE-1,1:SIZE-1,sourceZ)));
        HyAbs = max(HyAbs, abs(Hy(1:SIZE-1,1:SIZE-1,sourceZ)));
        HzAbs = max(HzAbs, abs(Hz(1:SIZE-1,1:SIZE-1,sourceZ)));
        
        %HthetaAbs(sourceX,:) = (HthetaAbs(sourceX+1,:) + HthetaAbs(sourceX-1,:))./2;
        %EAbs(sourceX,:) = (EphiAbs(sourceX+1,:) + EphiAbs(sourceX-1,:))./2;
        %HrAbs(sourceX,:) = (HrAbs(sourceX+1,:) + HrAbs(sourceX-1,:))./2;
        %HphiAbs(sourceX,:) = (HphiAbs(sourceX+1,:) + HphiAbs(sourceX-1,:))./2;
    end
    
    % Source.
    if q < floor(MaxTime)
        if SourceChoice == 1
        Ez(sourceX,sourceY,sourceZ) = Ez(sourceX,sourceY,sourceZ) + exp( -1*((q-td)/(PulseWidth/4))^2 ) * Sc/sqrt(3);
        elseif SourceChoice == 2
        Ez(sourceX,sourceY,sourceZ) = Ez(sourceX,sourceY,sourceZ) + sin(2*pi*f*(q)*dt) * Sc/sqrt(3);
        elseif SourceChoice == 3
        Ez(sourceX,sourceY,sourceZ) = Ez(sourceX,sourceY,sourceZ) + (1-2*(pi*fp*(q*dt-dr))^2)*exp(-1*(pi*fp*(q*dt-dr))^2) * Sc/sqrt(3);
        end
    end

%     i=2:SIZE;
%     j=i;
%     
%     % Calculation of Hx using update difference equation for Hx. This is time step q.
%     Hx(:,j) = (1-sigmy(:,j)*dt./(2*muHx(:,j)))./(1+sigmy(:,j)*dt./(2*muHx(:,j))).*Hx(:,j) + ((Ez(:,1:SIZE-1) - Ez(:,2:SIZE)) .* ((dt./(muHx(:,j)*delta))./(1+sigmy(:,j)*dt./(2*muHx(:,j)))));
%     
%     % Calculation of Hy using update difference equation for Hy. This is time step q.
%     Hy(i,:) = (1-sigmx(i,:)*dt./(2*muHy(i,:)))./(1+sigmx(i,:)*dt./(2*muHy(i,:))).*Hy(i,:) + ((Ez(2:SIZE,:) - Ez(1:SIZE-1,:)) .* ((dt./(muHy(i,:)*delta))./(1+sigmx(i,:)*dt./(2*muHy(i,:)))));
%     
%     % Calculation of Ezx using updated difference equation for Ez. This is time step q+1/2.
%     Ezx(:,:) = (1-sigx*dt./(2*eps))./(1+sigx*dt./(2*eps)).*Ezx(:,:) + ((Hy(2:SIZE+1,:) - Hy(1:SIZE,:)) .* ((dt./(eps*delta))./(1+sigx*dt./(2*eps))));
%     % Calculation of Ezy using updated difference equation for Ez. This is time step q+1/2.
%     Ezy(:,:) = (1-sigy*dt./(2*eps))./(1+sigy*dt./(2*eps)).*Ezy(:,:) + ((Hx(:,1:SIZE) - Hx(:,2:SIZE+1)) .* ((dt./(eps*delta))./(1+sigy*dt./(2*eps))));
%     
%     Ez = Ezx + Ezy;
%     EzSaved(:,q) = Ez(:,sourceY);
%     % Recording absolute value
%     if q > floor(MaxTime/Scalar)
%     EzAbs = max(EzAbs, abs(Ez(:,:)));
%     HxAbs = max(HxAbs, abs(Hx(:,:)));
%     HyAbs = max(HyAbs, abs(Hy(:,:)));
%     end
%     
%     % Source.
%     if q < floor(MaxTime)
%         if SourceChoice == 1
%         Ezy(sourceX,sourceY) = Ezy(sourceX,sourceY) + exp( -1*((q-td)/(PulseWidth/4))^2 ) * Sc;
%         elseif SourceChoice == 2
%         Ezy(sourceX,sourceY) = Ezy(sourceX,sourceY) + sin(2*pi*f*(q)*dt) * Sc;
%         elseif SourceChoice == 3
%         Ezy(sourceX,sourceY) = Ezy(sourceX,sourceY) + (1-2*(pi*fp*(q*dt-dr))^2)*exp(-1*(pi*fp*(q*dt-dr))^2) * Sc;
%         end
%     end
% 
    if (SaveFields == 1 && mod(q, SnapshotInterval) == 0)
        EzSnapshots(:,:,q/SnapshotInterval+1) = Ez(1+(0:SnapshotResolution:(SIZE-1)), 1+(0:SnapshotResolution:(SIZE-1)), Slice);
    end
    
end
toc
% Simulation animation.
if SaveFields == 1
    % Electric field snapshots.
    sizeS=size(EzSnapshots);
    for i=1:(MaxTime/SnapshotInterval)-1
        PlotData = EzSnapshots(:,:,i)/max(max(EzSnapshots(:,:,i)));
        %PlotData = PlotData/max(max(PlotData));
        figure (1)
        mesh(EzSnapshots(:,:,i));
        view(4, 4)
        xlim([0 sizeS(2)])
        ylim([0 sizeS(1)])
        zlim([-0.001 0.001])
        caxis([-0.0001 0.0001])
        xlabel('j-axis')
        ylabel('i-axis')        
        %colorbar

        figure (2)
        mesh(EzSnapshots(:,:,i));
        view(0, 90)
        xlim([0 sizeS(2)])
        ylim([0 sizeS(1)])
        zlim([-0.001 0.001])
        caxis([-0.0001 0.0001])
        xlabel ('j-axis')
        ylabel ('i-axis')
        %colorbar        
    end
end

 Partition = round(Ra/delta);
for i=1:(MaxTime/SnapshotInterval)-1
    figure (3)
    %hold off
    %plot([(sourceX+Partition)*delta/Ra (sourceX+Partition)*delta/Ra], [-1.1 1.1], 'Color', 'r');
    %hold on
    %plot([(sourceX-Partition)*delta/Ra (sourceX-Partition)*delta/Ra], [-1.1 1.1], 'Color', 'r');
    plot(EzSnapshots (:, sourceY, i))
    %xlim([0 (SIZE-1)*delta/Ra])
    ylim([-0.001 0.001])
    %xlabel('r/Ra')
    %ylabel('Electric field (Ez) at sourceY')
    
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
plot([(Partition)*delta/Ra (Partition)*delta/Ra], [-1.1 1.1], 'Color', 'r');
hold on
plot((0:SIZE/2-2)*delta/Ra,EAbs(SIZE/2:SIZE-2,sourceY)./HAbs(SIZE/2:SIZE-2,sourceY))
xlim([0 (SIZE/2-2)*delta/Ra])
%ylim([-1.1/imp0 1.1/imp0])
xlabel('r/Ra')
ylabel('Magnitude of wave impedance')

HxPhase = acos(Hx(1:SIZE-1,1:SIZE-1,sourceZ)./HxAbs);
HyPhase = acos(Hy(1:SIZE-1,1:SIZE-1,sourceZ)./HyAbs);
HzPhase = acos(Hz(1:SIZE-1,1:SIZE-1,sourceZ)./HzAbs);

ExPhase = acos(Ex(1:SIZE-1,1:SIZE-1,sourceZ)./ExAbs);
EyPhase = acos(Ey(1:SIZE-1,1:SIZE-1,sourceZ)./EyAbs);
EzPhase = acos(Ez(1:SIZE-1,1:SIZE-1,sourceZ)./EzAbs);

HxPhasor = HxAbs.*exp(1j.*HxPhase);
HyPhasor = HyAbs.*exp(1j.*HyPhase);
HzPhasor = HzAbs.*exp(1j.*HzPhase);

ExPhasor = ExAbs.*exp(1j.*ExPhase);
EyPhasor = EyAbs.*exp(1j.*EyPhase);
EzPhasor = EzAbs.*exp(1j.*EzPhase);

HPhasor = sqrt(HxPhasor.^2+HyPhasor.^2+HzPhasor.^2);
EPhasor = sqrt(ExPhasor.^2+EyPhasor.^2+EzPhasor.^2);

figure(5)
subplot(211)
hold off
plot([(Partition)*delta/Ra (Partition)*delta/Ra], [-1.1 1.1], 'Color', 'r');
hold on
plot((0:SIZE/2-2)*delta/Ra,1*real(EPhasor(SIZE/2:SIZE-2,sourceY)./HPhasor(SIZE/2:SIZE-2,sourceY)))
xlim([0 (SIZE/2-2)*delta/Ra])
%ylim([-1.1/imp0 1.1/imp0])
xlabel('r/Ra')
ylabel('Real part of wave impedance')

subplot(212)
hold off
plot([(Partition)*delta/Ra (Partition)*delta/Ra], [-1.1 1.1], 'Color', 'r');
hold on
plot((0:SIZE/2-2)*delta/Ra,imag(EPhasor(SIZE/2:SIZE-2,sourceY)./HPhasor(SIZE/2:SIZE-2,sourceY)))
xlim([0 (SIZE/2-2)*delta/Ra])
%ylim([-1.1/imp0 1.1/imp0])
xlabel('r/Ra')
ylabel('Imag part of wave impedance')

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