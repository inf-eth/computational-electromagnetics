clc
clear all

% Simulation parameters.
[pI pJ pPMLw pra prb pdt pdelta px0 py0] = Parameters;
I = pI; % No. of spatial steps in i direction.
J = pJ; % No. of spatial steps in j direction.
PMLw = pPMLw; % Width of PML layer.
% Cylinder dimensions.
ra = pra; % Inner radius of cylinder.
rb = prb; % Outer radius of cylinder.
MaxTime = 12*J; % No. of time steps
PulseWidth = 4*round(J/8); % Controls width of Gaussian Pulse
td = PulseWidth; % Temporal delay in pulse.
SaveFields = 1; % 0. No, 1. Yes.
SnapshotResolution = 1; % Snapshot resolution. 1 is best.
SnapshotInterval = 16; % Amount of time delay between snaps.
% Choice of source.
% 1. Gaussian 2. Sine wave 3. Ricker wavelet
SourceChoice = 2;
SourcePlane = 1; % Is the source a plane wave. 0. = Omni 1. Plane-wave.
SourceLocationX = I/2; % X Location of source. Only used for an omni-source.
SourceLocationY = PMLw+6; % Y Location of source.

% Constants.
c = 3e8;
pi = 3.141592654;
e0 = (1e-9)/(36*pi);
u0 = (1e-7)*4*pi;

dt = pdt;
delta = pdelta;
Sc = c * dt/delta

% Cylinder location.
x0 = px0;
y0 = py0;

l = PulseWidth*delta;
f = 2e9;%c/l
fmax = 1/(2*dt)
w = 2*pi*f;
k0 = w/c; % Free space wave number.
% Ricker wavelet parameters.
if SourceChoice == 3
    fp = f; % Peak frequency
    dr = PulseWidth*dt*2; % Delay
end

% Array sizes.
IEz = I;
JEz = J+2*PMLw;
IHx = I;
JHx = J+2*PMLw+1;
IHy = I;
JHy = J+2*PMLw;

% Initialization.
Ez = zeros(IEz, JEz, 3); % z-component of E-field
Dz = zeros(IEz, JEz, 3); % z-component of D
EzMask = zeros(IEz, JEz); % Ez mask for PEC.
Hx = zeros(IHx, JHx, 3); % i-component of H-field
Bx = zeros(IHx, JHx, 3); % i-component of B
HxMask = zeros(IHx, JHx); % Hx mask for PMC.
BxAve = zeros(IHy, JHy, 3); % Averaged i-component of B wrt to By.
Hy = zeros(IHy, JHy, 3); % j-component of H-field
By = zeros(IHy, JHy, 3); % j-component of B
HyMask = zeros(IHy, JHy); % Hy mask for PMC.
ByAve = zeros(IHx, JHx, 3); % Averaged j-component of B wrt Bx.

% Incident and Transmitted Fields.
Ezi = zeros(MaxTime, 1);
Ezt = zeros(MaxTime, 1);
Eztt = zeros(MaxTime, 1);
x1 = J/2+1; % Position of observation.

% Refractive Index calculations.
Y1 = J/2 + 5;
y1 = Y1*delta;
Y2 = J/2 + 6;
y2 = Y2*delta;
Ezy1 = zeros(MaxTime, 1);
Ezy2 = zeros(MaxTime, 1);

%einf = ones(IHx,JHx);
%einf(:,SlabLeft:SlabRight) = 1; % einf(Drude) or er in slab.
%uinf = ones(IHx,JHx);
%uinf(:,SlabLeft:SlabRight) = 1; % uinf(Drude) or ur in slab.
%wpesq = zeros(IHx,JHx);
%wpesq(:,SlabLeft:SlabRight) = 2*w^2; % DNG(Drude) value of wpe squared in slab.
%wpmsq = zeros(IHx,JHx);
%wpmsq(:,SlabLeft:SlabRight) = 2*w^2; % DNG(Drude) value of wpm squared in slab.
%ge = zeros(IHx,JHx);
%ge(:,SlabLeft:SlabRight) = w/32; % Electric collision frequency in slab.
%gm = zeros(IHx,JHx);
%gm(:,SlabLeft:SlabRight) = w/32; % Magnetic collision frequency in slab.

ax = zeros(IHx, JHx, 9);
ay = zeros(IHy, JHy, 9);
az = zeros(IEz, JEz, 5);

linecount = 0;
fprintf ( 1, 'Initialisation started! \n');
for i=1:IHx
    % Progress indicator.
    if mod(i,2) == 0
        fprintf(1, repmat('\b',1,linecount));
        linecount = fprintf(1, '%g %%', (i*100)/IHx);
    end
    for j=1:JHx
        ax(i,j,1) = (sinphi(i-1,j-1)^2)*(1/dt^2+wpmsquaredc(i-1,j-1,w)/4)+uphi(i-1,j-1)*(cosphi(i-1,j-1)^2)*(1/dt^2);
        ax(i,j,2) = (sinphi(i-1,j-1)^2)*(-2/dt^2+wpmsquaredc(i-1,j-1,w)/2)-uphi(i-1,j-1)*(cosphi(i-1,j-1)^2)*(2/dt^2);
        ax(i,j,3) = ax(i,j,1);
        ax(i,j,4) = (uphi(i-1,j-1)*(1/dt^2)-(1/dt^2+wpmsquaredc(i-1,j-1,w)/4))*sinphi(i-1,j-1)*cosphi(i-1,j-1);
        ax(i,j,5) = (uphi(i-1,j-1)*(-2/dt^2)-(-2/dt^2+wpmsquaredc(i-1,j-1,w)/2))*sinphi(i-1,j-1)*cosphi(i-1,j-1);
        ax(i,j,6) = ax(i,j,4);
        ax(i,j,7) = u0*uphi(i-1,j-1)*(-2/dt^2+wpmsquaredc(i-1,j-1,w)/2);
        ax(i,j,8) = u0*uphi(i-1,j-1)*(1/dt^2+wpmsquaredc(i-1,j-1,w)/4);
        ax(i,j,9) = ax(i,j,8);
        HxMask(i,j) = mask(i-1,j-0.5);
        if j < JHx
            ay(i,j,1) = (cosphi(i-1.5,j-1.5)^2)*(1/dt^2+wpmsquaredc(i-1.5,j-1.5,w)/4)+uphi(i-1.5,j-1.5)*(sinphi(i-1.5,j-1.5)^2)*(1/dt^2);
            ay(i,j,2) = (cosphi(i-1.5,j-1.5)^2)*(-2/dt^2+wpmsquaredc(i-1.5,j-1.5,w)/2)-uphi(i-1.5,j-1.5)*(sinphi(i-1.5,j-1.5)^2)*(2/dt^2);
            ay(i,j,3) = ay(i,j,1);
            ay(i,j,4) = (uphi(i-1.5,j-1.5)*(1/dt^2)-(1/dt^2+wpmsquaredc(i-1.5,j-1.5,w)/4))*sinphi(i-1.5,j-1.5)*cosphi(i-1.5,j-1.5);
            ay(i,j,5) = (uphi(i-1.5,j-1.5)*(-2/dt^2)-(-2/dt^2+wpmsquaredc(i-1.5,j-1.5,w)/2))*sinphi(i-1.5,j-1.5)*cosphi(i-1.5,j-1.5);
            ay(i,j,6) = ay(i,j,4);
            ay(i,j,7) = u0*uphi(i-1.5,j-1.5)*(-2/dt^2+wpmsquaredc(i-1.5,j-1.5,w)/2);
            ay(i,j,8) = u0*uphi(i-1.5,j-1.5)*(1/dt^2+wpmsquaredc(i-1.5,j-1.5,w)/4);
            ay(i,j,9) = ay(i,j,8);
            
            az0 = (4*dt^2)./(e0*(4*einf(i-1,j-1.5)+dt^2*wpsquared(i-1,j-1.5,w)+2*dt*einf(i-1,j-1.5)*gammae(i-1,j-1.5,w)));
            az(i,j,1) = (1/dt^2)*az0/A(i-1,j-1.5);
            az(i,j,2) = (1/(2*dt))*gammae(i-1,j-1.5,w)*az0/A(i-1,j-1.5);
            az(i,j,3) = (e0/dt^2)*einf(i-1,j-1.5)*az0;
            az(i,j,4) = (-1*e0/4)*wpsquared(i-1,j-1.5,w)*az0;
            az(i,j,5) = (1/(2*dt))*e0*einf(i-1,j-1.5)*gammae(i-1,j-1.5,w)*az0;
            
            EzMask(i,j) = mask(i-1,j-1.5);
            HyMask(i,j) = mask(i-1.5,j-1.5);
        end        
    end
end
fprintf(1, repmat('\b',1,linecount));
fprintf ( 1, 'Initialisation complete! \n');
mesh(ax(:,:,1))
% PML arrays.
PsiEzX = zeros(IEz, JEz);
PsiEzY = zeros(IEz, JEz);
PsiHyX = zeros(IHy, JHy);
PsiHxY = zeros(IHx, JHx);

% PML parameters.
kapp = 1;
a = 0.0004;
sig = 0.045;
kappe = 1.5;
kappm = 1.5;
% Electric.
kappex = kapp;
kappey = kapp;
aex = a;
aey = a;
sigex = 0;
sigey = sig;
bex = exp(-1*(aex/e0+sigex/(kappex*e0))*dt);
bey = exp(-1*(aey/e0+sigey/(kappey*e0))*dt);
Cex = (bex-1)*sigex/(sigex*kappex+kappe^2*aex);
Cey = (bey-1)*sigey/(sigey*kappey+kappe^2*aey);
% Magnetic.
kappmx = kapp;
kappmy = kapp;
amx = a;
amy = a;
sigmx = 0;
sigmy = u0/e0*sig;
bmx = exp(-1*(amx/u0+sigmx/(kappmx*u0))*dt);
bmy = exp(-1*(amy/u0+sigmy/(kappmy*u0))*dt);
Cmx = (bmx-1)*sigmx/(sigmx*kappmx+kappm^2*amx);
Cmy = (bmy-1)*sigmy/(sigmy*kappmy+kappm^2*amy);

if SaveFields == 1
    EzSnapshots = zeros(IEz/SnapshotResolution, (JEz)/SnapshotResolution, MaxTime/SnapshotInterval); % Data for plotting.
    frame = 1;
end

np = 1;
n0 = 2;
nf = 3;
linecount = 0;
% Outer loop for time-stepping.
fprintf ( 1, 'Dry run started! \n');
tic
% Test loop for incident field in free space.
for n = 0:MaxTime
    
    % Progress indicator.
    if mod(n,SnapshotInterval) == 0
        fprintf(1, repmat('\b',1,linecount));
        linecount = fprintf(1, '%g %%', (n*100)/MaxTime );
    end
    
    % ========================= Bx and Hx =============================
    % Hx Psi array.
    i=1:IHx;
    j=2:JHx-1;
    PsiHxY(i,j) = (Cmy/delta)*(-Ez(i,j,n0) + Ez(i,j-1,n0)) + bmy*PsiHxY(i,j);
    % Bx in normal space.
    j=(2+PMLw):(JHx-PMLw-1);
    Bx(i,j,nf) = Bx(i,j,n0) + (-Ez(i,j,n0) + Ez(i,j-1,n0)) * dt/delta;
    if PMLw > 0
        % Bx in lower PML layer.
        j=2:PMLw+1;
        Bx(i,j,nf) = Bx(i,j,n0) + dt*((1/kappmy)*(-Ez(i,j,n0) + Ez(i,j-1,n0)) * 1/delta + PsiHxY(i,j));
        % Bx in upper PML layer.
        j=JHx-PMLw:JHx-1;
        Bx(i,j,nf) = Bx(i,j,n0) + dt*((1/kappmy)*(-Ez(i,j,n0) + Ez(i,j-1,n0)) * 1/delta + PsiHxY(i,j));
    end
    Hx(:,:,nf) = Bx(:,:,nf)./u0;
    
    % ========================= By and Hy =============================
    % Hy Psi array.
    i=1:IHy-1;
    j=1:JHy;
    PsiHyX(i,j) = (Cmx/delta)*(Ez(i+1,j,n0)-Ez(i,j,n0)) + bmx*PsiHyX(i,j);
    PsiHyX(IHy,j) = (Cmx/delta)*(Ez(1,j,n0)-Ez(IHy,j,n0)) + bmx*PsiHyX(IHy,j);
    % By in normal space.
    j=(1+PMLw):JHy-PMLw;
    By(i,j,nf) = By(i,j,n0) + (Ez(i+1,j,n0) - Ez(i,j,n0)) * dt/delta;
    By(IHy,j,nf) = By(IHy,j,n0) + (Ez(1,j,n0) - Ez(IHy,j,n0)) * dt/delta; % PBC
    if PMLw > 0
        % By in lower PML layer.
        j=1:PMLw;
        By(i,j,nf) = By(i,j,n0) + dt*((1/kappmx)*(Ez(i+1,j,n0) - Ez(i,j,n0)) * 1/delta + PsiHyX(i,j));
        By(IHy,j,nf) = By(IHy,j,n0) + dt*((1/kappmx)*(Ez(1,j,n0) - Ez(IHy,j,n0)) * 1/delta + PsiHyX(IHy,j)); % PBC
        % By in upper PML layer.
        j=JHy-PMLw+1:JHy;
        By(i,j,nf) = By(i,j,n0) + dt*((1/kappmx)*(Ez(i+1,j,n0) - Ez(i,j,n0)) * 1/delta + PsiHyX(i,j));
        By(IHy,j,nf) = By(IHy,j,n0) + dt*((1/kappmx)*(Ez(1,j,n0) - Ez(IHy,j,n0)) * 1/delta + PsiHyX(IHy,j)); % PBC
    end
    Hy(:,:,nf) = By(:,:,nf)./u0;
    
    % ========================= Dz and Ez =============================
    % Psi arrays.
    i=2:IEz;
    j=1:JEz;
    PsiEzX(i,j) = (Cex/delta)*(Hy(i,j,nf)-Hy(i-1,j,nf)) + bex*PsiEzX(i,j);
    PsiEzX(1,j) = (Cex/delta)*(Hy(1,j,nf)-Hy(IEz,j,nf)) + bex*PsiEzX(1,j); % PBC
    PsiEzY(i,j) = (Cey/delta)*(-Hx(i,j+1,nf)+Hx(i,j,nf)) + bey*PsiEzY(i,j);
    PsiEzY(1,j) = (Cey/delta)*(-Hx(1,j+1,nf)+Hx(1,j,nf)) + bey*PsiEzY(1,j); % PBC
    % Dz in Normal Space.
    j=(1+PMLw):(JEz-PMLw);
    Dz(i,j,nf) = Dz(i,j,n0) + (Hy(i,j,nf)-Hy(i-1,j,nf)-Hx(i,j+1,nf)+Hx(i,j,nf)) * dt/delta;
    Dz(1,j,nf) = Dz(1,j,n0) + (Hy(1,j,nf)-Hy(IEz,j,nf)-Hx(1,j+1,nf)+Hx(1,j,nf)) * dt/delta; % PBC
    if PMLw > 0
        % Dz in lower PML layer.
        j=1:PMLw;
        Dz(i,j,nf) = Dz(i,j,n0) + dt*(((1/kappex)*(Hy(i,j,nf)-Hy(i-1,j,nf))+(1/kappey)*(-Hx(i,j+1,nf)+Hx(i,j,nf))) * 1/delta + PsiEzX(i,j) + PsiEzY(i,j));
        Dz(1,j,nf) = Dz(1,j,n0) + dt*(((1/kappex)*(Hy(1,j,nf)-Hy(IEz,j,nf))+(1/kappey)*(-Hx(1,j+1,nf)+Hx(1,j,nf))) * 1/delta + PsiEzX(1,j) + PsiEzY(1,j)); % PBC
        % Dz in upper PML layer.
        j=JEz-PMLw+1:JEz;
        Dz(i,j,nf) = Dz(i,j,n0) + dt*(((1/kappex)*(Hy(i,j,nf)-Hy(i-1,j,nf))+(1/kappey)*(-Hx(i,j+1,nf)+Hx(i,j,nf))) * 1/delta + PsiEzX(i,j) + PsiEzY(i,j));
        Dz(1,j,nf) = Dz(1,j,n0) + dt*(((1/kappex)*(Hy(1,j,nf)-Hy(IEz,j,nf))+(1/kappey)*(-Hx(1,j+1,nf)+Hx(1,j,nf))) * 1/delta + PsiEzX(1,j) + PsiEzY(1,j)); % PBC
    end
    Ez(:,:,nf) = Dz(:,:,nf)./e0;
    
    % ====================== Source ===================
    if SourcePlane == 1
        i = 1:IEz;
        j = SourceLocationY;        
    else
        i = SourceLocationX;
        j = SourceLocationY;
    end
    % Source.
    if SourceChoice == 1
        Ez(i,j,nf) = Ez(i,j,nf) + exp( -1*((n-td)/(PulseWidth/4))^2 ) * Sc;
    elseif SourceChoice == 2
        Ez(i,j,nf) = Ez(i,j,nf) + sin(2*pi*f*(n)*dt) * Sc;
    elseif SourceChoice == 3
        Ez(i,j,nf) = Ez(i,j,nf) + (1-2*(pi*fp*(n*dt-dr))^2)*exp(-1*(pi*fp*(n*dt-dr))^2) * Sc;
    end
    Dz(i,j,nf) = e0*Ez(i,j,nf);
    
    Ezi(n+1) = Ez(IEz/2,round((2*PMLw+J)/3)+1,nf);
    
    if (SaveFields == 1 && mod(n, SnapshotInterval) == 0)
        EzSnapshots(:,:,n/SnapshotInterval+1) = Ez(1+(0:SnapshotResolution:(IEz-1)), 1+(0:SnapshotResolution:(JEz-1)), nf);
    end

    
    np = mod(np, 3)+1;
    n0 = mod(n0, 3)+1;
    nf = mod(nf, 3)+1;
end
fprintf(1, repmat('\b',1,linecount));
fprintf ( 1, 'Dry run complete! \n');
toc
% Reinitialization of fields for actual simulation.
Ez = zeros(IEz, JEz, 3); % z-component of E-field
Dz = zeros(IEz, JEz, 3); % z-component of D
Hx = zeros(IHx, JHx, 3); % i-component of H-field
Bx = zeros(IHx, JHx, 3); % i-component of B
Hy = zeros(IHy, JHy, 3); % j-component of H-field
By = zeros(IHy, JHy, 3); % j-component of B

% Actual simulation with scatterer.
fprintf ( 1, 'Simulation started... \n');
np = 1;
n0 = 2;
nf = 3;
linecount = 0;
% Outer loop for time-stepping.
tic
% Test loop for incident field in free space.
for n = 0:MaxTime
    
    % Progress indicator.
    if mod(n,SnapshotInterval) == 0
        fprintf(1, repmat('\b',1,linecount));
        linecount = fprintf(1, '%g %%', (n*100)/MaxTime );
    end
    
    % ========================= Bx =============================
    % Hx Psi array.
    i=1:IHx;
    j=2:JHx-1;
    PsiHxY(i,j) = (Cmy/delta)*(-Ez(i,j,n0) + Ez(i,j-1,n0)) + bmy*PsiHxY(i,j);
    % Bx in normal space.
    j=(2+PMLw):(JHx-PMLw-1);
    Bx(i,j,nf) = Bx(i,j,n0) + (-Ez(i,j,n0) + Ez(i,j-1,n0)) * dt/delta;
    if PMLw > 0
        % Bx in lower PML layer.
        j=2:PMLw+1;
        Bx(i,j,nf) = Bx(i,j,n0) + dt*((1/kappmy)*(-Ez(i,j,n0) + Ez(i,j-1,n0)) * 1/delta + PsiHxY(i,j));
        % Bx in upper PML layer.
        j=JHx-PMLw:JHx-1;
        Bx(i,j,nf) = Bx(i,j,n0) + dt*((1/kappmy)*(-Ez(i,j,n0) + Ez(i,j-1,n0)) * 1/delta + PsiHxY(i,j));
    end
        
    % ========================= By =============================
    % Hy Psi array.
    i=1:IHy-1;
    j=1:JHy;
    PsiHyX(i,j) = (Cmx/delta)*(Ez(i+1,j,n0)-Ez(i,j,n0)) + bmx*PsiHyX(i,j);
    PsiHyX(IHy,j) = (Cmx/delta)*(Ez(1,j,n0)-Ez(IHy,j,n0)) + bmx*PsiHyX(IHy,j);
    % By in normal space.
    j=(1+PMLw):JHy-PMLw;
    By(i,j,nf) = By(i,j,n0) + (Ez(i+1,j,n0) - Ez(i,j,n0)) * dt/delta;
    By(IHy,j,nf) = By(IHy,j,n0) + (Ez(1,j,n0) - Ez(IHy,j,n0)) * dt/delta; % PBC
    if PMLw > 0
        % By in lower PML layer.
        i=1:IHy-1;
        j=1:PMLw;
        By(i,j,nf) = By(i,j,n0) + dt*((1/kappmx)*(Ez(i+1,j,n0) - Ez(i,j,n0)) * 1/delta + PsiHyX(i,j));
        By(IHy,j,nf) = By(IHy,j,n0) + dt*((1/kappmx)*(Ez(1,j,n0) - Ez(IHy,j,n0)) * 1/delta + PsiHyX(IHy,j)); % PBC
        % By in upper PML layer.
        i=1:IHy-1;
        j=JHy-PMLw+1:JHy;
        By(i,j,nf) = By(i,j,n0) + dt*((1/kappmx)*(Ez(i+1,j,n0) - Ez(i,j,n0)) * 1/delta + PsiHyX(i,j));
        By(IHy,j,nf) = By(IHy,j,n0) + dt*((1/kappmx)*(Ez(1,j,n0) - Ez(IHy,j,n0)) * 1/delta + PsiHyX(IHy,j)); % PBC
    end    
    
    % --- Sync ---
    
    % BxAve.
    i=3:IHy-3;
    j=(3+PMLw):JHy-PMLw-2;
    BxAve(i,j,nf) = (Bx(i,j,nf)+Bx(i+1,j,nf)+Bx(i,j+1,nf)+Bx(i+1,j+1,nf))/4;
    % ByAve.    
    i=3:IHx-2;
    j=(4+PMLw):(JHx-PMLw-3);
    ByAve(i,j,nf) = (By(i,j,nf)+By(i-1,j,nf)+By(i,j-1,nf)+By(i-1,j-1,nf))/4;
    
    % ========================= Hx =============================
    i=1:IHx;
    % Hx in normal space.
    j=(2+PMLw):(JHx-PMLw-1);
    %Hx(i,j,nf) = Bx(i,j,nf)./u0;
    Hx(i,j,nf) = (ax(i,j,1).*Bx(i,j,nf)+ax(i,j,2).*Bx(i,j,n0)+ax(i,j,3).*Bx(i,j,np)+ax(i,j,4).*ByAve(i,j,nf)+ax(i,j,5).*ByAve(i,j,n0)+ax(i,j,6).*ByAve(i,j,np)-ax(i,j,7).*Hx(i,j,n0)-ax(i,j,8).*Hx(i,j,np))./ax(i,j,9);
    if PMLw > 0
        % Hx in lower PML layer.
        j=2:PMLw+1;
        Hx(i,j,nf) = Bx(i,j,nf)./u0;
        % Hx in upper PML layer.
        j=JHx-PMLw:JHx-1;
        Hx(i,j,nf) = Bx(i,j,nf)./u0;
    end
    Hx(:,:,nf) = Hx(:,:,nf).*HxMask;
    Bx(:,:,nf) = Bx(:,:,nf).*HxMask;
    
    % ========================= Hy =============================
    % Hy in normal space.
    i=1:IHy;
    j=(1+PMLw):JHy-PMLw;
    %Hy(i,j,nf) = By(i,j,nf)./u0;
    Hy(i,j,nf) = (ay(i,j,1).*By(i,j,nf)+ay(i,j,2).*By(i,j,n0)+ay(i,j,3).*By(i,j,np)+ay(i,j,4).*BxAve(i,j,nf)+ay(i,j,5).*BxAve(i,j,n0)+ay(i,j,6).*BxAve(i,j,np)-ay(i,j,7).*Hy(i,j,n0)-ay(i,j,8).*Hy(i,j,np))./ay(i,j,9);
    if PMLw > 0
        % By in lower PML layer.
        i=1:IHy;
        j=1:PMLw;
        Hy(i,j,nf) = By(i,j,nf)./u0;
        % By in upper PML layer.
        i=1:IHy;
        j=JHy-PMLw+1:JHy;
        Hy(i,j,nf) = By(i,j,nf)./u0;
    end
    Hy(:,:,nf) = Hy(:,:,nf).*HyMask;
    By(:,:,nf) = By(:,:,nf).*HyMask;
    
    % ========================= Dz and Ez =============================
    % Psi arrays.
    i=2:IEz;
    j=1:JEz;
    PsiEzX(i,j) = (Cex/delta)*(Hy(i,j,nf)-Hy(i-1,j,nf)) + bex*PsiEzX(i,j);
    PsiEzX(1,j) = (Cex/delta)*(Hy(1,j,nf)-Hy(IEz,j,nf)) + bex*PsiEzX(1,j); % PBC
    PsiEzY(i,j) = (Cey/delta)*(-Hx(i,j+1,nf)+Hx(i,j,nf)) + bey*PsiEzY(i,j);
    PsiEzY(1,j) = (Cey/delta)*(-Hx(1,j+1,nf)+Hx(1,j,nf)) + bey*PsiEzY(1,j); % PBC
    % Dz in Normal Space.
    j=(1+PMLw):(JEz-PMLw);
    Dz(i,j,nf) = Dz(i,j,n0) + (Hy(i,j,nf)-Hy(i-1,j,nf)-Hx(i,j+1,nf)+Hx(i,j,nf)) * dt/delta;
    Dz(1,j,nf) = Dz(1,j,n0) + (Hy(1,j,nf)-Hy(IEz,j,nf)-Hx(1,j+1,nf)+Hx(1,j,nf)) * dt/delta; % PBC
    i=1:IEz;
    Ez(i,j,nf) = az(i,j,1).*(Dz(i,j,nf)-2*Dz(i,j,n0)+Dz(i,j,np))+az(i,j,2).*(Dz(i,j,nf)-Dz(i,j,np))+az(i,j,3).*(2*Ez(i,j,n0)-Ez(i,j,np))+az(i,j,4).*(2*Ez(i,j,n0)+Ez(i,j,np))+az(i,j,5).*Ez(i,j,np);
    if PMLw > 0
        % Dz in lower PML layer.
        i=2:IEz;
        j=1:PMLw;
        Dz(i,j,nf) = Dz(i,j,n0) + dt*(((1/kappex)*(Hy(i,j,nf)-Hy(i-1,j,nf))+(1/kappey)*(-Hx(i,j+1,nf)+Hx(i,j,nf))) * 1/delta + PsiEzX(i,j) + PsiEzY(i,j));
        Dz(1,j,nf) = Dz(1,j,n0) + dt*(((1/kappex)*(Hy(1,j,nf)-Hy(IEz,j,nf))+(1/kappey)*(-Hx(1,j+1,nf)+Hx(1,j,nf))) * 1/delta + PsiEzX(1,j) + PsiEzY(1,j)); % PBC
        i=1:IEz;
        Ez(i,j,nf) = Dz(i,j,nf)./e0;
        % Dz in upper PML layer.
        i=2:IEz;
        j=JEz-PMLw+1:JEz;
        Dz(i,j,nf) = Dz(i,j,n0) + dt*(((1/kappex)*(Hy(i,j,nf)-Hy(i-1,j,nf))+(1/kappey)*(-Hx(i,j+1,nf)+Hx(i,j,nf))) * 1/delta + PsiEzX(i,j) + PsiEzY(i,j));
        Dz(1,j,nf) = Dz(1,j,n0) + dt*(((1/kappex)*(Hy(1,j,nf)-Hy(IEz,j,nf))+(1/kappey)*(-Hx(1,j+1,nf)+Hx(1,j,nf))) * 1/delta + PsiEzX(1,j) + PsiEzY(1,j)); % PBC
        i=1:IEz;
        Ez(i,j,nf) = Dz(i,j,nf)./e0;
    end
            
    % ====================== Source ===================
    if SourcePlane == 1
        i = 1:IEz;
        j = SourceLocationY;        
    else
        i = SourceLocationX;
        j = SourceLocationY;
    end
    % Source.
    if SourceChoice == 1
        Ez(i,j,nf) = Ez(i,j,nf) + exp( -1*((n-td)/(PulseWidth/4))^2 ) * Sc;
    elseif SourceChoice == 2
        Ez(i,j,nf) = Ez(i,j,nf) + sin(2*pi*f*(n)*dt) * Sc;
    elseif SourceChoice == 3
        Ez(i,j,nf) = Ez(i,j,nf) + (1-2*(pi*fp*(n*dt-dr))^2)*exp(-1*(pi*fp*(n*dt-dr))^2) * Sc;
    end
    Dz(i,j,nf) = e0*Ez(i,j,nf);
    %Ez(:,:,nf) = Ez(:,:,nf).*EzMask;
    %Dz(:,:,nf) = Dz(:,:,nf).*EzMask;
    % Transmitted fields.
    Ezt(n+1) = Ez(IEz/2,round((2*PMLw+J)/3),nf);
    Eztt(n+1) = Ez(IEz/2,round((2*PMLw+J)/3)+10,nf);
    
    % Fields for calculation of refractive index.
    Ezy1(n+1) = Ez(IEz/2,Y1,nf);
    Ezy2(n+1) = Ez(IEz/2,Y2, nf);
    
    if (SaveFields == 1 && mod(n, SnapshotInterval) == 0)
        EzSnapshots(:,:,n/SnapshotInterval+1) = Ez(1+(0:SnapshotResolution:(IEz-1)), 1+(0:SnapshotResolution:(JEz-1)), nf);
    end
    
    np = mod(np, 3)+1;
    n0 = mod(n0, 3)+1;
    nf = mod(nf, 3)+1;
end
fprintf(1, repmat('\b',1,linecount));
fprintf ( 1, 'Simulation complete! \n');
toc
% Postprocessing.
Fs = 1/dt;                    % Sampling frequency
T = dt;                       % Sample time
L = length(Ezi);              % Length of signal
t = (0:L-1)*T;                % Time vector
fspan = 100;                  % Points to plot in frequency domain

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
ylim([-2 2])
axis 'auto i'
grid on
subplot(212)
plot(f(1:fspan), 1-TAU(1:fspan), 'LineWidth', 2.0, 'Color', 'b')
set(gca, 'FontSize', 10, 'FontWeight', 'b')
title('Reflection Coefficient', 'FontSize', 12, 'FontWeight', 'b')
xlabel('Frenuency (Hz)', 'FontSize', 11, 'FontWeight', 'b')
ylabel('1-|EZT(f)/EZI(f)|', 'FontSize', 11, 'FontWeight', 'b')
ylim([-2 2])
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

if SaveFields == 1
    % Electric field snapshots.
    sizeS=size(EzSnapshots);
    for i=1:(MaxTime/SnapshotInterval)-1

        figure (6)
        mesh ( EzSnapshots (:, :, i) );
        view (4, 4)
        xlim([0 sizeS(2)])
        ylim([0 sizeS(1)])
        zlim([-1 1])
        %caxis([-0.1 0.6])
        xlabel ('j-axis')
        ylabel ('i-axis')
        %colorbar

        figure (7)
        mesh ( EzSnapshots (:, :, i) );
        view (0, 90)
        xlim([0 sizeS(2)])
        ylim([0 sizeS(1)])
        zlim([-10 10])
        %caxis([-0.1 0.6])
        xlabel ('j-axis')
        ylabel ('i-axis')
        %colorbar

    end
end
