clc
clear all

% Simulation parameters.
I = 256; % No. of spatial steps in i direction.
J = 256; % No. of spatial steps in j direction.
PMLw = 64; % Width of PML layer.
SlabLeft = round(J/3+PMLw); % Location of left end of Slab.
SlabRight = round(2*J/3+PMLw); % Location of right end of Slab
MaxTime = 4*J; % No. of time steps
PulseWidth = round(J/8); % Controls width of Gaussian Pulse
td = PulseWidth; % Temporal delay in pulse.
SaveFields = 1; % 0. No, 1. Yes.
SnapshotResolution = 1; % Snapshot resolution. 1 is best.
SnapshotInterval = 16; % Amount of time delay between snaps.
% Choice of source.
% 1. Gaussian 2. Sine wave 3. Ricker wavelet
SourceChoice = 1;
SourcePlane = 1; % Is the source a plane wave. 0. = Omni 1. Plane-wave.
SourceLocationX = I/2; % X Location of source. Only used for an omni-source.
SourceLocationY = PMLw+6; % Y Location of source.

% Constants.
c = 3e8;
pi = 3.141592654;
e0 = (1e-9)/(36*pi);
u0 = (1e-7)*4*pi;

dt = 0.5e-11;
delta = 3e-3;
Sc = c * dt/delta

l = PulseWidth*delta;
f = c/l
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
Hx = zeros(IHx, JHx, 3); % i-component of H-field
Bx = zeros(IHx, JHx, 3); % i-component of B
Hy = zeros(IHy, JHy, 3); % j-component of H-field
By = zeros(IHy, JHy, 3); % j-component of B

% Incident and Transmitted Fields.
Ezi = zeros(MaxTime, 1);
Ezt = zeros(MaxTime, 1);
Eztt = zeros(MaxTime, 1);
x1 = SlabLeft+1; % Position of observation.

% Refractive Index calculations.
Y1 = SlabLeft + 5;
y1 = Y1*delta;
Y2 = SlabLeft + 6;
y2 = Y2*delta;
Ezy1 = zeros(MaxTime, 1);
Ezy2 = zeros(MaxTime, 1);

einf = ones(IHx,JHx);
einf(:,SlabLeft:SlabRight) = 1; % einf(Drude) or er in slab.
uinf = ones(IHx,JHx);
uinf(:,SlabLeft:SlabRight) = 1; % uinf(Drude) or ur in slab.
wpesq = zeros(IHx,JHx);
wpesq(:,SlabLeft:SlabRight) = 2*w^2; % DNG(Drude) value of wpe squared in slab.
wpmsq = zeros(IHx,JHx);
wpmsq(:,SlabLeft:SlabRight) = 2*w^2; % DNG(Drude) value of wpm squared in slab.
ge = zeros(IHx,JHx);
%ge(:,SlabLeft:SlabRight) = w/32; % Electric collision frequency in slab.
gm = zeros(IHx,JHx);
%gm(:,SlabLeft:SlabRight) = w/32; % Magnetic collision frequency in slab.

ae0 = (4*dt^2)./(e0*(4*einf+dt^2*wpesq+2*dt*einf.*ge));
ae = (1/dt^2)*ae0;
be = (1/(2*dt))*ge.*ae0;
ce = (e0/dt^2)*einf.*ae0;
de = (-1*e0/4).*wpesq.*ae0;
ee = (1/(2*dt))*e0*einf.*ge.*ae0;
am0 = (4*dt^2)./(u0*(4*uinf+dt^2*wpmsq+2*dt*uinf.*gm));
am = (1/dt^2)*am0;
bm = (1/(2*dt))*gm.*am0;
cm = (u0/dt^2)*uinf.*am0;
dm = (-1*u0/4).*wpmsq.*am0;
em = (1/(2*dt))*u0*uinf.*gm.*am0;

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
fprintf ( 1, '\rDry run started! \n');
tic
% Test loop for incident field in free space.
for n = 0:MaxTime
    
    % Progress indicator.
    if mod(n,2) == 0
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
    Hx(:,:,nf) = Bx(:,:,nf)./(u0*uinf);
    
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
    Hy(:,:,nf) = By(:,:,nf)./(u0*uinf(:,1:JHy));
    
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
    Ez(:,:,nf) = Dz(:,:,nf)./(e0*einf(:,1:JEz));
    
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
    
    Ezi(n+1) = Ez(IEz/2,SlabLeft+1,nf); % Incident field is left of slab.
    
    if (SaveFields == 1 && mod(n, SnapshotInterval) == 0)
        EzSnapshots(:,:,n/SnapshotInterval+1) = Ez(1+(0:SnapshotResolution:(IEz-1)), 1+(0:SnapshotResolution:(JEz-1)), nf);
    end

    
    np = mod(np, 3)+1;
    n0 = mod(n0, 3)+1;
    nf = mod(nf, 3)+1;
end
fprintf ( 1, '\rDry run complete! \n');
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
    if mod(n,2) == 0
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
    Hx(i,j,nf) = am(i,j).*(Bx(i,j,nf)-2*Bx(i,j,n0)+Bx(i,j,np))+bm(i,j).*(Bx(i,j,nf)-Bx(i,j,np))+cm(i,j).*(2*Hx(i,j,n0)-Hx(i,j,np))+dm(i,j).*(2*Hx(i,j,n0)+Hx(i,j,np))+em(i,j).*Hx(i,j,np);
    if PMLw > 0
        % Bx in lower PML layer.
        j=2:PMLw+1;
        Bx(i,j,nf) = Bx(i,j,n0) + dt*((1/kappmy)*(-Ez(i,j,n0) + Ez(i,j-1,n0)) * 1/delta + PsiHxY(i,j));
        Hx(i,j,nf) = Bx(i,j,nf)./(u0*uinf(i,j));
        % Bx in upper PML layer.
        j=JHx-PMLw:JHx-1;
        Bx(i,j,nf) = Bx(i,j,n0) + dt*((1/kappmy)*(-Ez(i,j,n0) + Ez(i,j-1,n0)) * 1/delta + PsiHxY(i,j));
        Hx(i,j,nf) = Bx(i,j,nf)./(u0*uinf(i,j));
    end
    
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
    i=1:IHy;
    Hy(i,j,nf) = am(i,j).*(By(i,j,nf)-2*By(i,j,n0)+By(i,j,np))+bm(i,j).*(By(i,j,nf)-By(i,j,np))+cm(i,j).*(2*Hy(i,j,n0)-Hy(i,j,np))+dm(i,j).*(2*Hy(i,j,n0)+Hy(i,j,np))+em(i,j).*Hy(i,j,np);
    if PMLw > 0
        % By in lower PML layer.
        i=1:IHy-1;
        j=1:PMLw;
        By(i,j,nf) = By(i,j,n0) + dt*((1/kappmx)*(Ez(i+1,j,n0) - Ez(i,j,n0)) * 1/delta + PsiHyX(i,j));
        By(IHy,j,nf) = By(IHy,j,n0) + dt*((1/kappmx)*(Ez(1,j,n0) - Ez(IHy,j,n0)) * 1/delta + PsiHyX(IHy,j)); % PBC
        i=1:IHy;
        Hy(i,j,nf) = By(i,j,nf)./(u0*uinf(i,j));
        % By in upper PML layer.
        i=1:IHy-1;
        j=JHy-PMLw+1:JHy;
        By(i,j,nf) = By(i,j,n0) + dt*((1/kappmx)*(Ez(i+1,j,n0) - Ez(i,j,n0)) * 1/delta + PsiHyX(i,j));
        By(IHy,j,nf) = By(IHy,j,n0) + dt*((1/kappmx)*(Ez(1,j,n0) - Ez(IHy,j,n0)) * 1/delta + PsiHyX(IHy,j)); % PBC
        i=1:IHy;
        Hy(i,j,nf) = By(i,j,nf)./(u0*uinf(i,j));
    end    
    
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
    Ez(i,j,nf) = ae(i,j).*(Dz(i,j,nf)-2*Dz(i,j,n0)+Dz(i,j,np))+be(i,j).*(Dz(i,j,nf)-Dz(i,j,np))+ce(i,j).*(2*Ez(i,j,n0)-Ez(i,j,np))+de(i,j).*(2*Ez(i,j,n0)+Ez(i,j,np))+ee(i,j).*Ez(i,j,np);
    if PMLw > 0
        % Dz in lower PML layer.
        i=2:IEz;
        j=1:PMLw;
        Dz(i,j,nf) = Dz(i,j,n0) + dt*(((1/kappex)*(Hy(i,j,nf)-Hy(i-1,j,nf))+(1/kappey)*(-Hx(i,j+1,nf)+Hx(i,j,nf))) * 1/delta + PsiEzX(i,j) + PsiEzY(i,j));
        Dz(1,j,nf) = Dz(1,j,n0) + dt*(((1/kappex)*(Hy(1,j,nf)-Hy(IEz,j,nf))+(1/kappey)*(-Hx(1,j+1,nf)+Hx(1,j,nf))) * 1/delta + PsiEzX(1,j) + PsiEzY(1,j)); % PBC
        i=1:IEz;
        Ez(i,j,nf) = Dz(i,j,nf)./(e0*einf(i,j));
        % Dz in upper PML layer.
        i=2:IEz;
        j=JEz-PMLw+1:JEz;
        Dz(i,j,nf) = Dz(i,j,n0) + dt*(((1/kappex)*(Hy(i,j,nf)-Hy(i-1,j,nf))+(1/kappey)*(-Hx(i,j+1,nf)+Hx(i,j,nf))) * 1/delta + PsiEzX(i,j) + PsiEzY(i,j));
        Dz(1,j,nf) = Dz(1,j,n0) + dt*(((1/kappex)*(Hy(1,j,nf)-Hy(IEz,j,nf))+(1/kappey)*(-Hx(1,j+1,nf)+Hx(1,j,nf))) * 1/delta + PsiEzX(1,j) + PsiEzY(1,j)); % PBC
        i=1:IEz;
        Ez(i,j,nf) = Dz(i,j,nf)./(e0*einf(i,j));
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
    
    % Transmitted fields.
    Ezt(n+1) = Ez(IEz/2,SlabLeft+1,nf);
    Eztt(n+1) = Ez(IEz/2,SlabRight+10,nf);
    
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
fprintf ( 1, '\rSimulation complete! \n');
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
