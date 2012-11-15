clc
clear all

% Simulation parameters.
SizeI = 512; % No. of spatial steps in x direction.
SizeJ = 512; % No. of spatial steps in y direction.
PMLw = 50; % Width of PML layer.
SlabLeft = round(SizeJ/3+PMLw); % Location of left end of Slab.
SlabRight = round(2*SizeJ/3+PMLw); % Location of right end of Slab
MaxTime = 4*SizeJ; % No. of time steps
PulseWidth = round(SizeJ/8); % Controls width of Gaussian Pulse
td = PulseWidth; % Temporal delay in pulse.
SnapshotResolution = 1; % Snapshot resolution. 1 is best.
SnapshotInterval = 4; % Amount of time delay between snaps.
% Choice of source.
% 1. Gaussian 2. Sine wave 3. Ricker wavelet
SourceChoice = 1;
SourcePlane = 1; % Is the source a plane wave. 0. = Omni 1. Plane-wave.
SourceLocationX = SizeI/2; % X Location of source. Only used for an omni-source.
SourceLocationY = PMLw+10; % Y Location of source.

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
    fp = f; % Peak frenuency
    dr = PulseWidth*dt*2; % Delay
end
% PML arrays.
PsiEzX = zeros(SizeI, SizeJ+2*PMLw);
PsiEzY = zeros(SizeI, SizeJ+2*PMLw);
PsiHyX = zeros(SizeI, SizeJ+2*PMLw);
PsiHxY = zeros(SizeI, SizeJ+2*PMLw+1);

% PML parameters.
kapp = 1;
a = 0.0004;
sig = 0.045;
% Electric.
kappex = kapp;
kappey = kapp;
aex = a;
aey = a;
sigex = 0;
sigey = sig;
bex = exp(-1*(aex/e0+sigex/(kappex*e0))*dt);
bey = exp(-1*(aey/e0+sigey/(kappey*e0))*dt);
Cex = (bex-1)*sigex/(sigex*kappex+kappex^2*aex);
Cey = (bey-1)*sigey/(sigey*kappey+kappey^2*aey);
% Magnetic.
kappmx = kapp;
kappmy = kapp;
amx = a;
amy = a;
sigmx = 0;
sigmy = u0/e0*sig;
bmx = exp(-1*(amx/u0+sigmx/(kappmx*u0))*dt);
bmy = exp(-1*(amy/u0+sigmy/(kappmy*u0))*dt);
Cmx = (bmx-1)*sigmx/(sigmx*kappmx+kappmx^2*amx);
Cmy = (bmy-1)*sigmy/(sigmy*kappmy+kappmy^2*amy);

% Initialization.
Ez = zeros(SizeI, SizeJ+2*PMLw, 3); % z-component of E-field
Dz = zeros(SizeI, SizeJ+2*PMLw, 3); % z-component of D
Hx = zeros(SizeI, SizeJ+2*PMLw+1, 3); % x-component of H-field
Bx = zeros(SizeI, SizeJ+2*PMLw+1, 3); % x-component of B
Hy = zeros(SizeI, SizeJ+2*PMLw, 3); % y-component of H-field
By = zeros(SizeI, SizeJ+2*PMLw, 3); % y-component of B

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

einf = ones(SizeI,SizeJ+2*PMLw+1);
einf(:,SlabLeft:SlabRight) = 1; % einf(Drude) or er in slab.
uinf = ones(SizeI,SizeJ+2*PMLw+1);
uinf(:,SlabLeft:SlabRight) = 1; % uinf(Drude) or ur in slab.
wpesn = zeros(SizeI,SizeJ+2*PMLw+1);
wpesn(:,SlabLeft:SlabRight) = 2*w^2; % DNG(Drude) value of wpe snuared in slab.
wpmsn = zeros(SizeI,SizeJ+2*PMLw+1);
wpmsn(:,SlabLeft:SlabRight) = 2*w^2; % DNG(Drude) value of wpm snuared in slab.
ge = zeros(SizeI,SizeJ+2*PMLw+1);
ge(:,SlabLeft:SlabRight) = w/32; % Electric collision frenuency in slab.
gm = zeros(SizeI,SizeJ+2*PMLw+1);
gm(:,SlabLeft:SlabRight) = w/32; % Magnetic collision frenuency in slab.

ae0 = (4*dt^2)./(e0*(4*einf+dt^2*wpesn+2*dt*einf.*ge));
ae = (1/dt^2)*ae0;
be = (1/(2*dt))*ge.*ae0;
ce = (e0/dt^2)*einf.*ae0;
de = (-1*e0/4).*wpesn.*ae0;
ee = (1/(2*dt))*e0*einf.*ge.*ae0;
am0 = (4*dt^2)./(u0*(4*uinf+dt^2*wpmsn+2*dt*uinf.*gm));
am = (1/dt^2)*am0;
bm = (1/(2*dt))*gm.*am0;
cm = (u0/dt^2)*uinf.*am0;
dm = (-1*u0/4).*wpmsn.*am0;
em = (1/(2*dt))*u0*uinf.*gm.*am0;

EzSnapshots = zeros(SizeI/SnapshotResolution, (SizeJ+2*PMLw)/SnapshotResolution, MaxTime/SnapshotInterval); % Data for plotting.
frame = 1;

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
    x=1:SizeI;
    y=2:SizeJ+2*PMLw;
    PsiHxY(x,y) = (Cmy/delta)*(-Ez(x,y,n0) + Ez(x,y-1,n0)) + bmy*PsiHxY(x,y);
    % Bx in normal space.
    y=(2+PMLw):((SizeJ+2*PMLw+1)-PMLw-1);
    Bx(x,y,nf) = Bx(x,y,n0) + (-Ez(x,y,n0) + Ez(x,y-1,n0)) * dt/delta;
    if PMLw > 0
        % Bx in lower PML layer.
        y=2:PMLw+1;
        Bx(x,y,nf) = Bx(x,y,n0) + dt*((1/kappmy)*(-Ez(x,y,n0) + Ez(x,y-1,n0)) * 1/delta + PsiHxY(x,y));
        % Bx in upper PML layer.
        y=(SizeJ+2*PMLw+1)-PMLw:(SizeJ+2*PMLw);
        Bx(x,y,nf) = Bx(x,y,n0) + dt*((1/kappmy)*(-Ez(x,y,n0) + Ez(x,y-1,n0)) * 1/delta + PsiHxY(x,y));
    end
    Hx(:,:,nf) = Bx(:,:,nf)./(u0*uinf);
    
    % ========================= By and Hy =============================
    % Hy Psi array.
    x=1:SizeI-1;
    y=1:SizeJ+2*PMLw;
    PsiHyX(x,y) = (Cmx/delta)*(Ez(x+1,y,n0)-Ez(x,y,n0)) + bmx*PsiHyX(x,y);
    PsiHyX(SizeI,y) = (Cmx/delta)*(Ez(1,y,n0)-Ez(SizeI,y,n0)) + bmx*PsiHyX(SizeI,y);
    % By in normal space.
    y=(1+PMLw):(SizeJ+2*PMLw)-PMLw;
    By(x,y,nf) = By(x,y,n0) + (Ez(x+1,y,n0) - Ez(x,y,n0)) * dt/delta;
    By(SizeI,y,nf) = By(SizeI,y,n0) + (Ez(1,y,n0) - Ez(SizeI,y,n0)) * dt/delta; % PBC
    if PMLw > 0
        % By in lower PML layer.
        y=1:PMLw;
        By(x,y,nf) = By(x,y,n0) + dt*((1/kappmx)*(Ez(x+1,y,n0) - Ez(x,y,n0)) * 1/delta + PsiHyX(x,y));
        By(SizeI,y,nf) = By(SizeI,y,n0) + dt*((1/kappmx)*(Ez(1,y,n0) - Ez(SizeI,y,n0)) * 1/delta + PsiHyX(SizeI,y)); % PBC
        % By in upper PML layer.
        y=(SizeJ+2*PMLw)-PMLw+1:(SizeJ+2*PMLw);
        By(x,y,nf) = By(x,y,n0) + dt*((1/kappmx)*(Ez(x+1,y,n0) - Ez(x,y,n0)) * 1/delta + PsiHyX(x,y));
        By(SizeI,y,nf) = By(SizeI,y,n0) + dt*((1/kappmx)*(Ez(1,y,n0) - Ez(SizeI,y,n0)) * 1/delta + PsiHyX(SizeI,y)); % PBC
    end
    Hy(:,:,nf) = By(:,:,nf)./(u0*uinf(:,1:SizeJ+2*PMLw));
    
    % ========================= Dz and Ez =============================
    % Psi arrays.
    x=2:SizeI;
    y=1:SizeJ+2*PMLw;
    PsiEzX(x,y) = (Cex/delta)*(Hy(x,y,nf)-Hy(x-1,y,nf)) + bex*PsiEzX(x,y);
    PsiEzX(1,y) = (Cex/delta)*(Hy(1,y,nf)-Hy(SizeI,y,nf)) + bex*PsiEzX(1,y); % PBC
    PsiEzY(x,y) = (Cey/delta)*(-Hx(x,y+1,nf)+Hx(x,y,nf)) + bey*PsiEzY(x,y);
    PsiEzY(1,y) = (Cey/delta)*(-Hx(1,y+1,nf)+Hx(1,y,nf)) + bey*PsiEzY(1,y); % PBC
    % Dz in Normal Space.
    y=(1+PMLw):((SizeJ+2*PMLw)-PMLw);
    Dz(x,y,nf) = Dz(x,y,n0) + (Hy(x,y,nf)-Hy(x-1,y,nf)-Hx(x,y+1,nf)+Hx(x,y,nf)) * dt/delta;
    Dz(1,y,nf) = Dz(1,y,n0) + (Hy(1,y,nf)-Hy(SizeI,y,nf)-Hx(1,y+1,nf)+Hx(1,y,nf)) * dt/delta; % PBC
    if PMLw > 0
        % Dz in lower PML layer.
        y=1:PMLw;
        Dz(x,y,nf) = Dz(x,y,n0) + dt*(((1/kappex)*(Hy(x,y,nf)-Hy(x-1,y,nf))+(1/kappey)*(-Hx(x,y+1,nf)+Hx(x,y,nf))) * 1/delta + PsiEzX(x,y) + PsiEzY(x,y));
        Dz(1,y,nf) = Dz(1,y,n0) + dt*(((1/kappex)*(Hy(1,y,nf)-Hy(SizeI,y,nf))+(1/kappey)*(-Hx(1,y+1,nf)+Hx(1,y,nf))) * 1/delta + PsiEzX(1,y) + PsiEzY(1,y)); % PBC
        % Dz in upper PML layer.
        y=(SizeJ+2*PMLw)-PMLw+1:(SizeJ+2*PMLw);
        Dz(x,y,nf) = Dz(x,y,n0) + dt*(((1/kappex)*(Hy(x,y,nf)-Hy(x-1,y,nf))+(1/kappey)*(-Hx(x,y+1,nf)+Hx(x,y,nf))) * 1/delta + PsiEzX(x,y) + PsiEzY(x,y));
        Dz(1,y,nf) = Dz(1,y,n0) + dt*(((1/kappex)*(Hy(1,y,nf)-Hy(SizeI,y,nf))+(1/kappey)*(-Hx(1,y+1,nf)+Hx(1,y,nf))) * 1/delta + PsiEzX(1,y) + PsiEzY(1,y)); % PBC
    end
    Ez(:,:,nf) = Dz(:,:,nf)./(e0*einf(:,1:SizeJ+2*PMLw));
    
    % ====================== Source ===================
    if SourcePlane == 1
        x = 1:SizeI;
        y = SourceLocationY;        
    else
        x = SourceLocationX;
        y = SourceLocationY;
    end
    % Source.
    if SourceChoice == 1
        Ez(x,y,nf) = Ez(x,y,nf) + exp( -1*((n-td)/(PulseWidth/4))^2 ) * Sc;
    elseif SourceChoice == 2
        Ez(x,y,nf) = Ez(x,y,nf) + sin(2*pi*f*(n)*dt) * Sc;
    elseif SourceChoice == 3
        Ez(x,y,nf) = Ez(x,y,nf) + (1-2*(pi*fp*(n*dt-dr))^2)*exp(-1*(pi*fp*(n*dt-dr))^2) * Sc;
    end
    Dz(x,y,nf) = e0*Ez(x,y,nf);
    
    Ezi(n+1) = Ez(SizeI/2,SlabLeft+1,nf); % Incident field is left of slab.
    
    if (mod(n, SnapshotInterval) == 0)
        EzSnapshots(:,:,n/SnapshotInterval+1) = Ez(1+(0:SnapshotResolution:(SizeI-1)), 1+(0:SnapshotResolution:((SizeJ+2*PMLw)-1)), nf);
    end

    
    np = mod(np, 3)+1;
    n0 = mod(n0, 3)+1;
    nf = mod(nf, 3)+1;
end
fprintf ( 1, '\rDry run complete! \n');
toc
% Reinitialization of fields for actual simulation.
Ez = zeros(SizeI, SizeJ+2*PMLw, 3); % z-component of E-field
Dz = zeros(SizeI, SizeJ+2*PMLw, 3); % z-component of D
Hx = zeros(SizeI, SizeJ+2*PMLw+1, 3); % x-component of H-field
Bx = zeros(SizeI, SizeJ+2*PMLw+1, 3); % x-component of B
Hy = zeros(SizeI, SizeJ+2*PMLw, 3); % y-component of H-field
By = zeros(SizeI, SizeJ+2*PMLw, 3); % y-component of B

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
    x=1:SizeI;
    y=2:SizeJ+2*PMLw;
    PsiHxY(x,y) = (Cmy/delta)*(-Ez(x,y,n0) + Ez(x,y-1,n0)) + bmy*PsiHxY(x,y);
    % Bx in normal space.
    y=(2+PMLw):((SizeJ+2*PMLw+1)-PMLw-1);
    Bx(x,y,nf) = Bx(x,y,n0) + (-Ez(x,y,n0) + Ez(x,y-1,n0)) * dt/delta;
    Hx(x,y,nf) = am(x,y).*(Bx(x,y,nf)-2*Bx(x,y,n0)+Bx(x,y,np))+bm(x,y).*(Bx(x,y,nf)-Bx(x,y,np))+cm(x,y).*(2*Hx(x,y,n0)-Hx(x,y,np))+dm(x,y).*(2*Hx(x,y,n0)+Hx(x,y,np))+em(x,y).*Hx(x,y,np);
    if PMLw > 0
        % Bx in lower PML layer.
        y=2:PMLw+1;
        Bx(x,y,nf) = Bx(x,y,n0) + dt*((1/kappmy)*(-Ez(x,y,n0) + Ez(x,y-1,n0)) * 1/delta + PsiHxY(x,y));
        Hx(x,y,nf) = Bx(x,y,nf)./(u0*uinf(x,y));
        % Bx in upper PML layer.
        y=(SizeJ+2*PMLw+1)-PMLw:(SizeJ+2*PMLw);
        Bx(x,y,nf) = Bx(x,y,n0) + dt*((1/kappmy)*(-Ez(x,y,n0) + Ez(x,y-1,n0)) * 1/delta + PsiHxY(x,y));
        Hx(x,y,nf) = Bx(x,y,nf)./(u0*uinf(x,y));
    end
    
    % ========================= By and Hy =============================
    % Hy Psi array.
    x=1:SizeI-1;
    y=1:SizeJ+2*PMLw;
    PsiHyX(x,y) = (Cmx/delta)*(Ez(x+1,y,n0)-Ez(x,y,n0)) + bmx*PsiHyX(x,y);
    PsiHyX(SizeI,y) = (Cmx/delta)*(Ez(1,y,n0)-Ez(SizeI,y,n0)) + bmx*PsiHyX(SizeI,y);
    % By in normal space.
    y=(1+PMLw):(SizeJ+2*PMLw)-PMLw;
    By(x,y,nf) = By(x,y,n0) + (Ez(x+1,y,n0) - Ez(x,y,n0)) * dt/delta;
    By(SizeI,y,nf) = By(SizeI,y,n0) + (Ez(1,y,n0) - Ez(SizeI,y,n0)) * dt/delta; % PBC
    x=1:SizeI;
    Hy(x,y,nf) = am(x,y).*(By(x,y,nf)-2*By(x,y,n0)+By(x,y,np))+bm(x,y).*(By(x,y,nf)-By(x,y,np))+cm(x,y).*(2*Hy(x,y,n0)-Hy(x,y,np))+dm(x,y).*(2*Hy(x,y,n0)+Hy(x,y,np))+em(x,y).*Hy(x,y,np);
    if PMLw > 0
        % By in lower PML layer.
        x=1:SizeI-1;
        y=1:PMLw;
        By(x,y,nf) = By(x,y,n0) + dt*((1/kappmx)*(Ez(x+1,y,n0) - Ez(x,y,n0)) * 1/delta + PsiHyX(x,y));
        By(SizeI,y,nf) = By(SizeI,y,n0) + dt*((1/kappmx)*(Ez(1,y,n0) - Ez(SizeI,y,n0)) * 1/delta + PsiHyX(SizeI,y)); % PBC
        x=1:SizeI;
        Hy(x,y,nf) = By(x,y,nf)./(u0*uinf(x,y));
        % By in upper PML layer.
        x=1:SizeI-1;
        y=(SizeJ+2*PMLw)-PMLw+1:(SizeJ+2*PMLw);
        By(x,y,nf) = By(x,y,n0) + dt*((1/kappmx)*(Ez(x+1,y,n0) - Ez(x,y,n0)) * 1/delta + PsiHyX(x,y));
        By(SizeI,y,nf) = By(SizeI,y,n0) + dt*((1/kappmx)*(Ez(1,y,n0) - Ez(SizeI,y,n0)) * 1/delta + PsiHyX(SizeI,y)); % PBC
        x=1:SizeI;
        Hy(x,y,nf) = By(x,y,nf)./(u0*uinf(x,y));
    end    
    
    % ========================= Dz and Ez =============================
    % Psi arrays.
    x=2:SizeI;
    y=1:SizeJ+2*PMLw;
    PsiEzX(x,y) = (Cex/delta)*(Hy(x,y,nf)-Hy(x-1,y,nf)) + bex*PsiEzX(x,y);
    PsiEzX(1,y) = (Cex/delta)*(Hy(1,y,nf)-Hy(SizeI,y,nf)) + bex*PsiEzX(1,y); % PBC
    PsiEzY(x,y) = (Cey/delta)*(-Hx(x,y+1,nf)+Hx(x,y,nf)) + bey*PsiEzY(x,y);
    PsiEzY(1,y) = (Cey/delta)*(-Hx(1,y+1,nf)+Hx(1,y,nf)) + bey*PsiEzY(1,y); % PBC
    % Dz in Normal Space.
    y=(1+PMLw):((SizeJ+2*PMLw)-PMLw);
    Dz(x,y,nf) = Dz(x,y,n0) + (Hy(x,y,nf)-Hy(x-1,y,nf)-Hx(x,y+1,nf)+Hx(x,y,nf)) * dt/delta;
    Dz(1,y,nf) = Dz(1,y,n0) + (Hy(1,y,nf)-Hy(SizeI,y,nf)-Hx(1,y+1,nf)+Hx(1,y,nf)) * dt/delta; % PBC
    x=1:SizeI;
    Ez(x,y,nf) = ae(x,y).*(Dz(x,y,nf)-2*Dz(x,y,n0)+Dz(x,y,np))+be(x,y).*(Dz(x,y,nf)-Dz(x,y,np))+ce(x,y).*(2*Ez(x,y,n0)-Ez(x,y,np))+de(x,y).*(2*Ez(x,y,n0)+Ez(x,y,np))+ee(x,y).*Ez(x,y,np);
    if PMLw > 0
        % Dz in lower PML layer.
        x=2:SizeI;
        y=1:PMLw;
        Dz(x,y,nf) = Dz(x,y,n0) + dt*(((1/kappex)*(Hy(x,y,nf)-Hy(x-1,y,nf))+(1/kappey)*(-Hx(x,y+1,nf)+Hx(x,y,nf))) * 1/delta + PsiEzX(x,y) + PsiEzY(x,y));
        Dz(1,y,nf) = Dz(1,y,n0) + dt*(((1/kappex)*(Hy(1,y,nf)-Hy(SizeI,y,nf))+(1/kappey)*(-Hx(1,y+1,nf)+Hx(1,y,nf))) * 1/delta + PsiEzX(1,y) + PsiEzY(1,y)); % PBC
        x=1:SizeI;
        Ez(x,y,nf) = Dz(x,y,nf)./(e0*einf(x,y));
        % Dz in upper PML layer.
        x=2:SizeI;
        y=(SizeJ+2*PMLw)-PMLw+1:(SizeJ+2*PMLw);
        Dz(x,y,nf) = Dz(x,y,n0) + dt*(((1/kappex)*(Hy(x,y,nf)-Hy(x-1,y,nf))+(1/kappey)*(-Hx(x,y+1,nf)+Hx(x,y,nf))) * 1/delta + PsiEzX(x,y) + PsiEzY(x,y));
        Dz(1,y,nf) = Dz(1,y,n0) + dt*(((1/kappex)*(Hy(1,y,nf)-Hy(SizeI,y,nf))+(1/kappey)*(-Hx(1,y+1,nf)+Hx(1,y,nf))) * 1/delta + PsiEzX(1,y) + PsiEzY(1,y)); % PBC
        x=1:SizeI;
        Ez(x,y,nf) = Dz(x,y,nf)./(e0*einf(x,y));
    end
            
    % ====================== Source ===================
    if SourcePlane == 1
        x = 1:SizeI;
        y = SourceLocationY;        
    else
        x = SourceLocationX;
        y = SourceLocationY;
    end
    % Source.
    if SourceChoice == 1
        Ez(x,y,nf) = Ez(x,y,nf) + exp( -1*((n-td)/(PulseWidth/4))^2 ) * Sc;
    elseif SourceChoice == 2
        Ez(x,y,nf) = Ez(x,y,nf) + sin(2*pi*f*(n)*dt) * Sc;
    elseif SourceChoice == 3
        Ez(x,y,nf) = Ez(x,y,nf) + (1-2*(pi*fp*(n*dt-dr))^2)*exp(-1*(pi*fp*(n*dt-dr))^2) * Sc;
    end
    Dz(x,y,nf) = e0*Ez(x,y,nf);
    
    % Transmitted fields.
    Ezt(n+1) = Ez(SizeI/2,SlabLeft+1,nf);
    Eztt(n+1) = Ez(SizeI/2,SlabRight+10,nf);
    
    % Fields for calculation of refractive index.
    Ezy1(n+1) = Ez(SizeI/2,Y1,nf);
    Ezy2(n+1) = Ez(SizeI/2,Y2, nf);
    
    if (mod(n, SnapshotInterval) == 0)
        EzSnapshots(:,:,n/SnapshotInterval+1) = Ez(1+(0:SnapshotResolution:(SizeI-1)), 1+(0:SnapshotResolution:((SizeJ+2*PMLw)-1)), nf);
    end
    
    np = mod(np, 3)+1;
    n0 = mod(n0, 3)+1;
    nf = mod(nf, 3)+1;
end
fprintf ( 1, '\rSimulation complete! \n');
toc
% Electric field snapshots.
for i=1:(MaxTime/SnapshotInterval)-1
    
    figure (6)
    mesh ( EzSnapshots (:, :, i) );
    view (4, 4)
    zlim ( [-1 1] )
    caxis([-0.1 0.6])
    xlabel ('y-axis')
    ylabel ('x-axis')
    %colorbar
    
    figure (7)
    mesh ( EzSnapshots (:, :, i) );
    view (0, 90)
    zlim ( [-10 10] )
    caxis([-0.1 0.6])
    xlabel ('y-axis')
    ylabel ('x-axis')
    %colorbar
    
end

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
axis 'auto x'
grid on
subplot(212)
plot(f(1:fspan), 1-TAU(1:fspan), 'LineWidth', 2.0, 'Color', 'b')
set(gca, 'FontSize', 10, 'FontWeight', 'b')
title('Reflection Coefficient', 'FontSize', 12, 'FontWeight', 'b')
xlabel('Frenuency (Hz)', 'FontSize', 11, 'FontWeight', 'b')
ylabel('1-|EZT(f)/EZI(f)|', 'FontSize', 11, 'FontWeight', 'b')
axis([-1 1 -2 2])
axis 'auto x'
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
