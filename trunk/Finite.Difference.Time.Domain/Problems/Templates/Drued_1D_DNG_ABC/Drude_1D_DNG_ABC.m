clc
clear all

% Simulation parameters.
SIZE = 128; % No. of spatial steps
SlabLeft = round(SIZE/3); % Location of left end of Slab.
SlabRight = round(2*SIZE/3); % Location of right end of Slab
MaxTime = SIZE; % No. of time steps
PulseWidth = round(SIZE/8); % Controls width of Gaussian Pulse
td = PulseWidth; % Temporal delay in pulse.
source = 10; % Location of source

% Constants.
c = 3e8;
pi = 3.141592654;
e0 = (1e-9)/(36*pi);
u0 = (1e-7)*4*pi;

dt = 1e-12;
dz = 3e-4;
Sc = c * dt/dz

l = PulseWidth*dz;
f = l/c;
w = 2*pi*f;

% Initialization.
Ex = zeros(SIZE, 3); % x-component of E-field
Dx = zeros(SIZE, 3); % x-component of D
Hy = zeros(SIZE, 3); % y-component of H-field
By = zeros(SIZE, 3); % y-component of B

wpsq = zeros(SIZE,1);
wpsq(SlabLeft:SlabRight) = 2*w^2;
wpmsq = zeros(SIZE,1);
wpmsq(SlabLeft:SlabRight) = 2*w^2;

a = 4./(e0*(4+wpsq*dt^2));
b = -1*e0*a;
c = -1*(wpsq*dt^2)./(4+(wpsq*dt^2));
am = 4./(u0*(4+wpmsq*dt^2));
bm = -1*u0*am;
cm = -1*(wpmsq*dt^2)./(4+(wpmsq*dt^2));

ExSnapshots = zeros(SIZE, MaxTime); % Data for plotting.
frame = 1;

n1 = 1;
n2 = 2;
% Outer loop for time-stepping.
tic
for q = 2:MaxTime
    
    % Storing past fields.
    Ex(:,3) = Ex(:,n2);
    Dx(:,3) = Dx(:,n2);
    Hy(:,3) = Hy(:,n2);
    By(:,3) = By(:,n2);
    
    % Calculation of Hy using update difference equation for Hy. This is time step q.
    By(1:SIZE-1,n2) = By(1:SIZE-1,n1) + ( ( Ex(1:SIZE-1,n1) - Ex(2:SIZE,n1) ) * dt/(dz) );
    Hy(:,n2) = am.*(By(:,n2)-2*By(:,n1)+By(:,3))+bm.*(-2*Hy(:,n1)+Hy(:,3))+cm.*(2*Hy(:,n1)+Hy(:,3));
    
    % ABC for H at SIZE.
    Hy(SIZE,n2) = Hy(SIZE-1,n1) + (Sc-1)/(Sc+1)*(Hy(SIZE-1,n2) - Hy(SIZE,n1) );
    By(SIZE,n2) = u0*Hy(SIZE,n2);

    % Calculation of Ex using updated difference equation for Ex. This is time step q+1/2.
    Dx(2:SIZE,n2) = Dx(2:SIZE, n1) + ( dt/(dz)*(Hy(1:SIZE-1, n2) - Hy(2:SIZE, n2)) );
    Ex(:,n2) = a.*(Dx(:,n2)-2*Dx(:,n1)+Dx(:,3))+b.*(-2*Ex(:,n1)+Ex(:,3))+c.*(2*Ex(:,n1)+Ex(:,3));
    
    % ABC for E at 1.
    Ex(1,n2) = Ex(2,n1) + (Sc-1)/(Sc+1)*(Ex(2,n2) - Ex(2,n1));
    Dx(1,n2) = e0*Ex(1,n2);
    
    % Activating a plane-wave source.
    Ex(source,n2) = Ex(source,n2) + exp( -1*((q-td)/(PulseWidth/4))^2 ) * Sc;
    Dx(source,n2) = e0*Ex(source,n2);

    ExSnapshots(:,frame) = Ex(:,n2);
    frame=frame+1;
    
    temp = n1;
    n1 = n2;
    n2 = temp;
end
toc
% Simulation animation.
for i=1:frame-1
    figure (1)
    plot ( ExSnapshots(:,i) )
    axis([0 SIZE -0.6 0.6])
    xlabel('Spatial step (k)')
    ylabel('Electric field (Ex)')
end