clc
clear all

% Simulation parameters.
SIZE = 64; % No. of spatial steps
MaxTime = 64; % No. of time steps
PulseWidth = 20; % Controls width of Gaussian Pulse
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

% Initialization.
Ex = zeros(SIZE, 3); % x-component of E-field
Dx = zeros(SIZE, 3); % x-component of D
Hy = zeros(SIZE, 3); % y-component of H-field
By = zeros(SIZE, 3); % y-component of B

ExSnapshots = zeros(SIZE, MaxTime); % Data for plotting.
frame = 1;

n1 = 1;
n2 = 2;
% Outer loop for time-stepping.
tic
for q = 2:MaxTime
    
    % Storing past fields.
    Ex(:,3) = Ex(:,n1);
    Dx(:,3) = Dx(:,n1);
    Hy(:,3) = Hy(:,n1);
    By(:,3) = By(:,n1);
    
    % Calculation of Hy using update difference equation for Hy. This is time step q.
    Hy(1:SIZE-1,n2) = Hy(1:SIZE-1,n1) + ( ( Ex(1:SIZE-1,n1) - Ex(2:SIZE,n1) ) * dt/(u0*dz) );
    % ABC for H at SIZE.
    Hy(SIZE,n2) = Hy(SIZE-1,n1) + (Sc-1)/(Sc+1)*(Hy(SIZE-1,n2) - Hy(SIZE,n1) );
    
    % Calculation of Ex using updated difference equation for Ex. This is time step q+1/2.
    Ex(2:SIZE,n2) = Ex(2:SIZE, n1) + ( dt/(e0*dz)*(Hy(1:SIZE-1, n2) - Hy(2:SIZE, n2)) );
    % ABC for E at 1.
    Ex(1,n2) = Ex(2,n1) + (Sc-1)/(Sc+1)*(Ex(2,n2) - Ex(2,n1));
    
    % Activating a plane-wave source.
    Ex(source,n2) = Ex(source,n2) + exp( -1*((q-td)/(PulseWidth/4))^2 ) * Sc;

    ExSnapshots(:,frame) = Ex(:,n2);
    frame=frame+1;
    
    temp = n1;
    n1 = n2;
    n2 = temp;
end
toc
% Simulation animation.
for i=1:frame-1
    figure (2)
    plot ( ExSnapshots(:,i) )
    axis([0 SIZE -0.6 0.6])
    xlabel('Spatial step (k)')
    ylabel('Electric field (Ex)')
end