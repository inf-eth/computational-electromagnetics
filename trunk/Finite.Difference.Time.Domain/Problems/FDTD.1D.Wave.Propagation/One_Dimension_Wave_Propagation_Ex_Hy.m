% Simulation parameters.
SIZE = 128; % No. of spatial steps
MaxTime = 128; % No. of time steps
PulseWidth = 32; % Width of Gaussian Pulse
imp0 = 377.0; % Impedence of free space
source = 2; % Location of source

% Constants.
c = 3e8;
pi = 3.141592654;
e0 = (1e-9)/(36*pi);
u0 = (1e-7)*4*pi;

dt = 1.0e-12;
dz = 3e-4;
Sc = c * dt/dz

% Initialization.
Ex = zeros ( SIZE, MaxTime ); % x-component of E-field
Hy = zeros ( SIZE, MaxTime ); % y-component of H-field
PLOT1(1) = 0; % Data for plotting.

% Outer loop for time-stepping.
tic
for q = 2:MaxTime
    % Calculation of Hy using update difference equation for Hy. This is time step q.
    Hy(1:SIZE-1,q) = Hy(1:SIZE-1,q-1) + ( ( Ex(1:SIZE-1,q-1) - Ex(2:SIZE,q-1) ) * dt/(u0*dz) );
    % ABC for H at SIZE.
    Hy(SIZE,q) = Hy(SIZE-1,q-1) + (Sc-1)/(Sc+1)*(Hy(SIZE-1,q) - Hy(SIZE,q-1) );
    
    % Calculation of Ex using updated difference equation for Ex. This is time step q+1/2.
    Ex(2:SIZE,q) = Ex(2:SIZE, q-1) + ( dt/(e0*dz)*(Hy(1:SIZE-1, q) - Hy(2:SIZE, q)) );
    % ABC for E at 1.
    Ex(1,q) = Ex(2,q-1) + (Sc-1)/(Sc+1)*(Ex(2,q) - Ex(2,q-1));
    
    % Activating a plane-wave source.
    Ex(source,q) = Ex(source,q) + exp( -1*((q-PulseWidth)^2)/PulseWidth ) * Sc;
end
toc
% Simulation animation.
for i=1:MaxTime
    figure (2)
    plot ( Ex(:,i) )
    axis([0 SIZE -0.2 1])
    xlabel('Spatial step (k)')
    ylabel('Electric field (Ex)')
end