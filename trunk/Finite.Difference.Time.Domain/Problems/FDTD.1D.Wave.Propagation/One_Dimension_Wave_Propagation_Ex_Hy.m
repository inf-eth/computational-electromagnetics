% Simulation parameters.
SIZE = 512; % No. of spatial steps
MaxTime = 512; % No. of time steps
PulseWidth = 512; % Width of Gaussian Pulse
imp0 = 377.0; % Impedence of free space

% Constants.
c = 3e8;
pi = 3.141592654;
e0 = (1e-9)/(36*pi);
u0 = (1e-7)*4*pi;

dt = 1e-12;
dz = 3e-4;
Sc = c * dt/dz

% Initialization.
Ex = zeros ( SIZE, MaxTime ); % z-component of E-field
Hy = zeros ( SIZE, MaxTime ); % y-component of H-field
PLOT1(1) = 0; % Data for plotting.

% Outer loop for time-stepping.
tic
for q = 2:MaxTime
    % Calculation of Hy using update difference equation for Hy. This is time step q.
    Hy(1:SIZE-1,q) = Hy(1:SIZE-1,q-1) + ( ( Ex(1:SIZE-1,q-1) - Ex(2:SIZE,q-1) ) * dt/(u0*dz) );
    % CBC for H at SIZE.
    Hy(SIZE,q) = Hy(SIZE,q-1) + ( ( Ex(SIZE,q-1) - Ex(1,q-1) ) * dt/(u0*dz) );
    
    % Calculation of Ez using updated difference equation for Ez. This is time step q+1/2.
    Ex(2:SIZE,q) = Ex(2:SIZE, q-1) + ( dt/(e0*dz)*(Hy(1:SIZE-1, q) - Hy(2:SIZE, q)) );
    % CBC for E at 1.
    Ex(1,q) = Ex(1,q-1) + ( dt/(e0*dz)*(Hy(SIZE, q) - Hy(1, q)) );
    
    % Activating a plane-wave source.
    Ex(100,q) = Ex(100,q) + exp( -1*(q-PulseWidth)^2/PulseWidth ) * Sc;
end
toc
% Simulation animation.
for i=1:MaxTime
    figure (2)
    plot ( Ex(:,i) )
end