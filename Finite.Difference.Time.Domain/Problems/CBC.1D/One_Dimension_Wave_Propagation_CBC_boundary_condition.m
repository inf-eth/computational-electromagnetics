% Simulation parameters.
SIZE = 1024; % No. of spatial steps
MaxTime = 2024; % No. of time steps
imp0 = 377.0; % Impedence of free space

% Constants.
c = 3e8;
pi = 3.141592654;
e0 = (1e-9)/(36*pi);
u0 = (1e-7)*4*pi;

dt = 0.8e-12;
dx = 3e-4;
Sc = c * dt/dx

% Initialization.
Ez = zeros ( SIZE, MaxTime ); % z-component of E-field
Hy = zeros ( SIZE, MaxTime ); % y-component of H-field
PLOT1(1) = 0; % Data for plotting.

% Outer loop for time-stepping.
tic
for q = 2:MaxTime
    % Calculation of Hy using update difference equation for Hy. This is time step q.
    Hy(1:SIZE-1,q) = Hy(1:SIZE-1,q-1) + ( ( Ez(2:SIZE,q-1) - Ez(1:SIZE-1,q-1) ) * Sc/imp0 );
    % CBC for H at SIZE.
    Hy(SIZE,q) = Hy(SIZE,q-1) + ( ( Ez(1,q-1) - Ez(SIZE,q-1) ) * Sc/imp0 );
    
    % Calculation of Ez using updated difference equation for Ez. This is time step q+1/2.
    Ez(2:SIZE,q) = Ez(2:SIZE, q-1) + ( Sc*imp0*(Hy(2:SIZE, q) - Hy(1:SIZE-1, q)) );
    % CBC for E at 1.
    Ez(1,q) = Ez(1,q-1) + ( Sc*imp0*(Hy(1, q) - Hy(SIZE, q)) );
    
    % Activating a plane-wave source.
    Ez(100,q) = Ez(100,q) + exp( -1*((q-256)^2)/1000 );
end
toc
% Simulation animation.
for i=1:MaxTime
    figure (2)
    plot ( Ez(:,i) )
end

% Plotting Electric field at differenct times. Can be used to view the
% electric field at discrete time steps.
% figure
% plot ( Ez(:,10) )
% figure
% plot ( Ez(:,30) )
% figure
% plot ( Ez(:,50) )
