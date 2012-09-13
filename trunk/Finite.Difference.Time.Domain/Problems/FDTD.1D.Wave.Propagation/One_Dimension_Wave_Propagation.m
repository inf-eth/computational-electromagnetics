% Simulation parameters.
SIZE = 200; % No. of spatial steps
MaxTime = 200; % No. of time steps
imp0 = 377.0; % Impedence of free space

% Initialization.
Ez = zeros ( SIZE, MaxTime ); % z-component of E-field
Hy = zeros ( SIZE, MaxTime ); % y-component of H-field
PLOT1(1) = 0; % Data for plotting.

% Outer loop for time-stepping.
for q = 2:MaxTime
    % Calculation of Hy using update difference equation for Hy. This is time step q.
    for m = 1:SIZE-1
        Hy(m,q) = Hy(m,q-1) + ( ( Ez(m+1,q-1) - Ez(m,q-1) ) / imp0 );
    end
    % Calculation of Ez using updated difference equation for Ez. This is time step q+1/2.
    for m = 2:SIZE
        Ez(m,q) = Ez(m,q-1) + ( imp0 * ( Hy(m,q) - Hy(m-1,q) ) );
    end
    % Activating a plane-wave source.
    Ez(1,q) = exp ( - 1 * ( (q-31)^2) / 100 );
end

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
