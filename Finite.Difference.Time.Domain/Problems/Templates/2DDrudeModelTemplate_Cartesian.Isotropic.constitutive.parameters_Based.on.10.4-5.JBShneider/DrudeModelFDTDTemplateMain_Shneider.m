% 2D FDTD template for Drude Model with isotropic constitutive parameters from 10.4-5 JB Shneider.

% * The relative permittivity and permeability returned from functions er() and ur() must be less than 1 or simulation will be unstable.
% ** er = ur = 1 simulates free space. Keep in mind, one or both these values can also be negative, e.g. to simulate DNG etc.

% Based on FDTD example 3.7 from Numerical techniques in Electromagnetics by Sadiku.
% Main m file for the FDTD simulation.

clc
clear all

% ============== Simulation related parameters ================
[Size XCenter YCenter delta ra rb DTp] = Parameters;
IHx = Size;
JHx = Size-1;
IHy = Size+1;
JHy = Size;
IEz = Size;
JEz = Size;
% Time indices for field calculation.
n0 = 1;
n1 = 2;
NNMax = 500;                   % Maximum time.
TimeResolutionFactor = 1;      % E field snapshots will be saved every x frames where x is resolution factor.
ResolutionFactor = 2;          % Resolution of plotted field is divided by this factor.
Js = 2;                         % J-position of the plane wave front.

% Different Constants.
Cl = 3e8;
f = 2.0e9;
pi = 3.141592654;
e0 = (1e-9) / (36*pi);
u0 = (1e-7) * 4 * pi;
einf = 1;
uinf = 1;
DT = DTp;
TwoPIFDeltaT = 2 * pi * f * DT;
NHW = 1/(2 * f * DT); % One half wave cycle.

% ====================== Data arrays =========================
wpsquaredEz = zeros( IEz, JEz );
wpmsquaredHx = zeros (IHx, JHx );
wpmsquaredHy = zeros ( IHy, JHy );

smaskEz = zeros ( IEz, JEz );

Hx = zeros ( IHx, JHx, 2 );
Hy = zeros ( IHy, JHy, 2 );

Ez = zeros ( IEz, JEz, 2 );

Jpz = zeros ( IEz, JEz, 2 );
Jmx = zeros ( IHx, JHx, 2 );
Jmy = zeros ( IHy, JHy, 2 );

EzSnapshots = zeros ( IEz/ResolutionFactor, JEz/ResolutionFactor, NNMax/TimeResolutionFactor ); % E field Snapshot storage.
% =============================================================

% ############ Initialization #############
fprintf ( 1, 'Initializing...' );
fprintf ( 1, '\nInitializing parametric arrays...' );
% Initializing er array.
for i=1:IHy     % IHy is size+1 or maximum I size.
    fprintf ( 1, '%g %% \n', ((i-1)*100)/IHy );
    for j=1:JHy
        
        % Ez related parameters.
        if i <= IEz
            smaskEz( i, j ) = s ( i, j-0.5 );
            wpsquaredEz(i, j) = wpsquared(i, j-0.5, 2*pi*f);
         end
        
        % Hx related parameters.
        if i <= IHx && j <= JHx
            wpmsquaredHx (i, j) = wpmsquared (i, j-0.5, 2*pi*f);
        end
        
        % Hy related parameters.
        wpmsquaredHy (i, j) = wpmsquared (i-0.5, j-1, 2*pi*f);
    end
end

figure (1)
mesh ( wpsquaredEz )
title ( 'wp squared' )
view (4, 4)
figure (2)
mesh ( wpmsquaredHx )
title ( 'wpm squared' )
view (4, 4)

fprintf ( 1, 'done.\n' );
% ############ Initialization Complete ##############
% ########### 2. Now running the Simulation #############
fprintf ( 1, 'Simulation started... \n' );
for n=0:NNMax-2
    fprintf ( 1, '%g %% \n', (n*100)/NNMax );

    % Hx and Hy.
    Hx ( :, :, n1 ) = Hx ( :, :, n0 ) - ( (DT/(delta*u0*uinf)) * ((Ez(:, 2:JHx+1, n0) - Ez(:, 1:JHx, n0)) + delta*Jmx(:, :, n0)) );
    Hy ( 2:IHy-1, :, n1 ) = Hy ( 2:IHy-1, :, n0 ) - ( (DT/(delta*u0*uinf)) * ((Ez(1:IHy-2, :, n0) - Ez(2:IHy-1, :, n0)) + delta*Jmy(2:IHy-1, :, n0)) );

    % Boundary conditions on Hy. Soft grid truncation.
    Hy ( 1, 2:JHy-1, n1 ) = (1/3) * ( Hy ( 2, 1:JHy-2, n0 ) + Hy ( 2, 2:JHy-1, n0 ) + Hy ( 2, 3:JHy, n0 ) );
    Hy ( IHy, 2:JHy-1, n1 ) = (1/3) * ( Hy ( IHy-1, 1:JHy-2, n0 ) + Hy ( IHy-1, 2:JHy-1, n0 ) + Hy ( IHy-1, 3:JHy, n0 ) );
    Hy ( 1, 1, n1 ) = (1/2) * ( Hy ( 2, 1, n0 ) + Hy ( 2, 2, n0 ) );
    Hy ( 1, JHy, n1 ) = (1/2) * ( Hy ( 2, JHy, n0 ) + Hy ( 2, JHy-1, n0 ) );
    Hy ( IHy, 1, n1 ) = (1/2) * ( Hy ( IHy-1, 1, n0 ) + Hy( IHy-1, 2, n0 ) );
    Hy ( IHy, JHy, n1 ) = (1/2) * ( Hy ( IHy-1, JHy, n0 ) + Hy ( IHy-1, JHy-1, n0 ) );
    
    % Magnetic polarization currents.
    Jmx ( :, :, n1 ) = Jmx ( :, :, n0 ) + (u0*DT)/2 .* wpmsquaredHx .* ( Hx(:, :, n1) + Hx (:, :, n0));
    Jmy ( 2:IHy-1, :, n1 ) = Jmy ( 2:IHy-1, :, n0 ) + (u0*DT)/2 .* wpmsquaredHy(2:IHy-1,:) .* ( Hy(2:IHy-1, :, n1) + Hy (2:IHy-1, :, n0));
    
    % Electric field.
    Ez ( :, 2:JEz-1, n1 ) = Ez ( :, 2:JEz-1, n0 ) + ( (DT/(delta*e0*einf)) * ( Hy ( 2:JEz+1, 2:JEz-1, n1 ) - Hy ( 1:JEz, 2:JEz-1, n1 ) + Hx ( :, 1:JEz-2,n1 ) - Hx ( :, 2:JEz-1, n1 ) - (delta/2)*(Jpz(:, 2:JEz-1, n1)+Jpz(:, 2:JEz-1, n0)) ));
    
    % Boundary conditions on Dz. Soft grid truncation.
    Ez ( 2:IEz-1, 1, n1 ) = (1/3) * ( Ez ( 1:IEz-2, 2, n0 ) + Ez ( 2:IEz-1, 2, n0 ) + Ez ( 3:IEz, 2, n0 ) );
    Ez ( 2:IEz-1, JEz, n1 ) = (1/3) * ( Ez ( 1:IEz-2, JEz-1, n0 ) + Ez ( 2:IEz-1, JEz-1, n0 ) + Ez ( 3:IEz, JEz-1, n0 ) );
    Ez ( 1, 1, n1 ) = (1/2) * ( Ez ( 1, 2, n0 ) + Ez ( 2, 2, n0 ) );
    Ez ( IEz, 1, n1 ) = (1/2) * ( Ez ( IEz, 2, n0 ) + Ez ( IEz-1, 2, n0 ) );
    Ez ( 1, JEz, n1 ) = (1/2) * ( Ez ( 1, JEz-1, n0 ) + Ez ( 2, JEz-1, n0 ) );
    Ez ( IEz, JEz, n1 ) = (1/2) * ( Ez ( IEz, JEz-1, n0 ) + Ez ( IEz-1, JEz-1, n0 ) );
    
    % Polarization current.
    Jpz ( :, :, n1 ) = Jpz ( :, :, n0 ) + (e0*DT)/2 .* wpsquaredEz .* ( Ez(:, :, n1) + Ez (:, :, n0));
    
    % Comment the if statement for a continous plane wave source. Otherwise only a single pulse will be incident.
%     if ( n < NHW )
        Ez ( :, Js, n1 ) = Ez ( :, Js, n1 ) + 1 * sin ( TwoPIFDeltaT * n );
%     end

    % This mask can be used to simulate PEC. The PEC locations can be defined in s.m file.
    %Ez ( :, :, n1 ) = smaskEz (:, :) .* Ez ( :, :, n1 );

    if ( mod(n, TimeResolutionFactor) == 0)
        EzSnapshots ( :, :, n/TimeResolutionFactor + 1 ) = Ez ( 1+(0:ResolutionFactor:(IEz-1)), 1+(0:ResolutionFactor:(JEz-1)), n1);
    end
    
    temp = n0;
    n0 = n1;
    n1 = temp;    
end
fprintf ( 1, '100 %% \n' );
fprintf ( 1, 'Calculations completed! \n' );

% Electric field snapshots.
for i=1:NNMax/TimeResolutionFactor-2
    
    figure (6)
    mesh ( EzSnapshots (:, :, i) );
    view (4, 4)
    zlim ( [-2 2] )
    
    figure (7)
    surf ( EzSnapshots (:, :, i) );
    view (0, 90)
    zlim ( [-10 10] )
    
end
fprintf ( 1, 'Simulation completed! \n' );