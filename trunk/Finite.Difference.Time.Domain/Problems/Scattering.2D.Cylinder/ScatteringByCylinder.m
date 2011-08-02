% FDTD example 3.7 from Numerical techniques in Electromagnetics by Sadiku.
% Main m file for the FDTD simulation.
clc
clear all
% Simulation related parameters.
IHx = 50;
JHx = 49;
IHy = 51;
JHy = 50;
IEz = 50;
JEz = 50;
NMax = 2; NNMax = 500; % Maximum time.
NHW = 40; % One half wave cycle.
Med = 2; % No of different media.
Js = 2; % J-position of the plane wave front.
% Different Constants.
delta = 3e-3;
Cl = 3e8;
f = 2.5e9;
pi = 3.141592654;
e0 = (1e-9) / (36*pi);
u0 = (1e-7) * 4 * pi;
DT = delta / ( 2 * Cl );
TwoPIFDeltaT = 2 * pi * f * DT;
% Data arrays.
CHx = zeros ( IHx, JHx ); % Conductance
CHy = zeros ( IHy, JHy ); % Conductance
REz = zeros ( IEz, JEz ); % Impedance
RaEz = zeros ( IEz, JEz ); % Scaling parameter dependent on material conductivity.
Hx = zeros ( IHx, JHx, NNMax );
Hy = zeros ( IHy, JHy, NNMax );
Ez = zeros ( IEz, JEz, NNMax );
AbsEz = zeros ( IEz, JEz ); % Absolute value of Ez for steady state solution

% ############ Initialization #############
fprintf ( 1, 'Initializing' );
fprintf ( 1, '.' );
% Initializing the eta and 1/eta arrays.
for i=1:IHx
    for j=1:JHx
        CHx ( i, j ) = DT / ( u0 * ur ( i, j-0.5 ) * delta );
    end
end

fprintf ( 1, '.' );
for i=1:IHy
    for j=1:JHy
        CHy ( i, j ) = DT / ( u0 * ur ( i-0.5, j-1 ) * delta );
    end
end
fprintf ( 1, '.' );
for i=1:IEz
    for j=1:JEz
        REz ( i, j ) = DT / ( e0 * er ( i, j-0.5 ) * delta );
        RaEz ( i, j ) = ( 1 - ( s(i, j-0.5) * DT )/( e0 * er( i, j-0.5 ) ) );
    end
end
fprintf ( 1, 'done.\n' );
% ############ Initialization Complete ##############
% ########### 2. Now running the Simulation #############
fprintf ( 1, 'Simulation started... \n', 0 );
for n=0:498
    fprintf ( 1, '%g %% \n', (n*100)/500 );
    % *** Calculation of Magnetic Field Components ***
    % * Calculation of Hx.
    Hx ( :, :, n+2 ) = Hx ( :, :, n+1 ) + ( CHx .* ( Ez ( :, 1:JHx, n+1 ) - Ez ( :, 2:JHx+1, n+1 ) ));
    % * Calculation of Hy.
    Hy ( 2:IHy-1, :, n+2 ) = Hy ( 2:IHy-1, :, n+1 ) + ( CHy( 2:IHy-1, : ) .* ( Ez ( 2:IHy-1, :, n+1 ) - Ez ( 1:IHy-2, :,n+1 ) ));
    % Boundary conditions on Hy. Soft grid truncation.
    Hy ( 1, 2:JHy-1, n+2 ) = (1/3) * ( Hy ( 2, 1:JHy-2, n+1 ) + Hy ( 2, 2:JHy-1, n+1 ) + Hy ( 2, 3:JHy, n+1 ) );
    Hy ( IHy, 2:JHy-1, n+2 ) = (1/3) * ( Hy ( IHy-1, 1:JHy-2, n+1 ) + Hy ( IHy-1, 2:JHy-1, n+1 ) + Hy ( IHy-1, 3:JHy, n+1 ) );
    Hy ( 1, 1, n+2 ) = (1/2) * ( Hy ( 2, 1, n+1 ) + Hy ( 2, 2, n+1 ) );
    Hy ( 1, JHy, n+2 ) = (1/2) * ( Hy ( 2, JHy, n+1 ) + Hy ( 2, JHy-1, n+1 ) );
    Hy ( IHy, 1, n+2 ) = (1/2) * ( Hy ( IHy-1, 1, n+1 ) + Hy( IHy-1, 2, n+1 ) );
    Hy ( IHy, JHy, n+2 ) = (1/2) * ( Hy ( IHy-1, JHy, n+1 ) + Hy ( IHy-1, JHy-1, n+1 ) );
    % ************************************************
    % *** Calculation of Electric Field Components ***
    % * Calculation of Ez.
    Ez ( :, 2:JEz-1, n+2 ) = ( RaEz ( :, 2:JEz-1 ) .* Ez ( :, 2:JEz-1, n+1 ) ) + ( REz ( :, 2:JEz-1 ) .* ( Hy ( 2:JEz+1, 2:JEz-1, n+2 ) - Hy ( 1:JEz, 2:JEz-1, n+2 ) + Hx ( :, 1:JEz-2,n+2 ) - Hx ( :, 2:JEz-1, n+2 ) ));
    % Boundary conditions on Ez. Soft grid truncation.
    Ez ( 2:IEz-1, 1, n+2 ) = (1/3) * ( Ez ( 1:IEz-2, 2, n+1 ) + Ez ( 2:IEz-1, 2, n+1 ) + Ez ( 3:IEz, 2, n+1 ) );
    Ez ( 2:IEz-1, JEz, n+2 ) = (1/3) * ( Ez ( 1:IEz-2, JEz-1, n+1 ) + Ez ( 2:IEz-1, JEz-1, n+1 ) + Ez ( 3:IEz, JEz-1, n+1 ) );
    Ez ( 1, 1, n+2 ) = (1/2) * ( Ez ( 1, 2, n+1 ) + Ez ( 2, 2, n+1 ) );
    Ez ( IEz, 1, n+2 ) = (1/2) * ( Ez ( IEz, 2, n+1 ) + Ez ( IEz-1, 2, n+1 ) );
    Ez ( 1, JEz, n+2 ) = (1/2) * ( Ez ( 1, JEz-1, n+1 ) + Ez ( 2, JEz-1, n+1 ) );
    Ez ( IEz, JEz, n+2 ) = (1/2) * ( Ez ( IEz, JEz-1, n+1 ) + Ez ( IEz-1, JEz-1, n+1 ) );
    % ************************************************
    % Applying a plane wave source at Js.
    Ez ( :, Js, n+2 ) = 1 * sin ( TwoPIFDeltaT * (n+2) );
    % 4. Retaining absolute values during the last half cycle.
    if ( n > 460 )
        A = abs ( Ez(:,:,n+2) ) >= AbsEz;
        B = abs ( Ez(:,:,n+2) ) < AbsEz;
        C = A .* abs ( Ez(:,:,n+2) );
        D = B .* AbsEz;
        AbsEz = C + D;
    end
end
fprintf ( 1, '100 %% \n' );
fprintf ( 1, 'Simulation complete! \n', 0 );
% Plotting Ez.
figure (1)
subplot ( 211 )
plot ( 5:43, AbsEz(25,5:43) )
title ( 'Ez ( 25, j )' )
subplot ( 212 )
plot ( 5:43, AbsEz(15,5:43) )
title ( 'Ez ( 15, j )' )


for i=1:NNMax
    figure (2)
    mesh ( Ez (:, :, i) );
    view (3, 34)
    zlim ( [-2 2] )
    
    figure (4)
    mesh ( Ez (:, :, i) );
    view (0, 90)
    zlim ( [-2 2] )
    
    %F(i)=getframe;
end