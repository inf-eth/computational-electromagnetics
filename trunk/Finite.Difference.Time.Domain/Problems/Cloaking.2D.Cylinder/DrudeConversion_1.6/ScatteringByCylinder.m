% Based on FDTD example 3.7 from Numerical techniques in Electromagnetics by Sadiku.
% Main m file for the FDTD simulation.
% Now this doesn't look at all like the original problem. Only thing
% similar are some of the constants and field orientation.

clc
clear all

% ============== Simulation related parameters ================
[Size XCenter YCenter delta ra rb] = Parameters;
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
ResolutionFactor = 1;          % Resolution of plotted field is divided by this factor.
Js = 2;                         % J-position of the plane wave front.
% Different Constants.
Cl = 3e8;
f = 2.0e9;
pi = 3.141592654;
e0 = (1e-9) / (36*pi);
u0 = (1e-7) * 4 * pi;
DT = delta / ( sqrt(2) * Cl );
TwoPIFDeltaT = 2 * pi * f * DT;
NHW = 1/(2 * f * DT); % One half wave cycle.
A = rb/(rb-ra);

% ====================== Data arrays =========================
uxxHx = zeros ( IHx, JHx );  % uxx for Hx
uxyHy = zeros ( IHy, JHy );  % uxy for Hy
uyxHx = zeros ( IHx, JHx );  % uyx for Hx
uyyHy = zeros ( IHy, JHy );  % uyy for Hy
ezzEz = zeros ( IEz, JEz );  % ezz for Ez

% ------------ Field-specific parameters ------------
ax = zeros ( IHx, JHx );
bx = zeros ( IHx, JHx );
cx = zeros ( IHx, JHx );
dx = zeros ( IHx, JHx );
ex = zeros ( IHx, JHx );
fx = zeros ( IHx, JHx );
gx = zeros ( IHx, JHx );
hx = zeros ( IHx, JHx );
lx = zeros ( IHx, JHx );

ay = zeros ( IHy, JHy );
by = zeros ( IHy, JHy );
cy = zeros ( IHy, JHy );
dy = zeros ( IHy, JHy );
ey = zeros ( IHy, JHy );
fy = zeros ( IHy, JHy );
gy = zeros ( IHy, JHy );
hy = zeros ( IHy, JHy );
ly = zeros ( IHy, JHy );

wpsquaredEz = zeros( IEz, JEz );

% **** Temp ****
erEz = zeros ( IEz, JEz );  % ezz for Ez
% *************

%RaEz = zeros ( IEz, JEz ); % Scaling parameter dependent on material conductivity.
cmaskEz = zeros ( IEz, JEz );       % cylinder mask.
cmaskHx = zeros ( IHx, JHx );
cmaskHy = zeros ( IHy, JHy );

smaskEz = zeros ( IEz, JEz );
smaskHx = zeros ( IHx, JHx );
smaskHy = zeros ( IHy, JHy );

Bx = zeros ( IHx, JHx, 3 );
By = zeros ( IHy, JHy, 3 );

Hx = zeros ( IHx, JHx, 3 );
Hy = zeros ( IHy, JHy, 3 );

Dz = zeros ( IEz, JEz, 3 );
Ez = zeros ( IEz, JEz, 3 );
EzSnapshots = zeros ( IEz/ResolutionFactor, JEz/ResolutionFactor, NNMax/TimeResolutionFactor ); % E field Snapshot storage.
% =============================================================

% ############ Initialization #############
fprintf ( 1, 'Initializing...' );
fprintf ( 1, '\nInitializing er array...' );
% Initializing er array.
for i=1:IEz
    for j=1:JEz
        ezzEz ( i, j ) = ezz ( i, j-0.5 );
%         erEz ( i, j ) = er( i, j-0.5 );
        erEz ( i, j ) = ezz( i, j-0.5 );
        smaskEz( i, j ) = s ( i, j-0.5 );
        wpsquaredEz(i, j) = wpsquared(i, j-0.5, 2*pi*f);
        cmaskEz( i, j ) = iscylinder (i, j-0.5);
        %er ( i, j-0.5 )
    end
end
fprintf ( 1, '\nInitializing ur (Hx) array...' );
% Initializing the ur arrays.
for i=1:IHx
    for j=1:JHx
        invurHx = inv ( ur ( i, j-0.5 ));
        uxxHx (i, j) = invurHx(1, 1);
        uyxHx (i, j) = invurHx(2, 1);
        cmaskHx (i, j) = iscylinder (i, j-0.5);        
%         smaskHx( i, j ) = s ( i, j-0.5 );
        ax (i, j) = 1;
        bx (i, j) = 1;
        cx (i, j) = 1;
        dx (i, j) = 1;
        ex (i, j) = 1;
        fx (i, j) = 1;
        gx (i, j) = 1;
        hx (i, j) = 1;
        lx (i, j) = 1;
    end
end
fprintf ( 1, '\nInitializing ur (Hy) array...' );
for i=1:IHy
    for j=1:JHy
        invurHy = inv ( ur ( i-0.5, j-1 ));
        uxyHy (i, j) = invurHy(1, 2);
        uyyHy (i, j) = invurHy(2, 2);
        cmaskHy(i, j) = iscylinder (i-0.5, j-1);
        %         smaskHy( i, j ) = s ( i-0.5, j-1 );
        ay (i, j) = 1;
        by (i, j) = 1;
        cy (i, j) = 1;
        dy (i, j) = 1;
        ey (i, j) = 1;
        fy (i, j) = 1;
        gy (i, j) = 1;
        hy (i, j) = 1;
        ly (i, j) = 1;
    end
end

figure (1)
mesh ( wpsquaredEz )
title ( 'wpsquaredEz' )
view (4, 4)
figure (2)
mesh ( erEz )
title ( 'ezz' )
view (4, 4)
% figure (5)
% mesh ( uxxHx )
% title ( 'uxxHx' )
% view (4, 4)
% figure (6)
% mesh ( uxyHy )
% title ( 'uxyHy' )
% view (4, 4)
% figure (7)
% mesh ( uyyHy )
% title ( 'uyyHy' )
% view (4, 4)

fprintf ( 1, 'done.\n' );
% ############ Initialization Complete ##############
% ########### 2. Now running the Simulation #############
fprintf ( 1, 'Simulation started... \n' );
for n=0:NNMax-2
    fprintf ( 1, '%g %% \n', (n*100)/NNMax );
    
    % Copying n1 fields as past fields. Will be used in second order time derivatives.
    Bx ( :, :, 3 ) = Bx ( :, :, n1 );
    By ( :, :, 3 ) = By ( :, :, n1 );
    Hx ( :, :, 3 ) = Hx ( :, :, n1 );
    Hy ( :, :, 3 ) = Hy ( :, :, n1 );
    Dz ( :, :, 3 ) = Dz ( :, :, n1 );
    Ez ( :, :, 3 ) = Ez ( :, :, n1 );
    
    % *** Calculation of Magnetic Field Components ***
    % * Calculation of Bx.
    Bx ( :, :, n1 ) = Bx ( :, :, n0 ) + ( (DT/delta) * ( Ez ( :, 1:JHx, n0 ) - Ez ( :, 2:JHx+1, n0 ) ));
    %Bx ( :, :, n1 ) = smaskHx ( :, : ) .* Bx ( :, :, n1 );
    
    % * Calculation of By.
    By ( 2:IHy-1, :, n1 ) = By ( 2:IHy-1, :, n0 ) + ( (DT/delta) * ( Ez ( 2:IHy-1, :, n0 ) - Ez ( 1:IHy-2, :,n0 ) ));
    %By ( :, :, n1 ) = smaskHy ( :, : ) .* By ( :, :, n1 );
   
    % Boundary conditions on By. Soft grid truncation.
    By ( 1, 2:JHy-1, n1 ) = (1/3) * ( By ( 2, 1:JHy-2, n0 ) + By ( 2, 2:JHy-1, n0 ) + By ( 2, 3:JHy, n0 ) );
    By ( IHy, 2:JHy-1, n1 ) = (1/3) * ( By ( IHy-1, 1:JHy-2, n0 ) + By ( IHy-1, 2:JHy-1, n0 ) + By ( IHy-1, 3:JHy, n0 ) );
    By ( 1, 1, n1 ) = (1/2) * ( By ( 2, 1, n0 ) + By ( 2, 2, n0 ) );
    By ( 1, JHy, n1 ) = (1/2) * ( By ( 2, JHy, n0 ) + By ( 2, JHy-1, n0 ) );
    By ( IHy, 1, n1 ) = (1/2) * ( By ( IHy-1, 1, n0 ) + By( IHy-1, 2, n0 ) );
    By ( IHy, JHy, n1 ) = (1/2) * ( By ( IHy-1, JHy, n0 ) + By ( IHy-1, JHy-1, n0 ) );
%     
%     By ( 1, :, n1 ) = By ( 2, :, n0 );
%     By ( IHy, :, n1 ) = By ( IHy-1, :, n0 );
    
    Hx ( :, :, n1 ) = (1/u0) * (uxxHx .* Bx ( :, :, n1 ) + uxyHy (1:IHy-1, 1:JHy-1 ) .* By (1:IHy-1, 1:JHy-1, n1 ));
%     ****** Major ****** Hx ( :, :, n1 ) = (1/u0) * (uxxHx .* Bx ( :, :, n1 ) + uxyHy (1:IHy-1, 1:JHy-1 ) .* By (1:IHy-1, 1:JHy-1, n1 ));
    %Hx ( :, :, n1 ) = smaskHx ( :, : ) .* Hx ( :, :, n1 );
    
    
    Hy (1:IHy-1, 1:JHy-1, n1) = (1/u0) * (uyxHx.*Bx ( :, :, n1 ) + uyyHy (1:IHy-1, 1:JHy-1) .* By (1:IHy-1, 1:JHy-1, n1 ));
    %Hy ( :, :, n1 ) = (1/u0) * By ( :, :, n1 );
%     ****** Major ******    Hy (1:IHy-1, 1:JHy-1, n1) = (1/u0) * (uyxHx.*Bx ( :, :, n1 ) + uyyHy (1:IHy-1, 1:JHy-1) .* By (1:IHy-1, 1:JHy-1, n1 ));
    
    % ABC for Hy;
%     Hy (1:IHy-1, JHy, n1) = Hy (1:IHy-1, JHy, n0);
    %Hy ( :, :, n1 ) = smaskHy ( :, : ) .* Hy ( :, :, n1 );
    
    Dz ( :, 2:JEz-1, n1 ) = (  Dz ( :, 2:JEz-1, n0 ) ) + ( (DT/delta) * ( Hy ( 2:JEz+1, 2:JEz-1, n1 ) - Hy ( 1:JEz, 2:JEz-1, n1 ) + Hx ( :, 1:JEz-2,n1 ) - Hx ( :, 2:JEz-1, n1 ) ));
    %Dz (:, :, n1) = Dz (:, :, n1) .* smask;
    
    % Boundary conditions on Dz. Soft grid truncation.
    Dz ( 2:IEz-1, 1, n1 ) = (1/3) * ( Dz ( 1:IEz-2, 2, n0 ) + Dz ( 2:IEz-1, 2, n0 ) + Dz ( 3:IEz, 2, n0 ) );
    Dz ( 2:IEz-1, JEz, n1 ) = (1/3) * ( Dz ( 1:IEz-2, JEz-1, n0 ) + Dz ( 2:IEz-1, JEz-1, n0 ) + Dz ( 3:IEz, JEz-1, n0 ) );
    Dz ( 1, 1, n1 ) = (1/2) * ( Dz ( 1, 2, n0 ) + Dz ( 2, 2, n0 ) );
    Dz ( IEz, 1, n1 ) = (1/2) * ( Dz ( IEz, 2, n0 ) + Dz ( IEz-1, 2, n0 ) );
    Dz ( 1, JEz, n1 ) = (1/2) * ( Dz ( 1, JEz-1, n0 ) + Dz ( 2, JEz-1, n0 ) );
    Dz ( IEz, JEz, n1 ) = (1/2) * ( Dz ( IEz, JEz-1, n0 ) + Dz ( IEz-1, JEz-1, n0 ) );
    
%     Dz ( :, 1, n1 ) = Dz ( :, 2, n0 );
%     Dz ( :, IEz, n1 ) = Dz ( :, IEz-1, n0 );
    % ************************************************

    
%     %     ****** Major ******   Ez ( :, :, n1 ) = (1/e0) * Dz ( :, :, n1 ) ./ (ezzEz);
    
    Ez ( :, :, n1 ) =  ( ~cmaskEz .* (1/e0) * (ezzEz) .* Dz ( :, :, n1 ) ) + ( cmaskEz .* ( (1/(e0*DT^2))*Dz ( :, :, n1 ) - (2/(e0*(DT^2)))*Dz ( :, :, n0) + (1/(e0*(DT^2)))*Dz( :, :, 3) + A*(2/(DT^2)-wpsquaredEz/2).*Ez(:, :, n0) - A*(1/(DT^2)+wpsquaredEz/4).*Ez (:, :, 3) ) ./ A*( 1/(DT^2) + wpsquaredEz/4) );
    % Ez calculation outside cylinder.
%     Ez ( :, :, n1 ) = ~cmaskEz .* (1/e0) * (ezzEz) .* Dz ( :, :, n1 );
    % Ez calculation inside cylinder using Drude dispersion model.
    % Ez ( :, :, n1 ) = cmaskEz .* ( (1/(e0*DT^2))*Dz ( :, :, n1 ) - (2/(e0*(DT^2)))*Dz ( :, :, n0) + (1/(e0*(DT^2)))*Dz( :, :, 3) + A*(2/(DT^2)-wpsquaredEz/2).*Ez(:, :, n0) - A*(1/(DT^2)+wpsquaredEz/4).*Ez (:, :, 3) ) ./ A*( 1/(DT^2) + wpsquaredEz/4);


%     if ( n < NHW )
    Ez ( :, Js, n1 ) = Ez ( :, Js, n1 ) + 1 * sin ( TwoPIFDeltaT * n );
    Dz ( :, Js, n1 ) = e0 * Ez ( :, Js, n1 );
%     end

    %Ez ( :, :, n1 ) = smaskEz (:, :) .* Ez ( :, :, n1 );
    %Dz ( :, :, n1 ) = smaskEz (:, :) .* Dz ( :, :, n1 );

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