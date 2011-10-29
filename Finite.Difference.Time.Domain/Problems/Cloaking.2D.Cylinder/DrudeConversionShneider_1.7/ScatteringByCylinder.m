% Based on FDTD example 3.7 from Numerical techniques in Electromagnetics by Sadiku.
% Main m file for the FDTD simulation.
% Now this doesn't look at all like the original problem. Only thing
% similar are some of the constants and field orientation.

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
NNMax = 1000;                   % Maximum time.
TimeResolutionFactor = 1;      % E field snapshots will be saved every x frames where x is time resolution factor.
ResolutionFactor = 5;          % Resolution of plotted field is divided by this factor.
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
A = (rb/(rb-ra))^2;

% ====================== Data arrays =========================
uxxHx = zeros ( IHx, JHx );  % uxx for Hx
uxyHy = zeros ( IHy, JHy );  % uxy for Hy
uyxHx = zeros ( IHx, JHx );  % uyx for Hx
uyyHy = zeros ( IHy, JHy );  % uyy for Hy
ezzEz = zeros ( IEz, JEz );  % ezz for Ez

wpsquaredEz = zeros( IEz, JEz );
wpmsquaredHx = zeros (IHx, JHx );
wpmsquaredHy = zeros ( IHy, JHy );

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

Jpz = zeros ( IEz, JEz, 3 );
Jmx = zeros ( IHx, JHx, 3 );
Jmy = zeros ( IHy, JHy, 3 );

EzSnapshots = zeros ( IEz/ResolutionFactor, JEz/ResolutionFactor, NNMax/TimeResolutionFactor ); % E field Snapshot storage.

% Space averaged B.
BxAve = zeros ( IHx, JHx, 3);
ByAve = zeros ( IHy, JHy, 3);
% =============================================================

% ############ Initialization #############
fprintf ( 1, 'Initializing...' );
fprintf ( 1, '\nInitializing parametric arrays...' );
% Initializing er array.
for i=1:IHy     % IHy is size+1 or maximum I size.
    for j=1:JHy
        
        % Ez related parameters.
        if i <= IEz
            ezzEz ( i, j ) = ezz ( i, j-0.5 );
    %         erEz ( i, j ) = er( i, j-0.5 );
            erEz ( i, j ) = er( i, j-0.5 );
            smaskEz( i, j ) = s ( i, j-0.5 );
            wpsquaredEz(i, j) = wpsquared(i, j-0.5, 2*pi*f);
            cmaskEz( i, j ) = iscylinder (i, j-0.5);
            %er ( i, j-0.5 )
        end
        
        % Hx related parameters.
        if i <= IHx && j <= JHx
            invurHx = inv ( ur ( i, j-0.5 ));
            uxxHx (i, j) = invurHx(1, 1);
            uyxHx (i, j) = invurHx(2, 1);
            cmaskHx (i, j) = iscylinder (i, j-0.5);
            wpmsquaredHx (i, j) = wpmsquareduxx (i, j-0.5, 2*pi*f);
    %         smaskHx( i, j ) = s ( i, j-0.5 );
        end
        
        % Hy related parameters.
        invurHy = inv ( ur ( i-0.5, j-1 ));
        uxyHy (i, j) = invurHy(1, 2);
        uyyHy (i, j) = invurHy(2, 2);
        cmaskHy(i, j) = iscylinder (i-0.5, j-1);
        wpmsquaredHy (i, j) = wpmsquareduyy (i-0.5, j-1, 2*pi*f);
        %         smaskHy( i, j ) = s ( i-0.5, j-1 );
       
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
figure (3)
mesh ( wpmsquaredHx )
title ( 'wpmsquaredHx' )
view (4, 4)

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
%     Bx ( :, :, n1 ) = Bx ( :, :, n0 ) + ( (DT/delta) * ( Ez ( :, 1:JHx, n0 ) - Ez ( :, 2:JHx+1, n0 ) ));
    %Bx ( :, :, n1 ) = smaskHx ( :, : ) .* Bx ( :, :, n1 );
    
    % Hx and Hy
    Hx ( :, :, n1 ) = Hx ( :, :, n0 ) - ( (DT/(delta*u0*uinf)) * ((Ez(:, 2:JHx+1, n0) - Ez(:, 1:JHx, n0)) + delta*Jmx(:, :, n0)) );
    Hy ( 2:IHy-1, :, n1 ) = Hy ( 2:IHy-1, :, n0 ) - ( (DT/(delta*u0*uinf)) * ((Ez(1:IHy-2, :, n0) - Ez(2:IHy-1, :, n0)) + delta*Jmy(2:IHy-1, :, n0)) );
    % * Calculation of By.
%     By ( 2:IHy-1, :, n1 ) = By ( 2:IHy-1, :, n0 ) + ( (DT/delta) * ( Ez ( 2:IHy-1, :, n0 ) - Ez ( 1:IHy-2, :,n0 ) ));
    %By ( :, :, n1 ) = smaskHy ( :, : ) .* By ( :, :, n1 );
   
    % Boundary conditions on By. Soft grid truncation.
    Hy ( 1, 2:JHy-1, n1 ) = (1/3) * ( Hy ( 2, 1:JHy-2, n0 ) + Hy ( 2, 2:JHy-1, n0 ) + Hy ( 2, 3:JHy, n0 ) );
    Hy ( IHy, 2:JHy-1, n1 ) = (1/3) * ( Hy ( IHy-1, 1:JHy-2, n0 ) + Hy ( IHy-1, 2:JHy-1, n0 ) + Hy ( IHy-1, 3:JHy, n0 ) );
    Hy ( 1, 1, n1 ) = (1/2) * ( Hy ( 2, 1, n0 ) + Hy ( 2, 2, n0 ) );
    Hy ( 1, JHy, n1 ) = (1/2) * ( Hy ( 2, JHy, n0 ) + Hy ( 2, JHy-1, n0 ) );
    Hy ( IHy, 1, n1 ) = (1/2) * ( Hy ( IHy-1, 1, n0 ) + Hy( IHy-1, 2, n0 ) );
    Hy ( IHy, JHy, n1 ) = (1/2) * ( Hy ( IHy-1, JHy, n0 ) + Hy ( IHy-1, JHy-1, n0 ) );
    
    Jmx ( :, :, n1 ) = Jmx ( :, :, n0 ) + (u0*DT)/2 .* wpmsquaredHx .* ( Hx(:, :, n1) + Hx (:, :, n0));
    Jmy ( 2:IHy-1, :, n1 ) = Jmy ( 2:IHy-1, :, n0 ) + (u0*DT)/2 .* wpmsquaredHy(2:IHy-1,:) .* ( Hy(2:IHy-1, :, n1) + Hy (2:IHy-1, :, n0));
%     
%     Hy ( 1, :, n1 ) = Hy ( 2, :, n0 );
%     Hy ( IHy, :, n1 ) = Hy ( IHy-1, :, n0 );

    % Space averaged B fields.
%     BxAve (2:IHx-1, 2:JHx-1, n1) = ( Bx(2:IHx-1, 2:JHx-1, n1) + Bx(2:IHx-1, 3:JHx, n1) + Bx(1:IHx-2, 2:JHx-1, n1) + Bx(1:IHx-2, 3:JHx, n1) )/4;
%     BxAve (2:IHx-1, 2:JHx-1, n0) = ( Bx(2:IHx-1, 2:JHx-1, n0) + Bx(2:IHx-1, 3:JHx, n0) + Bx(1:IHx-2, 2:JHx-1, n0) + Bx(1:IHx-2, 3:JHx, n0) )/4;
%     BxAve (2:IHx-1, 2:JHx-1, 3) = ( Bx(2:IHx-1, 2:JHx-1, 3) + Bx(2:IHx-1, 3:JHx, 3) + Bx(1:IHx-2, 2:JHx-1, 3) + Bx(1:IHx-2, 3:JHx, 3) )/4;
%     
%     ByAve (2:IHx-1, 2:JHx-1, n1) = ( Hy(2:IHx-1, 2:JHx-1, n1) + Hy(3:IHx, 3:JHx, n1) + Hy(2:IHx-1, 1:JHx-2, n1) + Hy(3:IHx, 1:JHx-2, n1) )/4;
%     ByAve (2:IHx-1, 2:JHx-1, n0) = ( Hy(2:IHx-1, 2:JHx-1, n0) + Hy(3:IHx, 3:JHx, n0) + Hy(2:IHx-1, 1:JHx-2, n0) + Hy(3:IHx, 1:JHx-2, n0) )/4;
%     ByAve (2:IHx-1, 2:JHx-1, 3) = ( Hy(2:IHx-1, 2:JHx-1, 3) + Hy(3:IHx, 3:JHx, 3) + Hy(2:IHx-1, 1:JHx-2, 3) + Hy(3:IHx, 1:JHx-2, 3) )/4;
    
    % Drude Model from paper.        
%     Hx ( :, :, n1 ) = (1/u0) * ~cmaskHx .* (uxxHx .* Bx ( :, :, n1 ) + uxyHy (1:IHy-1, 1:JHy-1 ) .* Hy (1:IHy-1, 1:JHy-1, n1 )) + cmaskHx .* ( ax.*Bx ( :, :, n1 ) + bx.*Bx ( :, :, n0) + cx.*Bx ( :, :, 3) + dx.*ByAve (1:IHy-1, 1:JHy-1, n1) + ex.*ByAve (1:IHy-1, 1:JHy-1, n0) + fx.*ByAve (1:IHy-1, 1:JHy-1, 3) - (gx.*Hx(:,:,n0) + hx.*Hx(:,:,3)) ) ./ lx;
%     Hx ( :, :, n1 ) = (1/u0) * (uxxHx .* Bx ( :, :, n1 ) + uxyHy (1:IHy-1, 1:JHy-1 ) .* Hy (1:IHy-1, 1:JHy-1, n1 ));
%     ****** Major ****** Hx ( :, :, n1 ) = (1/u0) * (uxxHx .* Bx ( :, :, n1 ) + uxyHy (1:IHy-1, 1:JHy-1 ) .* Hy (1:IHy-1, 1:JHy-1, n1 ));
    %Hx ( :, :, n1 ) = smaskHx ( :, : ) .* Hx ( :, :, n1 );
%     Hy ( 1:IHy-1, 1:JHy-1, n1 ) = (1/u0) * ~cmaskHy (1:IHy-1, 1:JHy-1) .* (uyxHx.*Bx ( :, :, n1 ) + uyyHy (1:IHy-1, 1:JHy-1) .* Hy (1:IHy-1, 1:JHy-1, n1 )) + cmaskHy (1:IHy-1, 1:JHy-1) .* ( ay(1:IHy-1, 1:JHy-1).*Hy ( 1:IHy-1, 1:JHy-1, n1 ) + by(1:IHy-1, 1:JHy-1).*Hy ( 1:IHy-1, 1:JHy-1, n0) + cy(1:IHy-1, 1:JHy-1).*Hy ( 1:IHy-1, 1:JHy-1, 3) + dy(1:IHy-1, 1:JHy-1).*BxAve (1:IHy-1, 1:JHy-1, n1) + ey(1:IHy-1, 1:JHy-1).*BxAve (1:IHy-1, 1:JHy-1, n0) + fy(1:IHy-1, 1:JHy-1).*BxAve (1:IHy-1, 1:JHy-1, 3) - (gy(1:IHy-1, 1:JHy-1).*Hy(1:IHy-1, 1:JHy-1,n0) + hy(1:IHy-1, 1:JHy-1).*Hy(1:IHy-1, 1:JHy-1,3)) ) ./ ly (1:IHy-1, 1:JHy-1);
    
%     Hy (1:IHy-1, 1:JHy-1, n1) = (1/u0) * (uyxHx.*Bx ( :, :, n1 ) + uyyHy (1:IHy-1, 1:JHy-1) .* Hy (1:IHy-1, 1:JHy-1, n1 ));
    %Hy ( :, :, n1 ) = (1/u0) * Hy ( :, :, n1 );
%     ****** Major ******    Hy (1:IHy-1, 1:JHy-1, n1) = (1/u0) * (uyxHx.*Bx ( :, :, n1 ) + uyyHy (1:IHy-1, 1:JHy-1) .* Hy (1:IHy-1, 1:JHy-1, n1 ));
    
    % ABC for Hy;
%     Hy (1:IHy-1, JHy, n1) = Hy (1:IHy-1, JHy, n0);
    %Hy ( :, :, n1 ) = smaskHy ( :, : ) .* Hy ( :, :, n1 );
    
%     Dz ( :, 2:JEz-1, n1 ) = (  Dz ( :, 2:JEz-1, n0 ) ) + ( (DT/delta) * ( Hy ( 2:JEz+1, 2:JEz-1, n1 ) - Hy ( 1:JEz, 2:JEz-1, n1 ) + Hx ( :, 1:JEz-2,n1 ) - Hx ( :, 2:JEz-1, n1 ) ));
    %Dz (:, :, n1) = Dz (:, :, n1) .* smask;
    
    Ez ( :, 2:JEz-1, n1 ) = Ez ( :, 2:JEz-1, n0 ) + ( (DT/(delta*e0*einf)) * ( Hy ( 2:JEz+1, 2:JEz-1, n1 ) - Hy ( 1:JEz, 2:JEz-1, n1 ) + Hx ( :, 1:JEz-2,n1 ) - Hx ( :, 2:JEz-1, n1 ) - (delta/2)*(Jpz(:, 2:JEz-1, n1)+Jpz(:, 2:JEz-1, n0)) ./ (erEz(:,2:JEz-1))));
    
    % Boundary conditions on Dz. Soft grid truncation.
    Ez ( 2:IEz-1, 1, n1 ) = (1/3) * ( Ez ( 1:IEz-2, 2, n0 ) + Ez ( 2:IEz-1, 2, n0 ) + Ez ( 3:IEz, 2, n0 ) );
    Ez ( 2:IEz-1, JEz, n1 ) = (1/3) * ( Ez ( 1:IEz-2, JEz-1, n0 ) + Ez ( 2:IEz-1, JEz-1, n0 ) + Ez ( 3:IEz, JEz-1, n0 ) );
    Ez ( 1, 1, n1 ) = (1/2) * ( Ez ( 1, 2, n0 ) + Ez ( 2, 2, n0 ) );
    Ez ( IEz, 1, n1 ) = (1/2) * ( Ez ( IEz, 2, n0 ) + Ez ( IEz-1, 2, n0 ) );
    Ez ( 1, JEz, n1 ) = (1/2) * ( Ez ( 1, JEz-1, n0 ) + Ez ( 2, JEz-1, n0 ) );
    Ez ( IEz, JEz, n1 ) = (1/2) * ( Ez ( IEz, JEz-1, n0 ) + Ez ( IEz-1, JEz-1, n0 ) );
    
    Jpz ( :, :, n1 ) = Jpz ( :, :, n0 ) + (e0*DT)/2 .* wpsquaredEz .* ( Ez(:, :, n1) + Ez (:, :, n0));
%     Ez ( :, 1, n1 ) = Ez ( :, 2, n0 );
%     Ez ( :, IEz, n1 ) = Ez ( :, IEz-1, n0 );
    % ************************************************

    
%     Ez ( :, :, n1 ) = (1/e0) * Ez ( :, :, n1 ) ./ (ezzEz);
%     %     ****** Major ******   Ez ( :, :, n1 ) = (1/e0) * Ez ( :, :, n1 ) ./ (ezzEz);
%     Ez ( :, :, n1 ) = (1/e0) * Ez ( :, :, n1 ) ./ (ezzEz);
%    Ez ( :, :, n1 ) =  ( ~cmaskEz .* (1/e0) * (ezzEz) .* Ez ( :, :, n1 ) ) + ( cmaskEz .* ( (1/(e0*DT^2))*Ez ( :, :, n1 ) - (2/(e0*(DT^2)))*Ez ( :, :, n0) + (1/(e0*(DT^2)))*Ez( :, :, 3) + A*(2/(DT^2)-wpsquaredEz/2).*Ez(:, :, n0) - A*(1/(DT^2)+wpsquaredEz/4).*Ez (:, :, 3) ) ./ A*( 1/(DT^2) + wpsquaredEz/4) );
%     ***** Drude ***** Ez ( :, :, n1 ) =  ( ~cmaskEz .* (1/e0) * (ezzEz) .* Ez ( :, :, n1 ) ) + ( cmaskEz .* ( (1/(e0*DT^2))*Ez ( :, :, n1 ) - (2/(e0*(DT^2)))*Ez ( :, :, n0) + (1/(e0*(DT^2)))*Ez( :, :, 3) + A*(2/(DT^2)-wpsquaredEz/2).*Ez(:, :, n0) - A*(1/(DT^2)+wpsquaredEz/4).*Ez (:, :, 3) ) ./ A*( 1/(DT^2) + wpsquaredEz/4) );
    % Ez calculation outside cylinder.
%     Ez ( :, :, n1 ) = ~cmaskEz .* (1/e0) * (ezzEz) .* Ez ( :, :, n1 );
    % Ez calculation inside cylinder using Drude dispersion model.
    % Ez ( :, :, n1 ) = cmaskEz .* ( (1/(e0*DT^2))*Ez ( :, :, n1 ) - (2/(e0*(DT^2)))*Ez ( :, :, n0) + (1/(e0*(DT^2)))*Ez( :, :, 3) + A*(2/(DT^2)-wpsquaredEz/2).*Ez(:, :, n0) - A*(1/(DT^2)+wpsquaredEz/4).*Ez (:, :, 3) ) ./ A*( 1/(DT^2) + wpsquaredEz/4);


%     if ( n < NHW )
    Ez ( :, Js, n1 ) = Ez ( :, Js, n1 ) + 1 * sin ( TwoPIFDeltaT * n );
%     Ez ( :, Js, n1 ) = e0 * Ez ( :, Js, n1 );
%     end

    %Ez ( :, :, n1 ) = smaskEz (:, :) .* Ez ( :, :, n1 );
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