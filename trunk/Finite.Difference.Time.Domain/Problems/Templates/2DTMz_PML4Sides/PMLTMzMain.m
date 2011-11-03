% 2D FDTD TMz simulation employing a PML. The solution region is bounded on all four sides with a PML and terminated by PEC.

clc
clear all

% ============== Simulation related parameters ================
[Size XCenter YCenter delta ra rb DTp PMLw] = Parameters;

IHx = Size+2*PMLw;
JHx = Size+2*PMLw-1;
IHy = Size+2*PMLw-1;
JHy = Size+2*PMLw;
IEz = Size+2*PMLw;
JEz = Size+2*PMLw;

% Time indices for field calculation.
n0 = 1;
n1 = 2;
NNMax = 200;                   % Maximum time.
TimeResolutionFactor = 1;      % E field snapshots will be saved every x frames where x is time resolution factor.
ResolutionFactor = 1;          % Resolution of plotted field is divided by this factor.
Js = 50;                       % J-position of the plane wave front.
Is = 50;                       % I-position of the plane wave front.

% Different Constants.
Cl = 299792458;
f = 2.0e9;
pi = 3.141592654;
e0 = (1e-9) / (36*pi);
u0 = (1e-7) * 4 * pi;
DT = DTp;
TwoPIFDeltaT = 2 * pi * f * DT;
NHW = 1/(2 * f * DT); % One half wave cycle.

% ====================== Data arrays =========================

% ------------ Field-specific parameters ------------
% Permeability.
urHx  = 1 + zeros (IHx, JHx);
urHy = 1 + zeros (IHy, JHy);
% Permittivity.
erEz = 1 + zeros (IEz, JEz);
% Magnetic conductance.
smHx  = zeros (IHx, JHx);
smHy =  zeros (IHy, JHy);
% Conductance.
sEz = zeros (IEz, JEz);
% PML conductance arrays.
%
%   Not yet declared.
%

% PEC mask array.
PECmaskEz = zeros ( IEz, JEz );
% ---------------------------------------------------

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
fprintf ( 1, '\nInitializing parametric arrays...' );
% Initializing er array.
for i=1:IHx
    fprintf ( 1, '%g %% \n', ((i-1)*100)/IHx );
    for j=1:JHx      
               
    end
end

figure (1)
mesh ( urHx )
title ( 'ur' )
view (4, 4)
figure (2)
mesh ( erEz )
title ( 'er' )
view (4, 4)
% figure (3)
% mesh ( uphiHx )
% title ( 'uphi' )
% view (4, 4)
% figure (4)
% mesh ( wpmsquaredHx )
% title ( 'wpm squared' )
% view (4, 4)
% figure (5)
% mesh ( wpsquaredEz )
% title ( 'wp squared' )
% view (4, 4)

fprintf ( 1, 'Initialization done.\n' );
% ############ Initialization Complete ##############

% ########### 2. Simulation #############
fprintf ( 1, 'Simulation started... \n' );
for n=0:NNMax-2
    fprintf ( 1, '%g %% \n', (n*100)/NNMax );
    
    % Copying n1 fields as past fields. Will be used in second order time derivatives.
    Bx (:, :, 3) = Bx (:, :, n1);
    By (:, :, 3) = By (:, :, n1);
    Hx (:, :, 3) = Hx (:, :, n1);
    Hy (:, :, 3) = Hy (:, :, n1);
    Dz (:, :, 3) = Dz (:, :, n1);
    Ez (:, :, 3) = Ez (:, :, n1);
    
    % Bx
    Bx (2:IHx-1, 1:JHx, n1) = Bx (2:IHx-1, 1:JHx, n0) + ( (DT/(delta)) * ( Ez (2:IEz-1, 1:JEz-1, n0) - Ez (2:IEz-1, 2:JEz, n0) )) + (DT*smHx (2:IHx-1, 1:JHx).*Hx (2:IHx-1, 1:JHx, n0));
    %Bx ( :, :, n1 ) = smaskHx ( :, : ) .* Bx ( :, :, n1 );
    
    % By
    By (1:IHy, 2:JHy-1, n1) = By (1:IHy, 2:JHy-1, n0) + ( (DT/(delta)) * ( Ez (2:IEz, 2:JEz-1, n0) - Ez (1:IEz-1, 2:JEz-1, n0) )) + (DT*smHy (1:IHy, 2:JHy-1).*Hy (1:IHy, 2:JHy-1, n0));
    
    % Magnetic fields.
    Hx (:, :, n1) = (1/u0) * Bx (:, :, n1) ./ (urHx);
    Hy (:, :, n1) = (1/u0) * By (:, :, n1) ./ (urHy);
%     By ( 2:IHy-1, :, n1 ) = By ( 2:IHy-1, :, n0 ) + ( (DT/delta) * ( Ez ( 2:IHy-1, :, n0 ) - Ez ( 1:IHy-2, :,n0 ) ));
    %By ( :, :, n1 ) = smaskHy ( :, : ) .* By ( :, :, n1 );
   
    % Boundary conditions on By. Soft grid truncation.
%     By ( 1, 2:JHy-1, n1 ) = (1/3) * ( By ( 2, 1:JHy-2, n0 ) + By ( 2, 2:JHy-1, n0 ) + By ( 2, 3:JHy, n0 ) );
%     By ( IHy, 2:JHy-1, n1 ) = (1/3) * ( By ( IHy-1, 1:JHy-2, n0 ) + By ( IHy-1, 2:JHy-1, n0 ) + By ( IHy-1, 3:JHy, n0 ) );
%     By ( 1, 1, n1 ) = (1/2) * ( By ( 2, 1, n0 ) + By ( 2, 2, n0 ) );
%     By ( 1, JHy, n1 ) = (1/2) * ( By ( 2, JHy, n0 ) + By ( 2, JHy-1, n0 ) );
%     By ( IHy, 1, n1 ) = (1/2) * ( By ( IHy-1, 1, n0 ) + By( IHy-1, 2, n0 ) );
%     By ( IHy, JHy, n1 ) = (1/2) * ( By ( IHy-1, JHy, n0 ) + By ( IHy-1, JHy-1, n0 ) );

    % Space averaged B fields.
%     BxAve (2:IHx-1, 2:JHx-1, n1) = ( Bx(3:IHx, 2:JHx-1, n1) + Bx(2:IHx-1, 3:JHx, n1) + Bx(3:IHx, 3:JHx, n1) + Bx(2:IHx-1, 2:JHx-1, n1) )/4;
%     BxAve (2:IHx-1, 2:JHx-1, n0) = ( Bx(3:IHx, 2:JHx-1, n0) + Bx(2:IHx-1, 3:JHx, n0) + Bx(3:IHx, 3:JHx, n0) + Bx(2:IHx-1, 2:JHx-1, n0) )/4;
%     BxAve (2:IHx-1, 2:JHx-1, 3) = ( Bx(3:IHx, 2:JHx-1, 3) + Bx(2:IHx-1, 3:JHx, 3) + Bx(3:IHx, 3:JHx, 3) + Bx(2:IHx-1, 2:JHx-1, 3) )/4;
%     
%     ByAve (2:IHx-1, 2:JHx-1, n1) = ( By(1:IHx-2, 2:JHx-1, n1) + By(2:IHx-1, 1:JHx-2, n1) + By(1:IHx-2, 1:JHx-2, n1) + By(2:IHx-1, 2:JHx-1, n1) )/4;
%     ByAve (2:IHx-1, 2:JHx-1, n0) = ( By(1:IHx-2, 2:JHx-1, n0) + By(2:IHx-1, 1:JHx-2, n0) + By(1:IHx-2, 1:JHx-2, n0) + By(2:IHx-1, 2:JHx-1, n0) )/4;
%     ByAve (2:IHx-1, 2:JHx-1, 3) = ( By(1:IHx-2, 2:JHx-1, 3) + By(2:IHx-1, 1:JHx-2, 3) + By(1:IHx-2, 1:JHx-2, 3) + By(2:IHx-1, 2:JHx-1, 3) )/4;

%     Hx ( :, :, n1 ) = ( ax.*Bx ( :, :, n1 ) + bx.*Bx ( :, :, n0) + cx.*Bx ( :, :, 3) + dx.*ByAve (1:IHy-1, 1:JHy-1, n1) + ex.*ByAve (1:IHy-1, 1:JHy-1, n0) + fx.*ByAve (1:IHy-1, 1:JHy-1, 3) - (gx.*Hx(:,:,n0) + hx.*Hx(:,:,3)) ) ./ lx;
% 
%     Hy ( 1:IHy-1, 1:JHy-1, n1 ) = ( ay(1:IHy-1, 1:JHy-1).*By ( 1:IHy-1, 1:JHy-1, n1 ) + by(1:IHy-1, 1:JHy-1).*By ( 1:IHy-1, 1:JHy-1, n0) + cy(1:IHy-1, 1:JHy-1).*By ( 1:IHy-1, 1:JHy-1, 3) + dy(1:IHy-1, 1:JHy-1).*BxAve (1:IHy-1, 1:JHy-1, n1) + ey(1:IHy-1, 1:JHy-1).*BxAve (1:IHy-1, 1:JHy-1, n0) + fy(1:IHy-1, 1:JHy-1).*BxAve (1:IHy-1, 1:JHy-1, 3) - (gy(1:IHy-1, 1:JHy-1).*Hy(1:IHy-1, 1:JHy-1,n0) + hy(1:IHy-1, 1:JHy-1).*Hy(1:IHy-1, 1:JHy-1,3)) ) ./ ly (1:IHy-1, 1:JHy-1);
%     
%     Dz ( :, 2:JEz-1, n1 ) = (  Dz ( :, 2:JEz-1, n0 ) ) + ( (DT/delta) * ( Hy ( 2:JEz+1, 2:JEz-1, n1 ) - Hy ( 1:JEz, 2:JEz-1, n1 ) + Hx ( :, 1:JEz-2,n1 ) - Hx ( :, 2:JEz-1, n1 ) ));
    
    Dz(2:IEz-1, 2:JEz-1, n1) = Dz(2:IEz-1, 2:JEz-1, n0) + ((DT/(delta))*( Hy(2:IHy, 2:JHy-1, n1) - Hy(1:IHy-1, 2:JHy-1, n1) - Hx(2:IHx-1, 2:JHx, n1) + Hx(2:IHx-1, 1:JHx-1, n1) )) + (DT*sEz(2:IEz-1, 2:JEz-1).*Ez(2:IEz-1, 2:JEz-1, n0));
    % Boundary conditions on Dz. Soft grid truncation.
%     Dz ( 2:IEz-1, 1, n1 ) = (1/3) * ( Dz ( 1:IEz-2, 2, n0 ) + Dz ( 2:IEz-1, 2, n0 ) + Dz ( 3:IEz, 2, n0 ) );
%     Dz ( 2:IEz-1, JEz, n1 ) = (1/3) * ( Dz ( 1:IEz-2, JEz-1, n0 ) + Dz ( 2:IEz-1, JEz-1, n0 ) + Dz ( 3:IEz, JEz-1, n0 ) );
%     Dz ( 1, 1, n1 ) = (1/2) * ( Dz ( 1, 2, n0 ) + Dz ( 2, 2, n0 ) );
%     Dz ( IEz, 1, n1 ) = (1/2) * ( Dz ( IEz, 2, n0 ) + Dz ( IEz-1, 2, n0 ) );
%     Dz ( 1, JEz, n1 ) = (1/2) * ( Dz ( 1, JEz-1, n0 ) + Dz ( 2, JEz-1, n0 ) );
%     Dz ( IEz, JEz, n1 ) = (1/2) * ( Dz ( IEz, JEz-1, n0 ) + Dz ( IEz-1, JEz-1, n0 ) );
%     
%     Ez ( :, :, n1 ) =  ( (1/(e0*(DT^2)))*Dz ( :, :, n1 ) - (2/(e0*(DT^2)))*Dz ( :, :, n0) + (1/(e0*(DT^2)))*Dz( :, :, 3) + AEz.*(2/(DT^2)-wpsquaredEz/2).*Ez(:, :, n0) - AEz.*(1/(DT^2)+wpsquaredEz/4).*Ez (:, :, 3) ) ./ (AEz.*( 1/(DT^2) + wpsquaredEz/4));

    Ez (:, :, n1) = (1/e0) * Dz (:, :, n1) ./ (erEz);
    % Comment out the if statement for a continuous source. Otherwise, a single pulse will be used.
%     if ( n < NHW )
    Ez (Is, 2:IEz-1, n1) = Ez (Is, 2:IEz-1, n1) + 1 * sin ( TwoPIFDeltaT * n );
%     Hy (2:IEz-1, Js-1, n1) = sqrt(e0/u0)*Ez (2:IEz-1, Js, n1);
%     Hx (2:IEz-1, Js, n1) = sqrt(e0/u0)*Ez (2:IEz-1, Js, n1);
    Dz (Is, 2:IEz-1, n1) = e0 * Ez (Is, 2:IEz-1, n1);
    
%     Ez (Js, Js+1, n1) = Ez (Js, Js, n1) + 1 * sin ( TwoPIFDeltaT * n );
%     Dz (Js, Js+1, n1) = e0*Ez (Js+1, Js, n1);
%     end
%     Ez (:, Js, n1) = Ez (:, Js, n1) + 1 * sin ( TwoPIFDeltaT * n );
%     Dz (:, Js, n1) = e0 * Ez (:, Js, n1);
%     Ez(2, 2:JEz-1, n1) = Ez(2, 2:JEz-1, n1) + sin ( TwoPIFDeltaT * n );

    % Uncomment this to zero out the field at PEC points. PEC points can be defined in s.m file.
%     Ez ( :, :, n1 ) = smaskEz (:, :) .* Ez ( :, :, n1 );
%     Dz ( :, :, n1 ) = smaskEz (:, :) .* Dz ( :, :, n1 );

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
    xlabel ('y-axis')
    ylabel ('x-axis')
    
    figure (7)
    surf ( EzSnapshots (:, :, i) );
    view (0, 90)
    zlim ( [-10 10] )
    xlabel ('y-axis')
    ylabel ('x-axis')
    
end
fprintf ( 1, 'Simulation completed! \n' );