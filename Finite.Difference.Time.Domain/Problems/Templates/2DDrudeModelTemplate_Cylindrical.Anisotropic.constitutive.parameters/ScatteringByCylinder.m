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
NNMax = 500;                   % Maximum time.
TimeResolutionFactor = 1;      % E field snapshots will be saved every x frames where x is resolution factor.
ResolutionFactor = 2;          % Resolution of plotted field is divided by this factor.
Js = 2;                         % J-position of the plane wave front.
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
uxxHx = zeros ( IHx, JHx );  % uxx for Hx
uxyHy = zeros ( IHy, JHy );  % uxy for Hy
uyxHx = zeros ( IHx, JHx );  % uyx for Hx
uyyHy = zeros ( IHy, JHy );  % uyy for Hy

urrHx  = zeros ( IHx, JHx );
uphiHx = zeros ( IHx, JHx );
ezzEz = zeros ( IEz, JEz );  % ezz for Ez
AEz = zeros ( IEz, JEz );

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
wpmsquaredHx = zeros (IHx, JHx );

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

% Space averaged B.
BxAve = zeros ( IHx, JHx, 3);
ByAve = zeros ( IHy, JHy, 3);
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
            ezzEz ( i, j ) = ezz ( i, j-0.5 );
    %         erEz ( i, j ) = er( i, j-0.5 );
            erEz ( i, j ) = ezz( i, j-0.5 );
            smaskEz( i, j ) = s ( i, j-0.5 );
            wpsquaredEz(i, j) = wpsquared(i, j-0.5, 2*pi*f);
            cmaskEz( i, j ) = iscylinder (i, j-0.5);
            AEz (i, j) = A (i, j-0.5);
            %er ( i, j-0.5 )
        end
        
        % Hx related parameters.
        if i <= IHx && j <= JHx
            invurHx = inv ( ur ( i, j-0.5 ));
            uxxHx (i, j) = invurHx(1, 1);
            uyxHx (i, j) = invurHx(2, 1);
            urrHx (i, j) = urr (i, j-0.5);
            uphiHx (i, j) = uphi (i, j-0.5);
            cmaskHx (i, j) = iscylinder (i, j-0.5);
            wpmsquaredHx (i, j) = wpmsquared (i, j-0.5, 2*pi*f);
    %         smaskHx( i, j ) = s ( i, j-0.5 );
            ax (i, j) = (sinphi(i, j-0.5))^2 * ( 1/(DT^2) + wpmsquared(i, j-0.5, 2*pi*f)/4) + uphi(i, j-0.5) * cosphi(i, j-0.5)^2 * (1/(DT^2));
            bx (i, j) = (sinphi(i, j-0.5))^2 * ( -2/(DT^2) + wpmsquared(i, j-0.5, 2*pi*f)/2) - uphi(i, j-0.5) * cosphi(i, j-0.5)^2 * (2/(DT^2));
            cx (i, j) = ax (i, j);
            dx (i, j) = ( uphi(i, j-0.5)/(DT^2) - ( 1/(DT^2) + (wpmsquared(i, j-0.5, 2*pi*f))/4 ) ) * sinphi(i, j-0.5) * cosphi(i, j-0.5);
            ex (i, j) = ( -2 * uphi(i, j-0.5)/(DT^2) - ( -2/(DT^2) + (wpmsquared(i, j-0.5, 2*pi*f))/2 ) ) * sinphi(i, j-0.5) * cosphi(i, j-0.5);
            fx (i, j) = dx (i, j);
            gx (i, j) = u0*uphi (i, j-0.5) * ( -2/(DT^2) + (wpmsquared(i, j-0.5, 2*pi*f))/2 );
            hx (i, j) = u0*uphi (i, j-0.5) * ( 1/(DT^2) + (wpmsquared(i, j-0.5, 2*pi*f))/4 );
            lx (i, j) = hx (i, j);
        end
        
        % Hy related parameters.
        invurHy = inv ( ur ( i-0.5, j-1 ));
        uxyHy (i, j) = invurHy(1, 2);
        uyyHy (i, j) = invurHy(2, 2);
        cmaskHy(i, j) = iscylinder (i-0.5, j-1);
        %         smaskHy( i, j ) = s ( i-0.5, j-1 );
        ay (i, j) = (cosphi(i-0.5, j-1))^2 * ( 1/(DT^2) + wpmsquared(i-0.5, j-1, 2*pi*f)/4) + uphi(i-0.5, j-1) * sinphi(i-0.5, j-1)^2 * (1/(DT^2));
        by (i, j) = (cosphi(i-0.5, j-1))^2 * ( -2/(DT^2) + wpmsquared(i-0.5, j-1, 2*pi*f)/2) - uphi(i-0.5, j-1) * sinphi(i-0.5, j-1)^2 * (2/(DT^2));
        cy (i, j) = ay (i, j);
        dy (i, j) = ( uphi(i-0.5, j-1)/(DT^2) - ( 1/(DT^2) + (wpmsquared(i-0.5, j-1, 2*pi*f))/4 ) ) * sinphi(i-0.5, j-1) * cosphi(i-0.5, j-1);
        ey (i, j) = ( -2 * uphi(i-0.5, j-1)/(DT^2) - ( -2/(DT^2) + (wpmsquared(i-0.5, j-1, 2*pi*f))/2 ) ) * sinphi(i-0.5, j-1) * cosphi(i-0.5, j-1);
        fy (i, j) = dy (i, j);
        gy (i, j) = u0*uphi (i-0.5, j-1) * ( -2/(DT^2) + (wpmsquared(i-0.5, j-1, 2*pi*f))/2 );
        hy (i, j) = u0*uphi (i-0.5, j-1) * ( 1/(DT^2) + (wpmsquared(i-0.5, j-1, 2*pi*f))/4 );
        ly (i, j) = hy (i, j);
        
    end
end

figure (1)
mesh ( urrHx )
title ( 'urrHx' )
view (4, 4)
figure (2)
mesh ( erEz )
title ( 'ezz' )
view (4, 4)
figure (3)
mesh ( uphiHx )
title ( 'uphiHx' )
view (4, 4)
figure (4)
mesh ( wpmsquaredHx )
title ( 'wpmsquaredHx' )
view (4, 4)
figure (5)
mesh ( wpsquaredEz )
title ( 'wpsquaredEz' )
view (4, 4)
% figure (6)
% mesh ( uxyHy )
% title ( 'uxyHy' )
% view (4, 4)
% figure (7)
% mesh ( uyyHy )
% title ( 'uyyHy' )
% view (4, 4)

fprintf ( 1, 'Initialization done.\n' );
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

    % Space averaged B fields.
%     BxAve (2:IHx-1, 2:JHx-1, n1) = ( Bx(2:IHx-1, 2:JHx-1, n1) + Bx(2:IHx-1, 3:JHx, n1) + Bx(1:IHx-2, 2:JHx-1, n1) + Bx(1:IHx-2, 3:JHx, n1) )/4;
%     BxAve (2:IHx-1, 2:JHx-1, n0) = ( Bx(2:IHx-1, 2:JHx-1, n0) + Bx(2:IHx-1, 3:JHx, n0) + Bx(1:IHx-2, 2:JHx-1, n0) + Bx(1:IHx-2, 3:JHx, n0) )/4;
%     BxAve (2:IHx-1, 2:JHx-1, 3) = ( Bx(2:IHx-1, 2:JHx-1, 3) + Bx(2:IHx-1, 3:JHx, 3) + Bx(1:IHx-2, 2:JHx-1, 3) + Bx(1:IHx-2, 3:JHx, 3) )/4;
%     
%     ByAve (2:IHx-1, 2:JHx-1, n1) = ( By(2:IHx-1, 2:JHx-1, n1) + By(3:IHx, 3:JHx, n1) + By(2:IHx-1, 1:JHx-2, n1) + By(3:IHx, 1:JHx-2, n1) )/4;
%     ByAve (2:IHx-1, 2:JHx-1, n0) = ( By(2:IHx-1, 2:JHx-1, n0) + By(3:IHx, 3:JHx, n0) + By(2:IHx-1, 1:JHx-2, n0) + By(3:IHx, 1:JHx-2, n0) )/4;
%     ByAve (2:IHx-1, 2:JHx-1, 3) = ( By(2:IHx-1, 2:JHx-1, 3) + By(3:IHx, 3:JHx, 3) + By(2:IHx-1, 1:JHx-2, 3) + By(3:IHx, 1:JHx-2, 3) )/4;
    
    BxAve (2:IHx-1, 2:JHx-1, n1) = ( Bx(3:IHx, 2:JHx-1, n1) + Bx(2:IHx-1, 3:JHx, n1) + Bx(3:IHx, 3:JHx, n1) + Bx(2:IHx-1, 2:JHx-1, n1) )/4;
    BxAve (2:IHx-1, 2:JHx-1, n0) = ( Bx(3:IHx, 2:JHx-1, n0) + Bx(2:IHx-1, 3:JHx, n0) + Bx(3:IHx, 3:JHx, n0) + Bx(2:IHx-1, 2:JHx-1, n0) )/4;
    BxAve (2:IHx-1, 2:JHx-1, 3) = ( Bx(3:IHx, 2:JHx-1, 3) + Bx(2:IHx-1, 3:JHx, 3) + Bx(3:IHx, 3:JHx, 3) + Bx(2:IHx-1, 2:JHx-1, 3) )/4;
    
    ByAve (2:IHx-1, 2:JHx-1, n1) = ( By(1:IHx-2, 2:JHx-1, n1) + By(2:IHx-1, 1:JHx-2, n1) + By(1:IHx-2, 1:JHx-2, n1) + By(2:IHx-1, 2:JHx-1, n1) )/4;
    ByAve (2:IHx-1, 2:JHx-1, n0) = ( By(1:IHx-2, 2:JHx-1, n0) + By(2:IHx-1, 1:JHx-2, n0) + By(1:IHx-2, 1:JHx-2, n0) + By(2:IHx-1, 2:JHx-1, n0) )/4;
    ByAve (2:IHx-1, 2:JHx-1, 3) = ( By(1:IHx-2, 2:JHx-1, 3) + By(2:IHx-1, 1:JHx-2, 3) + By(1:IHx-2, 1:JHx-2, 3) + By(2:IHx-1, 2:JHx-1, 3) )/4;
                % Drude Model from paper.        
%     Hx ( :, :, n1 ) = (1/u0) * Bx ( :, :, n1 );
    Hx ( :, :, n1 ) = ( ax.*Bx ( :, :, n1 ) + bx.*Bx ( :, :, n0) + cx.*Bx ( :, :, 3) + dx.*ByAve (1:IHy-1, 1:JHy-1, n1) + ex.*ByAve (1:IHy-1, 1:JHy-1, n0) + fx.*ByAve (1:IHy-1, 1:JHy-1, 3) - (gx.*Hx(:,:,n0) + hx.*Hx(:,:,3)) ) ./ lx;
%     Hy ( 1:IHy-1, 1:JHy-1, n1 ) = (1/u0) * By (1:IHy-1, 1:JHy-1, n1 );
    Hy ( 1:IHy-1, 1:JHy-1, n1 ) = ( ay(1:IHy-1, 1:JHy-1).*By ( 1:IHy-1, 1:JHy-1, n1 ) + by(1:IHy-1, 1:JHy-1).*By ( 1:IHy-1, 1:JHy-1, n0) + cy(1:IHy-1, 1:JHy-1).*By ( 1:IHy-1, 1:JHy-1, 3) + dy(1:IHy-1, 1:JHy-1).*BxAve (1:IHy-1, 1:JHy-1, n1) + ey(1:IHy-1, 1:JHy-1).*BxAve (1:IHy-1, 1:JHy-1, n0) + fy(1:IHy-1, 1:JHy-1).*BxAve (1:IHy-1, 1:JHy-1, 3) - (gy(1:IHy-1, 1:JHy-1).*Hy(1:IHy-1, 1:JHy-1,n0) + hy(1:IHy-1, 1:JHy-1).*Hy(1:IHy-1, 1:JHy-1,3)) ) ./ ly (1:IHy-1, 1:JHy-1);
%     Hx ( :, :, n1 ) = (1/u0) * ~cmaskHx .* (uxxHx .* Bx ( :, :, n1 ) + uxyHy (1:IHy-1, 1:JHy-1 ) .* By (1:IHy-1, 1:JHy-1, n1 )) + cmaskHx .* ( ax.*Bx ( :, :, n1 ) + bx.*Bx ( :, :, n0) + cx.*Bx ( :, :, 3) + dx.*ByAve (1:IHy-1, 1:JHy-1, n1) + ex.*ByAve (1:IHy-1, 1:JHy-1, n0) + fx.*ByAve (1:IHy-1, 1:JHy-1, 3) - (gx.*Hx(:,:,n0) + hx.*Hx(:,:,3)) ) ./ lx;
%     Hy ( 1:IHy-1, 1:JHy-1, n1 ) = (1/u0) * ~cmaskHy (1:IHy-1, 1:JHy-1) .* (uyxHx.*Bx ( :, :, n1 ) + uyyHy (1:IHy-1, 1:JHy-1) .* By (1:IHy-1, 1:JHy-1, n1 )) + cmaskHy (1:IHy-1, 1:JHy-1) .* ( ay(1:IHy-1, 1:JHy-1).*By ( 1:IHy-1, 1:JHy-1, n1 ) + by(1:IHy-1, 1:JHy-1).*By ( 1:IHy-1, 1:JHy-1, n0) + cy(1:IHy-1, 1:JHy-1).*By ( 1:IHy-1, 1:JHy-1, 3) + dy(1:IHy-1, 1:JHy-1).*BxAve (1:IHy-1, 1:JHy-1, n1) + ey(1:IHy-1, 1:JHy-1).*BxAve (1:IHy-1, 1:JHy-1, n0) + fy(1:IHy-1, 1:JHy-1).*BxAve (1:IHy-1, 1:JHy-1, 3) - (gy(1:IHy-1, 1:JHy-1).*Hy(1:IHy-1, 1:JHy-1,n0) + hy(1:IHy-1, 1:JHy-1).*Hy(1:IHy-1, 1:JHy-1,3)) ) ./ ly (1:IHy-1, 1:JHy-1);
%     Hx ( :, :, n1 ) = (1/u0) * (uxxHx .* Bx ( :, :, n1 ) + uxyHy (1:IHy-1, 1:JHy-1 ) .* By (1:IHy-1, 1:JHy-1, n1 ));
%     ****** Major ****** Hx ( :, :, n1 ) = (1/u0) * (uxxHx .* Bx ( :, :, n1 ) + uxyHy (1:IHy-1, 1:JHy-1 ) .* By (1:IHy-1, 1:JHy-1, n1 ));
    %Hx ( :, :, n1 ) = smaskHx ( :, : ) .* Hx ( :, :, n1 );
    
    
%     Hy (1:IHy-1, 1:JHy-1, n1) = (1/u0) * (uyxHx.*Bx ( :, :, n1 ) + uyyHy (1:IHy-1, 1:JHy-1) .* By (1:IHy-1, 1:JHy-1, n1 ));
    %Hy ( :, :, n1 ) = (1/u0) * By ( :, :, n1 );
%     ****** Major ******    Hy (1:IHy-1, 1:JHy-1, n1) = (1/u0) * (uyxHx.*Bx ( :, :, n1 ) + uyyHy (1:IHy-1, 1:JHy-1) .* By (1:IHy-1, 1:JHy-1, n1 ));
    
    % ABC for Hy;
%     Hy (1:IHy-1, JHy, n1) = Hy (1:IHy-1, JHy, n0);
    %Hy ( :, :, n1 ) = smaskHy ( :, : ) .* Hy ( :, :, n1 );
    
    Dz ( :, 2:JEz-1, n1 ) = (  Dz ( :, 2:JEz-1, n0 ) ) + ( (DT/delta) * ( Hy ( 2:JEz+1, 2:JEz-1, n1 ) - Hy ( 1:JEz, 2:JEz-1, n1 ) + Hx ( :, 1:JEz-2,n1 ) - Hx ( :, 2:JEz-1, n1 ) ));
%     Dz (:, :, n1) = Dz (:, :, n1) .* smask;
    
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

%     Ez ( :, :, n1 ) =  ( (1/(e0*(DT^2)))*Dz ( :, :, n1 ) - (2/(e0*(DT^2)))*Dz ( :, :, n0) + (1/(e0*(DT^2)))*Dz( :, :, 3) + AEz.*(2/(DT^2)-wpsquaredEz/2).*Ez(:, :, n0) - AEz.*(1/(DT^2)+wpsquaredEz/4).*Ez (:, :, 3) ) ./ (AEz.*( 1/(DT^2) + wpsquaredEz/4));
% 	Ez ( :, :, n1 ) = (1/e0) * Dz ( :, :, n1 ) ./ (ezzEz);
    Ez ( :, :, n1 ) =  ( (1/(e0*(DT^2)))*Dz ( :, :, n1 ) - (2/(e0*(DT^2)))*Dz ( :, :, n0) + (1/(e0*(DT^2)))*Dz( :, :, 3) + AEz.*(2/(DT^2)-wpsquaredEz/2).*Ez(:, :, n0) - AEz.*(1/(DT^2)+wpsquaredEz/4).*Ez (:, :, 3) ) ./ (AEz.*( 1/(DT^2) + wpsquaredEz/4));
% 	Ez ( :, :, n1 ) =  ( ~cmaskEz .* (1/e0) * (ezzEz) .* Dz ( :, :, n1 ) ) + ( cmaskEz .* ( (1/(e0*DT^2))*Dz ( :, :, n1 ) - (2/(e0*(DT^2)))*Dz ( :, :, n0) + (1/(e0*(DT^2)))*Dz( :, :, 3) + A*(2/(DT^2)-wpsquaredEz/2).*Ez(:, :, n0) - A*(1/(DT^2)+wpsquaredEz/4).*Ez (:, :, 3) ) ./ A*( 1/(DT^2) + wpsquaredEz/4) );
    
%     %     ****** Major ******   Ez ( :, :, n1 ) = (1/e0) * Dz ( :, :, n1 ) ./ (ezzEz);
%     Ez ( :, :, n1 ) = (1/e0) * Dz ( :, :, n1 ) ./ (ezzEz);
    
%     ***** Drude ***** Ez ( :, :, n1 ) =  ( ~cmaskEz .* (1/e0) * (ezzEz) .* Dz ( :, :, n1 ) ) + ( cmaskEz .* ( (1/(e0*DT^2))*Dz ( :, :, n1 ) - (2/(e0*(DT^2)))*Dz ( :, :, n0) + (1/(e0*(DT^2)))*Dz( :, :, 3) + A*(2/(DT^2)-wpsquaredEz/2).*Ez(:, :, n0) - A*(1/(DT^2)+wpsquaredEz/4).*Ez (:, :, 3) ) ./ A*( 1/(DT^2) + wpsquaredEz/4) );
    % Ez calculation outside cylinder.
%     Ez ( :, :, n1 ) = ~cmaskEz .* (1/e0) * (ezzEz) .* Dz ( :, :, n1 );
    % Ez calculation inside cylinder using Drude dispersion model.
    % Ez ( :, :, n1 ) = cmaskEz .* ( (1/(e0*DT^2))*Dz ( :, :, n1 ) - (2/(e0*(DT^2)))*Dz ( :, :, n0) + (1/(e0*(DT^2)))*Dz( :, :, 3) + A*(2/(DT^2)-wpsquaredEz/2).*Ez(:, :, n0) - A*(1/(DT^2)+wpsquaredEz/4).*Ez (:, :, 3) ) ./ A*( 1/(DT^2) + wpsquaredEz/4);


%     if ( n < NHW )
    Ez ( :, Js, n1 ) = Ez ( :, Js, n1 ) + 1 * sin ( TwoPIFDeltaT * n );
    Dz ( :, Js, n1 ) = e0 * Ez ( :, Js, n1 );
%     end

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
    
    figure (7)
    surf ( EzSnapshots (:, :, i) );
    view (0, 90)
    zlim ( [-10 10] )
    
end
fprintf ( 1, 'Simulation completed! \n' );