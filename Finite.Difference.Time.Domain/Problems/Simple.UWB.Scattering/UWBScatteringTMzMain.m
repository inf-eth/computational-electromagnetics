% 2D FDTD TMz simulation employing a PML. The solution region is bounded on two sides with a PML terminated by PEC. Useful for single plane wave simulations.

clc
clear all

% ============== Simulation related parameters ================
[ISize JSize XCenter YCenter delta ra rb DTp PMLw dtscalar SkinIW TissueIW TissueIIW PulseWidth TumourX TumourY TumourRadius] = Parameters;

IHx = ISize;
JHx = JSize+2*PMLw-1;
IHy = ISize+1;
JHy = JSize+2*PMLw;
IEz = ISize;
JEz = JSize+2*PMLw;

% Time indices for field calculation.
n0 = 1;
n1 = 2;
NNMax = 750;                   % Maximum time.
TimeResolutionFactor = 1;      % E field snapshots will be saved every x frames where x is time resolution factor.
ResolutionFactor = 2;          % Resolution of plotted field is divided by this factor.
Js = 100;                       % J-position of the plane wave front.
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

% s*DT/2*er;
Sc = zeros (IEz, JEz);
% sm*DT/2*ur;
ScmHx = zeros (IHx, JHx);
ScmHy = zeros (IHy, JHy);

% PML conductance arrays.
sex = zeros (IEz, JEz);     % sigma ex
sey = zeros (IEz, JEz);     % sigma ey
smx = zeros (IHy, JHy);     % sigma mx
smy = zeros (IHx, JHx);     % sigma my

Scsx = zeros (IEz, JEz);       % s*DT/2*er;
Scsy = zeros (IEz, JEz);       %
ScmsmxHy = zeros (IHy, JHy);    % sm*DT/2*ur;
ScmsmyHx = zeros (IHx, JHx);    %

% PEC mask array.
PECmaskEz = zeros ( IEz, JEz );
% ---------------------------------------------------

Bx = zeros ( IHx, JHx, 3 );
By = zeros ( IHy, JHy, 3 );

Hx = zeros ( IHx, JHx, 3 );
Hy = zeros ( IHy, JHy, 3 );

Dz = zeros ( IEz, JEz, 3 );
Ez = zeros ( IEz, JEz, 3 );

Dzx = zeros ( IEz, JEz, 3 );
Dzy = zeros ( IEz, JEz, 3 );
EzSnapshots = zeros ( IEz/ResolutionFactor, JEz/ResolutionFactor, NNMax/TimeResolutionFactor ); % E field Snapshot storage.

% =============================================================

% ############ Initialization #############
fprintf ( 1, 'Initializing...' );
fprintf ( 1, '\nInitializing parametric arrays...' );
% Initializing parametric arrays.
for i=1:IEz
    fprintf ( 1, '%g %% \n', ((i-1)*100)/IHx );
    for j=1:JEz      
          sEz (i, j) = sc (i, j);
          erEz (i, j) = er (i, j);
    end
end
% Initializing PML conductance arrays.
dpml = PMLw;
mpml = 250;          % Typical = 80;
semax = 2.6e7;      % Typical = 3.7e6;
for i=1:PMLw+1
    
%     sey(:, i+1) = semax*( (PMLw-i)/dpml )^mpml;
%     smy(:, i) = (1/1)*(u0/e0)*semax*( (PMLw-i+0.5)/dpml )^mpml;
    
    sey(:, i+1) = 1.7e10;
    smy(:, i) = 1.7e10;
    if i<PMLw+1
%         sey(:, JEz-PMLw+i) = semax*( (i)/dpml )^mpml;
%         smy(:, JHx-PMLw+i) = (1/1)*(u0/e0)*semax*( (i-0.5)/dpml )^mpml;  
        sey(:, JEz-PMLw+i) = 1.7e10;
        smy(:, JHx-PMLw+i) = 1.7e10;
    end
    
%     erEz(:, i+1) = 10;
%     erEz(:, JEz-PMLw+i-1) = 10;
end
% figure (1)
% mesh ( smy )
% title ( 'smy' )
% view (4, 4)
% figure (2)
% mesh ( sey )
% title ( 'sey' )
% view (4, 4)
figure (3)
surf ( erEz )
title ( 'erEz' )
view (4, 4)
% figure (4)
% mesh ( wpmsquaredHx )
% title ( 'wpm squared' )
% view (4, 4)
% figure (5)
% mesh ( wpsquaredEz )
% title ( 'wp squared' )
% view (4, 4)

% Normal space.
Sc = (DT*sEz)./(2*erEz);
ScmHx = (DT*smHx)./(2*urHx);
ScmHy = (DT*smHy)./(2*urHy);
% PML space.
Scsx = (DT*sex)./(2*erEz);
Scsy = (DT*sey)./(2*erEz);
ScmsmxHy = (DT*smx)./(2*urHy);
ScmsmyHx = (DT*smy)./(2*urHx);

fprintf ( 1, 'Initialization done.\n' );
% ############ Initialization Complete ##############

% ########### 2. Simulation #############
fprintf ( 1, 'Simulation started... \n' );
for n=0:NNMax-2
    fprintf ( 1, '%g %% \n', (n*100)/NNMax );
    
    % Copying n1 fields as past fields. Will be used in second order time derivatives.
    Bx(:, :, 3) = Bx(:, :, n1);
    By(:, :, 3) = By(:, :, n1);
    Hx(:, :, 3) = Hx(:, :, n1);
    Hy(:, :, 3) = Hy(:, :, n1);
    Dz(:, :, 3) = Dz(:, :, n1);
    Ez(:, :, 3) = Ez(:, :, n1);
    
    % Bx in normal space.
    Bx(:, (1+PMLw):(JHx-PMLw), n1) = ((1-ScmHx(:, (1+PMLw):(JHx-PMLw)))./(1+ScmHx(:, (1+PMLw):(JHx-PMLw)))) .* Bx(:, (1+PMLw):(JHx-PMLw), n0) + ( (DT/(delta))./(1+ScmHx(:, (1+PMLw):(JHx-PMLw))) .* ( Ez(:, (1+PMLw):(JEz-1-PMLw), n0) - Ez(:, (2+PMLw):(JEz-PMLw), n0) ));
    % Bx in lower PML region.
    Bx(:, 1:PMLw, n1) = ((1-ScmsmyHx(:, 1:PMLw))./(1+ScmsmyHx(:, 1:PMLw))) .* Bx(:, 1:PMLw, n0) + ( (DT/(delta))./(1+ScmsmyHx(:, 1:PMLw)) .* ( Ez(:, 1:PMLw, n0) - Ez(:, 2:PMLw+1, n0) ));
    % Bx in upper PML region.
    Bx(:, JHx-PMLw+1:JHx, n1) = ((1-ScmsmyHx(:, JHx-PMLw+1:JHx))./(1+ScmsmyHx(:, JHx-PMLw+1:JHx))) .* Bx(:, JHx-PMLw+1:JHx, n0) + ( (DT/(delta))./(1+ScmsmyHx(:, JHx-PMLw+1:JHx)) .* ( Ez(:, JEz-PMLw:JEz-1, n0) - Ez(:, JEz-PMLw+1:JEz, n0) ));
    
    % By in normal space.
    By(2:IHy-1, (2+PMLw):(JHy-PMLw-1), n1) = ((1-ScmHy(2:IHy-1, (2+PMLw):(JHy-PMLw-1)))./(1+ScmHy(2:IHy-1, (2+PMLw):(JHy-PMLw-1)))) .* By(2:IHy-1, (2+PMLw):(JHy-PMLw-1), n0) + ( (DT/(delta))./(1+ScmHy(2:IHy-1, (2+PMLw):(JHy-PMLw-1))) .* ( Ez (2:IEz, (2+PMLw):(JEz-PMLw-1), n0) - Ez (1:IEz-1, (2+PMLw):(JEz-PMLw-1), n0) ));
    % By in lower PML region.
    By(2:IHy-1, 2:PMLw+1, n1) = ((1-ScmsmxHy(2:IHy-1, 2:PMLw+1))./(1+ScmsmxHy(2:IHy-1, 2:PMLw+1))) .* By(2:IHy-1, 2:PMLw+1, n0) + ( (DT/(delta))./(1+ScmsmxHy(2:IHy-1, 2:PMLw+1)) .* ( Ez (2:IEz, 2:PMLw+1, n0) - Ez (1:IEz-1, 2:PMLw+1, n0) ));
    % By in upper PML region.
    By(2:IHy-1, JHy-PMLw:JHy-1, n1) = ((1-ScmsmxHy(2:IHy-1, JHy-PMLw:JHy-1))./(1+ScmsmxHy(2:IHy-1, JHy-PMLw:JHy-1))) .* By(2:IHy-1, JHy-PMLw:JHy-1, n0) + ( (DT/(delta))./(1+ScmsmxHy(2:IHy-1, JHy-PMLw:JHy-1)) .* ( Ez(2:IEz, JEz-PMLw:JEz-1, n0) - Ez(1:IEz-1, JEz-PMLw:JEz-1, n0) ));

    % Magnetic fields.
    Hx (:, :, n1) = (1/u0) * Bx (:, :, n1) ./ (urHx);
    Hy (:, :, n1) = (1/u0) * By (:, :, n1) ./ (urHy);
   
    % Boundary conditions on By. Soft grid truncation.
%     By ( 1, 2:JHy-1, n1 ) = (1/3) * ( By ( 2, 1:JHy-2, n0 ) + By ( 2, 2:JHy-1, n0 ) + By ( 2, 3:JHy, n0 ) );
%     By ( IHy, 2:JHy-1, n1 ) = (1/3) * ( By ( IHy-1, 1:JHy-2, n0 ) + By ( IHy-1, 2:JHy-1, n0 ) + By ( IHy-1, 3:JHy, n0 ) );
%     By ( 1, 1, n1 ) = (1/2) * ( By ( 2, 1, n0 ) + By ( 2, 2, n0 ) );
%     By ( 1, JHy, n1 ) = (1/2) * ( By ( 2, JHy, n0 ) + By ( 2, JHy-1, n0 ) );
%     By ( IHy, 1, n1 ) = (1/2) * ( By ( IHy-1, 1, n0 ) + By( IHy-1, 2, n0 ) );
%     By ( IHy, JHy, n1 ) = (1/2) * ( By ( IHy-1, JHy, n0 ) + By ( IHy-1, JHy-1, n0 ) );

    % Dz in normal space.
    Dz(:, (2+PMLw):(JEz-1-PMLw), n1) = ((1-Sc(:, (2+PMLw):(JEz-1-PMLw)))./(1+Sc(:, (2+PMLw):(JEz-1-PMLw)))) .* Dz(:, (2+PMLw):(JEz-1-PMLw), n0) + ( ((DT/delta)./(1+Sc(:, (2+PMLw):(JEz-1-PMLw)))) .* ( Hy(2:IHy, (2+PMLw):(JHy-1-PMLw), n1) - Hy(1:IHy-1, (2+PMLw):(JHy-1-PMLw), n1) - Hx(:, (2+PMLw):(JHx-PMLw), n1) + Hx(:, (1+PMLw):(JHx-1-PMLw), n1) ));
    
    % PML space.
%     Dz(:, 2:PMLw+1, n1) = ((1-Scsy(:, 2:PMLw+1))./(1+Scsy(:, 2:PMLw+1))) .* Dz(:, 2:PMLw+1, n0) + ( ((DT/delta)./(1+Scsy(:, 2:PMLw+1))) .* ( Hy(2:IHy, 2:PMLw+1, n1) - Hy(1:IHy-1, 2:PMLw+1, n1) - Hx(:, 2:PMLw+1, n1) + Hx(:, 1:PMLw, n1) ));
%     Dz(:, JEz-PMLw:JEz-1, n1) = ((1-Scsy(:, JEz-PMLw:JEz-1))./(1+Scsy(:, JEz-PMLw:JEz-1))) .* Dz(:, JEz-PMLw:JEz-1, n0) + ( ((DT/delta)./(1+Scsy(:, JEz-PMLw:JEz-1))) .* ( Hy(2:IHy, JHy-PMLw:JHy-1, n1) - Hy(1:IHy-1, JHy-PMLw:JHy-1, n1) - Hx(:, JHx-PMLw+1:JHx, n1) + Hx(:, JHx-PMLw-0:JHx-1, n1) ));
    
%     Dz(:, 2:PMLw+1, n1) = ((1-Scsy(:, 2:PMLw+1))./(1+Scsy(:, 2:PMLw+1))) .* Dz(:, 2:PMLw+1, n0) + ( ((DT/delta)./(1+Scsy(:, 2:PMLw+1))) .* ( - Hx(:, 2:PMLw+1, n1) + Hx(:, 1:PMLw, n1) ));
%     Dz(:, JEz-PMLw:JEz-1, n1) = ((1-Scsy(:, JEz-PMLw:JEz-1))./(1+Scsy(:, JEz-PMLw:JEz-1))) .* Dz(:, JEz-PMLw:JEz-1, n0) + ( ((DT/delta)./(1+Scsy(:, JEz-PMLw:JEz-1))) .* ( - Hx(:, JHx-PMLw+1:JHx, n1) + Hx(:, JHx-PMLw-0:JHx-1, n1) ));

    % PML space split Dz component calculations. Dzx and Dzy in lower PML region followed by Dzx and Dzy in upper PML regions.
%     Dzx(:, 2:PMLw+1, n1) = ((1-Scsx(:, 2:PMLw+1))./(1+Scsx(:, 2:PMLw+1))) .* Dzx(:, 2:PMLw+1, n0) + ( ((DT/delta)./(1+Scsx(:, 2:PMLw+1))) .* ( Hy(2:IHy, 2:PMLw+1, n1) - Hy(1:IHy-1, 2:PMLw+1, n1) ));
%     Dzy(:, 2:PMLw+1, n1) = ((1-Scsy(:, 2:PMLw+1))./(1+Scsy(:, 2:PMLw+1))) .* Dzy(:, 2:PMLw+1, n0) + ( ((DT/delta)./(1+Scsy(:, 2:PMLw+1))) .* ( - Hx(:, 2:PMLw+1, n1) + Hx(:, 1:PMLw+0, n1) ));
%     Dzx(:, JEz-PMLw:JEz-1, n1) = ((1-Scsx(:, JEz-PMLw:JEz-1))./(1+Scsx(:, JEz-PMLw:JEz-1))) .* Dzx(:, JEz-PMLw:JEz-1, n0) + ( ((DT/delta)./(1+Scsx(:, JEz-PMLw:JEz-1))) .* ( Hy(2:IHy, JHy-PMLw:JHy-1, n1) - Hy(1:IHy-1, JHy-PMLw:JHy-1, n1) ));
%     Dzy(:, JEz-PMLw:JEz-1, n1) = ((1-Scsy(:, JEz-PMLw:JEz-1))./(1+Scsy(:, JEz-PMLw:JEz-1))) .* Dzy(:, JEz-PMLw:JEz-1, n0) + ( ((DT/delta)./(1+Scsy(:, JEz-PMLw:JEz-1))) .* ( - Hx(:, JHx-PMLw+1:JHx, n1) + Hx(:, JHx-PMLw-0:JHx-1, n1) ));

    % Lower PML.
    Dzx(:, 2:PMLw+1, n1) = ((1-Scsx(:, 2:PMLw+1))./(1+Scsx(:, 2:PMLw+1))) .* Dzx(:, 2:PMLw+1, n0) + ( ((DT/delta)./(1+Scsx(:, 2:PMLw+1))) .* ( Hy(2:IHy, 2:PMLw+1, n1) - Hy(1:IHy-1, 2:PMLw+1, n1) ));
    Dzy(:, 2:PMLw+1, n1) = ((1-Scsy(:, 2:PMLw+1))./(1+Scsy(:, 2:PMLw+1))) .* Dzy(:, 2:PMLw+1, n0) + ( ((DT/delta)./(1+Scsy(:, 2:PMLw+1))) .* ( - Hx(:, 2:PMLw+1, n1) + Hx(:, 1:PMLw, n1) ));
    
    % Upper PML.
    Dzx(:, JEz-PMLw:JEz-1, n1) = ((1-Scsx(:, JEz-PMLw:JEz-1))./(1+Scsx(:, JEz-PMLw:JEz-1))) .* Dzx(:, JEz-PMLw:JEz-1, n0) + ( ((DT/delta)./(1+Scsx(:, JEz-PMLw:JEz-1))) .* ( Hy(2:IHy, JHy-PMLw:JHy-1, n1) - Hy(1:IHy-1, JHy-PMLw:JHy-1, n1) ));
    Dzy(:, JEz-PMLw:JEz-1, n1) = ((1-Scsy(:, JEz-PMLw:JEz-1))./(1+Scsy(:, JEz-PMLw:JEz-1))) .* Dzy(:, JEz-PMLw:JEz-1, n0) + ( ((DT/delta)./(1+Scsy(:, JEz-PMLw:JEz-1))) .* ( - Hx(:, JHx-PMLw+1:JHx, n1) + Hx(:, JHx-PMLw-0:JHx-1, n1) ));
    % Boundary conditions on Dz. Soft grid truncation.
%     Dz ( 2:IEz-1, 1, n1 ) = (1/3) * ( Dz ( 1:IEz-2, 2, n0 ) + Dz ( 2:IEz-1, 2, n0 ) + Dz ( 3:IEz, 2, n0 ) );
%     Dz ( 2:IEz-1, JEz, n1 ) = (1/3) * ( Dz ( 1:IEz-2, JEz-1, n0 ) + Dz ( 2:IEz-1, JEz-1, n0 ) + Dz ( 3:IEz, JEz-1, n0 ) );
%     Dz ( 1, 1, n1 ) = (1/2) * ( Dz ( 1, 2, n0 ) + Dz ( 2, 2, n0 ) );
%     Dz ( IEz, 1, n1 ) = (1/2) * ( Dz ( IEz, 2, n0 ) + Dz ( IEz-1, 2, n0 ) );
%     Dz ( 1, JEz, n1 ) = (1/2) * ( Dz ( 1, JEz-1, n0 ) + Dz ( 2, JEz-1, n0 ) );
%     Dz ( IEz, JEz, n1 ) = (1/2) * ( Dz ( IEz, JEz-1, n0 ) + Dz ( IEz-1, JEz-1, n0 ) );

    Dz(:, 2:PMLw+1, n1) = Dzx(:, 2:PMLw+1, n1) + Dzy(:, 2:PMLw+1, n1); 
    Dz(:, JEz-PMLw:JEz-1, n1) = Dzx(:, JEz-PMLw:JEz-1, n1) + Dzy(:, JEz-PMLw:JEz-1, n1);
%     Dz(:, :, n1) = Dz(:, :, n1)+Dzx(:, :, n1)+Dzy(:, :, n1);
    
    Ez (:, :, n1) = (1/e0) * Dz (:, :, n1) ./ (erEz);
    
    Ez (:, Js, n1) = Ez (:, Js, n1) + 1 * exp (- 1 * ((n-PulseWidth)/15)^2);% / dtscalar;
    Dz (:, Js, n1) = e0 * Ez (:, Js, n1);
    
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
    zlim ( [-1 1] )
    caxis([0 1])
    xlabel ('y-axis')
    ylabel ('x-axis')
    %     colorbar
    
    figure (7)
    surf ( EzSnapshots (:, :, i) );
    view (0, 90)
    zlim ( [-10 10] )
    caxis([0 1])
    xlabel ('y-axis')
    ylabel ('x-axis')
%     colorbar
    
end

figure (8)
% plot (reshape(EzSnapshots(IEz/2, 40, :), 500, 1))
plot ( squeeze(EzSnapshots(IEz/2, 40, :)))
axis ([0 length(squeeze(EzSnapshots(IEz/2, 40, :))) -1 1])

figure (9)
surf (squeeze (EzSnapshots (:, 40, 1:5:NNMax)));
view (0, 90)
zlim ( [-1 1] )
caxis([0 1])
title ('Scattering Image at j=40')
xlabel ('y-axis')
ylabel ('x-axis')

figure (10)
surf (squeeze (EzSnapshots (:, 150, 1:5:NNMax)));
view (0, 90)
zlim ( [-1 1] )
caxis([0 1])
title ('Scattering Image at j= 150')
xlabel ('y-axis')
ylabel ('x-axis')

fprintf ( 1, 'Simulation completed! \n' );