clc
clear all

% Simulation parameters.
SIZE = 1*2000; % No. of spatial steps
PotentialLeft = round(SIZE/2)-12; % Location of left end of potential barrier.
PotentialRight = round(SIZE/2)+12; % Location of right end of potential barrier.
V0 = 1200; % 600 eV barrier.
MaxTime = 12*SIZE; % No. of time steps
PulseWidth = round(SIZE/8); % Controls width of Gaussian Pulse
td = PulseWidth; % Temporal delay in pulse.
source = 1; % Location of source
SaveSnapshots = 1; % 0. No, 1. Yes.
SnapshotInterval = 32; % Amount of time delay between snaps.

% Choice of source.
% 1. Gaussian 2. Sine wave 3. Ricker wavelet
SourceChoice = 2;

% Constants.
h = 6.626068e-34;   % Planck's constant.
m = 5*9.10938188e-66; % Mass of electron.
e = 1.60217656e-19; % Charge on electron.
pi = 3.141592654;
c = 3e8;

dx = 0.01e-10;
dtc = h/(h^2/(m*dx^2)+V0/(e*2));
dt = dtc/2;
c1 = h*dt/(2*m*dx^2)
c2 = dt/h
E = 500/e;
Sc = c * dt/dx

l = PulseWidth*dx;
f = E/h
fmax = 1/(2*dt)
w = 2*pi*f;
k0 = w/c; % Free space wave number.
% Ricker wavelet parameters.
if SourceChoice == 3
    fp = f; % Peak frequency
    dr = PulseWidth*dt*2; % Delay
end

% Initialisation.
psiR = zeros(SIZE, 3); % real-component of wave function
psiI = zeros(SIZE, 3); % imaginary-component of wave function
V = zeros(SIZE); % Potential
V(PotentialLeft:PotentialRight) = V0/e;

if SaveSnapshots == 1
    psiISnapshots = zeros(SIZE', round(MaxTime/SnapshotInterval)); % Data for plotting.
    psiRSnapshots = zeros(SIZE', round(MaxTime/SnapshotInterval)); % Data for plotting.
    frame = 1;
end

n1 = 1;
n2 = 2;
linecount = 0;
% Pre-defined source.
% for x=1:1000
%    psiI(x,n1) = sin(x*dx*E/c-0*2*pi*f*dt);
%    psiR(x,n1) = cos(x*dx*E/c-0*2*pi*f*dt); 
% end
% figure(5)
% plot(psiI(:,n1))
% Outer loop for time-stepping.
tic
% Test loop for incident field in free space.
for q = 0:MaxTime
    
    % Progress indicator.
    if mod(q,SnapshotInterval) == 0
        fprintf(1, repmat('\b',1,linecount));
        linecount = fprintf(1, '%g %%', (q*100)/MaxTime );
    end
    
    l = 2:SIZE-1;
    % Calculation of psiI using update difference equation. This is time step q.
    psiI(l,n2) = c1*(psiR(l+1,n1) - 2*psiR(l,n1) + psiR(l-1,n1)) - c2*(V(l)'.*psiR(l,n1)) + psiI(l,n1);
    psiR(l,n2) = -1*c1*(psiI(l+1,n2) - 2*psiI(l,n2) + psiI(l-1,n2)) + c2*(V(l)'.*psiI(l,n2)) + psiR(l,n1);
    
    % Source.
    if SourceChoice == 1
    psiR(source,n2) = psiR(source,n2) + 0.5;%exp( -1*((q-td)/(PulseWidth/4))^2 );
    elseif SourceChoice == 2
        if q < MaxTime/5
            psiI(source,n2) = sin(2*pi*f*(q)*dt) * 1;
            psiR(source,n2) = cos(2*pi*f*(q)*dt) * 1;
        else
            psiI(source,n2) = 0;
            psiR(source,n2) = 0;
        end
    elseif SourceChoice == 3
    %Ex(source,n2) = Ex(source,n2) + (1-2*(pi*fp*(q*dt-dr))^2)*exp(-1*(pi*fp*(q*dt-dr))^2) * Sc;
    end
    
    if (SaveSnapshots == 1 && mod(q,SnapshotInterval) == 0)
        psiISnapshots(:,frame) = psiI(:,n2);
        psiRSnapshots(:,frame) = psiR(:,n2);
        frame=frame+1;
    end
    
    temp = n1;
    n1 = n2;
    n2 = temp;
end
fprintf(1, repmat('\b',1,linecount));
fprintf ( 1, 'Simulation complete! \n');

if SaveSnapshots == 1
    % Simulation animation.
    for i=1:frame-1
        figure (6)
        % Scatterer boundaries.
        hold off
        plot([PotentialLeft PotentialLeft], [-1 1], 'Color', 'r');
        hold on
        plot([PotentialRight PotentialRight], [-1 1], 'Color', 'r');
        plot(psiRSnapshots(:,i), 'LineWidth', 2.0, 'Color', 'b');
        plot(psiISnapshots(:,i), 'LineWidth', 2.0, 'Color', 'g');
        plot(sqrt(psiISnapshots(:,i).^2+psiRSnapshots(:,i).^2), 'LineWidth', 2.0, 'Color', 'm');
        set(gca, 'FontSize', 10, 'FontWeight', 'b')
        axis([0 SIZE -2 2])
        %xlim([0 SIZE])
        title('Time Domain Simulation', 'FontSize', 12, 'FontWeight', 'b')
        xlabel('Spatial step (k)', 'FontSize', 11, 'FontWeight', 'b')
        ylabel('Electric field (Ex)', 'FontSize', 11, 'FontWeight', 'b')
        grid on
    end
end