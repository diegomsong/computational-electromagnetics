% Simulation parameters.
SIZE = 512; % No. of spatial steps
MaxTime = SIZE*4; % No. of time steps
PulseWidth = SIZE/4; % Controls width of Gaussian Pulse
td = PulseWidth; % Temporal delay in pulse.
imp0 = 377.0; % Impedence of free space
source = 10; % Location of source

% Constants.
pi = 3.141592654;
e0 = (1e-9)/(36*pi);
u0 = (1e-7)*4*pi;
c = 1/sqrt(e0*u0);

dt = 1e-12;
dz = 3e-4;
Sc = c * dt/dz

% Material
Partition = SIZE/4;
er1 = 1;
er2 = 1;
ur1 = 1;
ur2 = 1;
sig1 = 0;
sig2 = 0.1;
sigm1 = 0;
sigm2 = 0;

c1 = 1/sqrt(er1*e0*ur1*u0);
c2 = 1/sqrt(er2*e0*ur2*u0);
Sc1 = c1 * dt/dz
Sc2 = c2 * dt/dz

eps = ones(SIZE, 1);
mu = ones(SIZE, 1);
sig = ones(SIZE, 1);
sigm = ones(SIZE, 1);

eps(1:Partition-1) = e0*er1;
eps(Partition) = e0*(er1+er2)/2;
eps(Partition+1:SIZE) = e0*er2;
mu(1:Partition-1) = u0*ur1;
mu(Partition:SIZE) = u0*ur2;

sig(1:Partition-1) = sig1;
sig(Partition) = (sig1+sig2)/2;
sig(Partition+1:SIZE) = sig2;
sigm(1:Partition-1) = sigm1;
sigm(Partition:SIZE) = sigm2;

% Initialization.
Ex = zeros(SIZE, MaxTime); % x-component of E-field
Hy = zeros(SIZE, MaxTime); % y-component of H-field
PLOT1(1) = 0; % Data for plotting.

% Outer loop for time-stepping.
tic
for q = 2:MaxTime
    % Calculation of Hy using update difference equation for Hy. This is time step q.
    Hy(1:SIZE-1,q) = (1-sigm(1:SIZE-1)*dt./(2*mu(1:SIZE-1)))./(1+sigm(1:SIZE-1)*dt./(2*mu(1:SIZE-1))).*Hy(1:SIZE-1,q-1) + ((Ex(1:SIZE-1,q-1) - Ex(2:SIZE,q-1)) .* ((dt./(mu(1:SIZE-1)*dz))./(1+sigm(1:SIZE-1)*dt./(2*mu(1:SIZE-1)))));
    % ABC for H at SIZE.
    Hy(SIZE,q) = Hy(SIZE-1,q-1) + (Sc-1)/(Sc+1)*(Hy(SIZE-1,q) - Hy(SIZE,q-1));
    
    % Calculation of Ex using updated difference equation for Ex. This is time step q+1/2.
    Ex(2:SIZE,q) = (1-sig(2:SIZE)*dt./(2*eps(2:SIZE)))./(1+sig(2:SIZE)*dt./(2*eps(2:SIZE))).*Ex(2:SIZE, q-1) + ((dt./(eps(2:SIZE)*dz))./(1+sig(2:SIZE)*dt./(2*eps(2:SIZE)))).*(Hy(1:SIZE-1, q) - Hy(2:SIZE, q));
    % ABC for E at 1.
    Ex(1,q) = Ex(2,q-1) + (Sc-1)/(Sc+1)*(Ex(2,q) - Ex(2,q-1));
    
    % Lossless code for reference.
    % Calculation of Hy using update difference equation for Hy. This is time step q.
    %Hy(1:SIZE-1,q) = Hy(1:SIZE-1,q-1) + ( ( Ex(1:SIZE-1,q-1) - Ex(2:SIZE,q-1) ) * dt/(u0*dz) );
    % ABC for H at SIZE.
    %Hy(SIZE,q) = Hy(SIZE-1,q-1) + (Sc-1)/(Sc+1)*(Hy(SIZE-1,q) - Hy(SIZE,q-1) );
    
    % Calculation of Ex using updated difference equation for Ex. This is time step q+1/2.
    %Ex(2:SIZE,q) = Ex(2:SIZE, q-1) + ( dt/(e0*dz)*(Hy(1:SIZE-1, q) - Hy(2:SIZE, q)) );
    % ABC for E at 1.
    %Ex(1,q) = Ex(2,q-1) + (Sc-1)/(Sc+1)*(Ex(2,q) - Ex(2,q-1));
    
    % Activating a plane-wave source.
    Ex(source,q) = Ex(source,q) + exp( -1*((q-td)/(PulseWidth/4))^2 ) * Sc;
end
toc
% Simulation animation.
for i=1:MaxTime
    figure (2)
    hold off
    plot([Partition-1 Partition-1], [-1 1], 'Color', 'r');
    hold on
    plot ( Ex(:,i) )
    axis([0 SIZE -0.6 0.6])
    xlabel('Spatial step (k)')
    ylabel('Electric field (Ex)')
end