% Simulation parameters.
Scalar = 8;
SIZE = 128*Scalar; % No. of spatial steps
MaxTime = SIZE*2*Scalar; % No. of time steps
PulseWidth = SIZE/Scalar; % Controls width of Gaussian Pulse
td = PulseWidth; % Temporal delay in pulse.
imp0 = 377.0; % Impedence of free space
source = 3; % Location of source
Ra = 150e-6;

% Constants.
pi = 3.141592654;
e0 = (1e-9)/(36*pi);
u0 = (1e-7)*4*pi;
c = 1/sqrt(e0*u0);

dt = .125e-14;
dz = 0.75e-6;
Sc = c * dt/dz

% Choice of source.
% 1. Gaussian 2. Sine wave 3. Ricker wavelet
SourceChoice = 2;

l = PulseWidth*dz;
f = 0.8727e12%c/(1*l)
fmax = 1/(2*dt)
w = 2*pi*f;
k0 = w/c; % Free space wave number.
% Ricker wavelet parameters.
if SourceChoice == 2
    fp = f; % Peak frequency
    dr = PulseWidth*dt*2; % Delay
end

% Material
Partition = floor(Ra/dz)%SIZE/4;
er1 = 1;
er2 = 1;
ur1 = 1;
ur2 = 1;
sig1 = 0;
sig2 = 300;
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

% Amplitude and phase calculations.
t1 = floor(MaxTime/Scalar)
t2 = floor(t1+1/f/4/dt)
ExAbs = zeros(SIZE, 1);
ExPhase = zeros(SIZE, 1);
HyAbs = zeros(SIZE, 1);
HyPhase = zeros(SIZE, 1);
ExPhasor = zeros(SIZE, 1);
HyPhasor = zeros(SIZE, 1);

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
    Ex(1,q) = Ex(2,q-1) + (Sc-1)/(Sc+1)*(Ex(2,q) - Ex(1,q-1));
    
    % Lossless code for reference.
    % Calculation of Hy using update difference equation for Hy. This is time step q.
    %Hy(1:SIZE-1,q) = Hy(1:SIZE-1,q-1) + ( ( Ex(1:SIZE-1,q-1) - Ex(2:SIZE,q-1) ) * dt/(u0*dz) );
    % ABC for H at SIZE.
    %Hy(SIZE,q) = Hy(SIZE-1,q-1) + (Sc-1)/(Sc+1)*(Hy(SIZE-1,q) - Hy(SIZE,q-1) );
    
    % Calculation of Ex using updated difference equation for Ex. This is time step q+1/2.
    %Ex(2:SIZE,q) = Ex(2:SIZE, q-1) + ( dt/(e0*dz)*(Hy(1:SIZE-1, q) - Hy(2:SIZE, q)) );
    % ABC for E at 1.
    %Ex(1,q) = Ex(2,q-1) + (Sc-1)/(Sc+1)*(Ex(2,q) - Ex(2,q-1));
    
    % Recording absolute value
    if q > floor(MaxTime/Scalar)
    ExAbs = max(ExAbs, abs(Ex(:,q)));
    HyAbs = max(HyAbs, abs(Hy(:,q)));
    end
    
    % Source.
    if SourceChoice == 1
    Ex(source,q) = Ex(source,q) + exp( -1*((q-td)/(PulseWidth/4))^2 ) * Sc;
    elseif SourceChoice == 2
    Ex(source,q) = Ex(source,q) + sin(2*pi*f*(q)*dt) * Sc;
    elseif SourceChoice == 3
    Ex(source,q) = Ex(source,q) + (1-2*(pi*fp*(q*dt-dr))^2)*exp(-1*(pi*fp*(q*dt-dr))^2) * Sc;
    end
    
    % Activating a plane-wave source.
    %Ex(source,q) = Ex(source,q) + exp( -1*((q-td)/(PulseWidth/4))^2 ) * Sc;
end
toc
% Simulation animation.
for i=1:Scalar*4:MaxTime
    figure (2)
    subplot(211)
    hold off
    plot([(Partition)*dz/Ra (Partition)*dz/Ra], [-1.1 1.1], 'Color', 'r');
    hold on
    plot((0:SIZE-1)*dz/Ra,Ex(:,i))
    xlim([0 (SIZE-1)*dz/Ra])
    ylim([-1.1 1.1])
    xlabel('r/Ra')
    ylabel('Electric field (Ex)')
    
    subplot(212)
    hold off
    plot([(Partition)*dz/Ra (Partition)*dz/Ra], [-1.1 1.1], 'Color', 'r');
    hold on
    plot((0:SIZE-1)*dz/Ra,Hy(:,i))
    xlim([0 (SIZE-1)*dz/Ra])
    ylim([-1.1/imp0 1.1/imp0])
    xlabel('r/Ra')
    ylabel('Magnetic field (Hy)')
end

figure(3)
subplot(211)
hold off
plot([(Partition)*dz/Ra (Partition)*dz/Ra], [-1.1 1.1], 'Color', 'r');
hold on
plot((0:SIZE-1)*dz/Ra,ExAbs)
xlim([0 (SIZE-1)*dz/Ra])
ylim([-1.1 1.1])
xlabel('r/Ra')
ylabel('Electric field (Ex)')

subplot(212)
hold off
plot([(Partition)*dz/Ra (Partition)*dz/Ra], [-1.1 1.1], 'Color', 'r');
hold on
plot((0:SIZE-1)*dz/Ra,HyAbs)
xlim([0 (SIZE-1)*dz/Ra])
ylim([-1.1/imp0 1.1/imp0])
xlabel('r/Ra')
ylabel('Magnetic field (Hy)')

figure(4)
hold off
plot([(Partition)*dz/Ra (Partition)*dz/Ra], [-1.1 1.1], 'Color', 'r');
hold on
plot((0:SIZE-1)*dz/Ra,ExAbs./HyAbs)
xlim([0 (SIZE-1)*dz/Ra])
%ylim([-1.1/imp0 1.1/imp0])
xlabel('r/Ra')
ylabel('Magnitude of wave impedance')

% ================== Postprocessing =====================
% Reference: Chapter 5 section 5.6 from Understanding FDTD.
% for i=1:SIZE
%     if abs(Ex(i,t1)) < 1e-12 && Ex(i,t2) > 0
%         ExPhase(i) = -pi/2;
%     elseif abs(Ex(i,t1)) < 1e-12 && Ex(i,t2) > 0
%         ExPhase(i) = -pi/2;
%     else
%         ExPhase(i) = atan((cos(2*pi*f*t2*dt)-Ex(i,t2)/Ex(i,t1))/(sin(2*pi*f*t2*dt)));
%     end
%     if abs(Hy(i,t1)) < 1e-12 && Hy(i,t2) > 0
%         HyPhase(i) = -pi/2;
%     elseif abs(Hy(i,t1)) < 1e-12 && Hy(i,t2) > 0
%         HyPhase(i) = -pi/2;
%     else
%         HyPhase(i) = atan((cos(2*pi*f*t2*dt)-Hy(i,t2)/Hy(i,t1))/(sin(2*pi*f*t2*dt)));
%     end
% end
% 
% for i=1:SIZE
%     if abs(Ex(i,t1)) >= abs(Ex(i,t2))
%         ExAbs(i) = Ex(i,t1)/cos(ExPhase(i));
%     else
%         ExAbs(i) = Ex(i,t2)/(cos(2*pi*f*t2*dt)*cos(ExPhase(i))-sin(2*pi*f*t2*dt)*sin(ExPhase(i)));
%     end
%     if abs(Hy(i,t1)) >= abs(Hy(i,t2))
%         HyAbs(i) = Hy(i,t1)/cos(HyPhase(i));
%     else
%         HyAbs(i) = Hy(i,t2)/(cos(2*pi*f*t2*dt)*cos(HyPhase(i))-sin(2*pi*f*t2*dt)*sin(HyPhase(i)));
%     end
% end
% 
% for i=1:SIZE
%     if ExAbs(i) < 0 && ExPhase(i) >= 0
%         ExAbs(i) = -1*ExAbs(i);
%         ExPhase(i) = ExPhase(i)-pi;
%     elseif ExAbs(i) < 0 && ExPhase(i) < 0
%         ExAbs(i) = -1*ExAbs(i);
%         ExPhase(i) = ExPhase(i)+pi;
%     end
%     if HyAbs(i) < 0 && HyPhase(i) >= 0
%         HyAbs(i) = -1*HyAbs(i);
%         HyPhase(i) = HyPhase(i)-pi;
%     elseif HyAbs(i) < 0 && HyPhase(i) < 0
%         HyAbs(i) = -1*HyAbs(i);
%         HyPhase(i) = HyPhase(i)+pi;
%     end
% end

ExPhase = acos(Ex(:,MaxTime)./ExAbs);
HyPhase = acos(Hy(:,MaxTime)./HyAbs);
ExPhasor = ExAbs.*exp(1j.*ExPhase);
HyPhasor = HyAbs.*exp(1j.*HyPhase);

figure(5)
subplot(211)
hold off
plot([(Partition)*dz/Ra (Partition)*dz/Ra], [-1.1 1.1], 'Color', 'r');
hold on
plot((0:SIZE-1)*dz/Ra,ExAbs)
xlim([0 (SIZE-1)*dz/Ra])
ylim([-1.1 1.1])
xlabel('r/Ra')
ylabel('Magnitude of Electric field (Ex)')

subplot(212)
hold off
plot([(Partition)*dz/Ra (Partition)*dz/Ra], [-1.1 1.1], 'Color', 'r');
hold on
plot((0:SIZE-1)*dz/Ra,ExPhase)
xlim([0 (SIZE-1)*dz/Ra])
ylim([-pi pi])
xlabel('r/Ra')
ylabel('Phase')

figure(6)
subplot(211)
hold off
plot([(Partition)*dz/Ra (Partition)*dz/Ra], [-1.1 1.1], 'Color', 'r');
hold on
plot((0:SIZE-1)*dz/Ra,HyAbs)
xlim([0 (SIZE-1)*dz/Ra])
ylim([-1.1/imp0 1.1/imp0])
xlabel('r/Ra')
ylabel('Magnitude of Magnetic field (Hy)')

subplot(212)
hold off
plot([(Partition)*dz/Ra (Partition)*dz/Ra], [-1.1 1.1], 'Color', 'r');
hold on
plot((0:SIZE-1)*dz/Ra,HyPhase)
xlim([0 (SIZE-1)*dz/Ra])
ylim([-pi pi])
xlabel('r/Ra')
ylabel('Phase')

figure(7)
subplot(211)
hold off
plot([(Partition)*dz/Ra (Partition)*dz/Ra], [-1.1 1.1], 'Color', 'r');
hold on
plot((0:SIZE-1)*dz/Ra,real(ExPhasor./HyPhasor))
xlim([0 (SIZE-1)*dz/Ra])
%ylim([-1.1/imp0 1.1/imp0])
xlabel('r/Ra')
ylabel('Real part of wave impedance')

subplot(212)
hold off
plot([(Partition)*dz/Ra (Partition)*dz/Ra], [-1.1 1.1], 'Color', 'r');
hold on
plot((0:SIZE-1)*dz/Ra,imag(ExPhasor./HyPhasor))
xlim([0 (SIZE-1)*dz/Ra])
%ylim([-1.1/imp0 1.1/imp0])
xlabel('r/Ra')
ylabel('Imag part of wave impedance')