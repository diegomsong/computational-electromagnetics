% Based on FDTD example 3.7 from Numerical techniques in Electromagnetics by Sadiku.
% Main m file for the FDTD simulation.
clc
clear all
% Simulation related parameters.
[Size XCenter YCenter delta ra rb] = Parameters;
IHx = Size;
JHx = Size-1;
IHy = Size+1;
JHy = Size;
IEz = Size;
JEz = Size;
NMax = 2; 
n0 = 1;
n1 = 2;
NNMax = 1000; % Maximum time.
ResolutionFactor = 5;  % E field snapshots will be saved every x frames where x is resolution factor.
%NHW = 40; % One half wave cycle.
Med = 2; % No of different media.
Js = 2; % J-position of the plane wave front.
% Different Constants.
%delta = 3e-3;
Cl = 3e8;
f = 2.0e9;
pi = 3.141592654;
e0 = (1e-9) / (36*pi);
u0 = (1e-7) * 4 * pi;
DT = delta / ( sqrt(2) * Cl );
TwoPIFDeltaT = 2 * pi * f * DT;
% Data arrays.
% CHx = zeros ( IHx, JHx ); % Conductance
% CHy = zeros ( IHy, JHy ); % Conductance
%REz = zeros ( IEz, JEz ); % Impedance
uxxHx = zeros ( IHx, JHx );  % uxx for Hx
uxyHy = zeros ( IHy, JHy );  % uxy for Hy
uyxHx = zeros ( IHx, JHx );  % uyx for Hx
uyyHy = zeros ( IHy, JHy );  % uyy for Hy
ezzEz = zeros ( IEz, JEz );  % ezz for Ez

% **** Temp ****
erEz = zeros ( IEz, JEz );  % ezz for Ez
% *************

RaEz = zeros ( IEz, JEz ); % Scaling parameter dependent on material conductivity.
smask = zeros ( IEz, JEz );

Bx = zeros ( IHx, JHx, 2 );
By = zeros ( IHy, JHy, 2 );

Hx = zeros ( IHx, JHx, 2 );
Hy = zeros ( IHy, JHy, 2 );

Dz = zeros ( IEz, JEz, 2 );
Ez = zeros ( IEz, JEz, 2 );
EzSnapshots = zeros ( IEz, JEz, NNMax/ResolutionFactor ); % E field Snapshot storage.

% ############ Initialization #############
fprintf ( 1, 'Initializing...' );
fprintf ( 1, '\nInitializing er array...' );
% Initializing er array.
for i=1:IEz
    for j=1:JEz
        ezzEz ( i, j ) = 1/(e0 * er ( i, j-0.5 ));
        erEz ( i, j ) = er( i, j-0.5 );
        %er ( i, j-0.5 )
    end
end
fprintf ( 1, '\nInitializing ur (Hx) array...' );
% Initializing the ur arrays.
for i=1:IHx
    for j=1:JHx
        invurHx = inv ( u0 * ur ( i, j-0.5 ));
        uxxHx (i, j) = invurHx(1, 1);
        uyxHx (i, j) = invurHx(2, 1);        
    end
end
fprintf ( 1, '\nInitializing ur (Hy) array...' );
for i=1:IHy
    for j=1:JHy
        invurHy = inv ( u0 * ur ( i-0.5, j-1 ));
        uxyHy (i, j) = invurHy(1, 2);
        uyyHy (i, j) = invurHy(2, 2);
    end
end
% Initializing the eta and 1/eta arrays.
% for i=1:IHx
%     for j=1:JHx
%         CHx ( i, j ) = DT / ( u0 * ur ( i, j-0.5 ) * delta );
%     end
% end
% 
% fprintf ( 1, '.' );
% for i=1:IHy
%     for j=1:JHy
%         CHy ( i, j ) = DT / ( u0 * ur ( i-0.5, j-1 ) * delta );
%     end
% end
% fprintf ( 1, '.' );
fprintf ( 1, '\nInitializing mask array...' );
for i=1:IEz
    for j=1:JEz
        REz ( i, j ) = DT / ( e0 * er ( i, j-0.5 ) * delta );
        %RaEz ( i, j ) = ( 1 - ( s(i, j-0.5) * DT )/( e0 * er( i, j-0.5 ) ) );
        RaEz ( i, j ) = 1;
        smask ( i, j ) = s ( i, j-0.5 );
%         if RaEz ( i, j ) < 0
%             RaEz ( i, j ) = -1 * RaEz ( i, j );
%         end
    end
end

figure (3)
mesh ( erEz )
title ( 'erEz' )
view (4, 4)
figure (5)
mesh ( u0 * uxxHx )
title ( 'u0 * uxxHx' )
view (4, 4)
figure (6)
mesh ( u0 * uxyHy )
title ( 'u0 * uxyHy' )
view (4, 4)
figure (7)
mesh ( u0 * uyyHy )
title ( 'u0 * uyyHy' )
view (4, 4)

%uxyHy == 0
% % max ( ezzEz )
% min ( ezzEz )
% average (ezzEz )
% REz
% RaEz
% uxxHx
% uxyHy
% uyxHx
% uyyHy

fprintf ( 1, 'done.\n' );
% ############ Initialization Complete ##############
% ########### 2. Now running the Simulation #############
fprintf ( 1, 'Simulation started... \n' );
for n=0:NNMax-2
    fprintf ( 1, '%g %% \n', (n*100)/NNMax );
    % *** Calculation of Magnetic Field Components ***
    % * Calculation of Bx.
    Bx ( :, :, n1 ) = Bx ( :, :, n0 ) + ( (DT/delta) * ( Ez ( :, 1:JHx, n0 ) - Ez ( :, 2:JHx+1, n0 ) ));

    % * Calculation of By.
    By ( 2:IHy-1, :, n1 ) = By ( 2:IHy-1, :, n0 ) + ( (DT/delta) * ( Ez ( 2:IHy-1, :, n0 ) - Ez ( 1:IHy-2, :,n0 ) ));
   
    % Boundary conditions on By. Soft grid truncation.
    By ( 1, 2:JHy-1, n1 ) = (1/3) * ( By ( 2, 1:JHy-2, n0 ) + By ( 2, 2:JHy-1, n0 ) + By ( 2, 3:JHy, n0 ) );
    By ( IHy, 2:JHy-1, n1 ) = (1/3) * ( By ( IHy-1, 1:JHy-2, n0 ) + By ( IHy-1, 2:JHy-1, n0 ) + By ( IHy-1, 3:JHy, n0 ) );
    By ( 1, 1, n1 ) = (1/2) * ( By ( 2, 1, n0 ) + By ( 2, 2, n0 ) );
    By ( 1, JHy, n1 ) = (1/2) * ( By ( 2, JHy, n0 ) + By ( 2, JHy-1, n0 ) );
    By ( IHy, 1, n1 ) = (1/2) * ( By ( IHy-1, 1, n0 ) + By( IHy-1, 2, n0 ) );
    By ( IHy, JHy, n1 ) = (1/2) * ( By ( IHy-1, JHy, n0 ) + By ( IHy-1, JHy-1, n0 ) );
    
    Hx ( :, :, n1 ) = uxxHx .* Bx ( :, :, n1 ) + uxyHy (1:IHy-1, 1:JHy-1 ) .* By (1:IHy-1, 1:JHy-1, n1 );
    Hy ( :, :, n1 ) = (1/u0) * By ( :, :, n1 );
    Hy (1:IHy-1, 1:JHy-1) = uyxHx.*Bx ( :, :, n1 ) + uyyHy (1:IHy-1, 1:JHy-1) .* By (1:IHy-1, 1:JHy-1, n1 );
    
    Dz ( :, 2:JEz-1, n1 ) = (  Dz ( :, 2:JEz-1, n0 ) ) + ( (DT/delta) * ( Hy ( 2:JEz+1, 2:JEz-1, n1 ) - Hy ( 1:JEz, 2:JEz-1, n1 ) + Hx ( :, 1:JEz-2,n1 ) - Hx ( :, 2:JEz-1, n1 ) ));
    %Dz (:, :, n1) = Dz (:, :, n1) .* smask;
    
    % Boundary conditions on Dz. Soft grid truncation.
    Dz ( 2:IEz-1, 1, n1 ) = (1/3) * ( Dz ( 1:IEz-2, 2, n0 ) + Dz ( 2:IEz-1, 2, n0 ) + Dz ( 3:IEz, 2, n0 ) );
    Dz ( 2:IEz-1, JEz, n1 ) = (1/3) * ( Dz ( 1:IEz-2, JEz-1, n0 ) + Dz ( 2:IEz-1, JEz-1, n0 ) + Dz ( 3:IEz, JEz-1, n0 ) );
    Dz ( 1, 1, n1 ) = (1/2) * ( Dz ( 1, 2, n0 ) + Dz ( 2, 2, n0 ) );
    Dz ( IEz, 1, n1 ) = (1/2) * ( Dz ( IEz, 2, n0 ) + Dz ( IEz-1, 2, n0 ) );
    Dz ( 1, JEz, n1 ) = (1/2) * ( Dz ( 1, JEz-1, n0 ) + Dz ( 2, JEz-1, n0 ) );
    Dz ( IEz, JEz, n1 ) = (1/2) * ( Dz ( IEz, JEz-1, n0 ) + Dz ( IEz-1, JEz-1, n0 ) );
    % ************************************************

    Ez ( :, :, n1 ) = (ezzEz) .* Dz ( :, :, n1 );
    % Applying a plane wave source at Js.
    % 1. Continuous source.
%     Ez ( :, Js, mod(n, 2)+2 ) = 1 * sin ( TwoPIFDeltaT * (mod(n, 2)+2) );
    % 2. Sinusoidal Source.
    Ez ( :, Js, n1 ) = 1 * sin ( TwoPIFDeltaT * n+2 );
    Dz ( :, Js, n1 ) = e0 * 1 * sin ( TwoPIFDeltaT * n+2 );

    if ( mod(n, ResolutionFactor) == 0)
        EzSnapshots ( :, :, n/ResolutionFactor + 1 ) = Ez ( :, :, n1);
    end
    temp = n0;
    n0 = n1;
    n1 = temp;    
end
fprintf ( 1, '100 %% \n' );
fprintf ( 1, 'Simulation complete! \n', 0 );

% Electric field snapshots.
for i=1:NNMax/ResolutionFactor-2
    
    figure (2)
    mesh ( EzSnapshots (:, :, i) );
    view (4, 4)
    zlim ( [-2 2] )
    
    figure (4)
    surf ( EzSnapshots (:, :, i) );
    view (0, 90)
    zlim ( [-10 10] )
    
end