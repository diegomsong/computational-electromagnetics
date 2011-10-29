function [Size XCenter YCenter delta ra rb DT] = Parameters
% Simulation related parameters.
Size = 100;
XCenter = (Size+1)/2;
YCenter = (Size-1)/2;
delta = 5e-3;

ra = 0.1;
rb = 0.2;
Cl = 3e8;
DT = delta / (4 * Cl);
