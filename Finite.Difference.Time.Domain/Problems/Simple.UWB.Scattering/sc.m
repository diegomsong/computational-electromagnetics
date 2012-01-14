% conductance
function rsc = sc (i, j)
[ISize JSize XCenter YCenter delta ra rb DT PMLw dtscalar Skin TissueIW TissueIIW PulseWidth] = Parameters;
rsc = 0;
abs = j-YCenter;
if (abs > 0 && abs <= Skin/delta)
    rsc = 0.0002;
end
if (abs > Skin/delta && abs <= TissueIW/delta)
    rsc = 0.02;
end
if (abs > TissueIW/delta && abs <= TissueIIW/delta)
    rsc = 0.036;
end
