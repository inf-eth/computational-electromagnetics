% relative permittivity
function rer = er (i, j)
[ISize JSize XCenter YCenter delta ra rb DT PMLw dtscalar Skin TissueIW TissueIIW PulseWidth] = Parameters;
rer = 1;
abs = j-YCenter;
if (abs > 0 && abs <= Skin/delta)
    rer = 36;
end
if (abs > Skin/delta && abs <= TissueIW/delta)
    rer = 4.108;
end
if (abs > TissueIW/delta && abs <= TissueIIW/delta)
    rer = 4.848;
end

