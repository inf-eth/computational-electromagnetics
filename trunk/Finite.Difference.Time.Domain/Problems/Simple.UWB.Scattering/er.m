% relative permittivity
function rer = er (i, j)
[ISize JSize XCenter YCenter delta ra rb DT PMLw dtscalar Skin TissueIW TissueIIW PulseWidth TumourX TumourY TumourRadius] = Parameters;
rer = 1;
abs = j-YCenter;
if (abs > 0 && abs <= Skin/delta)
    rer = 1;
end
if (abs > Skin/delta && abs <= (TissueIW+Skin)/delta)
    rer = 4.108;
end
if (abs > (TissueIW+Skin)/delta && abs <= (TissueIIW+TissueIW+Skin)/delta)
    rer = 4.848;
end
%tumour rer = 60.89;
if ( (i-TumourX/delta-XCenter)^2 + (j-TumourY/delta-YCenter)^2 < (TumourRadius/delta)^2 )
    rer = 60.89;
end

