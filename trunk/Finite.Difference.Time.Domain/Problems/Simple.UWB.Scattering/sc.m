% conductance
function rsc = sc (i, j)
[ISize JSize XCenter YCenter delta ra rb DT PMLw dtscalar Skin TissueIW TissueIIW PulseWidth TumourX TumourY TumourRadius] = Parameters;
rsc = 0;
abs = j-YCenter;
if (abs > 0 && abs <= Skin/delta)
    rsc = 0.0;
end
if (abs > Skin/delta && abs <= (TissueIW+Skin)/delta)
    rsc = 0.02;
end
if (abs > (TissueIW+Skin)/delta && abs <= (TissueIIW+TissueIW+Skin)/delta)
    rsc = 0.036;
end
% tumour rsc = 0.899;
%(x-a)2 + (y-b)2 = (r-ra)2
if ( (i-TumourX/delta-XCenter)^2 + (j-TumourY/delta-YCenter)^2 < (TumourRadius/delta)^2 )
    rsc = 0.899;
end